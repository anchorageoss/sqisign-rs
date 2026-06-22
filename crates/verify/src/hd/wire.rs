//! Phase 6 - SQIsignHD signature / public-key wire-format (de)serialization.
//!
//! Parses the raw bytes a caller would hand the verifier into the structured
//! fields [`hd_verify_l1`] consumes, then runs the full self-contained
//! verification. Deserialization only - the format is already the
//! determinant-recovery-compressed form (only three of the four response
//! scalars are present; the fourth, `d` or `c`, is recovered during stage 3
//! from the determinant relation `a·d - b·c ≡ k·q`).
//!
//! # Layout (Level 1, little-endian throughout)
//!
//! **Signature (108 bytes):**
//!
//! | field        | bytes | notes |
//! |--------------|-------|-------|
//! | `A_com`      | 64    | `Fp2::encode` (re ‖ im, each 32 B); the commitment curve `E_com` |
//! | `q`          | 17    | response degree, `q < 2¹³⁶` |
//! | `a`          | 9     | response scalar, reduced mod `2^r` (`r = 70`); canonical ⇒ top 2 bits zero |
//! | `b`          | 9     | response scalar mod `2^r` |
//! | `c_or_d`     | 9     | the transmitted one of `c`/`d` (the other is determinant-recovered); the selector is `a`'s parity |
//!
//! **Public key (64 bytes):** `A_pk` (64) - the curve coefficient, with the
//! public-key basis hints packed into its spare bits (below).
//!
//! The challenge `chal` is **not** transmitted: the verifier recomputes it as
//! `SHAKE256(j(E_com) ‖ j(E_pk) ‖ message)` (the SQIsignHD Fiat-Shamir; Phase 5
//! confirmed our SHAKE matches the reference byte-for-byte). Tampering the
//! message therefore changes the recomputed challenge and breaks verification.
//!
//! # Basis hints packed into spare bits (108, beating the paper's 109)
//!
//! Phase 6 transmitted two `2ᶠ`-torsion basis-selection hints per curve
//! (`hint_com_P/Q` on the signature, `hint_pk_P/Q` on the public key) as two
//! extra bytes. Each hint is a small table index (`0..20` ⇒ 5 bits). Rather
//! than send (or recompute) them, we **pack them into bits that the field
//! encodings already leave zero**:
//!
//! * A canonical `Fp` value is `< p = 5·2²⁴⁸ - 1 < 2²⁵¹`, so the **top 5 bits
//!   of each component's most-significant byte are always zero** (the top byte
//!   is `≤ 4`). `A_com` (and `A_pk`) is two `Fp` components, giving **two 5-bit
//!   slots** - exactly the two hints.
//!
//! Packing (encode): `byte[31] |= hint_P << 3` (top byte of `re`),
//! `byte[63] |= hint_Q << 3` (top byte of `im`). Unpacking (parse):
//! `hint = byte >> 3`; then mask `byte &= 0b0000_0111` to restore the canonical
//! `Fp` byte before [`Fp2::decode`]. The low 3 bits hold the real top-byte
//! value (`≤ 4`) and never collide with the hint in bits 3-7.
//!
//! This keeps the signature **108 bytes** and the public key **64 bytes** while
//! delivering the hints to the verifier *for free* - no recomputation on the
//! hot path, no extra bytes, no tag byte. (108 also beats the paper's 109-byte
//! figure, which counts an explicit `c`/`d` selector bit; here that selector is
//! implicit in `a`'s parity.) The hints remain a deterministic function of the
//! curve, so they are still recoverable independently via
//! [`crate::hd::canonical_hints_l1`], used as a cross-check, not on the verify path.
//!
//! Because a wrong hint selects a basis inconsistent with the response, any
//! tampering of the packed bits makes verification fail; the canonical hint is
//! the only value that verifies, so the packing introduces no malleability.
//!
//! # Strictness (the dim-2 SEC lessons)
//!
//! * exact length (no trailing bytes, no truncation) → [`HdReject::MalformedInput`];
//! * `Fp2::decode` (after masking the hint bits) rejects non-canonical / out-of-range coordinates;
//! * response scalars must be the canonical `[0, 2^r)` representative (the top
//!   two bits of the 9th byte must be zero);
//! * all failures return `Err`; nothing panics.

use crypto_bigint::U256;

use crate::{Fp2, Level1};

use crate::hd::hd_verify::{hd_challenge_from_curves, hd_challenge_len, HdReject};
use crate::hd::self_contained::{hd_verify_l1, HdSignatureL1};

/// `Fp2EncodedBytes` at Level 1 (`2 × 32`).
pub const FP2_BYTES: usize = 64;
/// Bytes for the response degree `q` (`q < 2¹³⁶`).
pub const Q_BYTES: usize = 17;
/// Bytes for a response scalar reduced mod `2^r` (`r = 70`).
pub const SCALAR_BYTES: usize = 9;
/// `r` (the response modulus exponent): scalars live in `[0, 2^r)`.
const R_BITS: u32 = 70;

/// Wire size of a Level-1 signature (hints packed into `A_com`; see module docs).
pub const SIG_WIRE_BYTES: usize = FP2_BYTES + Q_BYTES + 3 * SCALAR_BYTES;
/// Wire size of a Level-1 public key (hints packed into `A_pk`).
pub const PK_WIRE_BYTES: usize = FP2_BYTES;

// basis-hint packing into the spare top bits of an Fp2 encoding
//
// A canonical `Fp` value is `< p < 2^251`, so the most-significant byte of each
// 32-byte component is `≤ 4`: bits 3..=7 are always zero and hold a 5-bit hint.

/// Byte offset of the top byte of the real component within an `Fp2` encoding.
const HINT_BYTE_RE: usize = 31;
/// Byte offset of the top byte of the imaginary component within an `Fp2` encoding.
const HINT_BYTE_IM: usize = 63;
/// Hints occupy bits 3..=7 of those bytes.
const HINT_SHIFT: u32 = 3;
/// Mask for the canonical low bits of the top byte (`p`'s top byte is 4).
const FP_TOP_CANON_MASK: u8 = (1 << HINT_SHIFT) - 1; // 0b0000_0111
/// Largest hint that fits in the 5 spare bits (the NQR tables have 20 entries).
const MAX_PACKED_HINT: u32 = (1 << (8 - HINT_SHIFT)) - 1; // 31

/// Pack two basis hints into the spare top bits of a 64-byte `Fp2` encoding
/// in place. Returns `None` if a hint does not fit in 5 bits.
fn pack_hints(fp2_bytes: &mut [u8], hp: u32, hq: u32) -> Option<()> {
    if hp > MAX_PACKED_HINT || hq > MAX_PACKED_HINT {
        return None;
    }
    fp2_bytes[HINT_BYTE_RE] |= (hp as u8) << HINT_SHIFT;
    fp2_bytes[HINT_BYTE_IM] |= (hq as u8) << HINT_SHIFT;
    Some(())
}

/// Extract the two packed hints from a 64-byte `Fp2` encoding and return a copy
/// with the hint bits masked back to their canonical (zero) state, ready for
/// [`Fp2::decode`].
fn unpack_hints(fp2_bytes: &[u8]) -> (u32, u32, [u8; FP2_BYTES]) {
    let hp = (fp2_bytes[HINT_BYTE_RE] >> HINT_SHIFT) as u32;
    let hq = (fp2_bytes[HINT_BYTE_IM] >> HINT_SHIFT) as u32;
    let mut clean = [0u8; FP2_BYTES];
    clean.copy_from_slice(&fp2_bytes[..FP2_BYTES]);
    clean[HINT_BYTE_RE] &= FP_TOP_CANON_MASK;
    clean[HINT_BYTE_IM] &= FP_TOP_CANON_MASK;
    (hp, hq, clean)
}

/// A parsed signature (structured fields, no challenge - that is recomputed).
#[derive(Clone, Debug)]
pub struct ParsedSignature {
    pub a_com: Fp2<Level1>,
    pub a: i128,
    pub b: i128,
    pub c_or_d: i128,
    pub q: U256,
    pub hint_com_p: u32,
    pub hint_com_q: u32,
}

/// A parsed public key.
#[derive(Clone, Debug)]
pub struct ParsedPublicKey {
    pub a_pk: Fp2<Level1>,
    pub hint_pk_p: u32,
    pub hint_pk_q: u32,
}

/// Read a 17-byte little-endian value into `U256` (`q < 2¹³⁶`, so words 3 and
/// the high bytes of word 2 are zero).
fn read_q(b: &[u8]) -> U256 {
    let w0 = u64::from_le_bytes(b[0..8].try_into().unwrap());
    let w1 = u64::from_le_bytes(b[8..16].try_into().unwrap());
    let w2 = b[16] as u64;
    U256::from_words([w0, w1, w2, 0])
}

/// Encode `q` (`< 2¹³⁶`) into 17 little-endian bytes. Returns `None` if `q`
/// does not fit (a programming error on the encode side).
fn write_q(q: &U256, out: &mut [u8]) -> Option<()> {
    let w = q.to_words();
    if w[3] != 0 || w[2] > 0xFF {
        return None;
    }
    out[0..8].copy_from_slice(&w[0].to_le_bytes());
    out[8..16].copy_from_slice(&w[1].to_le_bytes());
    out[16] = w[2] as u8;
    Some(())
}

/// Read a canonical response scalar: 9 little-endian bytes holding a value in
/// `[0, 2^r)` (`r = 70`). The top two bits of the 9th byte must be zero.
fn read_scalar(b: &[u8]) -> Result<i128, HdReject> {
    // r = 70 = 8*8 + 6, so byte 8 holds bits 64..69; bits 70,71 must be zero.
    if b[8] >= (1u8 << (R_BITS as usize - 64)) {
        return Err(HdReject::MalformedInput);
    }
    let mut v = 0u128;
    for (i, &byte) in b.iter().enumerate().take(SCALAR_BYTES) {
        v |= (byte as u128) << (8 * i);
    }
    // v < 2^70 < 2^127, so the cast is exact and non-negative.
    Ok(v as i128)
}

/// Encode a response scalar reduced to its canonical `[0, 2^r)` representative.
fn write_scalar(s: i128, out: &mut [u8]) {
    let v = s.rem_euclid(1i128 << R_BITS) as u128; // [0, 2^r)
    for (i, b) in out.iter_mut().enumerate().take(SCALAR_BYTES) {
        *b = (v >> (8 * i)) as u8;
    }
}

/// Parse a Level-1 signature from its wire bytes.
pub fn parse_signature(bytes: &[u8]) -> Result<ParsedSignature, HdReject> {
    if bytes.len() != SIG_WIRE_BYTES {
        return Err(HdReject::MalformedInput);
    }
    let mut pos = 0;
    // The commitment basis hints are packed into A_com's spare top bits; unpack
    // them first, then decode the masked (canonical) A_com.
    let (hint_com_p, hint_com_q, a_com_clean) = unpack_hints(&bytes[pos..pos + FP2_BYTES]);
    let a_com = Fp2::<Level1>::decode(&a_com_clean).ok_or(HdReject::MalformedInput)?;
    pos += FP2_BYTES;
    let q = read_q(&bytes[pos..pos + Q_BYTES]);
    pos += Q_BYTES;
    let a = read_scalar(&bytes[pos..pos + SCALAR_BYTES])?;
    pos += SCALAR_BYTES;
    let b = read_scalar(&bytes[pos..pos + SCALAR_BYTES])?;
    pos += SCALAR_BYTES;
    let c_or_d = read_scalar(&bytes[pos..pos + SCALAR_BYTES])?;
    pos += SCALAR_BYTES;
    debug_assert_eq!(pos, SIG_WIRE_BYTES);
    Ok(ParsedSignature {
        a_com,
        a,
        b,
        c_or_d,
        q,
        hint_com_p,
        hint_com_q,
    })
}

/// Parse a Level-1 public key from its wire bytes.
pub fn parse_public_key(bytes: &[u8]) -> Result<ParsedPublicKey, HdReject> {
    if bytes.len() != PK_WIRE_BYTES {
        return Err(HdReject::MalformedInput);
    }
    // The public-key basis hints are packed into A_pk's spare top bits.
    let (hint_pk_p, hint_pk_q, a_pk_clean) = unpack_hints(&bytes[0..FP2_BYTES]);
    let a_pk = Fp2::<Level1>::decode(&a_pk_clean).ok_or(HdReject::MalformedInput)?;
    Ok(ParsedPublicKey {
        a_pk,
        hint_pk_p,
        hint_pk_q,
    })
}

/// Serialize a signature to its `SIG_WIRE_BYTES` wire form. Scalars are reduced
/// to their canonical `[0, 2^r)` representative. The two commitment basis hints
/// are packed into `A_com`'s spare top bits (see the module docs). Returns
/// `None` if `q ≥ 2¹³⁶` or a hint does not fit in 5 bits.
pub fn encode_signature(
    a_com: &Fp2<Level1>,
    a: i128,
    b: i128,
    c_or_d: i128,
    q: &U256,
    hint_com_p: u32,
    hint_com_q: u32,
) -> Option<[u8; SIG_WIRE_BYTES]> {
    let mut out = [0u8; SIG_WIRE_BYTES];
    let mut pos = 0;
    out[pos..pos + FP2_BYTES].copy_from_slice(&a_com.encode());
    pack_hints(&mut out[pos..pos + FP2_BYTES], hint_com_p, hint_com_q)?;
    pos += FP2_BYTES;
    write_q(q, &mut out[pos..pos + Q_BYTES])?;
    pos += Q_BYTES;
    write_scalar(a, &mut out[pos..pos + SCALAR_BYTES]);
    pos += SCALAR_BYTES;
    write_scalar(b, &mut out[pos..pos + SCALAR_BYTES]);
    pos += SCALAR_BYTES;
    write_scalar(c_or_d, &mut out[pos..pos + SCALAR_BYTES]);
    pos += SCALAR_BYTES;
    debug_assert_eq!(pos, SIG_WIRE_BYTES);
    Some(out)
}

/// Serialize a public key to its `PK_WIRE_BYTES` wire form: the curve
/// coefficient with the two basis hints packed into `A_pk`'s spare top bits.
/// Returns `None` if a hint does not fit in 5 bits.
pub fn encode_public_key(
    a_pk: &Fp2<Level1>,
    hint_pk_p: u32,
    hint_pk_q: u32,
) -> Option<[u8; PK_WIRE_BYTES]> {
    let mut out = [0u8; PK_WIRE_BYTES];
    out[0..FP2_BYTES].copy_from_slice(&a_pk.encode());
    pack_hints(&mut out[0..FP2_BYTES], hint_pk_p, hint_pk_q)?;
    Some(out)
}

/// The raw-bytes verification entry point: parse `(signature, public_key)`,
/// recompute the challenge from the curves + message, and run the full
/// self-contained Level-1 verification. Returns `Ok(())` on accept, or the
/// rejection reason ([`HdReject::MalformedInput`] for unparseable bytes).
pub fn hd_verify_bytes_l1(
    signature: &[u8],
    public_key: &[u8],
    message: &[u8],
) -> Result<(), HdReject> {
    let sig = parse_signature(signature)?;
    let pk = parse_public_key(public_key)?;
    hd_verify_l1_parsed(&sig, &pk, message)
}

/// Verify already-parsed `(signature, public_key)` structures against a message.
///
/// This is the post-parse half of [`hd_verify_bytes_l1`]: it recomputes the
/// Fiat-Shamir challenge from the parsed curves and runs the self-contained
/// verification. It lets a higher-level dispatcher (e.g. the unified
/// `AnySignature` autodetect in `sqisign-rs`) reuse the parsed structures it
/// already holds without re-serializing to bytes.
pub fn hd_verify_l1_parsed(
    sig: &ParsedSignature,
    pk: &ParsedPublicKey,
    message: &[u8],
) -> Result<(), HdReject> {
    // Recompute the Fiat-Shamir challenge from the parsed curves + message.
    let n = hd_challenge_len::<Level1>();
    let mut chal = [0u8; 64];
    if !hd_challenge_from_curves(&sig.a_com, &pk.a_pk, message, &mut chal[..n]) {
        return Err(HdReject::BadCurve);
    }
    let mut chal_limbs = [0u64; 4];
    for (i, limb) in chal_limbs.iter_mut().enumerate() {
        *limb = u64::from_le_bytes(chal[i * 8..i * 8 + 8].try_into().unwrap());
    }

    let hdsig = HdSignatureL1 {
        a_pk: pk.a_pk.clone(),
        a_com: sig.a_com.clone(),
        hint_pk_p: pk.hint_pk_p,
        hint_pk_q: pk.hint_pk_q,
        hint_com_p: sig.hint_com_p,
        hint_com_q: sig.hint_com_q,
        message,
        chal_limbs: &chal_limbs,
        claimed_chal: &chal[..n],
        resp_a: sig.a,
        resp_b: sig.b,
        resp_c_or_d: sig.c_or_d,
        q: sig.q,
    };
    hd_verify_l1(&hdsig)
}

/// Convenience: `true` iff the raw-bytes signature verifies.
#[inline]
pub fn hd_verify_bytes_l1_bool(signature: &[u8], public_key: &[u8], message: &[u8]) -> bool {
    hd_verify_bytes_l1(signature, public_key, message).is_ok()
}
