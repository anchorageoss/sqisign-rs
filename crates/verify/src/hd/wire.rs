//! The compact wire format and the byte-level verification entry.
//!
//! A signature is one bit string, least significant bit first, padded with
//! zeros to a whole byte:
//!
//! ```text
//! Re(A_com) | Im(A_com) | (q − 3)/4 | a | b | c_or_d | hint_com_P | hint_com_Q
//!  ⌈log2 p⌉   ⌈log2 p⌉     e − 2      r   r     r          8            8      bits
//! ```
//!
//! `2 ⌈log2 p⌉ + (e − 2) + 3 r + 16` bits: 1107, 1690 and 2223 bits, that
//! is 139, 212 and 278 bytes at levels I, III and V. `q ≡ 3 (mod 4)` for
//! every good degree, so its two low bits are not sent. Each scalar is the
//! representative in `[0, 2^r)`; the selector between `c` and `d` is `a`'s
//! parity. The hints are the library's table indices.
//!
//! A public key is `A_pk | hint_pk_P | hint_pk_Q` with `A_pk` the canonical
//! `F_p^2` encoding in whole bytes: 84, 130 and 170 bytes.
//!
//! The decoder rejects a wrong length, a coordinate not below `p` and a
//! nonzero padding bit; the challenge is recomputed, never transmitted.

use crate::fp::{Fp2, FpBackend};

use super::hd_verify::{hd_challenge, hd_challenge_len, HdReject, MAX_CHAL_BYTES};
use super::params::{
    pk_wire_bytes, sig_wire_bits, sig_wire_bytes, HdLevel, MAX_FP2_BYTES, MAX_PK_WIRE_BYTES,
    MAX_SIG_WIRE_BYTES,
};
use super::self_contained::{hd_verify_with, HdSignature, HdVerifyOptions, HdVerifyReport};
use super::uint::{NL, U6};
use super::wide::{W, W_LIMBS};

/// A decoded compact signature.
#[derive(Clone, Debug)]
pub struct ParsedSignature<L: FpBackend> {
    pub a_com: Fp2<L>,
    pub a: W,
    pub b: W,
    pub c_or_d: W,
    pub q: U6,
    pub hint_com_p: u32,
    pub hint_com_q: u32,
}

/// A decoded compact public key.
#[derive(Clone, Debug)]
pub struct ParsedPublicKey<L: FpBackend> {
    pub a_pk: Fp2<L>,
    pub hint_pk_p: u32,
    pub hint_pk_q: u32,
}

/// `n` bits of `src` from bit `*pos`, least significant first, into the
/// low bits of `dst` (little-endian bytes, the rest zero); advances `*pos`.
fn get_bits(src: &[u8], pos: &mut usize, n: usize, dst: &mut [u8]) {
    dst.iter_mut().for_each(|b| *b = 0);
    let (mut byte, shift) = (*pos / 8, (*pos % 8) as u32);
    let full = n / 8;
    for d in dst.iter_mut().take(full) {
        let lo = src[byte] >> shift;
        let hi = if shift == 0 {
            0
        } else {
            src[byte + 1] << (8 - shift)
        };
        *d = lo | hi;
        byte += 1;
    }
    let rem = (n % 8) as u32;
    if rem != 0 {
        let lo = src[byte] >> shift;
        let hi = if shift + rem > 8 {
            src[byte + 1] << (8 - shift)
        } else {
            0
        };
        dst[full] = (lo | hi) & ((1u8 << rem) - 1);
    }
    *pos += n;
}

/// The low `n` bits of `src` (little-endian bytes) into `dst` at bit
/// `*pos`, least significant first; `dst` is zero from `*pos` on. Advances
/// `*pos`.
fn put_bits(dst: &mut [u8], pos: &mut usize, n: usize, src: &[u8]) {
    let (mut byte, shift) = (*pos / 8, (*pos % 8) as u32);
    let full = n / 8;
    for s in src.iter().take(full) {
        dst[byte] |= s << shift;
        if shift != 0 {
            dst[byte + 1] |= s >> (8 - shift);
        }
        byte += 1;
    }
    let rem = (n % 8) as u32;
    if rem != 0 {
        let s = src[full] & ((1u8 << rem) - 1);
        dst[byte] |= s << shift;
        if shift + rem > 8 {
            dst[byte + 1] |= s >> (8 - shift);
        }
    }
    *pos += n;
}

/// Decode a signature (strict).
pub fn parse_signature<L: HdLevel>(bytes: &[u8]) -> Result<ParsedSignature<L>, HdReject> {
    if bytes.len() != sig_wire_bytes::<L>() {
        return Err(HdReject::MalformedInput);
    }
    let fp = L::FP2_BYTES / 2;
    let p_bits = L::P_BITS as usize;
    let mut pos = 0usize;

    // A_com: two integers of ⌈log2 p⌉ bits, into the canonical encoding;
    // the field decoder rejects a coordinate not below p.
    let mut enc = [0u8; MAX_FP2_BYTES];
    get_bits(bytes, &mut pos, p_bits, &mut enc[..fp]);
    get_bits(bytes, &mut pos, p_bits, &mut enc[fp..2 * fp]);
    let a_com = Fp2::<L>::decode(&enc[..2 * fp]).ok_or(HdReject::MalformedInput)?;

    // (q − 3)/4 in e − 2 bits: q = 4 x + 3 < 2^e.
    let mut buf = [0u8; 8 * NL];
    get_bits(bytes, &mut pos, (L::E_EMBED - 2) as usize, &mut buf);
    let q = U6::from_le_bytes(&buf)
        .shl(2)
        .wrapping_add(&U6::from_u64(3));

    let mut scalar = || {
        let mut buf = [0u8; 8 * W_LIMBS];
        get_bits(bytes, &mut pos, L::R as usize, &mut buf);
        W::from_le_bytes(&buf)
    };
    let a = scalar();
    let b = scalar();
    let c_or_d = scalar();

    let mut h = [0u8; 2];
    get_bits(bytes, &mut pos, 16, &mut h);

    // the padding: every bit past the fields is zero
    debug_assert_eq!(pos, sig_wire_bits::<L>());
    let last = bytes[bytes.len() - 1];
    if pos % 8 != 0 && last >> (pos % 8) != 0 {
        return Err(HdReject::MalformedInput);
    }
    Ok(ParsedSignature {
        a_com,
        a,
        b,
        c_or_d,
        q,
        hint_com_p: h[0] as u32,
        hint_com_q: h[1] as u32,
    })
}

/// Decode a public key (strict).
pub fn parse_public_key<L: HdLevel>(bytes: &[u8]) -> Result<ParsedPublicKey<L>, HdReject> {
    if bytes.len() != pk_wire_bytes::<L>() {
        return Err(HdReject::MalformedInput);
    }
    let a_pk = Fp2::<L>::decode(&bytes[..L::FP2_BYTES]).ok_or(HdReject::MalformedInput)?;
    Ok(ParsedPublicKey {
        a_pk,
        hint_pk_p: bytes[L::FP2_BYTES] as u32,
        hint_pk_q: bytes[L::FP2_BYTES + 1] as u32,
    })
}

/// Encode a signature into `out` (at least [`sig_wire_bytes`] long);
/// returns the length. `None` if `q` is not below `2^e` or not congruent
/// to 3 modulo 4 (no good degree is either), or a hint exceeds a byte.
/// Scalars are reduced modulo `2^r`.
#[allow(clippy::too_many_arguments)]
pub fn encode_signature<L: HdLevel>(
    a_com: &Fp2<L>,
    a: W,
    b: W,
    c_or_d: W,
    q: &U6,
    hint_com_p: u32,
    hint_com_q: u32,
    out: &mut [u8],
) -> Option<usize> {
    let n = sig_wire_bytes::<L>();
    if out.len() < n || hint_com_p > 255 || hint_com_q > 255 {
        return None;
    }
    if q.bits() > L::E_EMBED || q.0[0] & 3 != 3 {
        return None;
    }
    let fp = L::FP2_BYTES / 2;
    let p_bits = L::P_BITS as usize;
    out[..n].iter_mut().for_each(|b| *b = 0);
    let mut pos = 0usize;
    let enc = a_com.encode();
    let enc: &[u8] = enc.as_ref();
    put_bits(out, &mut pos, p_bits, &enc[..fp]);
    put_bits(out, &mut pos, p_bits, &enc[fp..2 * fp]);
    put_bits(
        out,
        &mut pos,
        (L::E_EMBED - 2) as usize,
        &q.shr(2).to_le_bytes(),
    );
    let mask = W::mask(L::R);
    for s in [a, b, c_or_d] {
        put_bits(out, &mut pos, L::R as usize, &(s & mask).to_le_bytes());
    }
    put_bits(out, &mut pos, 16, &[hint_com_p as u8, hint_com_q as u8]);
    debug_assert_eq!(pos, sig_wire_bits::<L>());
    Some(n)
}

/// Encode a public key into `out` (at least [`pk_wire_bytes`] long).
pub fn encode_public_key<L: HdLevel>(
    a_pk: &Fp2<L>,
    hint_pk_p: u32,
    hint_pk_q: u32,
    out: &mut [u8],
) -> Option<usize> {
    let n = pk_wire_bytes::<L>();
    if out.len() < n || hint_pk_p > 255 || hint_pk_q > 255 {
        return None;
    }
    out[..L::FP2_BYTES].copy_from_slice(a_pk.encode().as_ref());
    out[L::FP2_BYTES] = hint_pk_p as u8;
    out[L::FP2_BYTES + 1] = hint_pk_q as u8;
    Some(n)
}

/// The challenge of a parsed signature under a parsed key, as limbs.
pub fn challenge_limbs<L: HdLevel>(
    sig: &ParsedSignature<L>,
    pk: &ParsedPublicKey<L>,
    message: &[u8],
) -> Option<[u64; 4]> {
    let mut pkb = [0u8; MAX_PK_WIRE_BYTES];
    let n = encode_public_key::<L>(&pk.a_pk, pk.hint_pk_p, pk.hint_pk_q, &mut pkb)?;
    let len = hd_challenge_len::<L>();
    let mut chal = [0u8; MAX_CHAL_BYTES];
    hd_challenge::<L>(&pkb[..n], &sig.a_com, message, &mut chal[..len]);
    let mut limbs = [0u64; 4];
    for (i, limb) in limbs.iter_mut().enumerate() {
        *limb = u64::from_le_bytes(chal[i * 8..i * 8 + 8].try_into().unwrap());
    }
    Some(limbs)
}

/// Verify a parsed signature with the default options (the four-point
/// verifier).
pub fn hd_verify_parsed<L: HdLevel>(
    sig: &ParsedSignature<L>,
    pk: &ParsedPublicKey<L>,
    message: &[u8],
) -> Result<(), HdReject> {
    hd_verify_parsed_with::<L>(sig, pk, message, &HdVerifyOptions::default()).0
}

/// Verify a parsed signature with options: the verdict and the run's report
/// (see [`hd_verify_with`]).
pub fn hd_verify_parsed_with<L: HdLevel>(
    sig: &ParsedSignature<L>,
    pk: &ParsedPublicKey<L>,
    message: &[u8],
    opts: &HdVerifyOptions,
) -> (Result<(), HdReject>, HdVerifyReport) {
    let chal = match challenge_limbs::<L>(sig, pk, message) {
        Some(c) => c,
        None => return (Err(HdReject::BadCurve), HdVerifyReport::default()),
    };
    let hdsig = HdSignature {
        a_pk: pk.a_pk.clone(),
        a_com: sig.a_com.clone(),
        hint_pk_p: pk.hint_pk_p,
        hint_pk_q: pk.hint_pk_q,
        hint_com_p: sig.hint_com_p,
        hint_com_q: sig.hint_com_q,
        chal_limbs: &chal,
        resp_a: sig.a,
        resp_b: sig.b,
        resp_c_or_d: sig.c_or_d,
        q: sig.q,
    };
    hd_verify_with::<L>(&hdsig, opts)
}

/// [`hd_verify_bytes`] without the canonical-hint check. Measurement only
/// (feature `instrument`); not a verifier.
#[cfg(feature = "instrument")]
pub fn hd_verify_bytes_without_hint_check<L: HdLevel>(
    signature: &[u8],
    public_key: &[u8],
    message: &[u8],
) -> Result<(), HdReject> {
    let sig = parse_signature::<L>(signature)?;
    let pk = parse_public_key::<L>(public_key)?;
    let chal = challenge_limbs::<L>(&sig, &pk, message).ok_or(HdReject::BadCurve)?;
    let hdsig = HdSignature {
        a_pk: pk.a_pk.clone(),
        a_com: sig.a_com.clone(),
        hint_pk_p: pk.hint_pk_p,
        hint_pk_q: pk.hint_pk_q,
        hint_com_p: sig.hint_com_p,
        hint_com_q: sig.hint_com_q,
        chal_limbs: &chal,
        resp_a: sig.a,
        resp_b: sig.b,
        resp_c_or_d: sig.c_or_d,
        q: sig.q,
    };
    crate::hd::self_contained::hd_verify_without_hint_check::<L>(
        &hdsig,
        &HdVerifyOptions::default(),
    )
}

/// Verify from bytes with the default options (the four-point verifier).
pub fn hd_verify_bytes<L: HdLevel>(
    signature: &[u8],
    public_key: &[u8],
    message: &[u8],
) -> Result<(), HdReject> {
    hd_verify_bytes_with::<L>(signature, public_key, message, &HdVerifyOptions::default()).0
}

/// Verify from bytes with options: the verdict and the run's report (see
/// [`hd_verify_with`]). A malformed input is rejected before any half-chain,
/// with an empty report.
pub fn hd_verify_bytes_with<L: HdLevel>(
    signature: &[u8],
    public_key: &[u8],
    message: &[u8],
    opts: &HdVerifyOptions,
) -> (Result<(), HdReject>, HdVerifyReport) {
    let sig = match parse_signature::<L>(signature) {
        Ok(s) => s,
        Err(e) => return (Err(e), HdVerifyReport::default()),
    };
    let pk = match parse_public_key::<L>(public_key) {
        Ok(p) => p,
        Err(e) => return (Err(e), HdVerifyReport::default()),
    };
    hd_verify_parsed_with(&sig, &pk, message, opts)
}

/// The wire size of a signature at this level.
pub const fn signature_bytes<L: HdLevel>() -> usize {
    sig_wire_bytes::<L>()
}

/// The wire size of a public key at this level.
pub const fn public_key_bytes<L: HdLevel>() -> usize {
    pk_wire_bytes::<L>()
}

const _: () = assert!(MAX_SIG_WIRE_BYTES >= sig_wire_bytes::<crate::params::P324_3>());
const _: () = assert!(MAX_SIG_WIRE_BYTES >= sig_wire_bytes::<crate::params::P500_27>());
