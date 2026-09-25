//! The compact wire format and the byte-level verification entry.
//!
//! ```text
//! signature  = A_com | q | a | b | c_or_d | hint_com_P | hint_com_Q
//! public key = A_pk | hint_pk_P | hint_pk_Q
//! ```
//!
//! `A` is the canonical `F_p^2` encoding (`2 FP` bytes), `q` the response
//! degree in `⌈e/8⌉` bytes little-endian (`q < 2^e`), each scalar the
//! canonical representative in `[0, 2^r)` in `⌈r/8⌉` bytes (the bits above
//! `r` must be zero), and the hints one byte each (the library's table
//! indices; values past the table are the `x + i` fallback). The selector
//! between `c` and `d` is `a`'s parity. Sizes at level I: 142 and 84 bytes
//! (`docs/COMPACT_R3.md`).
//!
//! The decoder is strict about lengths, field elements and the scalars'
//! range; the challenge is recomputed, never transmitted.

use crate::fp::{Fp2, FpBackend};

use super::hd_verify::{hd_challenge, hd_challenge_len, HdReject, MAX_CHAL_BYTES};
use super::params::{
    pk_wire_bytes, sig_wire_bytes, HdLevel, MAX_PK_WIRE_BYTES, MAX_SIG_WIRE_BYTES,
};
use super::self_contained::{hd_verify, HdSignature};
use super::uint::U4;

/// A decoded compact signature.
#[derive(Clone, Debug)]
pub struct ParsedSignature<L: FpBackend> {
    pub a_com: Fp2<L>,
    pub a: i128,
    pub b: i128,
    pub c_or_d: i128,
    pub q: U4,
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

fn read_q<L: HdLevel>(b: &[u8]) -> Result<U4, HdReject> {
    let q = U4::from_le_bytes(&b[..L::Q_BYTES]);
    if q.bits() > L::E_EMBED {
        return Err(HdReject::MalformedInput);
    }
    Ok(q)
}

fn write_q<L: HdLevel>(q: &U4, out: &mut [u8]) -> Option<()> {
    let bytes = q.to_le_bytes();
    if bytes[L::Q_BYTES..].iter().any(|&x| x != 0) {
        return None;
    }
    out[..L::Q_BYTES].copy_from_slice(&bytes[..L::Q_BYTES]);
    Some(())
}

fn read_scalar<L: HdLevel>(b: &[u8]) -> Result<i128, HdReject> {
    let mut v = 0u128;
    for (i, &byte) in b.iter().enumerate().take(L::SCALAR_BYTES) {
        v |= (byte as u128) << (8 * i);
    }
    if v >> L::R != 0 {
        return Err(HdReject::MalformedInput);
    }
    Ok(v as i128)
}

fn write_scalar<L: HdLevel>(s: i128, out: &mut [u8]) {
    let v = s.rem_euclid(1i128 << L::R) as u128;
    for (i, b) in out.iter_mut().enumerate().take(L::SCALAR_BYTES) {
        *b = (v >> (8 * i)) as u8;
    }
}

/// Decode a signature (strict).
pub fn parse_signature<L: HdLevel>(bytes: &[u8]) -> Result<ParsedSignature<L>, HdReject> {
    if bytes.len() != sig_wire_bytes::<L>() {
        return Err(HdReject::MalformedInput);
    }
    let mut pos = 0;
    let a_com =
        Fp2::<L>::decode(&bytes[pos..pos + L::FP2_BYTES]).ok_or(HdReject::MalformedInput)?;
    pos += L::FP2_BYTES;
    let q = read_q::<L>(&bytes[pos..pos + L::Q_BYTES])?;
    pos += L::Q_BYTES;
    let a = read_scalar::<L>(&bytes[pos..pos + L::SCALAR_BYTES])?;
    pos += L::SCALAR_BYTES;
    let b = read_scalar::<L>(&bytes[pos..pos + L::SCALAR_BYTES])?;
    pos += L::SCALAR_BYTES;
    let c_or_d = read_scalar::<L>(&bytes[pos..pos + L::SCALAR_BYTES])?;
    pos += L::SCALAR_BYTES;
    let hint_com_p = bytes[pos] as u32;
    let hint_com_q = bytes[pos + 1] as u32;
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
/// returns the length, `None` if `q` does not fit or a hint exceeds a byte.
pub fn encode_signature<L: HdLevel>(
    a_com: &Fp2<L>,
    a: i128,
    b: i128,
    c_or_d: i128,
    q: &U4,
    hint_com_p: u32,
    hint_com_q: u32,
    out: &mut [u8],
) -> Option<usize> {
    let n = sig_wire_bytes::<L>();
    if out.len() < n || hint_com_p > 255 || hint_com_q > 255 {
        return None;
    }
    let mut pos = 0;
    out[pos..pos + L::FP2_BYTES].copy_from_slice(a_com.encode().as_ref());
    pos += L::FP2_BYTES;
    write_q::<L>(q, &mut out[pos..pos + L::Q_BYTES])?;
    pos += L::Q_BYTES;
    for s in [a, b, c_or_d] {
        write_scalar::<L>(s, &mut out[pos..pos + L::SCALAR_BYTES]);
        pos += L::SCALAR_BYTES;
    }
    out[pos] = hint_com_p as u8;
    out[pos + 1] = hint_com_q as u8;
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

/// Verify a parsed signature.
pub fn hd_verify_parsed<L: HdLevel>(
    sig: &ParsedSignature<L>,
    pk: &ParsedPublicKey<L>,
    message: &[u8],
) -> Result<(), HdReject> {
    let chal = challenge_limbs::<L>(sig, pk, message).ok_or(HdReject::BadCurve)?;
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
    hd_verify::<L>(&hdsig)
}

/// Verify from bytes.
pub fn hd_verify_bytes<L: HdLevel>(
    signature: &[u8],
    public_key: &[u8],
    message: &[u8],
) -> Result<(), HdReject> {
    let sig = parse_signature::<L>(signature)?;
    let pk = parse_public_key::<L>(public_key)?;
    hd_verify_parsed(&sig, &pk, message)
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
