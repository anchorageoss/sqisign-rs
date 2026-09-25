//! The compact format's challenge and the recovery of the dropped response
//! scalar.
//!
//! The challenge is the round-3 hash shape applied to the compact key:
//! `SHAKE256("SQI" || pk_bytes || A_com || msg)`, the first `λ/8` bytes as a
//! little-endian scalar (the specification's `hash_to_challenge` with the
//! compact public key in place of the round-3 one). It is not transmitted:
//! the verifier recomputes it from the transmitted commitment curve.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use crate::fp::{Fp2, FpBackend};

use super::params::HdLevel;

/// The largest challenge of the implemented levels, in bytes.
pub const MAX_CHAL_BYTES: usize = 32;

/// Bytes of the challenge at this level, `λ/8`.
#[inline]
pub fn hd_challenge_len<L: HdLevel>() -> usize {
    L::CHALLENGE_BYTES
}

/// `SHAKE256("SQI" || pk_bytes || A_com || msg)`, `out.len()` bytes.
#[inline]
pub fn hd_challenge<L: FpBackend>(pk_bytes: &[u8], a_com: &Fp2<L>, message: &[u8], out: &mut [u8]) {
    let mut hasher = Shake256::default();
    hasher.update(b"SQI");
    hasher.update(pk_bytes);
    hasher.update(a_com.encode().as_ref());
    hasher.update(message);
    hasher.finalize_xof().read(out);
}

/// Why a compact verification rejects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HdReject {
    /// A curve coefficient that is not a valid Montgomery curve.
    BadCurve,
    /// The challenge does not recompute (unused on the wire path, kept for
    /// callers that carry a challenge).
    ChallengeMismatch,
    /// The challenge isogeny or its basis could not be recovered.
    ChallengeRecovery,
    /// The response images could not be recovered (the pairing discrete
    /// logarithm failed: the commitment and challenge bases are inconsistent).
    ResponseRecovery,
    /// `2^e − q` is not a sum of two squares the verifier can find.
    NormEquation,
    /// The gluing parameter `m = v_2(a_2)` is out of range.
    GluingBound,
    /// The canonical four-torsion of the commitment side could not be formed.
    CanonicalCom,
    /// The canonical four-torsion of the challenge side could not be formed
    /// (the response images are not a basis of the `2^r`-torsion).
    CanonicalChal,
    /// A half-chain failed (gluing or a product along the chain).
    ChainFailed,
    /// The two half-chains do not meet.
    MiddleCodomainMismatch,
    /// The embedded isogeny is not the response.
    HdImageMismatch,
    /// Wrong length or a non-canonical field element or scalar.
    MalformedInput,
}

#[inline]
fn mask_r(r: u32) -> u128 {
    if r >= 128 {
        u128::MAX
    } else {
        (1u128 << r) - 1
    }
}

#[inline]
fn red_i(x: i128, r: u32) -> u128 {
    (x as u128) & mask_r(r)
}

#[inline]
fn mul_r(x: u128, y: u128, r: u32) -> u128 {
    x.wrapping_mul(y) & mask_r(r)
}

/// `a^-1 mod 2^r` for odd `a` (Newton).
#[inline]
fn inv_2r(a: u128, r: u32) -> u128 {
    let m = mask_r(r);
    let a = a & m;
    debug_assert!(a & 1 == 1, "inverse mod 2^r requires odd input");
    let mut inv = 1u128;
    let mut prec = 1u32;
    while prec < r {
        inv = inv.wrapping_mul(2u128.wrapping_sub(a.wrapping_mul(inv))) & m;
        prec *= 2;
    }
    inv & m
}

/// The dropped scalar from `a d − b c ≡ k q (mod 2^r)`: `(c, d)` with the
/// transmitted one of them being `c_or_d` (`c` when `a` is odd, else `d`).
#[inline]
pub fn recover_response_cd(
    a: i128,
    b: i128,
    c_or_d: i128,
    q: u128,
    k: u128,
    r: u32,
) -> (u128, u128) {
    let a_r = red_i(a, r);
    let b_r = red_i(b, r);
    let cod = red_i(c_or_d, r);
    let q_r = q & mask_r(r);
    let k_r = k & mask_r(r);
    let kq = mul_r(k_r, q_r, r);
    if a_r & 1 == 1 {
        let c = cod;
        let bc = mul_r(b_r, c, r);
        let d = mul_r(inv_2r(a_r, r), (kq.wrapping_add(bc)) & mask_r(r), r);
        (c, d)
    } else {
        let d = cod;
        let ad = mul_r(a_r, d, r);
        let c = mul_r(inv_2r(b_r, r), ad.wrapping_sub(kq) & mask_r(r), r);
        (c, d)
    }
}
