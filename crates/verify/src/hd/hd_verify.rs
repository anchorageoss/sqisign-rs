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
use super::wide::W;

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
    /// A basis hint is not the canonical one of its curve.
    NonCanonicalHint,
    /// A pairing does not have exact order `2^r`: `e(P, Q)` on the challenge
    /// curve, or `e(R, S)` on the commitment curve (its logarithm `k` is
    /// even).
    PairingOrder,
    /// The response degree is not good: `q ≥ 2^e`, or `2^e − q` is not a
    /// probable prime congruent to 1 modulo 4.
    DegreeNotGood,
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

/// The dropped scalar from `a d − b c ≡ k q (mod 2^r)`: `(c, d)` with the
/// transmitted one of them being `c_or_d` (`c` when `a` is odd, else `d`).
/// The caller has rejected `a` and `b` both even, where neither inverse
/// exists.
#[inline]
pub fn recover_response_cd(a: W, b: W, c_or_d: W, q: W, k: W, r: u32) -> (W, W) {
    let m = W::mask(r);
    let (a_r, b_r, cod) = (a & m, b & m, c_or_d & m);
    let kq = (k & m).wrapping_mul(q & m) & m;
    if a_r.is_odd() {
        let c = cod;
        let bc = b_r.wrapping_mul(c) & m;
        let d = a_r.inv_odd().wrapping_mul(kq.wrapping_add(bc)) & m;
        (c, d)
    } else {
        let d = cod;
        let ad = a_r.wrapping_mul(d) & m;
        let c = b_r.inv_odd().wrapping_mul(ad.wrapping_sub(kq)) & m;
        (c, d)
    }
}
