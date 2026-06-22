//! SQIsignHD FastVerify orchestration and the SHAKE256 hash-to-challenge.
//!
//! This ties together the six verification stages into a single entry point.
//! The genuinely new piece - not covered by the sage oracle - is the
//! **hash-to-challenge** ([`hd_challenge`]), ported from the C reference
//! (`Signature/src/sqisignhd/ref/sqisignhdx/sign.c`, `hash_to_challenge`):
//!
//! ```text
//!   chal = LE_int( SHAKE256_n( encode(j(E_com)) || encode(j(E_pk)) || message ) )
//! ```
//!
//! i.e. a **single** SHAKE256 pass (no iteration/grinding), the commitment
//! curve's j-invariant **first**, then the public key's, then the message,
//! squeezing `n` = the field byte length and reading it as a little-endian
//! integer. (Contrast the dim-2 SQIsign challenge in `crate::hash`,
//! which hashes `pk` first, iterates, and masks.) The `Fp2::encode` byte format
//! and `sha3::Shake256` are exactly those the dim-2 verifier uses to pass the
//! 300 NIST KATs, so they match `fp2_encode` / `fips202.c` byte-for-byte.
//!
//! The orchestration recomputes the challenge from the signature's curves and
//! the message and checks it against the signature's challenge (the binding of
//! message to signature), then runs the dimension-4 chain and the
//! middle-codomain check.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use crate::ec::EcCurve;
use crate::{Fp2, FpBackend};
use typenum::Unsigned as _;

use crate::hd::chain::{middle_codomain_matches, run_half_chain};
use crate::hd::point::ThetaPointDim4;

/// Maximum challenge byte length across levels (Level 5 field is < 64 bytes).
pub const MAX_CHAL_BYTES: usize = 64;

/// The HD challenge byte length for level `L`: the field word-length the C
/// reference squeezes (`NWORDS_FIELD * 8`), which equals `FpEncodedBytes` at
/// Level 1 (32 bytes). Level-1 is the validated target of this phase.
#[inline]
pub fn hd_challenge_len<L: FpBackend>() -> usize {
    <L as crate::SecurityLevel>::FpEncodedBytes::USIZE
}

/// Compute the SQIsignHD challenge into `out` (sized by the caller to
/// [`hd_challenge_len`]): `SHAKE256(encode(j_com) || encode(j_pk) || message)`.
#[inline]
pub fn hd_challenge<L: FpBackend>(
    j_com: &Fp2<L>,
    j_pk: &Fp2<L>,
    message: &[u8],
    out: &mut [u8],
) {
    let mut hasher = Shake256::default();
    hasher.update(j_com.encode().as_ref());
    hasher.update(j_pk.encode().as_ref());
    hasher.update(message);
    let mut reader = hasher.finalize_xof();
    reader.read(out);
}

/// Compute the challenge from the commitment and public-key Montgomery
/// `A`-coefficients (stage 1 curve recovery + the hash). Returns `false` if
/// either coefficient is not a valid Montgomery curve.
#[inline]
pub fn hd_challenge_from_curves<L: FpBackend>(
    a_com: &Fp2<L>,
    a_pk: &Fp2<L>,
    message: &[u8],
    out: &mut [u8],
) -> bool {
    let (e_com, e_pk) = match (EcCurve::from_a(a_com), EcCurve::from_a(a_pk)) {
        (Some(c), Some(p)) => (c, p),
        _ => return false,
    };
    hd_challenge(&e_com.j_inv(), &e_pk.j_inv(), message, out);
    true
}

// Stage 3 (response recovery): recover (c, d) from (a, b, c_or_d, q, k) mod 2^r.
// Modular arithmetic mod 2^r (r <= 70 < 128) using u128: the low r bits of a
// wrapping product/sum are exactly the value mod 2^r.

#[inline]
fn mask_r(r: u32) -> u128 {
    if r >= 128 {
        u128::MAX
    } else {
        (1u128 << r) - 1
    }
}

/// `x mod 2^r` for a possibly-negative `x` (two's-complement low bits).
#[inline]
fn red_i(x: i128, r: u32) -> u128 {
    (x as u128) & mask_r(r)
}

/// `x * y mod 2^r`.
#[inline]
fn mul_r(x: u128, y: u128, r: u32) -> u128 {
    x.wrapping_mul(y) & mask_r(r)
}

/// Inverse of an odd `a` modulo `2^r`, by Newton/Hensel doubling.
#[inline]
fn inv_2r(a: u128, r: u32) -> u128 {
    let m = mask_r(r);
    let a = a & m;
    debug_assert!(a & 1 == 1, "inverse mod 2^r requires odd input");
    let mut inv = 1u128; // correct mod 2
    let mut prec = 1u32;
    while prec < r {
        // inv <- inv * (2 - a*inv) mod 2^r
        inv = inv.wrapping_mul(2u128.wrapping_sub(a.wrapping_mul(inv))) & m;
        prec *= 2;
    }
    inv & m
}

/// Recover the response coefficients `(c, d)` mod `2^r` from the signature
/// scalars `a, b`, the stored coefficient `c_or_d`, the response degree `q`,
/// and the discrete log `k`. Mirrors `image_response` in the sage reference:
/// if `a` is odd then `c = c_or_d`, `d = a⁻¹(k·q + b·c)`; otherwise
/// `d = c_or_d`, `c = b⁻¹(a·d - k·q)`. (`a·d - b·c ≡ k·q (mod 2^r)`.)
#[inline]
pub fn recover_response_cd(a: i128, b: i128, c_or_d: i128, q: u128, k: u128, r: u32) -> (u128, u128) {
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

// Orchestration

/// Inputs to [`hd_verify`].
///
/// The challenge binding (curves + message) is checked from real signature
/// fields; the dimension-4 chain runs from the per-step kernels. Wire parsing
/// (Phase 6) populates this struct from signature bytes; future phases replace
/// the supplied kernels with ones derived from stages 1-3.
pub struct HdVerifyInputs<'a, L: FpBackend> {
    /// Public-key curve Montgomery `A`.
    pub a_pk: Fp2<L>,
    /// Commitment curve Montgomery `A`.
    pub a_com: Fp2<L>,
    /// The signed message.
    pub message: &'a [u8],
    /// The challenge from the signature, as little-endian bytes
    /// ([`hd_challenge_len`] of them).
    pub claimed_chal: &'a [u8],
    /// Per-step kernels of the F1 half-chain (gluing then plain).
    pub f1_kernels: &'a [[ThetaPointDim4<L>; 4]],
    /// Per-step kernels of the F2_dual half-chain.
    pub f2_dual_kernels: &'a [[ThetaPointDim4<L>; 4]],
}

/// Why a verification was rejected (for diagnostics; `hd_verify` returns a bool).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HdReject {
    /// A curve coefficient was not a valid Montgomery curve.
    BadCurve,
    /// The recomputed challenge did not match the signature's challenge.
    ChallengeMismatch,
    /// A half-chain codomain was not computable.
    ChainFailed,
    /// `F1`'s codomain did not match the Hadamard of `F2_dual`'s.
    MiddleCodomainMismatch,
    /// The stage-6 HD-image condition `F(T) = (±a₁P, ±a₂P, *, 0)` failed.
    HdImageMismatch,
    /// The signature or public-key bytes could not be parsed (bad length,
    /// out-of-range field element, or non-canonical scalar).
    MalformedInput,
}

/// Run the verification and return `Ok(())` on accept or the rejection reason.
///
/// Checks performed (all real computation): (0/1) recompute the challenge from
/// the curves + message and compare to the signature's challenge - this binds
/// the message; (4) run both dimension-4 half-chains from their kernels;
/// (5) the projective middle-codomain match. Stage 6 (the HD-image condition)
/// is not yet computed here (it needs the full chain *evaluation*); see NOTES.
pub fn hd_verify_checked<L: FpBackend>(inp: &HdVerifyInputs<L>) -> Result<(), HdReject> {
    // Stage 0/1: challenge binding.
    let n = hd_challenge_len::<L>();
    let mut chal = [0u8; MAX_CHAL_BYTES];
    if !hd_challenge_from_curves(&inp.a_com, &inp.a_pk, inp.message, &mut chal[..n]) {
        return Err(HdReject::BadCurve);
    }
    if inp.claimed_chal.len() != n || chal[..n] != *inp.claimed_chal {
        return Err(HdReject::ChallengeMismatch);
    }

    // Stage 4: the dimension-4 half-chains.
    let f1_last = run_half_chain(inp.f1_kernels).ok_or(HdReject::ChainFailed)?;
    let f2_last = run_half_chain(inp.f2_dual_kernels).ok_or(HdReject::ChainFailed)?;

    // Stage 5: the middle-codomain match.
    if !middle_codomain_matches(&f1_last, &f2_last) {
        return Err(HdReject::MiddleCodomainMismatch);
    }
    Ok(())
}

/// Convenience wrapper: `true` iff the signature verifies.
#[inline]
pub fn hd_verify<L: FpBackend>(inp: &HdVerifyInputs<L>) -> bool {
    hd_verify_checked(inp).is_ok()
}
