//! Phase 5b.2 - challenge isogeny recovery (stage 2).
//!
//! Recovers the challenge isogeny `φ_chal : E_pk → E_chal` of degree `2^λ`
//! from the public-key curve and the challenge scalar, and computes the data
//! that the response stage (5b.3) consumes: the codomain `E_chal`, the image
//! Weil pairing `w_chal`, and the rescaled image basis `(P_chal_resc,
//! Q_chal_resc)`. This mirrors `Verify.py::recover_chal`.
//!
//! # The sage reference (Level 1: `e = 248`, `λ = 128`, `r = 70`)
//!
//! ```text
//! rescale1 = 2^(e-λ-r) = 2^50;  rescale2 = 2^r = 2^70;  rescale3 = 2^λ = 2^128
//! B_pk_rplamb = (2^50·P_pk, 2^50·Q_pk)                       # order 2^198
//! B_pk_lamb   = (2^70·B_pk_rplamb[0], 2^70·B_pk_rplamb[1])   # order 2^128
//! φ_chal, E_chal = isogeny_from_scalar_x_only(E_pk, 2^128, chal, B_pk_lamb)
//!                                          # ker = B_pk_lamb[0] + chal·B_pk_lamb[1]
//! φP, φQ, w_chal = evaluate_isogeny_x_only_with_image_pairing(
//!                      φ_chal, B_pk_rplamb[0], B_pk_rplamb[1], 2^198, 2^128)
//! P_chal_resc = φP + chal·φQ ;  Q_chal_resc = 2^128·φQ
//! ```
//!
//! # What is reused from the dim-2 crate (`sqisign-verify`)
//!
//! The challenge isogeny is computed exactly as the dim-2 verifier computes its
//! own challenge curve (`verify::compute_challenge_curve`): an `ec_ladder3pt`
//! kernel `P + [chal]Q` followed by `ec_eval_even` (which also pushes the
//! rescaled basis through the isogeny). The pairing is the biextension
//! [`weil`]; the rescaled image basis is assembled with the Jacobian full-point
//! arithmetic ([`jac_add`], [`jac_dbl`]). The torsion basis itself comes from
//! Phase 5b.1 ([`crate::hd::hd_torsion_basis_l1`]).
//!
//! # What is HD-specific (new here)
//!
//! The image points are lifted from their x-coordinates with the **FESTA**
//! square root (`utilities/fast_sqrt.py::sqrt_Fp2`), *not* the canonical-even
//! `Fp2::sqrt`. The two roots differ by a sign, and `recover_chal`'s
//! `P_chal_resc`/`Q_chal_resc` are pinned to the FESTA convention, so matching
//! them byte-for-byte requires reproducing that exact (non-canonical) lift.

use crate::ec::isogeny::ec_eval_even;
use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::point::{ec_dbl_iter_basis, ec_ladder3pt};
use crate::ec::{EcBasis, EcCurve, EcIsogEven, JacPoint};
use crate::{Fp2, FpBackend, Level1};

use crate::hd::hd_torsion_basis_l1;

/// Level-1 SQIsignHD parameters used by stage 2.
const RESCALE1_BITS: usize = 50; // e - λ - r  = 248 - 128 - 70
const RESCALE2_BITS: usize = 70; // r
const LAMBDA_BITS: usize = 128; // λ (challenge isogeny degree exponent)
const PAIRING_E: u32 = 198; // e - rescale1_bits = order of B_pk_rplamb (2^198)

/// `(p - 3) / 4 = 5·2^246 - 1`, little-endian u64 limbs (FESTA `sqrt_Fp2`).
const P_MINUS_3_DIV_4: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    90_071_992_547_409_919,
];
/// `(p - 1) / 2 = 5·2^247 - 1`, little-endian u64 limbs (FESTA `sqrt_Fp2`).
const P_MINUS_1_DIV_2: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    180_143_985_094_819_839,
];

/// The recovered challenge data (stage 2 outputs).
///
/// `e_chal` is normalised (`C = 1`, so `e_chal.a` is the affine Montgomery
/// coefficient). `p_chal_resc`/`q_chal_resc` are full Jacobian points on
/// `e_chal`; use [`crate::hd::jac_to_affine`] to compare against affine references.
pub struct ChallengeRecovery {
    pub e_chal: EcCurve<Level1>,
    pub w_chal: Fp2<Level1>,
    pub p_chal_resc: JacPoint<Level1>,
    pub q_chal_resc: JacPoint<Level1>,
}

/// FESTA square root (`utilities/fast_sqrt.py::sqrt_Fp2`).
///
/// A valid (non-canonical) square root of `a`, matching the sign convention
/// the HD `curve_point` lift is pinned to. Branches on field data - fine for a
/// verifier operating on public values.
fn festa_sqrt<L: FpBackend>(a: &Fp2<L>) -> Fp2<L> {
    let a1 = a.pow_vartime(&P_MINUS_3_DIV_4);
    let x0 = a1.mul(a);
    let alpha = a1.mul(&x0);
    let neg_one = Fp2::<L>::one().neg();
    if bool::from(alpha.ct_equal(&neg_one)) {
        Fp2::<L>::i_element().mul(&x0)
    } else {
        let b = Fp2::<L>::one().add(&alpha).pow_vartime(&P_MINUS_1_DIV_2);
        b.mul(&x0)
    }
}

/// Lift an affine x-coordinate to a full point on `y² = x³ + A x² + x` using
/// the FESTA square root (`kummer_line.py::curve_point`).
fn curve_point_lift<L: FpBackend>(x: &Fp2<L>, a: &Fp2<L>) -> JacPoint<L> {
    // y² = x·(x² + A·x + 1)
    let y2 = x.sqr().add(&a.mul(x)).add(&Fp2::one()).mul(x);
    let y = festa_sqrt(&y2);
    JacPoint::new(x.clone(), y, Fp2::one())
}

/// Full-point scalar multiplication `[k]·P` (MSB-first double-and-add).
///
/// `k` is little-endian u64 limbs; all `64·len` bits are processed (leading
/// zeros are harmless). [`jac_dbl`]/[`jac_add`] handle the identity, so the
/// accumulator may legitimately pass through `O`.
pub(crate) fn jac_scalar_mul<L: FpBackend>(
    p: &JacPoint<L>,
    k: &[u64],
    curve: &EcCurve<L>,
) -> JacPoint<L> {
    let mut acc = JacPoint::identity();
    for i in (0..k.len() * 64).rev() {
        acc = jac_dbl(&acc, curve);
        if (k[i >> 6] >> (i & 63)) & 1 == 1 {
            acc = jac_add(&acc, p, curve);
        }
    }
    acc
}

/// Iterated full-point doubling `[2ⁿ]·P`.
pub(crate) fn jac_dbl_iter<L: FpBackend>(
    p: &JacPoint<L>,
    n: usize,
    curve: &EcCurve<L>,
) -> JacPoint<L> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, curve);
    }
    acc
}

/// Recover the challenge isogeny `φ_chal : E_pk → E_chal` and the stage-2 data
/// (`Verify.py::recover_chal`) for Level 1.
///
/// `a_pk` is the public-key Montgomery coefficient, `(hp, hq)` its torsion-basis
/// hints, and `chal` the challenge scalar as little-endian u64 limbs. Returns
/// `None` if any curve/kernel is degenerate.
pub fn recover_challenge_l1(
    a_pk: &Fp2<Level1>,
    hp: u32,
    hq: u32,
    chal: &[u64],
) -> Option<ChallengeRecovery> {
    let mut e_pk = EcCurve::from_a(a_pk)?;
    e_pk.normalize_a24(); // required by ec_ladder3pt; e_pk.a stays affine (C = 1)

    // 2^248-torsion basis (P_pk, Q_pk) from the hint (Phase 5b.1), as full
    // points so the x-only difference x(P_pk - Q_pk) is unambiguous.
    let (p_pk, q_pk) = hd_torsion_basis_l1(a_pk, hp, hq)?;
    let pmq = jac_add(&p_pk, &q_pk.neg(), &e_pk);
    let base = EcBasis::new(p_pk.to_xz(), q_pk.to_xz(), pmq.to_xz());

    // B_pk_rplamb (order 2^198) and B_pk_lamb (order 2^128 = the kernel basis).
    let basis_rplamb = ec_dbl_iter_basis(&base, RESCALE1_BITS, &mut e_pk);
    let basis_lamb = ec_dbl_iter_basis(&basis_rplamb, RESCALE2_BITS, &mut e_pk);

    // Kernel K = B_pk_lamb[0] + [chal]·B_pk_lamb[1], then the 2^λ-isogeny.
    let kernel = ec_ladder3pt(chal, &basis_lamb.p, &basis_lamb.q, &basis_lamb.pmq, &e_pk)?;
    let phi = EcIsogEven {
        curve: e_pk.clone(),
        kernel,
        length: LAMBDA_BITS as u32,
    };

    // Push B_pk_rplamb through φ: images are x(φP), x(φQ) on E_chal.
    let mut images = [basis_rplamb.p.clone(), basis_rplamb.q.clone()];
    let mut e_chal = e_pk.clone();
    ec_eval_even(&mut e_chal, &phi, &mut images)?;
    e_chal.normalize(); // C = 1 → e_chal.a is the affine coefficient
    e_chal.normalize_a24();

    // Lift the images to full points via the FESTA square root.
    let mut xp = images[0].clone();
    let mut xq = images[1].clone();
    xp.normalize();
    xq.normalize();
    let im_p = curve_point_lift(&xp.x, &e_chal.a);
    let mut im_q = curve_point_lift(&xq.x, &e_chal.a);

    // Sign-correct φQ via the Weil-pairing compatibility e(φP,φQ) = e(P,Q)^deg.
    // The comparison is convention-independent (any global pairing convention
    // raised to the same power on both curves preserves the equality), so the
    // flip matches the oracle regardless of how `weil` normalises the value.
    let im_pmq = jac_add(&im_p, &im_q.neg(), &e_chal);
    let pair_e1 = weil(
        PAIRING_E,
        &im_p.to_xz(),
        &im_q.to_xz(),
        &im_pmq.to_xz(),
        &mut e_chal,
    );
    let pair_e0 = weil(
        PAIRING_E,
        &basis_rplamb.p,
        &basis_rplamb.q,
        &basis_rplamb.pmq,
        &mut e_pk,
    );
    let mut w_chal = pair_e0;
    for _ in 0..LAMBDA_BITS {
        w_chal = w_chal.sqr(); // w_chal = e(P,Q)^(2^λ)
    }
    if !bool::from(w_chal.ct_equal(&pair_e1)) {
        im_q = im_q.neg();
    }

    // Rescaled image basis on E_chal.
    let p_chal_resc = jac_add(&im_p, &jac_scalar_mul(&im_q, chal, &e_chal), &e_chal);
    let q_chal_resc = jac_dbl_iter(&im_q, LAMBDA_BITS, &e_chal);

    Some(ChallengeRecovery {
        e_chal,
        w_chal,
        p_chal_resc,
        q_chal_resc,
    })
}
