//! Stage 2 of the compact verification: the challenge isogeny and the
//! rescaled basis of `E_chall` (the library's `recover_chal`).
//!
//! From the hinted basis `(P, Q)` of `E_pk[2^f]`: `B_rpl = 2^(f − λ − r)
//! (P, Q)` of order `2^(λ + r)`, `B_λ = 2^r B_rpl` of order `2^λ`; the
//! challenge isogeny `φ` has kernel `B_λ.P + [c] B_λ.Q` and degree `2^λ`
//! (computed from a generator of order `2^(λ+2)`, whose top four-torsion
//! fixes the codomain's Montgomery model, as the round-3 chain does);
//! `φ(B_rpl)` is lifted to full points with the sign of `φ(Q)` fixed by
//! `e(φP, φQ) = e(P, Q)^(2^λ)`; then `P_resc = φP + [c] φQ`, `Q_resc = 2^λ
//! φQ`, of order `2^r` on `E_chall`, and `w_chall = e(P, Q)^(2^λ)` at
//! `2^(λ + r)`.

use crate::ec::isogeny::iso_isogeny_2chain_eval;
use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::point::{ec_dbl_iter_basis, ec_ladder3pt};
use crate::ec::{EcBasis, EcCurve, JacPoint};
use crate::fp::{Fp2, FpBackend};

use super::basis::hd_torsion_basis;
use super::params::HdLevel;

/// The challenge curve, its rescaled response basis and the pairing value
/// the response recovery divides by.
pub struct ChallengeRecovery<L: FpBackend> {
    pub e_chal: EcCurve<L>,
    pub w_chal: Fp2<L>,
    pub p_chal_resc: JacPoint<L>,
    pub q_chal_resc: JacPoint<L>,
}

/// The library's `sqrt_Fp2` (the FESTA formula: `a1 = a^((p−3)/4)`, `x0 =
/// a1 a`, `α = a1 x0`; `i x0` if `α = −1`, else `(1 + α)^((p−1)/2) x0`),
/// which its `KummerPoint.curve_point` uses to lift the challenge images;
/// a specific root, not the even-normalised one of the basis.
fn festa_sqrt<L: FpBackend>(a: &Fp2<L>) -> Fp2<L> {
    let (e34, e12) = festa_exponents::<L>();
    let a1 = a.pow_vartime(&e34);
    let x0 = a1.mul(a);
    let alpha = a1.mul(&x0);
    let neg_one = Fp2::<L>::one().neg();
    if bool::from(alpha.ct_equal(&neg_one)) {
        Fp2::<L>::i_element().mul(&x0)
    } else {
        Fp2::<L>::one().add(&alpha).pow_vartime(&e12).mul(&x0)
    }
}

/// `((p − 3)/4, (p − 1)/2)` as little-endian limbs, from the prime's bytes.
fn festa_exponents<L: FpBackend>() -> ([u64; 11], [u64; 11]) {
    let bytes = L::prime_le_bytes();
    let mut p = [0u64; 11];
    for (i, b) in bytes.iter().enumerate() {
        p[i / 8] |= (*b as u64) << (8 * (i % 8));
    }
    let sub_small = |mut x: [u64; 11], k: u64| {
        let (d, borrow) = x[0].overflowing_sub(k);
        x[0] = d;
        let mut i = 1;
        let mut b = borrow;
        while b && i < 11 {
            let (d, nb) = x[i].overflowing_sub(1);
            x[i] = d;
            b = nb;
            i += 1;
        }
        x
    };
    let shr = |mut x: [u64; 11], k: u32| {
        for _ in 0..k {
            let mut carry = 0u64;
            for w in x.iter_mut().rev() {
                let next = *w & 1;
                *w = (*w >> 1) | (carry << 63);
                carry = next;
            }
        }
        x
    };
    (shr(sub_small(p, 3), 2), shr(sub_small(p, 1), 1))
}

/// A full point above `x` on `y^2 = x^3 + A x^2 + x`, the library's lift.
fn curve_point_lift<L: FpBackend>(x: &Fp2<L>, a: &Fp2<L>) -> JacPoint<L> {
    let y2 = x.sqr().add(&a.mul(x)).add(&Fp2::one()).mul(x);
    JacPoint::new(x.clone(), festa_sqrt(&y2), Fp2::one())
}

/// `[k] P` on Jacobian points, double-and-add, variable time (public data).
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

/// `[2^n] P` on Jacobian points.
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

/// Recover the challenge from the public key and the challenge scalar
/// (`chal` little-endian limbs, `λ` bits).
pub fn recover_challenge<L: HdLevel>(
    a_pk: &Fp2<L>,
    hp: u32,
    hq: u32,
    chal: &[u64],
) -> Option<ChallengeRecovery<L>> {
    let f = L::TWO_ADIC_EXPONENT;
    let rescale1 = (f - L::LAMBDA - L::R) as usize;
    let rescale2 = L::R as usize;
    let pairing_e = f - rescale1 as u32; // λ + r

    let mut e_pk = EcCurve::from_a(a_pk)?;
    e_pk.normalize_a24();

    let (p_pk, q_pk) = hd_torsion_basis::<L>(a_pk, hp, hq)?;
    let pmq = jac_add(&p_pk, &q_pk.neg(), &e_pk);
    let base = EcBasis::new(p_pk.to_xz(), q_pk.to_xz(), pmq.to_xz());

    let basis_rplamb = ec_dbl_iter_basis(&base, rescale1, &mut e_pk);
    // the kernel generator with the two extra torsion bits the chain uses
    // to fix the codomain's model: order 2^(λ + 2), [4] of it is P_λ + [c] Q_λ
    let basis_lamb2 = ec_dbl_iter_basis(&basis_rplamb, rescale2 - 2, &mut e_pk);

    let kernel = ec_ladder3pt(
        chal,
        &basis_lamb2.p,
        &basis_lamb2.q,
        &basis_lamb2.pmq,
        &e_pk,
    )?;
    let mut images = [basis_rplamb.p.clone(), basis_rplamb.q.clone()];
    let mut e_chal = iso_isogeny_2chain_eval(&e_pk, &kernel, L::LAMBDA + 2, &mut images)?;
    e_chal.normalize();
    e_chal.normalize_a24();

    let mut xp = images[0].clone();
    let mut xq = images[1].clone();
    xp.normalize();
    xq.normalize();
    let im_p = curve_point_lift(&xp.x, &e_chal.a);
    let mut im_q = curve_point_lift(&xq.x, &e_chal.a);

    let im_pmq = jac_add(&im_p, &im_q.neg(), &e_chal);
    let pair_e1 = weil(
        pairing_e,
        &im_p.to_xz(),
        &im_q.to_xz(),
        &im_pmq.to_xz(),
        &mut e_chal,
    );
    let pair_e0 = weil(
        pairing_e,
        &basis_rplamb.p,
        &basis_rplamb.q,
        &basis_rplamb.pmq,
        &mut e_pk,
    );
    let mut w_chal = pair_e0;
    for _ in 0..L::LAMBDA {
        w_chal = w_chal.sqr();
    }
    if !bool::from(w_chal.ct_equal(&pair_e1)) {
        im_q = im_q.neg();
    }

    let p_chal_resc = jac_add(&im_p, &jac_scalar_mul(&im_q, chal, &e_chal), &e_chal);
    let q_chal_resc = jac_dbl_iter(&im_q, L::LAMBDA as usize, &e_chal);

    Some(ChallengeRecovery {
        e_chal,
        w_chal,
        p_chal_resc,
        q_chal_resc,
    })
}
