//!
//! Provides routines to generate canonical bases for the `2^f`-torsion
//! subgroup, including hint-based generation for fast recomputation
//! during verification, y-coordinate recovery, and on-curve checks.

use super::point::{xdbl_a24, xdbl_e0};
use super::{EcBasis, EcCurve, EcPoint, JacPoint};
use crate::fp::{Fp, Fp2, FpBackend};
use subtle::Choice;

/// Recover the y-coordinate of a point on the Montgomery curve
/// `y^2 = x^3 + Ax^2 + x` with `C = 1`.
///
/// Takes the x-coordinate `px` and the curve. Computes `y = sqrt(x^3 + Ax^2 + x)`
/// and returns `(y, is_on_curve)` indicating whether `px` is on the curve.
#[inline]
pub fn ec_recover_y<L: FpBackend>(px: &Fp2<L>, curve: &EcCurve<L>) -> (Fp2<L>, Choice) {
    let t0 = px.sqr();
    let y = t0.mul(&curve.a); // Ax^2
    let y = y.add(px); // Ax^2 + x
    let t0 = t0.mul(px);
    let mut y = y.add(&t0); // x^3 + Ax^2 + x
    let valid = y.sqrt_verify();
    (y, valid)
}

/// Compute a deterministic difference point `P - Q` from x-only projective
/// coordinates, using Proposition 3 of <https://eprint.iacr.org/2017/518>.
#[inline]
pub fn difference_point<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> EcPoint<L> {
    let t0 = p.x.mul(&q.x);
    let t1 = p.z.mul(&q.z);
    let bxx = t0.sub(&t1);
    let bxx = bxx.sqr();
    let bxx = bxx.mul(&curve.c); // C*(P.x*Q.x - P.z*Q.z)^2

    let bxz = t0.add(&t1);
    let t0 = p.x.mul(&q.z);
    let t1 = p.z.mul(&q.x);
    let bzz_sum = t0.add(&t1);
    let bxz = bxz.mul(&bzz_sum); // (P.x*Q.x + P.z*Q.z)(P.x*Q.z + P.z*Q.x)

    let bzz = t0.sub(&t1);
    let bzz = bzz.sqr();
    let bzz = bzz.mul(&curve.c); // C*(P.x*Q.z - P.z*Q.x)^2

    let bxz = bxz.mul(&curve.c); // C*(...)*(...)
    let t0_t1 = t0.mul(&t1);
    let two_a_pzqz = t0_t1.mul(&curve.a);
    let two_a_pzqz = two_a_pzqz.add(&two_a_pzqz);
    let bxz = bxz.add(&two_a_pzqz);

    // Normalize by C * conj(C)^2 * conj(P.z)^2 * conj(Q.z)^2
    // to ensure the denominator is a fourth power in Fp
    let norm = curve.c.conjugate();
    let norm = norm.sqr();
    let norm = norm.mul(&curve.c);

    let pz_bar = p.z.conjugate();
    let pz_bar_sq = pz_bar.sqr();
    let norm = norm.mul(&pz_bar_sq);

    let qz_bar = q.z.conjugate();
    let qz_bar_sq = qz_bar.sqr();
    let norm = norm.mul(&qz_bar_sq);

    let bxx = bxx.mul(&norm);
    let bxz = bxz.mul(&norm);
    let bzz = bzz.mul(&norm);

    // Solve quadratic: discriminant = Bxz^2 - Bxx*Bzz
    let disc = bxz.sqr();
    let prod = bxx.mul(&bzz);
    let disc = disc.sub(&prod);
    let disc = disc.sqrt();

    let pq_x = bxz.add(&disc);
    let pq_z = bzz;

    EcPoint::new(pq_x, pq_z)
}

/// Like [`difference_point`], but selects between the two quadratic roots.
///
/// The quadratic formula yields two candidate x-coordinates for P ± Q.
/// When `negate_disc` is false, returns `(Bxz + √Δ : Bzz)` (same as
/// `difference_point`). When true, returns `(Bxz - √Δ : Bzz)`.
#[inline]
pub fn difference_point_with_hint<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    curve: &EcCurve<L>,
    negate_disc: bool,
) -> Option<EcPoint<L>> {
    let t0 = p.x.mul(&q.x);
    let t1 = p.z.mul(&q.z);
    let bxx = t0.sub(&t1);
    let bxx = bxx.sqr();
    let bxx = bxx.mul(&curve.c);

    let bxz = t0.add(&t1);
    let t0 = p.x.mul(&q.z);
    let t1 = p.z.mul(&q.x);
    let bzz_sum = t0.add(&t1);
    let bxz = bxz.mul(&bzz_sum);

    let bzz = t0.sub(&t1);
    let bzz = bzz.sqr();
    let bzz = bzz.mul(&curve.c);

    let bxz = bxz.mul(&curve.c);
    let t0_t1 = t0.mul(&t1);
    let two_a_pzqz = t0_t1.mul(&curve.a);
    let two_a_pzqz = two_a_pzqz.add(&two_a_pzqz);
    let bxz = bxz.add(&two_a_pzqz);

    let norm = curve.c.conjugate();
    let norm = norm.sqr();
    let norm = norm.mul(&curve.c);

    let pz_bar = p.z.conjugate();
    let pz_bar_sq = pz_bar.sqr();
    let norm = norm.mul(&pz_bar_sq);

    let qz_bar = q.z.conjugate();
    let qz_bar_sq = qz_bar.sqr();
    let norm = norm.mul(&qz_bar_sq);

    let bxx = bxx.mul(&norm);
    let bxz = bxz.mul(&norm);
    let bzz = bzz.mul(&norm);

    let disc = bxz.sqr();
    let prod = bxx.mul(&bzz);
    let disc = disc.sub(&prod);
    let disc = disc.sqrt();

    // When disc == 0 both sign choices yield the same point, so only
    // the canonical encoding (negate_disc = false) is accepted.
    if negate_disc && bool::from(disc.ct_is_zero()) {
        return None;
    }

    let pq_x = if negate_disc {
        bxz.sub(&disc)
    } else {
        bxz.add(&disc)
    };

    Some(EcPoint::new(pq_x, bzz))
}

/// Lift an x-only basis `{P, Q, P-Q}` to Jacobian coordinates, assuming
/// `P.z = 1` and `E.C = 1` (normalized curve and point).
///
/// Returns `Choice(1)` if P's x-coordinate is on the curve.
#[inline]
pub fn lift_basis_normalized<L: FpBackend>(
    basis: &mut EcBasis<L>,
    curve: &EcCurve<L>,
) -> (JacPoint<L>, JacPoint<L>, Choice) {
    let (py, ret) = ec_recover_y(&basis.p.x, curve);

    let jp = JacPoint::new(basis.p.x.clone(), py.clone(), Fp2::one());

    // Okeya-Sakurai y-recovery for Q
    let v1 = jp.x.mul(&basis.q.z);
    let v2 = basis.q.x.add(&v1);
    let v3 = basis.q.x.sub(&v1);
    let v3 = v3.sqr();
    let v3 = v3.mul(&basis.pmq.x);

    let v1_2a = curve.a.add(&curve.a);
    let v1_2a = v1_2a.mul(&basis.q.z);
    let v2 = v2.add(&v1_2a);

    let v4 = jp.x.mul(&basis.q.x);
    let v4 = v4.add(&basis.q.z);
    let v2 = v2.mul(&v4);

    let v1_sub = v1_2a.mul(&basis.q.z);
    let v2 = v2.sub(&v1_sub);
    let v2 = v2.mul(&basis.pmq.z);

    let qy_num = v3.sub(&v2);

    let v1_denom = py.add(&py);
    let v1_denom = v1_denom.mul(&basis.q.z);
    let v1_denom = v1_denom.mul(&basis.pmq.z);

    let qx = basis.q.x.mul(&v1_denom);
    let qz = basis.q.z.mul(&v1_denom);

    // Transform to Jacobian coordinates
    let qz_sq = qz.sqr();
    let qy_jac = qy_num.mul(&qz_sq);
    let qx_jac = qx.mul(&qz);

    let jq = JacPoint::new(qx_jac, qy_jac, qz);

    (jp, jq, ret)
}

/// Lift an x-only basis to Jacobian coordinates (general case).
///
/// Normalizes the curve and P.z first, then delegates to
/// `lift_basis_normalized`.
#[inline]
pub fn lift_basis<L: FpBackend>(
    basis: &mut EcBasis<L>,
    curve: &mut EcCurve<L>,
) -> (JacPoint<L>, JacPoint<L>, Choice) {
    // Batch-invert P.z and C
    let mut inverses = [basis.p.z.clone(), curve.c.clone()];
    let mut t1 = [Fp2::<L>::zero(), Fp2::zero()];
    let mut t2 = [Fp2::<L>::zero(), Fp2::zero()];
    Fp2::batched_inv(&mut inverses, &mut t1, &mut t2);

    basis.p.x = basis.p.x.mul(&inverses[0]);
    basis.p.z = Fp2::one();
    curve.a = curve.a.mul(&inverses[1]);
    curve.c = Fp2::one();

    lift_basis_normalized(basis, curve)
}

/// Check if an x-coordinate is on the curve `y^2 = x^3 + Ax^2 + x`.
/// Assumes `C = 1`.
#[inline]
pub fn is_on_curve<L: FpBackend>(x: &Fp2<L>, curve: &EcCurve<L>) -> Choice {
    let t0 = x.add(&curve.a); // x + A
    let t0 = t0.mul(x); // x^2 + Ax
    let t0 = t0.add_one(); // x^2 + Ax + 1
    let t0 = t0.mul(x); // x^3 + Ax^2 + x
    t0.is_square()
}

/// Clear the odd cofactor and excess powers of 2 to get a point of
/// order exactly `2^f`.
#[inline]
fn clear_cofactor_for_maximal_even_order<L: FpBackend>(
    p: &EcPoint<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    cofactor: &[u64],
    cofactor_bitlen: usize,
    torsion_even_power: u32,
) -> EcPoint<L> {
    // Clear the odd cofactor
    let mut result = super::point::ec_mul(p, cofactor, cofactor_bitlen, curve);

    // Clear excess power of two to get order 2^f
    for _ in 0..(torsion_even_power - f) {
        result = xdbl_a24(&result, &curve.a24, curve.is_a24_computed_and_normalized);
    }
    result
}

/// Generate the canonical `2^f`-torsion basis on E0 (A = 0).
///
/// Uses the precomputed basis points `BASIS_E0_PX` and `BASIS_E0_QX`.
#[inline]
fn ec_basis_e0_2f<L: FpBackend>(
    curve: &EcCurve<L>,
    f: u32,
    px_bytes: &[u8],
    qx_bytes: &[u8],
    torsion_even_power: u32,
) -> EcBasis<L> {
    let px = Fp2::<L>::decode(px_bytes).expect("invariant: precomputed BASIS_E0_PX must decode");
    let qx = Fp2::<L>::decode(qx_bytes).expect("invariant: precomputed BASIS_E0_QX must decode");

    let mut p = EcPoint::new(px, Fp2::one());
    let mut q = EcPoint::new(qx, Fp2::one());

    for _ in 0..(torsion_even_power - f) {
        p = xdbl_e0(&p);
        q = xdbl_e0(&q);
    }

    let pmq = difference_point(&p, &q, curve);
    EcBasis::new(p, q, pmq)
}

/// Find a non-quadratic-residue factor for entangled basis generation.
///
/// Finds the smallest `b >= start` such that `1 + b^2` is a NQR in Fp
/// and `-A/(1 + i*b)` is a valid x-coordinate. Returns `(x, hint)` where
/// `hint` encodes `b` as a u8 (0 if b >= 128, signaling fallback).
#[inline]
fn find_nqr_factor<L: FpBackend>(curve: &EcCurve<L>, start: u8) -> Result<(Fp2<L>, u16), ()> {
    let mut n = start as u16;
    loop {
        let mut qr_b = true;
        while qr_b {
            if n >= 1024 {
                return Err(());
            }
            let val = (n as u64) * (n as u64) + 1;
            let tmp = Fp::<L>::from_small(val);
            qr_b = bool::from(tmp.is_square());
            n += 1;
        }

        let b = Fp::<L>::from_small((n - 1) as u64);
        let z = Fp2 {
            re: Fp::one(),
            im: b.clone(),
        };
        let t0 = Fp2 {
            re: Fp::zero(),
            im: b,
        };

        let a_sq = curve.a.sqr();
        let lhs = a_sq.mul(&t0);
        let z_sq = z.sqr();
        let disc = lhs.sub(&z_sq);
        if !bool::from(disc.is_square()) {
            let z_inv = z.inv();
            let x = curve.a.mul(&z_inv).neg();
            let hint = if n <= 128 { n - 1 } else { 0 };
            return Ok((x, hint));
        }
    }
}

/// Find x(P) = n*A for the smallest n >= start such that n*A is on the curve.
///
/// Only called when A is a NQR. Returns `(x, hint)`.
#[inline]
fn find_na_x_coord<L: FpBackend>(curve: &EcCurve<L>, start: u8) -> Result<(Fp2<L>, u16), ()> {
    let mut n = start as u16;
    let mut x = if n == 1 {
        curve.a.clone()
    } else {
        curve.a.mul_small(n as u32)
    };

    while !bool::from(is_on_curve(&x, curve)) {
        if n >= 1024 {
            return Err(());
        }
        x = x.add(&curve.a);
        n += 1;
    }

    let hint = if n < 128 { n } else { 0 };
    Ok((x, hint))
}

/// Generate a `2^f`-torsion basis `<P, Q>` on the given curve, where Q
/// is above the Montgomery point `(0:0)`.
///
/// Stores a hint byte for fast recomputation via `ec_curve_to_basis_2f_from_hint`.
/// Returns `(basis, hint)`.
///
/// The `precomp` parameter provides the Level-specific constants.
#[inline]
pub fn ec_curve_to_basis_2f_to_hint<L: FpBackend>(
    curve: &mut EcCurve<L>,
    f: u32,
    px_bytes: &[u8],
    qx_bytes: &[u8],
    cofactor: &[u64],
    cofactor_bitlen: usize,
    torsion_even_power: u32,
) -> Result<(EcBasis<L>, u8), ()> {
    curve.normalize_curve_and_a24();

    if bool::from(curve.a.ct_is_zero()) {
        let basis = ec_basis_e0_2f(curve, f, px_bytes, qx_bytes, torsion_even_power);
        return Ok((basis, 0));
    }

    let hint_a = bool::from(curve.a.is_square());

    let (px, hint) = if !hint_a {
        find_na_x_coord(curve, 1)?
    } else {
        find_nqr_factor(curve, 1)?
    };

    let p = EcPoint::new(px.clone(), Fp2::one());

    let qx = curve.a.add(&px).neg();
    let q = EcPoint::new(qx, Fp2::one());

    let p = clear_cofactor_for_maximal_even_order(
        &p,
        curve,
        f,
        cofactor,
        cofactor_bitlen,
        torsion_even_power,
    );
    let q = clear_cofactor_for_maximal_even_order(
        &q,
        curve,
        f,
        cofactor,
        cofactor_bitlen,
        torsion_even_power,
    );

    // Compute P-Q, then set output so Q is above (0,0)
    let basis_q = difference_point(&p, &q, curve);

    let hint_byte = ((hint as u8) << 1) | (hint_a as u8);
    Ok((EcBasis::new(p, basis_q, q), hint_byte))
}

/// Reconstruct a `2^f`-torsion basis `<P, Q>` from a hint byte, where Q
/// is above the Montgomery point `(0:0)`.
///
/// Returns `(basis, ok)` where `ok` is always 1 (failure only in debug
/// builds of the C reference).
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn ec_curve_to_basis_2f_from_hint<L: FpBackend>(
    curve: &mut EcCurve<L>,
    f: u32,
    hint: u8,
    px_bytes: &[u8],
    qx_bytes: &[u8],
    cofactor: &[u64],
    cofactor_bitlen: usize,
    torsion_even_power: u32,
) -> Result<(EcBasis<L>, i32), ()> {
    curve.normalize_curve_and_a24();

    if bool::from(curve.a.ct_is_zero()) {
        let basis = ec_basis_e0_2f(curve, f, px_bytes, qx_bytes, torsion_even_power);
        return Ok((basis, 1));
    }

    let hint_a = (hint & 1) != 0;
    let hint_p = hint >> 1;

    let px = if hint_p == 0 {
        if !hint_a {
            find_na_x_coord(curve, 128)?.0
        } else {
            find_nqr_factor(curve, 128)?.0
        }
    } else if !hint_a {
        curve.a.mul_small(hint_p as u32)
    } else {
        let z = Fp2 {
            re: Fp::one(),
            im: Fp::<L>::from_small(hint_p as u64),
        };
        let z_inv = z.inv();
        curve.a.mul(&z_inv).neg()
    };

    let p = EcPoint::new(px.clone(), Fp2::one());

    let qx = curve.a.add(&px).neg();
    let q = EcPoint::new(qx, Fp2::one());

    let p = clear_cofactor_for_maximal_even_order(
        &p,
        curve,
        f,
        cofactor,
        cofactor_bitlen,
        torsion_even_power,
    );
    let q = clear_cofactor_for_maximal_even_order(
        &q,
        curve,
        f,
        cofactor,
        cofactor_bitlen,
        torsion_even_power,
    );

    let basis_q = difference_point(&p, &q, curve);
    Ok((EcBasis::new(p, basis_q, q), 1))
}
