//! The SQIsignHD library's torsion basis from hints (`torsion_basis_2f_from_hint`
//! in its `Theta_dim4` package): `x(P)` is the `hp`-th entry of the
//! non-residue table (or `hp + i` past the table), `x(Q)` the `hq`-th entry of
//! the second table times `α = (−A + √(A² − 4)) / 2`, both cleared of the odd
//! cofactor, then lifted to full points. This convention differs from the
//! round-3 specification's `TorsionBasisFromHint`; it is kept so that the
//! library's Sage verifier is the oracle for the compact format.

use crate::ec::basis::is_on_curve;
use crate::ec::point::ec_mul;
use crate::ec::{EcCurve, EcPoint, JacPoint};
use crate::fp::{Fp2, FpBackend};

use super::params::HdLevel;

// The library's square roots: `sqrt_Fp2_det` (the root with an even real
// part, or an even imaginary part when the real part is zero) for the
// basis, the difference point and `α`. The round-3 crate's `Fp2::sqrt`
// returns a different root; the compact format follows the library so the
// hinted basis is the oracle's, point for point.

/// `sqrt_Fp2_det`.
#[inline]
pub(crate) fn sqrt_det<L: FpBackend>(a: &Fp2<L>) -> Fp2<L> {
    a.sqrt_canonical_even()
}

/// The library's `difference_point`: `x(P − Q)` from `x(P)`, `x(Q)` (affine)
/// on `y^2 = x^3 + A x^2 + x`, as a projective `(X : Z)`.
fn hd_difference_point<L: FpBackend>(xp: &Fp2<L>, xq: &Fp2<L>, a: &Fp2<L>) -> EcPoint<L> {
    let zpq = xp.sub(xq);
    let t2 = xp.mul(xq);
    let t3 = t2.sub(&Fp2::one());
    let t0 = zpq.mul(&t3);
    let zpq = zpq.sqr();
    let t0 = t0.sqr();
    let t1 = t2.add(&Fp2::one());
    let t3 = xp.add(xq);
    let t1 = t1.mul(&t3);
    let t2 = t2.mul(a);
    let t2 = t2.add(&t2);
    let t1 = t1.add(&t2);
    let t2 = t1.sqr();
    let t0 = t2.sub(&t0);
    let t0 = sqrt_det(&t0);
    EcPoint::new(t0.add(&t1), zpq)
}

/// The library's `recover_y`: `sqrt_Fp2_det(x (x^2 + A x + 1))`.
fn hd_recover_y<L: FpBackend>(x: &Fp2<L>, a: &Fp2<L>) -> Fp2<L> {
    let y2 = x.sqr().add(&a.mul(x)).add(&Fp2::one()).mul(x);
    sqrt_det(&y2)
}

/// The library's `lift_basis`: `P = (x_P, recover_y(x_P))` and `Q` by
/// Okeya–Sakurai from `x(P)`, `(X : Z)(Q)` and `(X : Z)(P − Q)`.
fn hd_lift_basis<L: FpBackend>(
    xp: &Fp2<L>,
    kq: &EcPoint<L>,
    kpq: &EcPoint<L>,
    a: &Fp2<L>,
) -> (JacPoint<L>, JacPoint<L>) {
    let yp = hd_recover_y(xp, a);
    let v1 = xp.mul(&kq.z);
    let v2 = kq.x.add(&v1);
    let v3 = kq.x.sub(&v1);
    let v3 = v3.sqr();
    let v3 = v3.mul(&kpq.x);
    let v1 = a.add(a);
    let v1 = v1.mul(&kq.z);
    let v2 = v2.add(&v1);
    let v4 = xp.mul(&kq.x);
    let v4 = v4.add(&kq.z);
    let v2 = v2.mul(&v4);
    let v1 = v1.mul(&kq.z);
    let v2 = v2.sub(&v1);
    let v2 = v2.mul(&kpq.z);
    let yq = v3.sub(&v2);
    let v1 = yp.add(&yp);
    let v1 = v1.mul(&kq.z);
    let v1 = v1.mul(&kpq.z);
    let xq = kq.x.mul(&v1);
    let zq = kq.z.mul(&v1);
    // (xq : yq : zq) in projective (X : Y : Z) with x = X/Z, y = Y/Z, to
    // Jacobian (X Z : Y Z^2 : Z)
    let zq2 = zq.sqr();
    let p = JacPoint::new(xp.clone(), yp, Fp2::one());
    let q = JacPoint::new(xq.mul(&zq), yq.mul(&zq2), zq);
    (p, q)
}

#[inline]
fn fp2_from_le<L: FpBackend>(bytes: &[u8]) -> Fp2<L> {
    Fp2::<L>::decode(bytes).expect("non-residue table entry must be a canonical field element")
}

/// Affine `(x, y)` of a Jacobian point (not the identity).
#[inline]
pub fn jac_to_affine<L: FpBackend>(j: &JacPoint<L>) -> (Fp2<L>, Fp2<L>) {
    let z2 = j.z.sqr();
    let z3 = z2.mul(&j.z);
    (j.x.mul(&z2.inv()), j.y.mul(&z3.inv()))
}

/// The basis of `E_A[2^f]` (`f` the full 2-adic exponent) the hints
/// `(hp, hq)` select, as Jacobian points, exactly as the library's
/// `torsion_basis_2f_from_hint`; `None` on an invalid curve.
pub fn torsion_basis_2f_from_hint<L: FpBackend>(
    curve_a: &Fp2<L>,
    hp: u32,
    hq: u32,
    nqr_table: &[Fp2<L>],
    z_nqr_table: &[Fp2<L>],
    cofactor: &[u64],
    cofactor_bits: usize,
) -> Option<(JacPoint<L>, JacPoint<L>)> {
    let mut curve = EcCurve::from_a(curve_a)?;
    let i = Fp2::<L>::i_element();

    let xp = if (hp as usize) < nqr_table.len() {
        nqr_table[hp as usize].clone()
    } else {
        Fp2::<L>::from_small(hp as u64).add(&i)
    };

    let disc = curve_a.sqr().sub(&Fp2::<L>::from_small(4));
    let alpha = sqrt_det(&disc).sub(curve_a).half();

    let xq_base = if (hq as usize) < z_nqr_table.len() {
        z_nqr_table[hq as usize].clone()
    } else {
        Fp2::<L>::from_small(hq as u64).add(&i)
    };
    let xq = xq_base.mul(&alpha);

    let mut kp = ec_mul(
        &EcPoint::new(xp, Fp2::one()),
        cofactor,
        cofactor_bits,
        &mut curve,
    );
    let mut kq = ec_mul(
        &EcPoint::new(xq, Fp2::one()),
        cofactor,
        cofactor_bits,
        &mut curve,
    );
    kp.normalize();
    kq.normalize();

    let kpq = hd_difference_point(&kp.x, &kq.x, curve_a);
    Some(hd_lift_basis(&kp.x, &kq, &kpq, curve_a))
}

/// The library's canonical hints for a curve: the first table entries that
/// give points on the curve. `None` if no entry does (a crafted curve).
pub fn canonical_hints<L: HdLevel>(curve_a: &Fp2<L>) -> Option<(u32, u32)> {
    let curve = EcCurve::from_a(curve_a)?;
    let nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::nqr_table()[k]));
    let z_nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::z_nqr_table()[k]));

    let hp = (0..20).find(|&k| bool::from(is_on_curve(&nqr[k], &curve)))? as u32;

    let alpha = sqrt_det(&curve_a.sqr().sub(&Fp2::<L>::from_small(4)))
        .sub(curve_a)
        .half();
    let hq = (0..20).find(|&k| bool::from(is_on_curve(&z_nqr[k].mul(&alpha), &curve)))? as u32;

    Some((hp, hq))
}

/// The hinted basis of `E_A[2^f]` at this level's tables and cofactor.
pub fn hd_torsion_basis<L: HdLevel>(
    curve_a: &Fp2<L>,
    hp: u32,
    hq: u32,
) -> Option<(JacPoint<L>, JacPoint<L>)> {
    let nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::nqr_table()[k]));
    let z_nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::z_nqr_table()[k]));
    let (cof, bits) = L::cofactor();
    torsion_basis_2f_from_hint(curve_a, hp, hq, &nqr, &z_nqr, &cof, bits)
}
