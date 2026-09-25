//! Canonical Montgomery models (spec Algorithm 4.96, `NormalizeMontgomery`).
//! A curve has six Montgomery coefficients; given a point of order 4 not
//! above `(0, 0)`, the reference computes all six and keeps the largest under
//! the `Fp2` ordering, together with the change of coordinates to it.

use super::{EcBasis, EcChangeCoordMatrix, EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};
use subtle::Choice;

/// `(X : Z) -> (a X + b Z : c X + d Z)`.
#[inline]
pub fn ec_apply_isomorphism<L: FpBackend>(p: &mut EcPoint<L>, m: &EcChangeCoordMatrix<L>) {
    let x = p.x.mul(&m.a).add(&p.z.mul(&m.b));
    let z = p.x.mul(&m.c).add(&p.z.mul(&m.d));
    p.x = x;
    p.z = z;
}

/// The Montgomery coefficient `(2 (x^4 + z^4) : x^4 - z^4)` attached to the
/// theta null point `(x : z)`.
#[inline]
fn compute_montgomery_coefficient<L: FpBackend>(th: &EcPoint<L>) -> EcPoint<L> {
    let xx = th.x.sqr().sqr();
    let zz = th.z.sqr().sqr();
    let sum = xx.add(&zz);
    EcPoint::new(sum.add(&sum), xx.sub(&zz))
}

/// Index of the largest element under [`Fp2::less_than`], constant time.
#[inline]
fn find_max_coefficient<L: FpBackend>(list: &[Fp2<L>]) -> (Fp2<L>, u8) {
    let mut max = list[0].clone();
    let mut max_index = 0u8;
    for (i, item) in list.iter().enumerate().skip(1) {
        let ctl = max.less_than(item);
        max = Fp2::select(&max, item, ctl);
        let mask = 0u8.wrapping_sub(ctl.unwrap_u8());
        max_index ^= (max_index ^ i as u8) & mask;
    }
    (max, max_index)
}

/// From a theta null point `(x : z)` of a curve (a point of order 4 not
/// above `(0, 0)`), the canonical Montgomery model `(A_max : 1)` with `A24`
/// normalised, and the change of coordinates onto it. `None` if a coefficient
/// denominator vanishes, i.e. the theta point was invalid.
#[inline]
pub fn ec_theta_to_montgomery<L: FpBackend>(
    th: &EcPoint<L>,
) -> Option<(EcCurve<L>, EcChangeCoordMatrix<L>)> {
    let ma = th.x.neg();
    let mb = th.z.neg();
    let ia = th.x.mul_by_i(Choice::from(0));
    let ib = th.z.mul_by_i(Choice::from(0));

    let mont1 = compute_montgomery_coefficient(th);
    let th2 = EcPoint::new(th.x.add(&th.z), th.x.sub(&th.z));
    let mont2 = compute_montgomery_coefficient(&th2);
    let th3 = EcPoint::new(ia.add(&th.z), th.x.add(&ib));
    let mont3 = compute_montgomery_coefficient(&th3);

    let mut dens = [mont1.z.clone(), mont2.z.clone(), mont3.z.clone()];
    let mut t1 = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero()];
    let mut t2 = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero()];
    Fp2::batched_inv(&mut dens, &mut t1, &mut t2);
    if bool::from(dens[0].ct_is_zero()) {
        return None;
    }
    let a1 = mont1.x.mul(&dens[0]);
    let a2 = mont2.x.mul(&dens[1]);
    let a3 = mont3.x.mul(&dens[2]);
    let coefs = [
        a1.clone(),
        a2.clone(),
        a3.clone(),
        a1.neg(),
        a2.neg(),
        a3.neg(),
    ];
    let (a_max, max) = find_max_coefficient(&coefs);

    let mut curve = EcCurve::<L>::init();
    curve.a = a_max;
    curve.normalize_curve_and_a24();

    // 0: [ b,  a;  b, -a]   1: [-a,  b;  b, -a]   2: [ia, ib; -b,  a]
    // 3: [ b,  a; -b,  a]   4: [-a,  b; -b,  a]   5: [ia, ib;  b, -a]
    let is = |set: &[u8]| Choice::from(set.contains(&max) as u8);
    let ma_sel = Fp2::select(&th.z, &ma, is(&[1, 4]));
    let a = Fp2::select(&ma_sel, &ia, is(&[2, 5]));
    let mb_sel = Fp2::select(&th.x, &th.z, is(&[1, 4]));
    let b = Fp2::select(&mb_sel, &ib, is(&[2, 5]));
    let c = Fp2::select(&th.z, &mb, is(&[2, 3, 4]));
    let d = Fp2::select(&th.x, &ma, is(&[0, 1, 5]));
    Some((curve, EcChangeCoordMatrix { a, b, c, d }))
}

/// `NormalizeMontgomery`: replace `curve` by its canonical model using a
/// point of order 4 not above `(0, 0)`, moving `basis` along. `None` if the
/// point was invalid.
#[inline]
pub fn ec_normalize_montgomery<L: FpBackend>(
    curve: &mut EcCurve<L>,
    four_torsion: &EcPoint<L>,
    basis: Option<&mut EcBasis<L>>,
) -> Option<()> {
    let (new_curve, m) = ec_theta_to_montgomery(four_torsion)?;
    *curve = new_curve;
    if let Some(b) = basis {
        ec_apply_isomorphism(&mut b.p, &m);
        ec_apply_isomorphism(&mut b.q, &m);
        ec_apply_isomorphism(&mut b.pmq, &m);
    }
    Some(())
}
