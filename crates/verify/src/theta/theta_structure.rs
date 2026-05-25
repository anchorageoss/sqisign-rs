use crate::fp::FpBackend;

use super::{ThetaPoint, ThetaStructure};

/// Hadamard transform on a theta point.
///
/// (x, y, z, t) -> (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)
#[inline]
pub fn hadamard<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    let t1 = p.x.add(&p.y);
    let t2 = p.x.sub(&p.y);
    let t3 = p.z.add(&p.t);
    let t4 = p.z.sub(&p.t);

    ThetaPoint {
        x: t1.add(&t3),
        y: t2.add(&t4),
        z: t1.sub(&t3),
        t: t2.sub(&t4),
    }
}

/// Coordinate-wise squaring of a theta point.
///
/// (x, y, z, t) -> (x^2, y^2, z^2, t^2)
#[inline]
pub fn pointwise_square<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    ThetaPoint {
        x: p.x.sqr(),
        y: p.y.sqr(),
        z: p.z.sqr(),
        t: p.t.sqr(),
    }
}

/// Square all coordinates then apply Hadamard.
///
/// (x, y, z, t) -> hadamard(x^2, y^2, z^2, t^2)
#[inline]
pub fn to_squared_theta<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    hadamard(&pointwise_square(p))
}

/// Precompute the 8 Fp2 constants needed for efficient theta doubling
/// and (2,2)-isogeny computation.
///
/// Given null_point = (x, y, z, t) and its squared theta dual
/// (XX, YY, ZZ, TT) = to_squared_theta(null_point), computes:
///   - cap_{XYZ0, YZT0, XZT0, XYT0}: triple products from the dual
///   - {xyz0, yzt0, xzt0, xyt0}: triple products from the null point
#[inline]
pub fn theta_precomputation<L: FpBackend>(a: &mut ThetaStructure<L>) {
    if a.precomputation {
        return;
    }

    let a_dual = to_squared_theta(&a.null_point);

    let t1 = a_dual.x.mul(&a_dual.y);
    let t2 = a_dual.z.mul(&a_dual.t);
    a.cap_xyz0 = t1.mul(&a_dual.z);
    a.cap_xyt0 = t1.mul(&a_dual.t);
    a.cap_yzt0 = t2.mul(&a_dual.y);
    a.cap_xzt0 = t2.mul(&a_dual.x);

    let t1 = a.null_point.x.mul(&a.null_point.y);
    let t2 = a.null_point.z.mul(&a.null_point.t);
    a.xyz0 = t1.mul(&a.null_point.z);
    a.xyt0 = t1.mul(&a.null_point.t);
    a.yzt0 = t2.mul(&a.null_point.y);
    a.xzt0 = t2.mul(&a.null_point.x);

    a.precomputation = true;
}

/// Double a theta point on the given theta structure.
///
/// Assumes no coordinate of the null point is zero. Performs
/// precomputation lazily if not yet done.
#[inline]
pub fn double_point<L: FpBackend>(p: &ThetaPoint<L>, a: &mut ThetaStructure<L>) -> ThetaPoint<L> {
    let mut out = to_squared_theta(p);

    out.x = out.x.sqr();
    out.y = out.y.sqr();
    out.z = out.z.sqr();
    out.t = out.t.sqr();

    if !a.precomputation {
        theta_precomputation(a);
    }

    out.x = out.x.mul(&a.cap_yzt0);
    out.y = out.y.mul(&a.cap_xzt0);
    out.z = out.z.mul(&a.cap_xyt0);
    out.t = out.t.mul(&a.cap_xyz0);

    out = hadamard(&out);

    out.x = out.x.mul(&a.yzt0);
    out.y = out.y.mul(&a.xzt0);
    out.z = out.z.mul(&a.xyt0);
    out.t = out.t.mul(&a.xyz0);

    out
}

/// Iterated theta doubling: compute [2^exp] P.
#[inline]
pub fn double_iter<L: FpBackend>(
    p: &ThetaPoint<L>,
    a: &mut ThetaStructure<L>,
    exp: u32,
) -> ThetaPoint<L> {
    if exp == 0 {
        p.clone()
    } else {
        let mut out = double_point(p, a);
        for _ in 1..exp {
            out = double_point(&out, a);
        }
        out
    }
}

/// Check if a theta point is a product theta point.
///
/// A theta point (x : y : z : t) is a product point when x*t == y*z.
#[inline]
pub fn is_product_theta_point<L: FpBackend>(p: &ThetaPoint<L>) -> subtle::Choice {
    let t1 = p.x.mul(&p.t);
    let t2 = p.y.mul(&p.z);
    t1.ct_equal(&t2)
}
