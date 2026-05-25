use crate::fp::FpBackend;

use super::theta_structure::{hadamard, to_squared_theta};
use super::{ThetaIsogeny, ThetaPoint, ThetaStructure};

/// Compute a (2,2)-isogeny from 8-torsion kernel generators.
///
/// Given a theta structure `a` and two points `t1_8`, `t2_8` of order 8,
/// computes the isogeny with kernel `[4](t1_8, t2_8)`.
///
/// `hadamard_bool_1` controls whether the domain uses standard or dual
/// coordinates; `hadamard_bool_2` controls the same for the codomain.
///
/// When `verify` is true, checks compatibility conditions on the kernel images.
/// Returns `None` on failure.
#[inline]
pub fn theta_isogeny_compute<L: FpBackend>(
    a: &ThetaStructure<L>,
    t1_8: &ThetaPoint<L>,
    t2_8: &ThetaPoint<L>,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
    verify: bool,
) -> Option<ThetaIsogeny<L>> {
    let mut out = ThetaIsogeny {
        t1_8: t1_8.clone(),
        t2_8: t2_8.clone(),
        hadamard_bool_1,
        hadamard_bool_2,
        domain: a.clone(),
        precomputation: ThetaPoint::default(),
        codomain: ThetaStructure::default(),
    };
    out.codomain.precomputation = false;

    let (tt1, tt2) = if hadamard_bool_1 {
        let h1 = hadamard(t1_8);
        let h2 = hadamard(t2_8);
        (to_squared_theta(&h1), to_squared_theta(&h2))
    } else {
        (to_squared_theta(t1_8), to_squared_theta(t2_8))
    };

    // Projective factors must be nonzero
    if bool::from(
        tt2.x.ct_is_zero()
            | tt2.y.ct_is_zero()
            | tt2.z.ct_is_zero()
            | tt2.t.ct_is_zero()
            | tt1.x.ct_is_zero()
            | tt1.y.ct_is_zero(),
    ) {
        return None;
    }

    let t1 = tt1.x.mul(&tt2.y);
    let t2 = tt1.y.mul(&tt2.x);

    out.codomain.null_point.x = tt2.x.mul(&t1);
    out.codomain.null_point.y = tt2.y.mul(&t2);
    out.codomain.null_point.z = tt2.z.mul(&t1);
    out.codomain.null_point.t = tt2.t.mul(&t2);

    let t3 = tt2.z.mul(&tt2.t);
    out.precomputation.x = t3.mul(&tt1.y);
    out.precomputation.y = t3.mul(&tt1.x);
    out.precomputation.z = out.codomain.null_point.t.clone();
    out.precomputation.t = out.codomain.null_point.z.clone();

    if verify {
        let v1 = tt1.x.mul(&out.precomputation.x);
        let v2 = tt1.y.mul(&out.precomputation.y);
        if !bool::from(v1.ct_equal(&v2)) {
            return None;
        }
        let v1 = tt1.z.mul(&out.precomputation.z);
        let v2 = tt1.t.mul(&out.precomputation.t);
        if !bool::from(v1.ct_equal(&v2)) {
            return None;
        }
        let v1 = tt2.x.mul(&out.precomputation.x);
        let v2 = tt2.z.mul(&out.precomputation.z);
        if !bool::from(v1.ct_equal(&v2)) {
            return None;
        }
        let v1 = tt2.y.mul(&out.precomputation.y);
        let v2 = tt2.t.mul(&out.precomputation.t);
        if !bool::from(v1.ct_equal(&v2)) {
            return None;
        }
    }

    if hadamard_bool_2 {
        out.codomain.null_point = hadamard(&out.codomain.null_point);
    }

    Some(out)
}

/// Compute a (2,2)-isogeny from 4-torsion kernel generators.
///
/// When only the 4-torsion above the kernel is known (not the 8-torsion),
/// computes the codomain using square roots of the null point's
/// squared-theta coordinate products (sqrt(AA*BB) and sqrt(AA*CC)).
#[inline]
pub fn theta_isogeny_compute_4<L: FpBackend>(
    a: &ThetaStructure<L>,
    t1_4: &ThetaPoint<L>,
    t2_4: &ThetaPoint<L>,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) -> ThetaIsogeny<L> {
    let mut out = ThetaIsogeny {
        t1_8: t1_4.clone(),
        t2_8: t2_4.clone(),
        hadamard_bool_1,
        hadamard_bool_2,
        domain: a.clone(),
        precomputation: ThetaPoint::default(),
        codomain: ThetaStructure::default(),
    };
    out.codomain.precomputation = false;

    let (tt1, tt2) = if hadamard_bool_1 {
        let h1 = hadamard(t1_4);
        let sq1 = to_squared_theta(&h1);
        let h2 = hadamard(&a.null_point);
        let sq2 = to_squared_theta(&h2);
        (sq1, sq2)
    } else {
        (to_squared_theta(t1_4), to_squared_theta(&a.null_point))
    };

    // sqrt(AA*BB) and sqrt(AA*CC)
    let sqaabb = tt2.x.mul(&tt2.y).sqrt();
    let sqaacc = tt2.x.mul(&tt2.z).sqrt();

    // Codomain null point
    out.codomain.null_point.y = sqaabb.mul(&sqaacc);
    out.precomputation.t = out.codomain.null_point.y.mul(&tt1.z);
    out.codomain.null_point.y = out.codomain.null_point.y.mul(&tt1.x);

    out.codomain.null_point.t = tt1.z.mul(&sqaabb).mul(&tt2.x);

    let t1x_aa = tt1.x.mul(&tt2.x);
    out.codomain.null_point.z = t1x_aa.mul(&tt2.z);
    out.codomain.null_point.x = t1x_aa.mul(&sqaacc);

    // Precomputation
    let t1x_dd = tt1.x.mul(&tt2.t);
    out.precomputation.z = t1x_dd.mul(&tt2.y);
    let t1x_dd_cc = t1x_dd.mul(&tt2.z);
    out.precomputation.y = t1x_dd_cc.mul(&sqaabb);
    out.precomputation.x = t1x_dd_cc.mul(&tt2.y);
    out.precomputation.z = out.precomputation.z.mul(&sqaacc);
    out.precomputation.t = out.precomputation.t.mul(&tt2.y);

    if hadamard_bool_2 {
        out.codomain.null_point = hadamard(&out.codomain.null_point);
    }

    out
}

/// Compute a (2,2)-isogeny from 2-torsion kernel generators.
///
/// When only the kernel is known (not the 4- or 8-torsion above it),
/// uses the null point to extract the necessary square roots.
#[inline]
pub fn theta_isogeny_compute_2<L: FpBackend>(
    a: &ThetaStructure<L>,
    t1_2: &ThetaPoint<L>,
    t2_2: &ThetaPoint<L>,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) -> ThetaIsogeny<L> {
    let mut out = ThetaIsogeny {
        t1_8: t1_2.clone(),
        t2_8: t2_2.clone(),
        hadamard_bool_1,
        hadamard_bool_2,
        domain: a.clone(),
        precomputation: ThetaPoint::default(),
        codomain: ThetaStructure::default(),
    };
    out.codomain.precomputation = false;

    let tt2 = if hadamard_bool_1 {
        let h = hadamard(&a.null_point);
        to_squared_theta(&h)
    } else {
        to_squared_theta(&a.null_point)
    };

    // Codomain: (AA, sqrt(AA*BB), sqrt(AA*CC), sqrt(AA*DD))
    out.codomain.null_point.x = tt2.x.clone();
    out.codomain.null_point.y = tt2.x.mul(&tt2.y).sqrt();
    out.codomain.null_point.z = tt2.x.mul(&tt2.z).sqrt();
    out.codomain.null_point.t = tt2.x.mul(&tt2.t).sqrt();

    // Precomputation
    let cc_dd = tt2.z.mul(&tt2.t);
    out.precomputation.y = cc_dd.mul(&out.codomain.null_point.y);
    out.precomputation.x = cc_dd.mul(&tt2.y);
    out.precomputation.z = tt2.t.mul(&out.codomain.null_point.z).mul(&tt2.y);
    out.precomputation.t = tt2.z.mul(&out.codomain.null_point.t).mul(&tt2.y);

    if hadamard_bool_2 {
        out.codomain.null_point = hadamard(&out.codomain.null_point);
    }

    out
}

/// Evaluate a (2,2) theta isogeny on a theta point.
#[inline]
pub fn theta_isogeny_eval<L: FpBackend>(phi: &ThetaIsogeny<L>, p: &ThetaPoint<L>) -> ThetaPoint<L> {
    let mut out = if phi.hadamard_bool_1 {
        let h = hadamard(p);
        to_squared_theta(&h)
    } else {
        to_squared_theta(p)
    };

    out.x = out.x.mul(&phi.precomputation.x);
    out.y = out.y.mul(&phi.precomputation.y);
    out.z = out.z.mul(&phi.precomputation.z);
    out.t = out.t.mul(&phi.precomputation.t);

    if phi.hadamard_bool_2 {
        out = hadamard(&out);
    }

    out
}
