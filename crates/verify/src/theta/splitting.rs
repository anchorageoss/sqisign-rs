//! The splitting step: from the dual theta null point of the last codomain,
//! recognise a product of elliptic curves and recover the two Montgomery
//! curves in canonical form together with the images of the pushed points.

use super::structure::is_product_theta_point;
use super::{ChainMode, ThetaCoupleCurve, ThetaCouplePoint, ThetaPoint};
use crate::ec::normalize::{ec_apply_isomorphism, ec_theta_to_montgomery};
use crate::ec::EcPoint;
use crate::fp::FpBackend;

/// Maximum number of points a chain can push (mirrors `chain::MAX_POINTS`).
pub const MAX_POINTS: usize = 8;

/// A split codomain with the images of the pushed points.
pub type SplitProduct<L> = (
    ThetaCoupleCurve<L>,
    [Option<ThetaCouplePoint<L>>; MAX_POINTS],
);

/// Change to the product theta structure: `(x + z : y + t : y - t : x - z)`.
#[inline]
fn to_product_structure<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    ThetaPoint {
        x: p.x.add(&p.z),
        y: p.y.add(&p.t),
        z: p.y.sub(&p.t),
        t: p.x.sub(&p.z),
    }
}

/// A product theta point `(x : y : z : t)` with `x t = y z` as the pair of
/// Kummer points `((x : z), (x : y))`. `None` if it is not a product point.
#[inline]
fn tensor_product_to_couple_point<L: FpBackend>(p: &ThetaPoint<L>) -> Option<ThetaCouplePoint<L>> {
    if !bool::from(is_product_theta_point(p)) {
        return None;
    }
    Some(ThetaCouplePoint::new(
        EcPoint::new(p.x.clone(), p.z.clone()),
        EcPoint::new(p.x.clone(), p.y.clone()),
    ))
}

/// Whether `(x : z)` is a valid elliptic theta null point: both coordinates
/// non-zero and `x^4 != z^4`.
#[inline]
fn valid_theta_null_point<L: FpBackend>(th: &EcPoint<L>) -> bool {
    let xx = th.x.sqr().sqr();
    let zz = th.z.sqr().sqr();
    !bool::from(th.x.ct_is_zero() | th.z.ct_is_zero() | xx.ct_equal(&zz))
}

/// Split the final codomain. `dual_null_point` is the dual theta null point
/// of the last structure; `pts` are the pushed points. Which factors are
/// brought to canonical Montgomery form depends on `mode`. Returns the
/// product and the images, or `None` if the codomain is not a product of
/// elliptic curves.
pub fn splitting_to_elliptic_product<L: FpBackend>(
    dual_null_point: &ThetaPoint<L>,
    pts: &[ThetaPoint<L>],
    mode: ChainMode,
) -> Option<SplitProduct<L>> {
    let null = to_product_structure(dual_null_point);
    let zero =
        null.x.ct_is_zero() | null.y.ct_is_zero() | null.z.ct_is_zero() | null.t.ct_is_zero();
    if bool::from(zero) {
        return None;
    }
    let th1th2 = tensor_product_to_couple_point(&null)?;

    let mut images: [Option<ThetaCouplePoint<L>>; MAX_POINTS] = Default::default();
    for (slot, p) in images.iter_mut().zip(pts.iter()) {
        *slot = Some(tensor_product_to_couple_point(&to_product_structure(p))?);
    }

    if mode == ChainMode::Verify && !valid_theta_null_point(&th1th2.p2) {
        return None;
    }

    let mut e12 = ThetaCoupleCurve {
        e1: crate::ec::EcCurve::init(),
        e2: crate::ec::EcCurve::init(),
    };
    if matches!(mode, ChainMode::Both | ChainMode::E1 | ChainMode::Verify) {
        let (e1, m) = ec_theta_to_montgomery(&th1th2.p1)?;
        e12.e1 = e1;
        for img in images.iter_mut().flatten() {
            ec_apply_isomorphism(&mut img.p1, &m);
        }
    }
    if matches!(mode, ChainMode::Both | ChainMode::E2) {
        let (e2, m) = ec_theta_to_montgomery(&th1th2.p2)?;
        e12.e2 = e2;
        for img in images.iter_mut().flatten() {
            ec_apply_isomorphism(&mut img.p2, &m);
        }
    }
    Some((e12, images))
}
