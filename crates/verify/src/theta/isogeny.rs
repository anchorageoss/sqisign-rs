//! Generic (2, 2)-isogeny steps between theta structures (spec Algorithms
//! 4.115 to 4.117).

use super::structure::{
    hadamard, hadamard_assign, pointwise_product_assign, to_squared_theta, to_squared_theta_assign,
};
use super::{ThetaPoint, ThetaStructure};
use crate::fp::FpBackend;

/// A generic (2, 2)-isogeny, carried as its codomain structure.
#[derive(Clone, Debug)]
pub struct ThetaIsogeny<L: FpBackend> {
    /// The codomain theta structure.
    pub codomain: ThetaStructure<L>,
}

/// The isogeny with kernel `[4] <T1_8, T2_8>`. With `dual_domain` the input
/// points are in dual coordinates. With `verify`, rejects kernels whose
/// 4-torsion is not isotropic or whose intermediate values vanish. `None`
/// on rejection.
pub fn theta_isogeny_compute<L: FpBackend>(
    t1_8: &ThetaPoint<L>,
    t2_8: &ThetaPoint<L>,
    dual_domain: bool,
    verify: bool,
) -> Option<ThetaIsogeny<L>> {
    let (tt1, tt2) = if dual_domain {
        (
            to_squared_theta(&hadamard(t1_8)),
            to_squared_theta(&hadamard(t2_8)),
        )
    } else {
        (to_squared_theta(t1_8), to_squared_theta(t2_8))
    };
    // TT1 = (Ax, Bx, Cy, Dy), TT2 = (Az, Bw, Cz, Dw); zeros mean an unexpected splitting
    if verify {
        let zero = tt2.x.ct_is_zero()
            | tt2.y.ct_is_zero()
            | tt2.z.ct_is_zero()
            | tt2.t.ct_is_zero()
            | tt1.x.ct_is_zero()
            | tt1.y.ct_is_zero();
        if bool::from(zero) {
            return None;
        }
    }
    let t1 = tt1.y.mul(&tt2.t); // Bx Dw
    let t2 = tt1.x.mul(&tt2.z); // Ax Cz
    let inv_dual_null_point = ThetaPoint {
        x: t1.mul(&tt2.z), // A^-1 = Cz Bx Dw
        y: t2.mul(&tt2.t), // B^-1 = Ax Cz Dw
        z: t1.mul(&tt2.x), // C^-1 = Az Bx Dw
        t: t2.mul(&tt2.y), // D^-1 = Ax Cz Bw
    };
    let dbl_data = ThetaPoint {
        x: tt2.x.clone(),
        y: tt1.x.clone(),
        z: tt1.y.clone(),
        t: tt2.y.clone(),
    };
    if verify {
        // isotropy of the 4-torsion: Cy C^-1 = Dy D^-1
        let a = tt1.z.mul(&inv_dual_null_point.z);
        let b = tt1.t.mul(&inv_dual_null_point.t);
        if !bool::from(a.ct_equal(&b)) {
            return None;
        }
    }
    Some(ThetaIsogeny {
        codomain: ThetaStructure {
            inv_dual_null_point,
            dbl_data,
            inv_sqr_null_point: ThetaPoint::zero(),
            precomputation: false,
        },
    })
}

impl<L: FpBackend> ThetaIsogeny<L> {
    /// `P <- image of P` in place (P28), the same arithmetic as
    /// [`Self::eval`] with no point constructed and moved.
    #[inline]
    pub fn eval_assign(&self, p: &mut ThetaPoint<L>, dual_domain: bool) {
        if dual_domain {
            hadamard_assign(p);
        }
        to_squared_theta_assign(p);
        pointwise_product_assign(p, &self.codomain.inv_dual_null_point);
    }

    /// Image of `P`, given in dual coordinates when `dual_domain`.
    #[inline]
    pub fn eval(&self, p: &ThetaPoint<L>, dual_domain: bool) -> ThetaPoint<L> {
        let s = if dual_domain {
            to_squared_theta(&hadamard(p))
        } else {
            to_squared_theta(p)
        };
        let n = &self.codomain.inv_dual_null_point;
        ThetaPoint {
            x: s.x.mul(&n.x),
            y: s.y.mul(&n.y),
            z: s.z.mul(&n.z),
            t: s.t.mul(&n.t),
        }
    }
}
