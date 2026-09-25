//! Coordinate-wise operations on theta points, the Hadamard transform, and
//! doubling on a theta structure (spec Algorithms 4.108 to 4.114).

use super::{ThetaPoint, ThetaStructure};
use crate::fp::FpBackend;
use subtle::Choice;

/// `(x + y + z + t : x - y + z - t : x + y - z - t : x - y - z + t)`.
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

/// Coordinate-wise projective inverse `(yzt : xzt : xyt : xyz)`, six
/// multiplications and no inversion.
#[inline]
pub fn invert_point<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    let a = p.x.mul(&p.y);
    let d = p.z.mul(&p.t);
    ThetaPoint {
        x: p.y.mul(&d),
        y: p.x.mul(&d),
        z: p.t.mul(&a),
        t: p.z.mul(&a),
    }
}

/// The Hadamard transform in place (P28): two temporaries, the four
/// coordinates rewritten where they are.
#[inline]
pub fn hadamard_assign<L: FpBackend>(p: &mut ThetaPoint<L>) {
    // (x, y, z, t) -> (x + y, x - y, z + t, z - t) -> (t1 + t3, t2 + t4, t1 - t3, t2 - t4)
    let mut t1 = p.x.clone();
    t1.add_assign(&p.y); // t1 = x + y
    p.x.sub_assign(&p.y); // x <- t2 = x - y
    let mut t3 = p.z.clone();
    t3.add_assign(&p.t); // t3 = z + t
    p.z.sub_assign(&p.t); // z <- t4 = z - t
                          // then y <- t2 + t4 and t <- t2 - t4
    p.y = p.x.clone();
    p.y.add_assign(&p.z);
    p.t = p.x.clone();
    p.t.sub_assign(&p.z);
    // x <- t1 + t3, z <- t1 - t3
    p.x = t1.clone();
    p.x.add_assign(&t3);
    p.z = t1;
    p.z.sub_assign(&t3);
}

/// Coordinate-wise square in place (P28).
#[inline]
pub fn pointwise_square_assign<L: FpBackend>(p: &mut ThetaPoint<L>) {
    p.x.sqr_assign();
    p.y.sqr_assign();
    p.z.sqr_assign();
    p.t.sqr_assign();
}

/// Coordinate-wise product in place (P28): `p <- p * q`.
#[inline]
pub fn pointwise_product_assign<L: FpBackend>(p: &mut ThetaPoint<L>, q: &ThetaPoint<L>) {
    p.x.mul_assign(&q.x);
    p.y.mul_assign(&q.y);
    p.z.mul_assign(&q.z);
    p.t.mul_assign(&q.t);
}

/// [`to_squared_theta`] in place (P28).
#[inline]
pub fn to_squared_theta_assign<L: FpBackend>(p: &mut ThetaPoint<L>) {
    pointwise_square_assign(p);
    hadamard_assign(p);
}

/// Coordinate-wise square.
#[inline]
pub fn pointwise_square<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    ThetaPoint {
        x: p.x.sqr(),
        y: p.y.sqr(),
        z: p.z.sqr(),
        t: p.t.sqr(),
    }
}

/// Coordinate-wise product.
#[inline]
pub fn pointwise_product<L: FpBackend>(p: &ThetaPoint<L>, q: &ThetaPoint<L>) -> ThetaPoint<L> {
    ThetaPoint {
        x: p.x.mul(&q.x),
        y: p.y.mul(&q.y),
        z: p.z.mul(&q.z),
        t: p.t.mul(&q.t),
    }
}

/// Hadamard transform of the coordinate-wise square: the squared theta
/// coordinates in the dual structure.
#[inline]
pub fn to_squared_theta<L: FpBackend>(p: &ThetaPoint<L>) -> ThetaPoint<L> {
    hadamard(&pointwise_square(p))
}

/// Whether `P = (x : y : z : t)` lies on a product theta structure, i.e.
/// `x t = y z`.
#[inline]
pub fn is_product_theta_point<L: FpBackend>(p: &ThetaPoint<L>) -> Choice {
    p.x.mul(&p.t).ct_equal(&p.y.mul(&p.z))
}

impl<L: FpBackend> ThetaStructure<L> {
    /// The dual theta null point, rebuilt from the isogeny data: with
    /// `dbl_data = (Az, Ax, Bx, Bw)`, `A = Az Ax Bw`, `B = Az Bx Bw`, and
    /// `C`, `D` are the `t` and `z` coordinates of the inverse dual null point.
    #[inline]
    pub fn dual_null_point(&self) -> ThetaPoint<L> {
        let tmp = self.dbl_data.x.mul(&self.dbl_data.t);
        ThetaPoint {
            x: tmp.mul(&self.dbl_data.y),
            y: tmp.mul(&self.dbl_data.z),
            z: self.inv_dual_null_point.t.clone(),
            t: self.inv_dual_null_point.z.clone(),
        }
    }

    /// Compute the doubling constant `inv_sqr_null_point` once.
    #[inline]
    pub fn precompute(&mut self) {
        if self.precomputation {
            return;
        }
        let dual = self.dual_null_point();
        self.inv_sqr_null_point = invert_point(&to_squared_theta(&dual));
        self.precomputation = true;
    }

    /// `[2] P` on this structure (spec Algorithm 4.114).
    #[inline]
    pub fn double(&mut self, p: &ThetaPoint<L>) -> ThetaPoint<L> {
        self.precompute();
        let out = pointwise_square(&to_squared_theta(p));
        let out = pointwise_product(&out, &self.inv_sqr_null_point);
        let out = hadamard(&out);
        pointwise_product(&out, &self.inv_dual_null_point)
    }

    /// `P <- [2] P` in place (P28): the same formulas as [`Self::double`],
    /// every coordinate rewritten where it is.
    #[inline]
    pub fn double_assign(&mut self, p: &mut ThetaPoint<L>) {
        self.precompute();
        to_squared_theta_assign(p);
        pointwise_square_assign(p);
        pointwise_product_assign(p, &self.inv_sqr_null_point);
        hadamard_assign(p);
        pointwise_product_assign(p, &self.inv_dual_null_point);
    }

    /// `[2^exp] P`.
    #[inline]
    pub fn double_iter(&mut self, p: &ThetaPoint<L>, exp: u16) -> ThetaPoint<L> {
        let mut out = p.clone();
        for _ in 0..exp {
            self.double_assign(&mut out);
        }
        out
    }

    /// A structure in regular (not dual) coordinates from its theta null
    /// point, as used right after the gluing step.
    #[inline]
    pub fn from_null_point(null: &ThetaPoint<L>) -> Self {
        Self {
            inv_dual_null_point: invert_point(null),
            dbl_data: ThetaPoint::zero(),
            inv_sqr_null_point: invert_point(&to_squared_theta(null)),
            precomputation: true,
        }
    }
}
