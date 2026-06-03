//!
//! Provides addition, doubling (both Montgomery and Weierstrass forms),
//! negation, and conversions between Jacobian, XZ, and Weierstrass
//! representations.

use super::{AddComponents, EcCurve, EcPoint, JacPoint};
use crate::fp::{Fp, Fp2, FpBackend};
use subtle::Choice;

/// Constant-time select for Jacobian points.
#[inline]
fn select_jac_point<L: FpBackend>(p1: &JacPoint<L>, p2: &JacPoint<L>, ctl: Choice) -> JacPoint<L> {
    JacPoint {
        x: Fp2::select(&p1.x, &p2.x, ctl),
        y: Fp2::select(&p1.y, &p2.y, ctl),
        z: Fp2::select(&p1.z, &p2.z, ctl),
    }
}

impl<L: FpBackend> JacPoint<L> {
    /// Test if two Jacobian points are equal.
    #[inline]
    pub fn ct_equal(&self, other: &Self) -> Choice {
        let t0 = other.z.sqr();
        let t2 = self.x.mul(&t0); // x1*z2^2
        let t1 = self.z.sqr();
        let t3 = other.x.mul(&t1); // x2*z1^2
        let dx = t2.sub(&t3);

        let t0 = t0.mul(&other.z);
        let t0 = self.y.mul(&t0); // y1*z2^3
        let t1 = t1.mul(&self.z);
        let t1 = other.y.mul(&t1); // y2*z1^3
        let dy = t0.sub(&t1);

        dy.ct_is_zero() & dx.ct_is_zero()
    }

    /// Convert from Jacobian to x-only (X:Z) by dropping Y and squaring Z.
    #[inline]
    pub fn to_xz(&self) -> EcPoint<L> {
        let x = self.x.clone();
        let z = self.z.sqr();

        // If this is the identity (0:1:0), we have (0:0) but want (1:0)
        let c1 = x.ct_is_zero();
        let c2 = z.ct_is_zero();
        let one = Fp2::<L>::one();
        let x = Fp2::select(&x, &one, c1 & c2);

        EcPoint { x, z }
    }

    /// Negate: `(X : Y : Z) -> (X : -Y : Z)`.
    #[inline]
    pub fn neg(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.neg(),
            z: self.z.clone(),
        }
    }

    /// Convert from Montgomery Jacobian to Weierstrass.
    ///
    /// Returns `(Q_ws, t, ao3)` where:
    /// - `Q_ws` is the point in Weierstrass Jacobian coords
    /// - `t = a * Z^4` (the modified Jacobian extra coordinate)
    /// - `ao3 = A/3`
    #[inline]
    pub fn to_ws(&self, curve: &EcCurve<L>) -> (JacPoint<L>, Fp2<L>, Fp2<L>) {
        let mut ao3 = Fp2::<L>::zero();
        let t;

        if !bool::from(curve.a.ct_is_zero()) {
            ao3.re = curve.a.re.div3();
            ao3.im = curve.a.im.div3();

            let z2 = self.z.sqr();
            let qx = ao3.mul(&z2);
            let qx = qx.add(&self.x);

            let z4 = z2.sqr();
            let a = ao3.mul(&curve.a);
            let one = Fp::<L>::one();
            let a_re = one.sub(&a.re);
            let a_im = a.im.neg();
            let a = Fp2 { re: a_re, im: a_im };
            t = z4.mul(&a);

            let q = JacPoint {
                x: qx,
                y: self.y.clone(),
                z: self.z.clone(),
            };
            (q, t, ao3)
        } else {
            let z2 = self.z.sqr();
            t = z2.sqr();

            let q = JacPoint {
                x: self.x.clone(),
                y: self.y.clone(),
                z: self.z.clone(),
            };
            (q, t, ao3)
        }
    }

    /// Convert from Weierstrass Jacobian back to Montgomery Jacobian.
    #[inline]
    pub fn from_ws(p: &JacPoint<L>, ao3: &Fp2<L>, curve: &EcCurve<L>) -> JacPoint<L> {
        let x = if !bool::from(curve.a.ct_is_zero()) {
            let t = p.z.sqr();
            let t = t.mul(ao3);
            p.x.sub(&t)
        } else {
            p.x.clone()
        };

        JacPoint {
            x,
            y: p.y.clone(),
            z: p.z.clone(),
        }
    }
}

/// Doubling on a Montgomery curve in Jacobian coordinates.
///
/// Cost: 6M + 6S.
#[inline]
pub fn jac_dbl<L: FpBackend>(p: &JacPoint<L>, ac: &EcCurve<L>) -> JacPoint<L> {
    let is_identity = p.x.ct_is_zero() & p.z.ct_is_zero();

    let t0 = p.x.sqr(); // x1^2
    let t1 = t0.add(&t0);
    let t0 = t0.add(&t1); // 3x1^2
    let t1 = p.z.sqr(); // z1^2
    let t2 = p.x.mul(&ac.a);
    let t2 = t2.add(&t2); // 2Ax1
    let t2 = t1.add(&t2); // 2Ax1+z1^2
    let t2 = t1.mul(&t2); // z1^2(2Ax1+z1^2)
    let t2 = t0.add(&t2); // alpha = 3x1^2 + z1^2(2Ax1+z1^2)

    let qz = p.y.mul(&p.z);
    let qz = qz.add(&qz); // z2 = 2y1z1

    let t0 = qz.sqr();
    let t0 = t0.mul(&ac.a); // 4Ay1^2z1^2
    let t1 = p.y.sqr();
    let t1 = t1.add(&t1); // 2y1^2
    let t3 = p.x.add(&p.x); // 2x1
    let t3 = t1.mul(&t3); // 4x1y1^2
    let qx = t2.sqr(); // alpha^2
    let qx = qx.sub(&t0); // alpha^2 - 4Ay1^2z1^2
    let qx = qx.sub(&t3);
    let qx = qx.sub(&t3); // alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    let qy = t3.sub(&qx); // 4x1y1^2 - x2
    let qy = qy.mul(&t2); // alpha(4x1y1^2 - x2)
    let t1 = t1.sqr(); // 4y1^4
    let qy = qy.sub(&t1);
    let qy = qy.sub(&t1); // alpha(4x1y1^2 - x2) - 8y1^4

    let qx = Fp2::select(&qx, &p.x, is_identity);
    let qz = Fp2::select(&qz, &p.z, is_identity);

    JacPoint {
        x: qx,
        y: qy,
        z: qz,
    }
}

/// Doubling on a Weierstrass curve in modified Jacobian coordinates.
///
/// Takes `(X:Y:Z)` and the extra coordinate `t = a*Z^4`, returns
/// `(Q, u)` where `u` is the updated `t` for the result point.
///
/// Cost: 3M + 5S.
#[inline]
pub fn jac_dbl_ws<L: FpBackend>(p: &JacPoint<L>, t: &Fp2<L>) -> (JacPoint<L>, Fp2<L>) {
    let is_identity = p.x.ct_is_zero() & p.z.ct_is_zero();

    let xx = p.x.sqr();
    let c = p.y.sqr();
    let c = c.add(&c); // A = 2*Y^2
    let cc = c.sqr(); // AA = A^2
    let r = cc.add(&cc); // R = 2*AA
    let s = p.x.add(&c);
    let s = s.sqr();
    let s = s.sub(&xx);
    let s = s.sub(&cc); // S = (X+A)^2-XX-AA
    let m = xx.add(&xx);
    let m = m.add(&xx);
    let m = m.add(t); // M = 3*XX+T1
    let qx = m.sqr();
    let qx = qx.sub(&s);
    let qx = qx.sub(&s); // X3 = M^2-2*S
    let qz = p.y.mul(&p.z);
    let qz = qz.add(&qz); // Z3 = 2*Y*Z
    let qy = s.sub(&qx);
    let qy = qy.mul(&m);
    let qy = qy.sub(&r); // Y3 = M*(S-X3)-R
    let u = t.mul(&r);
    let u = u.add(&u); // T3 = 2*R*T1

    let qx = Fp2::select(&qx, &p.x, is_identity);
    let qz = Fp2::select(&qz, &p.z, is_identity);

    (
        JacPoint {
            x: qx,
            y: qy,
            z: qz,
        },
        u,
    )
}

/// Complete addition on a Montgomery curve in Jacobian coordinates.
///
/// Handles all edge cases (identity, doubling, inverses) in constant time.
/// Cost: 17M + 6S + 13a.
#[inline]
pub fn jac_add<L: FpBackend>(p: &JacPoint<L>, q: &JacPoint<L>, ac: &EcCurve<L>) -> JacPoint<L> {
    let ctl1 = p.z.ct_is_zero();
    let ctl2 = q.z.ct_is_zero();

    let t0 = p.z.sqr(); // z1^2
    let t1 = q.z.sqr(); // z2^2

    // dy and dx for the ordinary case
    let v1 = t1.mul(&q.z); // z2^3
    let t2 = t0.mul(&p.z); // z1^3
    let v1 = v1.mul(&p.y); // y1*z2^3
    let t2 = t2.mul(&q.y); // y2*z1^3
    let dy_ord = t2.sub(&v1);
    let u2 = t0.mul(&q.x); // x2*z1^2
    let u1 = t1.mul(&p.x); // x1*z2^2
    let dx_ord = u2.sub(&u1);

    // dy and dx for the doubling case
    let t1_dbl = p.y.add(&p.y); // dx_dbl = 2y1
    let t2_dbl = {
        let two_a = ac.a.add(&ac.a);
        let t2 = two_a.mul(&p.x); // 2Ax1
        let t2 = t2.add(&t0); // 2Ax1 + z1^2
        let t2 = t2.mul(&t0); // z1^2 * (2Ax1 + z1^2)
        let x_sq = p.x.sqr();
        let t2 = t2.add(&x_sq);
        let t2 = t2.add(&x_sq);
        let t2 = t2.add(&x_sq); // 3*x1^2 + z1^2*(2Ax1 + z1^2)
        t2.mul(&q.z) // dy_dbl = z2 * (3*x1^2 + ...)
    };

    // Switch to double variables if dx == 0 and dy == 0
    let ctl = dx_ord.ct_is_zero() & dy_ord.ct_is_zero();
    let dx = Fp2::select(&dx_ord, &t1_dbl, ctl);
    let dy = Fp2::select(&dy_ord, &t2_dbl, ctl);

    let t0 = p.z.mul(&q.z); // z1*z2
    let t1 = t0.sqr(); // (z1*z2)^2
    let t2 = dx.sqr(); // dx^2
    let t3 = dy.sqr(); // dy^2

    // x3 = dy^2 - dx^2 * (A*(z1*z2)^2 + u1 + u2)
    let rx = ac.a.mul(&t1);
    let rx = rx.add(&u1);
    let rx = rx.add(&u2);
    let rx = rx.mul(&t2);
    let rx = t3.sub(&rx);

    // y3 = dy * (u1 * dx^2 - x3) - v1 * dx^3
    let ry = u1.mul(&t2);
    let ry = ry.sub(&rx);
    let ry = ry.mul(&dy);
    let t3 = t2.mul(&dx); // dx^3
    let t3 = t3.mul(&v1); // v1 * dx^3
    let ry = ry.sub(&t3);

    // z3 = dx * z1 * z2
    let rz = dx.mul(&t0);

    let mut r = JacPoint {
        x: rx,
        y: ry,
        z: rz,
    };

    // If P.z == 0, return Q; if Q.z == 0, return P
    r = select_jac_point(&r, q, ctl1);
    r = select_jac_point(&r, p, ctl2);

    r
}

/// Compute the addition components (u, v, w) such that
/// `P+Q = (u-v : w)` and `P-Q = (u+v : w)` in (X:Z) coordinates.
#[inline]
pub fn jac_to_xz_add_components<L: FpBackend>(
    p: &JacPoint<L>,
    q: &JacPoint<L>,
    ac: &EcCurve<L>,
) -> AddComponents<L> {
    let t0 = p.z.sqr(); // z1^2
    let t1 = q.z.sqr(); // z2^2
    let t2 = p.x.mul(&t1); // x1*z2^2
    let t3 = t0.mul(&q.x); // z1^2*x2
    let t4 = p.y.mul(&q.z);
    let t4 = t4.mul(&t1); // y1*z2^3
    let t5 = p.z.mul(&q.y);
    let t5 = t5.mul(&t0); // z1^3*y2
    let t0 = t0.mul(&t1); // (z1*z2)^2
    let t6 = t4.mul(&t5); // (z1*z2)^3*y1*y2
    let v = t6.add(&t6); // 2*(z1*z2)^3*y1*y2
    let t4 = t4.sqr(); // y1^2*z2^6
    let t5 = t5.sqr(); // z1^6*y2^2
    let t4 = t4.add(&t5);
    let t5 = t2.add(&t3); // x1*z2^2 + z1^2*x2
    let t6 = t3.add(&t3); // 2*z1^2*x2
    let t6 = t5.sub(&t6); // lambda = x1*z2^2 - z1^2*x2
    let t6 = t6.sqr(); // lambda^2
    let t1 = ac.a.mul(&t0); // A*(z1*z2)^2
    let t1 = t5.add(&t1); // gamma
    let t1 = t1.mul(&t6); // gamma*lambda^2
    let u = t4.sub(&t1);
    let w = t6.mul(&t0); // (z1*z2)^2 * lambda^2

    AddComponents { u, v, w }
}
