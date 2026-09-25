//! Build-time arithmetic on `E0: y^2 = x^3 + x` in Jacobian coordinates over
//! any `FieldOps`. Variable time; generator only.

use crate::field::FieldOps;
use num_bigint::BigUint;

/// An affine point, `None` being the point at infinity.
pub type Affine<E> = Option<(E, E)>;

/// Jacobian point `(X : Y : Z)` with `x = X/Z^2`, `y = Y/Z^3`; `None` is infinity.
#[derive(Clone)]
pub struct Jac<E> {
    pub x: E,
    pub y: E,
    pub z: E,
}

pub struct Curve<'a, F: FieldOps> {
    pub f: &'a F,
}

impl<'a, F: FieldOps> Curve<'a, F> {
    pub fn new(f: &'a F) -> Self {
        Self { f }
    }

    pub fn lift(&self, p: &Affine<F::E>) -> Option<Jac<F::E>> {
        p.as_ref().map(|(x, y)| Jac {
            x: x.clone(),
            y: y.clone(),
            z: self.f.one(),
        })
    }

    pub fn to_affine(&self, p: &Option<Jac<F::E>>) -> Affine<F::E> {
        let p = p.as_ref()?;
        let zi = self.f.inv(&p.z);
        let zi2 = self.f.sqr(&zi);
        let zi3 = self.f.mul(&zi2, &zi);
        Some((self.f.mul(&p.x, &zi2), self.f.mul(&p.y, &zi3)))
    }

    /// Does `(x, y)` satisfy `y^2 = x^3 + x`?
    pub fn on_curve(&self, x: &F::E, y: &F::E) -> bool {
        let lhs = self.f.sqr(y);
        let x3 = self.f.mul(&self.f.sqr(x), x);
        let rhs = self.f.add(&x3, x);
        self.f.eq(&lhs, &rhs)
    }

    pub fn neg(&self, p: &Option<Jac<F::E>>) -> Option<Jac<F::E>> {
        p.as_ref().map(|p| Jac {
            x: p.x.clone(),
            y: self.f.neg(&p.y),
            z: p.z.clone(),
        })
    }

    /// Doubling on `y^2 = x^3 + a x` with `a = 1` (dbl-2007-bl).
    pub fn double(&self, p: &Option<Jac<F::E>>) -> Option<Jac<F::E>> {
        let f = self.f;
        let p = p.as_ref()?;
        if f.is_zero(&p.y) {
            return None;
        }
        let xx = f.sqr(&p.x);
        let yy = f.sqr(&p.y);
        let yyyy = f.sqr(&yy);
        let zz = f.sqr(&p.z);
        // S = 2((X + YY)^2 - XX - YYYY)
        let t = f.sqr(&f.add(&p.x, &yy));
        let s = f.double(&f.sub(&f.sub(&t, &xx), &yyyy));
        // M = 3 XX + a ZZ^2, a = 1
        let m = f.add(&f.add(&f.double(&xx), &xx), &f.sqr(&zz));
        // X3 = M^2 - 2S
        let x3 = f.sub(&f.sqr(&m), &f.double(&s));
        // Y3 = M (S - X3) - 8 YYYY
        let eight_yyyy = f.double(&f.double(&f.double(&yyyy)));
        let y3 = f.sub(&f.mul(&m, &f.sub(&s, &x3)), &eight_yyyy);
        // Z3 = (Y + Z)^2 - YY - ZZ
        let z3 = f.sub(&f.sub(&f.sqr(&f.add(&p.y, &p.z)), &yy), &zz);
        Some(Jac {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    /// Addition (add-2007-bl), falling back to doubling / infinity when the
    /// inputs coincide up to sign.
    pub fn add(&self, p: &Option<Jac<F::E>>, q: &Option<Jac<F::E>>) -> Option<Jac<F::E>> {
        let f = self.f;
        let (p, q) = match (p, q) {
            (None, _) => return q.clone(),
            (_, None) => return p.clone(),
            (Some(p), Some(q)) => (p, q),
        };
        let z1z1 = f.sqr(&p.z);
        let z2z2 = f.sqr(&q.z);
        let u1 = f.mul(&p.x, &z2z2);
        let u2 = f.mul(&q.x, &z1z1);
        let s1 = f.mul(&f.mul(&p.y, &q.z), &z2z2);
        let s2 = f.mul(&f.mul(&q.y, &p.z), &z1z1);
        let h = f.sub(&u2, &u1);
        let r = f.double(&f.sub(&s2, &s1));
        if f.is_zero(&h) {
            return if f.is_zero(&r) {
                self.double(&Some(p.clone()))
            } else {
                None
            };
        }
        let i = f.sqr(&f.double(&h));
        let j = f.mul(&h, &i);
        let v = f.mul(&u1, &i);
        let x3 = f.sub(&f.sub(&f.sqr(&r), &j), &f.double(&v));
        let y3 = f.sub(&f.mul(&r, &f.sub(&v, &x3)), &f.double(&f.mul(&s1, &j)));
        let z3 = f.mul(&f.sub(&f.sub(&f.sqr(&f.add(&p.z, &q.z)), &z1z1), &z2z2), &h);
        Some(Jac {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    pub fn sub(&self, p: &Option<Jac<F::E>>, q: &Option<Jac<F::E>>) -> Option<Jac<F::E>> {
        self.add(p, &self.neg(q))
    }

    pub fn eq(&self, p: &Option<Jac<F::E>>, q: &Option<Jac<F::E>>) -> bool {
        let f = self.f;
        match (p, q) {
            (None, None) => true,
            (None, _) | (_, None) => false,
            (Some(p), Some(q)) => {
                let z1z1 = f.sqr(&p.z);
                let z2z2 = f.sqr(&q.z);
                let lhs_x = f.mul(&p.x, &z2z2);
                let rhs_x = f.mul(&q.x, &z1z1);
                let lhs_y = f.mul(&f.mul(&p.y, &q.z), &z2z2);
                let rhs_y = f.mul(&f.mul(&q.y, &p.z), &z1z1);
                f.eq(&lhs_x, &rhs_x) && f.eq(&lhs_y, &rhs_y)
            }
        }
    }

    /// `[k] P` by double-and-add.
    pub fn mul(&self, p: &Option<Jac<F::E>>, k: &BigUint) -> Option<Jac<F::E>> {
        let mut acc: Option<Jac<F::E>> = None;
        for i in (0..k.bits()).rev() {
            acc = self.double(&acc);
            if k.bit(i) {
                acc = self.add(&acc, p);
            }
        }
        acc
    }

    /// `[2^n] P`.
    pub fn mul_pow2(&self, p: &Option<Jac<F::E>>, n: u64) -> Option<Jac<F::E>> {
        let mut acc = p.clone();
        for _ in 0..n {
            acc = self.double(&acc);
        }
        acc
    }
}
