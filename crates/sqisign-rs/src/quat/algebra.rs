//! Arithmetic of quaternion algebra elements.

use super::{QuatAlg, QuatAlgElem, QuatLattice, Vec4};
use crate::mp::Ibz;

impl<const N: usize> Vec4<N> {
    /// Product of two coordinate vectors in the basis `(1, i, j, ij)` with
    /// `i^2 = -1`, `j^2 = -p`.
    pub fn coord_mul(&self, b: &Self, alg: &QuatAlg<N>) -> Self {
        let a = &self.0;
        let b = &b.0;
        let mut s0 = Ibz::zero();
        s0 = s0.sub(&a[2].mul(&b[2]));
        s0 = s0.sub(&a[3].mul(&b[3]));
        s0 = s0.mul(&alg.p);
        s0 = s0.add(&a[0].mul(&b[0]));
        s0 = s0.sub(&a[1].mul(&b[1]));
        let mut s1 = Ibz::zero();
        s1 = s1.add(&a[2].mul(&b[3]));
        s1 = s1.sub(&a[3].mul(&b[2]));
        s1 = s1.mul(&alg.p);
        s1 = s1.add(&a[0].mul(&b[1]));
        s1 = s1.add(&a[1].mul(&b[0]));
        let mut s2 = Ibz::zero();
        s2 = s2.add(&a[0].mul(&b[2]));
        s2 = s2.add(&a[2].mul(&b[0]));
        s2 = s2.sub(&a[1].mul(&b[3]));
        s2 = s2.add(&a[3].mul(&b[1]));
        let mut s3 = Ibz::zero();
        s3 = s3.add(&a[0].mul(&b[3]));
        s3 = s3.add(&a[3].mul(&b[0]));
        s3 = s3.sub(&a[2].mul(&b[1]));
        s3 = s3.add(&a[1].mul(&b[2]));
        Self([s0, s1, s2, s3])
    }
}

impl<const N: usize> QuatAlgElem<N> {
    /// Bring two elements to the same denominator.
    pub fn equal_denom(a: &Self, b: &Self) -> (Self, Self) {
        let gcd = a.denom.gcd(&b.denom);
        let (da, _) = a.denom.div(&gcd);
        let (db, _) = b.denom.div(&gcd);
        let mut ra = *a;
        let mut rb = *b;
        for i in 0..4 {
            ra.coord.0[i] = a.coord.0[i].mul(&db);
            rb.coord.0[i] = b.coord.0[i].mul(&da);
        }
        let d = da.mul(&db).mul(&gcd);
        ra.denom = d;
        rb.denom = d;
        (ra, rb)
    }

    /// `a + b`.
    pub fn add(&self, b: &Self) -> Self {
        let (ra, rb) = Self::equal_denom(self, b);
        Self {
            denom: ra.denom,
            coord: ra.coord.add(&rb.coord),
        }
    }

    /// `a - b`.
    pub fn sub(&self, b: &Self) -> Self {
        let (ra, rb) = Self::equal_denom(self, b);
        Self {
            denom: ra.denom,
            coord: ra.coord.sub(&rb.coord),
        }
    }

    /// `a * b` (in this order).
    pub fn mul(&self, b: &Self, alg: &QuatAlg<N>) -> Self {
        Self {
            denom: self.denom.mul(&b.denom),
            coord: self.coord.coord_mul(&b.coord, alg),
        }
    }

    /// Reduced norm as `(numerator, denominator)`, both positive, in lowest
    /// terms.
    pub fn norm(&self, alg: &QuatAlg<N>) -> (Ibz<N>, Ibz<N>) {
        let c = &self.coord.0;
        let mut norm = c[3].mul(&c[3]).add(&c[2].mul(&c[2]));
        norm = norm.mul(&alg.p);
        norm = norm.add(&c[0].mul(&c[0]));
        norm = norm.add(&c[1].mul(&c[1]));
        let nd = self.denom.mul(&self.denom);
        let g = norm.gcd(&nd);
        let (num, _) = norm.div(&g);
        let (den, _) = nd.div(&g);
        (num.abs(), den.abs())
    }

    /// Standard involution `x -> conj(x)`.
    pub fn conj(&self) -> Self {
        Self {
            denom: self.denom,
            coord: Vec4([
                self.coord.0[0],
                self.coord.0[1].neg(),
                self.coord.0[2].neg(),
                self.coord.0[3].neg(),
            ]),
        }
    }

    /// `scalar * elem`, not normalised.
    pub fn scalar_mul(&self, scalar: &Ibz<N>) -> Self {
        Self {
            denom: self.denom,
            coord: Vec4([
                self.coord.0[0].mul(scalar),
                self.coord.0[1].mul(scalar),
                self.coord.0[2].mul(scalar),
                self.coord.0[3].mul(scalar),
            ]),
        }
    }

    /// Divide out the content of the coordinates and the denominator, and
    /// make the denominator positive.
    pub fn normalize(&mut self) {
        let mut gcd = self.coord.content();
        gcd = gcd.gcd(&self.denom);
        let (d, _) = self.denom.div(&gcd);
        self.denom = d;
        let _ = self.coord.scalar_div_inplace(&gcd);
        let sign = Ibz::set(if Ibz::zero() > self.denom { -1 } else { 1 }, 2);
        self.coord = self.coord.scalar_mul(&sign);
        self.denom = self.denom.mul(&sign);
    }

    /// Whether the element is zero.
    pub fn is_zero(&self) -> bool {
        self.coord.is_zero()
    }

    /// Whether two representations denote the same element.
    pub fn equals(&self, b: &Self) -> bool {
        self.sub(b).is_zero()
    }

    /// Write `x` in `order` as `content * primitive` with `primitive` a
    /// primitive coordinate vector of `order`'s basis. `x` must lie in
    /// `order`, which must have a triangular basis.
    pub fn make_primitive(&self, order: &QuatLattice<N>) -> (Vec4<N>, Ibz<N>) {
        let (ok, mut prim) = order.contains(self);
        debug_assert!(ok);
        let content = prim.content();
        if !content.is_zero() {
            for v in prim.0.iter_mut() {
                let (q, _) = v.div(&content);
                *v = q;
            }
        }
        (prim, content)
    }
}
