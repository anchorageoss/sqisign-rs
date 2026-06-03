use super::{EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};
use subtle::Choice;

impl<L: FpBackend> EcCurve<L> {
    /// Initialize a curve with `(A:C) = (0:1)` and uncached A24.
    #[inline]
    pub fn init() -> Self {
        Self::default()
    }

    /// Initialize a curve from a Montgomery coefficient A with `C = 1`.
    /// Returns `None` if `A^2 - 4 = 0` (i.e. A = +-2).
    #[inline]
    pub fn from_a(a: &Fp2<L>) -> Option<Self> {
        if !Self::verify_a(a) {
            return None;
        }
        Some(Self {
            a: a.clone(),
            ..Self::default()
        })
    }

    /// Check that A is a valid Montgomery coefficient (A^2 - 4 != 0).
    #[inline]
    pub fn verify_a(a: &Fp2<L>) -> bool {
        let mut t = Fp2::<L>::one();
        t.re = t.re.add(&t.re); // t = 2
        if bool::from(a.ct_equal(&t)) {
            return false;
        }
        t.re = t.re.neg(); // t = -2
        if bool::from(a.ct_equal(&t)) {
            return false;
        }
        true
    }

    /// Reduce `(A : C)` to `(A/C : 1)` in place.
    #[inline]
    pub fn normalize(&mut self) {
        let c_inv = self.c.inv();
        self.a = self.a.mul(&c_inv);
        self.c = Fp2::one();
    }

    /// Compute `(A+2C : 4C)` from the curve coefficients.
    #[inline]
    pub fn ac_to_a24(&self) -> EcPoint<L> {
        if self.is_a24_computed_and_normalized {
            return self.a24.clone();
        }
        let z = self.c.add(&self.c);
        let x = self.a.add(&z);
        let z = z.add(&z);
        EcPoint { x, z }
    }

    /// Given `(A+2C : 4C)`, recover `(A : C)`.
    #[inline]
    pub fn from_a24(a24: &EcPoint<L>) -> Self {
        // (A:C) = (4*(A+2C) - 2*4C : 4C) = (4*a24.x - 2*a24.z : a24.z)
        let mut a = a24.x.add(&a24.x);
        a = a.sub(&a24.z);
        a = a.add(&a);
        let c = a24.z.clone();
        Self {
            a,
            c,
            a24: EcPoint::identity(),
            is_a24_computed_and_normalized: false,
        }
    }

    /// Compute and cache `A24 = ((A+2C)/(4C) : 1)` if not already done.
    #[inline]
    pub fn normalize_a24(&mut self) {
        if !self.is_a24_computed_and_normalized {
            self.a24 = self.ac_to_a24();
            let z_inv = self.a24.z.inv();
            self.a24.x = self.a24.x.mul(&z_inv);
            self.a24.z = Fp2::one();
            self.is_a24_computed_and_normalized = true;
        }
    }

    /// Normalize both `(A:C)` and A24 in place.
    #[inline]
    pub fn normalize_curve_and_a24(&mut self) {
        if !bool::from(self.c.ct_is_one()) {
            self.normalize();
        }

        if !self.is_a24_computed_and_normalized {
            // A24 = ((A + 2) / 4 : 1) since C = 1 after normalization
            let mut a24x = self.a.add_one();
            a24x = a24x.add_one();
            // Preserve the imaginary part
            let mut result = a24x;
            result.im = self.a.im.clone();
            result = result.half();
            result = result.half();
            self.a24.x = result;
            self.a24.z = Fp2::one();
            self.is_a24_computed_and_normalized = true;
        }
    }

    /// Compute the j-invariant `j(E)` from the curve `(A:C)`.
    #[inline]
    pub fn j_inv(&self) -> Fp2<L> {
        let mut t0;
        let mut t1;

        t1 = self.c.sqr();
        let mut j_inv = self.a.sqr();
        t0 = t1.add(&t1);
        t0 = j_inv.sub(&t0);
        t0 = t0.sub(&t1);
        j_inv = t0.sub(&t1);
        t1 = t1.sqr();
        j_inv = j_inv.mul(&t1);
        t0 = t0.add(&t0);
        t0 = t0.add(&t0);
        t1 = t0.sqr();
        t0 = t0.mul(&t1);
        t0 = t0.add(&t0);
        t0 = t0.add(&t0);
        j_inv = j_inv.inv();
        j_inv = t0.mul(&j_inv);
        j_inv
    }

    /// Recover a y-coordinate from x on the curve `y^2 = x^3 + (A/C)x^2 + x`.
    /// Assumes `C = 1`. Returns `(y, is_on_curve)` where `is_on_curve`
    /// indicates whether x is on the curve (y is a valid square root).
    #[inline]
    pub fn recover_y(&self, px: &Fp2<L>) -> (Fp2<L>, Choice) {
        let t0 = px.sqr();
        let mut y = t0.mul(&self.a); // Ax^2
        y = y.add(px); // Ax^2 + x
        let t0 = t0.mul(px);
        y = y.add(&t0); // x^3 + Ax^2 + x
        let valid = y.sqrt_verify();
        (y, valid)
    }
}
