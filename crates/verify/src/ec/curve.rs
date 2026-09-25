//! Curve-level helpers: validity of `A`, normalisation of `(A : C)` and of
//! `A24`, and the j-invariant.

use super::{EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};

impl<L: FpBackend> EcCurve<L> {
    /// `E0: (A : C) = (0 : 1)`.
    #[inline]
    pub fn init() -> Self {
        Self::default()
    }

    /// A curve from its Montgomery coefficient with `C = 1`, or `None` when
    /// `A = +-2` (singular).
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

    /// `A^2 - 4 != 0`.
    #[inline]
    pub fn verify_a(a: &Fp2<L>) -> bool {
        let two = Fp2::<L>::from_small(2);
        let is_two = a.ct_equal(&two);
        let is_minus_two = a.ct_equal(&two.neg());
        !bool::from(is_two | is_minus_two)
    }

    /// `(A : C) -> (A/C : 1)`.
    #[inline]
    pub fn normalize(&mut self) {
        let c_inv = self.c.inv();
        self.a = self.a.mul(&c_inv);
        self.c = Fp2::one();
    }

    /// `(A + 2C : 4C)` from `(A : C)`, or the cached normalised `A24`.
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

    /// The curve `(A : C) = (4 A24 - 2 C24 : C24)` from `(A24 : C24)`.
    #[inline]
    pub fn from_a24(a24: &EcPoint<L>) -> Self {
        let mut a = a24.x.add(&a24.x);
        a = a.sub(&a24.z);
        a = a.add(&a);
        Self {
            a,
            c: a24.z.clone(),
            a24: EcPoint::identity(),
            is_a24_computed_and_normalized: false,
        }
    }

    /// Cache `A24 = ((A + 2C)/4C : 1)`.
    #[inline]
    pub fn normalize_a24(&mut self) {
        if !self.is_a24_computed_and_normalized {
            let a24 = self.ac_to_a24();
            let z_inv = a24.z.inv();
            self.a24 = EcPoint {
                x: a24.x.mul(&z_inv),
                z: Fp2::one(),
            };
            self.is_a24_computed_and_normalized = true;
        }
    }

    /// Normalise `(A : C)` to `(A/C : 1)` and cache `A24 = ((A + 2)/4 : 1)`.
    #[inline]
    pub fn normalize_curve_and_a24(&mut self) {
        if !bool::from(self.c.ct_is_one()) {
            self.normalize();
        }
        if !self.is_a24_computed_and_normalized {
            // (A + 2) / 4 with C = 1
            let a24x = self.a.add_one().add_one().half().half();
            self.a24 = EcPoint {
                x: a24x,
                z: Fp2::one(),
            };
            self.is_a24_computed_and_normalized = true;
        }
    }

    /// The j-invariant `256 (A^2 - 3C^2)^3 / (C^4 (A^2 - 4C^2))`.
    #[inline]
    pub fn j_inv(&self) -> Fp2<L> {
        let c2 = self.c.sqr();
        let a2 = self.a.sqr();
        // t0 = A^2 - 3C^2
        let t0 = a2.sub(&c2.add(&c2)).sub(&c2);
        // denominator: (A^2 - 4C^2) C^4
        let den = t0.sub(&c2).mul(&c2.sqr());
        // numerator: 256 (A^2 - 3C^2)^3 = 4 (4 t0)^3
        let t = t0.add(&t0).add(&t0.add(&t0)); // 4 t0
        let num = t.mul(&t.sqr());
        let num = num.add(&num).add(&num.add(&num)); // 4 (4 t0)^3
        num.mul(&den.inv())
    }
}
