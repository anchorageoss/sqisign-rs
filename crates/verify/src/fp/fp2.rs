//!
//! `Fp2 = Fp[i] / (i^2 + 1)`. All operations are expressed in terms of
//! the per-level `Fp` operations from [`super::FpBackend`], so `Fp2`
//! itself is fully generic over the level. Multiplication uses
//! Karatsuba (three `Fp` multiplications). Square root follows the
//! algorithm from ePrint 2024/1563.

use super::{Fp, Fp2, FpBackend};
use hybrid_array::Array;
use subtle::{Choice, ConstantTimeEq, ConstantTimeLess};

impl<L: FpBackend> Default for Fp2<L> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<L: FpBackend> Fp2<L> {
    /// Construct from raw u64 limb slices already in Montgomery form.
    #[inline]
    pub fn from_limbs(re: &[u64], im: &[u64]) -> Self {
        Self {
            re: Fp::<L>::from_limbs(re),
            im: Fp::<L>::from_limbs(im),
        }
    }

    /// `0 + 0i`.
    #[inline]
    pub fn zero() -> Self {
        Self {
            re: Fp::<L>::zero(),
            im: Fp::<L>::zero(),
        }
    }

    /// `1 + 0i`.
    #[inline]
    pub fn one() -> Self {
        Self {
            re: Fp::<L>::one(),
            im: Fp::<L>::zero(),
        }
    }

    /// `0 + 1i` (the algebraic generator `i`).
    #[inline]
    pub fn i_element() -> Self {
        Self {
            re: Fp::<L>::zero(),
            im: Fp::<L>::one(),
        }
    }

    /// Construct `val + 0i` from a small integer.
    #[inline]
    pub fn from_small(val: u64) -> Self {
        Self {
            re: Fp::<L>::from_small(val),
            im: Fp::<L>::zero(),
        }
    }

    /// Constant-time zero check.
    #[inline]
    pub fn ct_is_zero(&self) -> Choice {
        self.re.ct_is_zero() & self.im.ct_is_zero()
    }

    /// Constant-time one check.
    #[inline]
    pub fn ct_is_one(&self) -> Choice {
        self.re.ct_equal(&Fp::<L>::one()) & self.im.ct_is_zero()
    }

    /// Constant-time equality.
    #[inline]
    pub fn ct_equal(&self, other: &Self) -> Choice {
        self.re.ct_equal(&other.re) & self.im.ct_equal(&other.im)
    }

    /// `self + rhs`.
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        Self {
            re: self.re.add(&rhs.re),
            im: self.im.add(&rhs.im),
        }
    }

    /// `self + 1` (specialized add that avoids constructing `Fp2::one()`
    /// for the right-hand side).
    #[inline]
    pub fn add_one(&self) -> Self {
        Self {
            re: self.re.add(&Fp::<L>::one()),
            im: self.im.clone(),
        }
    }

    /// `self - rhs`.
    #[inline]
    pub fn sub(&self, rhs: &Self) -> Self {
        Self {
            re: self.re.sub(&rhs.re),
            im: self.im.sub(&rhs.im),
        }
    }

    /// `-self`.
    #[inline]
    pub fn neg(&self) -> Self {
        Self {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }

    /// `self * rhs`: the backend's [`FpBackend::fp2_mul`], Karatsuba by
    /// default (three `Fp` multiplications plus five adds/subs), fused
    /// products in the assembly backend.
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = Self::zero();
        L::fp2_mul_pair(&mut out, self, rhs);
        out
    }

    /// `self^2`: the backend's [`FpBackend::fp2_sqr`], by default the
    /// factored identity `re = (y.re + y.im)(y.re - y.im)`, `im = 2 y.re
    /// y.im`, two `Fp` multiplications instead of three.
    #[inline]
    pub fn sqr(&self) -> Self {
        let mut out = Self::zero();
        L::fp2_sqr(
            &mut out.re.limbs,
            &mut out.im.limbs,
            &self.re.limbs,
            &self.im.limbs,
        );
        out
    }

    /// `self <- self * rhs` in place (P28).
    #[inline]
    pub fn mul_assign(&mut self, rhs: &Self) {
        L::fp2_mul_assign(self, rhs);
    }

    /// `self <- self^2` in place (P28).
    #[inline]
    pub fn sqr_assign(&mut self) {
        L::fp2_sqr_assign(self);
    }

    /// `self <- self + rhs` in place (P28).
    #[inline]
    pub fn add_assign(&mut self, rhs: &Self) {
        self.re.add_assign(&rhs.re);
        self.im.add_assign(&rhs.im);
    }

    /// `self <- self - rhs` in place (P28).
    #[inline]
    pub fn sub_assign(&mut self, rhs: &Self) {
        self.re.sub_assign(&rhs.re);
        self.im.sub_assign(&rhs.im);
    }

    /// Multiply by a small (32-bit) integer.
    #[inline]
    pub fn mul_small(&self, val: u32) -> Self {
        Self {
            re: self.re.mul_small(val),
            im: self.im.mul_small(val),
        }
    }

    /// `self / 2`.
    #[inline]
    pub fn half(&self) -> Self {
        Self {
            re: self.re.half(),
            im: self.im.half(),
        }
    }

    /// Conjugate `a + bi -> a - bi`.
    #[inline]
    pub fn conjugate(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: self.im.neg(),
        }
    }

    /// Multiplicative inverse: `(a + bi)^-1 = (a - bi) / (a^2 + b^2)`.
    #[inline]
    pub fn inv(&self) -> Self {
        // t0 = re^2 + im^2
        let t_re2 = self.re.sqr();
        let t_im2 = self.im.sqr();
        let norm = t_re2.add(&t_im2);
        let n_inv = norm.inv();
        let new_re = self.re.mul(&n_inv);
        let new_im = self.im.mul(&n_inv).neg();
        Self {
            re: new_re,
            im: new_im,
        }
    }

    /// Returns `Choice(1)` if `self` is a square in 𝔽p². Equivalent to
    /// checking whether `re^2 + im^2` is a square in Fp.
    #[inline]
    pub fn is_square(&self) -> Choice {
        let t_re2 = self.re.sqr();
        let t_im2 = self.im.sqr();
        let norm = t_re2.add(&t_im2);
        norm.is_square()
    }

    /// Square root in 𝔽p² by the algorithm of Aardal et al. (ePrint
    /// 2024/1563), computed exactly as the SQIsign round-3 reference does:
    /// the result is one of the two roots, chosen deterministically by the
    /// algorithm, with no further sign normalisation. Only meaningful when
    /// `self` is a square.
    #[inline]
    pub fn sqrt(&self) -> Self {
        // x0 = delta = sqrt(a.re^2 + a.im^2)
        let re2 = self.re.sqr();
        let im2 = self.im.sqr();
        let mut x0 = re2.add(&im2);
        x0 = x0.sqrt();
        // If a.im = 0, restore x0 = a.re to avoid x0 + a.re collapsing
        // to zero when delta = -a.re.
        let im_is_zero = self.im.ct_is_zero();
        x0 = Fp::<L>::select(&x0, &self.re, im_is_zero);
        // x0 = delta + a.re;   t0 = 2 * x0
        x0 = x0.add(&self.re);
        let t0_first = x0.add(&x0);

        // x1 = t0^((p-3)/4)
        let mut x1 = t0_first.exp3div4();

        // x0 = x0 * x1;   x1 = x1 * a.im;   t1 = (2*x0)^2
        x0 = x0.mul(&x1);
        x1 = x1.mul(&self.im);
        let two_x0 = x0.add(&x0);
        let t1 = two_x0.sqr();

        // If t1 == t0 return x0 + x1*i, else x1 - x0*i
        let f = t0_first.sub(&t1).ct_is_zero();
        let re = Fp::<L>::select(&x1, &x0, f);
        let im = Fp::<L>::select(&x0.neg(), &x1, f);
        Self { re, im }
    }

    /// The SQIsign round-2 reference's square root: [`Self::sqrt`] followed
    /// by a sign normalisation so that the real part is even when non-zero,
    /// and otherwise the imaginary part is even. Kept for byte-exact
    /// comparison with code built on the round-2 reference (PRISM_v2); the
    /// round-3 protocol layers use [`Self::sqrt`].
    #[inline]
    pub fn sqrt_canonical_even(&self) -> Self {
        let s = self.sqrt();
        let re_is_zero = s.re.ct_is_zero();
        let re_is_odd = lsb_choice(s.re.encode()[0]);
        let im_is_odd = lsb_choice(s.im.encode()[0]);
        let negate = re_is_odd | (re_is_zero & im_is_odd);
        Self {
            re: Fp::<L>::select(&s.re, &s.re.neg(), negate),
            im: Fp::<L>::select(&s.im, &s.im.neg(), negate),
        }
    }

    /// Frobenius `a + bi -> a - bi` (the same map as [`Self::conjugate`],
    /// under the name the round-3 reference uses).
    #[inline]
    pub fn frob(&self) -> Self {
        self.conjugate()
    }

    /// `i * self` when `ctl` is set, `-i * self` when it is clear. Mirrors
    /// the round-3 reference's `fp2_mul_by_i`, which always multiplies by a
    /// fourth root of unity and uses `ctl` only to pick its sign.
    #[inline]
    pub fn mul_by_i(&self, ctl: Choice) -> Self {
        Self {
            re: Fp::<L>::select(&self.im, &self.im.neg(), ctl),
            im: Fp::<L>::select(&self.re.neg(), &self.re, ctl),
        }
    }

    /// Constant-time `self < other` on canonical encodings, compared as the
    /// little-endian integer `encode(re) || encode(im)`: the imaginary part
    /// is the most significant. This is the round-3 reference's
    /// `fp2_less_than` ordering (used for canonical curve models). It is not
    /// the real-part-first ordering PRISM's Section 5 defines; that one is
    /// built where it is needed.
    #[inline]
    pub fn less_than(&self, other: &Self) -> Choice {
        let a = self.encode();
        let b = other.encode();
        let mut less = Choice::from(0);
        let mut equal_so_far = Choice::from(1);
        for (x, y) in a.iter().rev().zip(b.iter().rev()) {
            less |= equal_so_far & x.ct_lt(y);
            equal_so_far &= x.ct_eq(y);
        }
        less
    }

    /// Replace `self` with `sqrt(self)` and return `Choice(1)` if
    /// the original value was indeed a square (verified by squaring
    /// the result and comparing). The replacement happens
    /// unconditionally.
    #[inline]
    pub fn sqrt_verify(&mut self) -> Choice {
        let original = self.clone();
        let s = self.sqrt();
        let check = s.sqr();
        *self = s;
        original.ct_equal(&check)
    }

    /// Montgomery batch inversion: computes the inverse of each
    /// element in `x` in place, using a single 𝔽p² inversion plus
    /// `3 * (len - 1)` 𝔽p² multiplications.
    ///
    /// Caller must provide scratch slices `t1` and `t2` of the same
    /// length as `x`. This avoids heap allocation, keeping the crate
    /// `no_std`-compatible.
    ///
    #[inline]
    pub fn batched_inv(x: &mut [Self], t1: &mut [Self], t2: &mut [Self]) {
        let len = x.len();
        debug_assert_eq!(t1.len(), len);
        debug_assert_eq!(t2.len(), len);
        if len == 0 {
            return;
        }
        // t1[i] = x[0] * x[1] * ... * x[i]
        t1[0] = x[0].clone();
        for i in 1..len {
            t1[i] = t1[i - 1].mul(&x[i]);
        }
        // inverse = 1 / (x[0] * ... * x[len-1])
        let inverse = t1[len - 1].inv();
        // t2[0] = inverse;
        // t2[i] = t2[i-1] * x[len-i]   (so t2[i] = 1 / (x[0] * ... * x[len-1-i]))
        t2[0] = inverse;
        for i in 1..len {
            t2[i] = t2[i - 1].mul(&x[len - i]);
        }
        // x[0] = t2[len-1] = 1 / x[0]
        x[0] = t2[len - 1].clone();
        // x[i] = t1[i-1] * t2[len-i-1] = (x[0]*..*x[i-1]) * (1 / (x[0]*..*x[i])) = 1/x[i]
        for i in 1..len {
            x[i] = t1[i - 1].mul(&t2[len - i - 1]);
        }
    }

    /// Variable-time square-and-multiply exponentiation. `exp` is a
    /// little-endian array of 64-bit limbs. Result is `self^exp` in GF(p²).
    ///
    /// **Not constant-time.** The control flow depends on the bits of
    /// `exp`; use only with non-secret exponents.
    #[inline]
    pub fn pow_vartime(&self, exp: &[u64]) -> Self {
        let mut acc = self.clone();
        let mut out = Self::one();
        for &word in exp {
            for i in 0..64 {
                let bit = (word >> i) & 1;
                if bit == 1 {
                    out = out.mul(&acc);
                }
                acc = acc.sqr();
            }
        }
        out
    }

    /// Serialize: `FpEncodedBytes` of `re` followed by `FpEncodedBytes` of `im`.
    #[inline]
    pub fn encode(&self) -> Array<u8, L::Fp2EncodedBytes> {
        let mut out = Array::<u8, L::Fp2EncodedBytes>::default();
        let n = L::FpEncodedBytes::USIZE;
        let re_bytes = self.re.encode();
        let im_bytes = self.im.encode();
        out[..n].copy_from_slice(re_bytes.as_ref());
        out[n..2 * n].copy_from_slice(im_bytes.as_ref());
        out
    }

    /// Deserialize from canonical bytes (`re` first, then `im`). Returns
    /// `None` if either component is out of range.
    #[inline]
    pub fn decode(bytes: &[u8]) -> Option<Self> {
        let n = L::FpEncodedBytes::USIZE;
        if bytes.len() != 2 * n {
            return None;
        }
        let re = Fp::<L>::decode(&bytes[..n])?;
        let im = Fp::<L>::decode(&bytes[n..])?;
        Some(Self { re, im })
    }

    /// Constant-time conditional swap.
    #[inline]
    pub fn cswap(&mut self, other: &mut Self, ctl: Choice) {
        self.re.cswap(&mut other.re, ctl);
        self.im.cswap(&mut other.im, ctl);
    }

    /// Constant-time conditional select.
    #[inline]
    pub fn select(a0: &Self, a1: &Self, ctl: Choice) -> Self {
        Self {
            re: Fp::<L>::select(&a0.re, &a1.re, ctl),
            im: Fp::<L>::select(&a0.im, &a1.im, ctl),
        }
    }
}

impl<L: FpBackend> ConstantTimeEq for Fp2<L> {
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.ct_equal(other)
    }
}

#[inline]
fn lsb_choice(b: u8) -> Choice {
    Choice::from(b & 1)
}

// Re-export the `USIZE` const from `hybrid_array::ArraySize` so that
// `L::FpEncodedBytes::USIZE` resolves above. (typenum::Unsigned is the
// underlying trait.)
use typenum::Unsigned as _;
