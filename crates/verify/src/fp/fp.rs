//!
//! Methods forward to the per-level [`super::FpBackend`] implementation.
//! Downstream code writes `where L: FpBackend` and gets identical Rust
//! APIs on top of whatever limb layout each level chooses.

use super::{Fp, FpBackend};
use hybrid_array::Array;
use subtle::{Choice, ConstantTimeEq};

impl<L: FpBackend> Default for Fp<L> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<L: FpBackend> Fp<L> {
    /// The zero element.
    #[inline]
    pub fn zero() -> Self {
        let mut limbs = Array::<u64, L::FpLimbs>::default();
        L::set_zero(&mut limbs);
        Self { limbs }
    }

    /// The multiplicative identity (internal Montgomery form).
    #[inline]
    pub fn one() -> Self {
        let mut limbs = Array::<u64, L::FpLimbs>::default();
        L::set_one(&mut limbs);
        Self { limbs }
    }

    /// Construct from a raw u64 limb slice already in Montgomery form.
    #[inline]
    pub fn from_limbs(limbs: &[u64]) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        debug_assert_eq!(limbs.len(), out.len());
        out.as_mut_slice().copy_from_slice(limbs);
        Self { limbs: out }
    }

    /// Construct from a small integer, in internal Montgomery form.
    #[inline]
    pub fn from_small(val: u64) -> Self {
        let mut limbs = Array::<u64, L::FpLimbs>::default();
        L::set_small(&mut limbs, val);
        Self { limbs }
    }

    /// Constant-time equality test.
    #[inline]
    pub fn ct_equal(&self, other: &Self) -> Choice {
        L::is_equal(&self.limbs, &other.limbs)
    }

    /// Constant-time zero test.
    #[inline]
    pub fn ct_is_zero(&self) -> Choice {
        L::is_zero(&self.limbs)
    }

    /// `self + rhs mod 2p`.
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::add(&mut out, &self.limbs, &rhs.limbs);
        Self { limbs: out }
    }

    /// `self - rhs mod 2p`.
    #[inline]
    pub fn sub(&self, rhs: &Self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::sub(&mut out, &self.limbs, &rhs.limbs);
        Self { limbs: out }
    }

    /// `-self mod 2p`.
    #[inline]
    pub fn neg(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::neg(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// `self * rhs mod 2p` (Montgomery multiplication).
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::mul(&mut out, &self.limbs, &rhs.limbs);
        Self { limbs: out }
    }

    /// `self^2 mod 2p` (specialized squaring).
    #[inline]
    pub fn sqr(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::sqr(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// `self * val mod 2p` for a small integer `val`.
    #[inline]
    pub fn mul_small(&self, val: u32) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::mul_small(&mut out, &self.limbs, val);
        Self { limbs: out }
    }

    /// `1 / self mod p`. Returns zero if `self == 0` (no panic on
    /// zero input).
    #[inline]
    pub fn inv(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::inv(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// Square root in Fp. The result is well-defined up to sign;
    /// callers that need to verify squareness should use
    /// [`Self::is_square`] separately.
    #[inline]
    pub fn sqrt(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::sqrt(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// Returns `Choice(1)` if `self` is a quadratic residue (or zero).
    #[inline]
    pub fn is_square(&self) -> Choice {
        L::is_square(&self.limbs)
    }

    /// `self / 2 mod p`.
    #[inline]
    pub fn half(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::half(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// `self / 3 mod p`.
    #[inline]
    pub fn div3(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::div3(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// `self^((p-3)/4) mod p` (square root progenitor).
    #[inline]
    pub fn exp3div4(&self) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::exp3div4(&mut out, &self.limbs);
        Self { limbs: out }
    }

    /// Serialize to canonical little-endian bytes. The returned array
    /// has length `L::FpEncodedBytes`.
    #[inline]
    pub fn encode(&self) -> Array<u8, L::FpEncodedBytes> {
        let mut out = Array::<u8, L::FpEncodedBytes>::default();
        L::encode(&mut out, &self.limbs);
        out
    }

    /// Deserialize from canonical little-endian bytes. Returns `None`
    /// if `bytes` does not represent an integer in `[0, p)`.
    #[inline]
    pub fn decode(bytes: &[u8]) -> Option<Self> {
        let mut limbs = Array::<u64, L::FpLimbs>::default();
        let ok = L::decode(&mut limbs, bytes);
        if bool::from(ok) {
            Some(Self { limbs })
        } else {
            None
        }
    }

    /// Decode-with-reduce: deserialize a possibly-longer little-endian
    /// byte string by reducing it modulo `p`. Always succeeds.
    #[inline]
    pub fn decode_reduce(bytes: &[u8]) -> Self {
        let mut limbs = Array::<u64, L::FpLimbs>::default();
        L::decode_reduce(&mut limbs, bytes);
        Self { limbs }
    }

    /// Constant-time conditional swap. If `ctl` is set, swap
    /// `self` and `other`; otherwise no-op.
    #[inline]
    pub fn cswap(&mut self, other: &mut Self, ctl: Choice) {
        L::cswap(&mut self.limbs, &mut other.limbs, ctl);
    }

    /// Constant-time conditional select. Returns `a0` if `ctl` is
    /// clear, `a1` if `ctl` is set.
    #[inline]
    pub fn select(a0: &Self, a1: &Self, ctl: Choice) -> Self {
        let mut out = Array::<u64, L::FpLimbs>::default();
        L::select(&mut out, &a0.limbs, &a1.limbs, ctl);
        Self { limbs: out }
    }
}

impl<L: FpBackend> ConstantTimeEq for Fp<L> {
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.ct_equal(other)
    }
}
