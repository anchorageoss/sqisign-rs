//! Fixed-precision signed integers, a port of the SQIsign round-3 reference's
//! `src/mp` library (spec Section 4.1).
//!
//! An [`Ibz`] is `N` 64-bit limbs in two's complement together with a
//! public `bitlen` bound: the value lies strictly between `-2^(bitlen-1)`
//! and `2^(bitlen-1)`. Every operation reads and writes only the limbs the
//! bounds of its operands cover, so running time depends on the declared
//! bounds (which are public by construction) and not on the values, except
//! in the functions documented as variable time (division, gcd, primality,
//! anything with `vartime` in its name). The container size `N` is fixed
//! per parameter set: `5 * ceil(log2 p / 64)` limbs, as in the reference,
//! so no honest computation ever reaches the top of the container. Where
//! the reference documents a bound, the port asserts it in debug builds.
//!
//! No heap, no floats, no external big-integer crate.

// The code below is a limb-by-limb port of C; index loops, explicit range
// checks and long argument lists mirror the reference on purpose.
#![allow(
    clippy::needless_range_loop,
    clippy::manual_range_contains,
    clippy::manual_div_ceil,
    clippy::too_many_arguments,
    clippy::manual_clamp
)]

use core::cmp::Ordering;

mod arith;
pub mod ct;
mod div;
mod gcd;
mod modq;
mod prime;
mod rand;
mod sqrt;
mod string;

pub use ct::*;
pub(crate) use div::div_qlen as div_qlen_pub;
pub use rand::{DefaultDomain, Rng, ShakeRng, SplitRng};

/// Largest container this module supports; scratch buffers are sized from
/// it so that no operation needs `2 * N` as an array length.
pub const IBZ_MAX_LIMBS: usize = 64;

/// Bits per limb.
pub const LIMB_BITS: i32 = 64;

/// Number of limbs needed to hold `bits` bits.
#[inline]
pub const fn nlimbs(bits: i32) -> usize {
    ((bits + LIMB_BITS - 1) / LIMB_BITS) as usize
}

/// Round `bits` up to a whole limb.
#[inline]
pub const fn limb_align(bits: i32) -> i32 {
    (nlimbs(bits) as i32) * LIMB_BITS
}

/// A fixed-precision signed integer of at most `N` limbs.
#[derive(Clone, Copy)]
pub struct Ibz<const N: usize> {
    /// Bound on the two's complement bit length, sign bit included:
    /// `-2^(bitlen-1) < value < 2^(bitlen-1)`. Always at least 1.
    pub(crate) bitlen: i32,
    /// Little-endian limbs in two's complement. Only the first
    /// `nlimbs(bitlen)` are meaningful.
    pub(crate) limbs: [u64; N],
}

impl<const N: usize> Default for Ibz<N> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const N: usize> PartialEq for Ibz<N> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<const N: usize> Eq for Ibz<N> {}

impl<const N: usize> PartialOrd for Ibz<N> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for Ibz<N> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        Ibz::cmp_value(self, other)
    }
}

impl<const N: usize> core::fmt::Debug for Ibz<N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut buf = [0u8; 2 * IBZ_MAX_LIMBS * 16 + 2];
        let s = self.to_hex(&mut buf);
        write!(f, "Ibz({s}, bound {})", self.bitlen)
    }
}

impl<const N: usize> Ibz<N> {
    /// Full width of the container in bits: no value may exceed it.
    pub const MAX_BITS: i32 = (N as i32) * LIMB_BITS;

    const CHECK: () = assert!(N <= IBZ_MAX_LIMBS && N >= 2, "unsupported Ibz limb count");

    /// The integer 0, bound 1.
    #[inline]
    pub const fn zero() -> Self {
        #[allow(clippy::let_unit_value)]
        let _ = Self::CHECK;
        Self {
            bitlen: 1,
            limbs: [0; N],
        }
    }

    /// The integer 1, bound 2.
    #[inline]
    pub const fn one() -> Self {
        let mut z = Self::zero();
        z.limbs[0] = 1;
        z.bitlen = 2;
        z
    }

    /// The integer 2, bound 3.
    #[inline]
    pub const fn two() -> Self {
        let mut z = Self::zero();
        z.limbs[0] = 2;
        z.bitlen = 3;
        z
    }

    /// The integer 3, bound 3.
    #[inline]
    pub const fn three() -> Self {
        let mut z = Self::zero();
        z.limbs[0] = 3;
        z.bitlen = 3;
        z
    }

    /// Number of active limbs.
    #[inline]
    pub(crate) fn n(&self) -> usize {
        nlimbs(self.bitlen)
    }

    /// The current bound.
    #[inline]
    pub fn get_bound(&self) -> i32 {
        debug_assert!(self.bitlen > 0);
        self.bitlen
    }

    /// Interpret the lowest `bitlen` bits of `x` as a signed two's complement
    /// integer and set the bound to `bitlen` (clamped to `[1, 32]`). The
    /// value `-2^(bitlen-1)` is not representable and folds to zero.
    pub fn set(x: i64, bitlen: i32) -> Self {
        debug_assert!((0..=32).contains(&bitlen));
        let bitlen = if bitlen == 0 { 1 } else { bitlen.min(32) };
        let mut a = Self::zero();
        a.bitlen = bitlen;
        a.limbs[0] = x as u64;
        a.normalise();
        a
    }

    /// From a small unsigned value with a tight bound.
    #[inline]
    pub fn from_u32(x: u32) -> Self {
        let bits = 33 - (x.leading_zeros() as i32).min(32);
        Self::set(x as i64, bits.min(32).max(1))
    }

    /// From a little-endian non-negative unsigned limb array of `bitlen`
    /// value bits; the result has bound `bitlen + 1`. A `source` shorter
    /// than `bitlen` bits is zero-extended, a longer one is truncated to
    /// `bitlen` bits.
    ///
    /// # Panics
    ///
    /// If `bitlen` is negative or `bitlen + 1` does not fit the `N` limbs
    /// (P19 B2: an explicit check, not an out-of-bounds index).
    pub fn from_bits(source: &[u64], bitlen: i32) -> Self {
        assert!(
            bitlen >= 0 && bitlen < Self::MAX_BITS,
            "Ibz::from_bits: {bitlen} bits do not fit {N} limbs"
        );
        let n = nlimbs(bitlen);
        let mut t = Self::zero();
        for (d, s) in t.limbs[..n].iter_mut().zip(source.iter()) {
            *d = *s;
        }
        t.bitlen = bitlen + 1;
        let n1 = nlimbs(t.bitlen);
        if n1 > n {
            t.limbs[n1 - 1] = 0;
        } else {
            let j = bitlen - LIMB_BITS * (n1 as i32 - 1);
            t.limbs[n1 - 1] &= (1u64 << j) - 1;
        }
        t
    }

    /// From a little-endian unsigned limb array, all of it.
    ///
    /// # Panics
    ///
    /// If `source` has `N` limbs or more (the bound `64 len + 1` must fit).
    #[inline]
    pub fn from_digits(source: &[u64]) -> Self {
        Self::from_bits(source, (source.len() as i32) * LIMB_BITS)
    }

    /// From little-endian bytes of a non-negative integer.
    ///
    /// # Panics
    ///
    /// If `8 bytes.len() + 1` bits do not fit the `N` limbs.
    pub fn from_le_bytes(bytes: &[u8]) -> Self {
        assert!(
            8 * bytes.len() < Self::MAX_BITS as usize,
            "Ibz::from_le_bytes: {} bytes do not fit {N} limbs",
            bytes.len()
        );
        let mut limbs = [0u64; IBZ_MAX_LIMBS];
        for (i, b) in bytes.iter().enumerate() {
            limbs[i / 8] |= (*b as u64) << (8 * (i % 8));
        }
        let n = bytes.len().div_ceil(8);
        Self::from_bits(&limbs[..n], 8 * bytes.len() as i32)
    }

    /// The active limbs of a non-negative integer, little-endian, exactly
    /// enough for its bit size.
    pub fn to_digits(&self, target: &mut [u64]) {
        debug_assert!(self.is_positive());
        let n = nlimbs(self.bitsize());
        for t in target.iter_mut() {
            *t = 0;
        }
        target[..n].copy_from_slice(&self.limbs[..n]);
    }

    /// Low 32 bits as a signed integer.
    #[inline]
    pub fn get(&self) -> i32 {
        (self.limbs[0] & 0xFFFF_FFFF) as u32 as i32
    }

    /// Little-endian bytes of a non-negative value, `len` of them.
    pub fn to_le_bytes(&self, out: &mut [u8]) {
        debug_assert!(self.is_positive());
        let n = self.n();
        for (i, b) in out.iter_mut().enumerate() {
            let limb = i / 8;
            *b = if limb < n {
                (self.limbs[limb] >> (8 * (i % 8))) as u8
            } else {
                0
            };
        }
    }

    /// Sets bits `i+1..` of `a` to bit `i`.
    #[inline]
    fn propagate_bit(a: &mut u64, i: i32) {
        let v = 1u64 << i;
        *a = (*a & (v - 1)).wrapping_sub(*a & v);
    }

    /// Clear the bits above the bound in the top limb and sign-extend from
    /// the bound; `-2^(bitlen-1)` becomes zero. Variable time in the low
    /// limbs (the scan for zero).
    pub(crate) fn normalise(&mut self) {
        debug_assert!(self.bitlen > 0);
        let n = self.n();
        let b = self.bitlen - LIMB_BITS * (n as i32 - 1);
        debug_assert!(b > 0 && b <= LIMB_BITS);
        self.limbs[n - 1] &= u64::MAX >> (LIMB_BITS - b);
        let mut i = 0;
        while i < n - 1 && self.limbs[i] == 0 {
            i += 1;
        }
        if i == n {
            self.limbs[0] = 0;
            self.bitlen = 1;
            return;
        }
        Self::propagate_bit(&mut self.limbs[n - 1], b - 1);
    }

    /// Swap two integers.
    #[inline]
    pub fn swap(a: &mut Self, b: &mut Self) {
        core::mem::swap(a, b);
    }

    /// Swap when `cond` is all ones, leave alone when it is zero; touches
    /// every limb of the container.
    #[inline]
    pub fn cswap(a: &mut Self, b: &mut Self, cond: u64) {
        for i in 0..N {
            let t = cond & (a.limbs[i] ^ b.limbs[i]);
            a.limbs[i] ^= t;
            b.limbs[i] ^= t;
        }
        let t1 = (cond as i32) & (a.bitlen ^ b.bitlen);
        a.bitlen ^= t1;
        b.bitlen ^= t1;
    }

    /// Change the bound. Increasing sign-extends into the new limbs;
    /// decreasing truncates and re-normalises (variable time in the value).
    pub fn set_bound(&mut self, bitlen: i32) {
        debug_assert!(self.bitlen > 0);
        debug_assert!(bitlen >= 0 && bitlen <= Self::MAX_BITS);
        let bitlen = bitlen.clamp(0, Self::MAX_BITS);
        if bitlen == self.bitlen {
            return;
        }
        if bitlen > self.bitlen {
            let n0 = self.n();
            let n1 = nlimbs(bitlen);
            let bit = self.limbs[n0 - 1] >> (LIMB_BITS - 1);
            for i in n0..n1 {
                self.limbs[i] = 0u64.wrapping_sub(bit);
            }
            self.bitlen = bitlen;
        } else {
            self.bitlen = bitlen + (bitlen == 0) as i32;
            self.normalise();
        }
    }

    /// Constant-time [`Ibz::set_bound`]: the decrease branch masks and
    /// sign-extends the top limb without scanning the value.
    pub fn set_bound_ct(&mut self, bitlen: i32) {
        debug_assert!(self.bitlen > 0);
        debug_assert!(bitlen >= 0 && bitlen <= Self::MAX_BITS);
        let bitlen = bitlen.clamp(0, Self::MAX_BITS);
        if bitlen >= self.bitlen {
            let n0 = self.n();
            let n1 = nlimbs(bitlen);
            let bit = self.limbs[n0 - 1] >> (LIMB_BITS - 1);
            for i in n0..n1 {
                self.limbs[i] = 0u64.wrapping_sub(bit);
            }
            self.bitlen = bitlen;
        } else {
            self.bitlen = bitlen + (bitlen == 0) as i32;
            let n = self.n();
            let b = self.bitlen - LIMB_BITS * (n as i32 - 1);
            let mask_p = if b >= LIMB_BITS {
                u64::MAX
            } else {
                (1u64 << b) - 1
            };
            let signmask = 0u64.wrapping_sub((self.limbs[n - 1] >> (b - 1)) & 1);
            self.limbs[n - 1] = (self.limbs[n - 1] & mask_p) | (signmask & !mask_p);
        }
    }

    /// `self` with the bound set (non constant time).
    #[inline]
    pub fn with_bound(mut self, bitlen: i32) -> Self {
        self.set_bound(bitlen);
        self
    }

    /// `self` with the bound set in constant time.
    #[inline]
    pub fn with_bound_ct(mut self, bitlen: i32) -> Self {
        self.set_bound_ct(bitlen);
        self
    }

    /// Whether the value is non-negative (zero counts as positive).
    #[inline]
    pub fn is_positive(&self) -> bool {
        (self.limbs[self.n() - 1] >> (LIMB_BITS - 1)) == 0
    }

    /// All-ones mask when negative, zero otherwise.
    #[inline]
    pub(crate) fn neg_mask(&self) -> u64 {
        0u64.wrapping_sub(self.limbs[self.n() - 1] >> (LIMB_BITS - 1))
    }

    /// Whether the value is zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.limbs[..self.n()].iter().all(|&l| l == 0)
    }

    /// Whether the value is one.
    #[inline]
    pub fn is_one(&self) -> bool {
        self.limbs[0] == 1 && self.limbs[1..self.n()].iter().all(|&l| l == 0)
    }

    /// Whether the value is even.
    #[inline]
    pub fn is_even(&self) -> bool {
        self.limbs[0] & 1 == 0
    }

    /// Whether the value is odd.
    #[inline]
    pub fn is_odd(&self) -> bool {
        self.limbs[0] & 1 == 1
    }

    /// Three-way comparison of values.
    pub fn cmp_value(a: &Self, b: &Self) -> Ordering {
        let la = a.n();
        let lb = b.n();
        let a_hi = a.neg_mask();
        let b_hi = b.neg_mask();
        if a_hi == 0 && b_hi != 0 {
            return Ordering::Greater;
        }
        if a_hi != 0 && b_hi == 0 {
            return Ordering::Less;
        }
        for i in (0..la.max(lb)).rev() {
            let al = if i < la { a.limbs[i] } else { a_hi };
            let bl = if i < lb { b.limbs[i] } else { b_hi };
            if al != bl {
                return al.cmp(&bl);
            }
        }
        Ordering::Equal
    }

    /// Compare with a small signed integer.
    #[inline]
    pub fn cmp_i32(&self, y: i32) -> Ordering {
        debug_assert!(y != i32::MIN);
        Self::cmp_value(self, &Self::set(y as i64, 32))
    }

    /// Bit size `ceil(log2(|a| + 1))`, zero for zero.
    pub fn bitsize(&self) -> i32 {
        let n = self.n();
        let sign = self.neg_mask();
        let mut i = n as i32 - 1;
        while i >= 0 && self.limbs[i as usize] == sign {
            i -= 1;
        }
        if i < 0 {
            return (sign as i32).wrapping_neg();
        }
        let l = i * LIMB_BITS;
        let mut d = self.limbs[i as usize] ^ sign;
        if sign != 0 {
            let mut carry = 1u64;
            for j in (0..i).rev() {
                if self.limbs[j as usize] != 0 {
                    carry = 0;
                    break;
                }
            }
            d = d.wrapping_add(carry);
        }
        if d == 0 {
            return l + LIMB_BITS + 1;
        }
        l + LIMB_BITS - d.leading_zeros() as i32
    }

    /// Constant-time [`Ibz::bitsize`]: scans every active limb.
    pub fn bitsize_ct(&self) -> i32 {
        let nwords = self.n();
        let abs_a = self.cneg(self.neg_mask());
        let mut top_word = 0u64;
        let mut top_idx = 0u64;
        let mut found = 0u32;
        for i in (0..nwords).rev() {
            let w = abs_a.limbs[i];
            let nz = ct::is_nonzero(w);
            let take = ct::barrier(0u64.wrapping_sub((nz & (1 - found)) as u64));
            top_word = (top_word & !take) | (w & take);
            top_idx = (top_idx & !take) | ((i as u64) & take);
            found |= nz;
        }
        (top_idx as i32) * LIMB_BITS + ct::limb_bit_length(top_word) as i32
    }

    /// Size in base `base`, `ceil(log_base(|a| + 1))`.
    pub fn size_in_base(&self, base: u32) -> i32 {
        if base < 2 {
            return 0;
        }
        let n = self.bitsize();
        if n == 0 {
            return 0;
        }
        let k = 32 - base.leading_zeros() as i32 - 1;
        if base & (base - 1) == 0 {
            (n + k - 1) / k
        } else {
            (n + k - 1) / k + 1
        }
    }

    /// Two-adic valuation: position of the lowest set bit (zero for zero).
    /// Variable time (stops at the first non-zero limb).
    pub fn two_adic(&self) -> i32 {
        for (i, &w) in self.limbs[..self.n()].iter().enumerate() {
            if w != 0 {
                return (i as i32) * LIMB_BITS + w.trailing_zeros() as i32;
            }
        }
        0
    }

    /// Bit `i` of the two's complement representation.
    #[inline]
    pub(crate) fn bit_at(&self, i: i32) -> u32 {
        ((self.limbs[(i / LIMB_BITS) as usize] >> (i % LIMB_BITS)) & 1) as u32
    }
}

impl<const N: usize> zeroize::Zeroize for Ibz<N> {
    fn zeroize(&mut self) {
        self.limbs.zeroize();
        self.bitlen = 1;
    }
}
