//! Basic arithmetic: addition, subtraction, multiplication, negation,
//! shifts, powers and the bounded variants used by the lattice code.

use super::{ct, nlimbs, Ibz, IBZ_MAX_LIMBS, LIMB_BITS};

/// `a + b + carry_in` with carry out, both carries 0 or 1. Branch-free.
#[inline(always)]
pub(crate) fn addc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s1, c1) = a.overflowing_add(b);
    let (s2, c2) = s1.overflowing_add(carry);
    (s2, (c1 as u64) | (c2 as u64))
}

/// `a - b - borrow_in` with borrow out, both borrows 0 or 1. Branch-free.
#[inline(always)]
pub(crate) fn subb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d1, b1) = a.overflowing_sub(b);
    let (d2, b2) = d1.overflowing_sub(borrow);
    (d2, (b1 as u64) | (b2 as u64))
}

impl<const N: usize> Ibz<N> {
    /// `a + b`; the bound grows by one bit.
    pub fn add(&self, b: &Self) -> Self {
        debug_assert!(self.bitlen > 0 && b.bitlen > 0);
        let la = self.n();
        let lb = b.n();
        let a_hi = self.neg_mask();
        let b_hi = b.neg_mask();
        let n = (self.bitlen.max(b.bitlen) + 1).min(Self::MAX_BITS);
        let mut sum = Self::zero();
        sum.bitlen = n;
        let mut carry = 0u64;
        for i in 0..nlimbs(n) {
            let al = if i < la { self.limbs[i] } else { a_hi };
            let bl = if i < lb { b.limbs[i] } else { b_hi };
            let (s, c) = addc(al, bl, carry);
            sum.limbs[i] = s;
            carry = c;
        }
        sum
    }

    /// `a - b`; the bound grows by one bit.
    pub fn sub(&self, b: &Self) -> Self {
        debug_assert!(self.bitlen > 0 && b.bitlen > 0);
        let la = self.n();
        let lb = b.n();
        let a_hi = self.neg_mask();
        let b_hi = b.neg_mask();
        let n = (self.bitlen.max(b.bitlen) + 1).min(Self::MAX_BITS);
        let mut diff = Self::zero();
        diff.bitlen = n;
        let mut borrow = 0u64;
        for i in 0..nlimbs(n) {
            let al = if i < la { self.limbs[i] } else { a_hi };
            let bl = if i < lb { b.limbs[i] } else { b_hi };
            let (d, c) = subb(al, bl, borrow);
            diff.limbs[i] = d;
            borrow = c;
        }
        diff
    }

    /// `a * b`; the bound is the sum of the bounds minus one, capped at the
    /// container. Constant time in the values (the loops run over the
    /// bounds).
    pub fn mul(&self, b: &Self) -> Self {
        debug_assert!(self.bitlen > 0 && b.bitlen > 0);
        let m1_neg = self.neg_mask();
        let m2_neg = b.neg_mask();
        let m1 = self.cneg(m1_neg);
        let m2 = b.cneg(m2_neg);
        let mut res = Self::zero();
        res.bitlen = (self.bitlen + b.bitlen - 1).min(Self::MAX_BITS);
        let n1 = m1.n();
        let n2 = m2.n();
        for i in 0..n1 {
            let mut carry = 0u64;
            let n3 = n2.min(N - i);
            for j in 0..n3 {
                let prod = (m1.limbs[i] as u128) * (m2.limbs[j] as u128);
                let sum = prod + res.limbs[i + j] as u128 + carry as u128;
                res.limbs[i + j] = sum as u64;
                carry = (sum >> LIMB_BITS) as u64;
            }
            if i + n2 < N {
                res.limbs[i + n2] = carry;
            }
        }
        res.cneg(m1_neg ^ m2_neg)
    }

    /// `-a`, same bound.
    pub fn neg(&self) -> Self {
        self.cneg(u64::MAX)
    }

    /// `-a` when `mask` is all ones, `a` when it is zero. Same bound.
    pub fn cneg(&self, mask: u64) -> Self {
        debug_assert!(self.bitlen > 0);
        let n = self.n();
        let mut out = Self::zero();
        let mut carry = (mask & 1) as u128;
        for i in 0..n {
            let t = ((self.limbs[i] ^ mask) as u128) + carry;
            out.limbs[i] = t as u64;
            carry = t >> LIMB_BITS;
        }
        out.bitlen = self.bitlen;
        out
    }

    /// `|a|`.
    #[inline]
    pub fn abs(&self) -> Self {
        self.cneg(self.neg_mask())
    }

    /// `floor(a / 2^exp)`, rounding towards negative infinity.
    pub fn div_2exp(&self, exp: u32) -> Self {
        debug_assert!(self.bitlen > 0);
        let mut q = *self;
        if exp == 0 {
            return q;
        }
        debug_assert!(exp as i32 <= Self::MAX_BITS);
        let iexp = (exp as i32).min(Self::MAX_BITS);
        let nwords = self.n();
        multiple_shiftr(&mut q.limbs[..nwords], iexp as u32);
        q.bitlen = if iexp < self.bitlen {
            self.bitlen - iexp
        } else {
            1
        };
        q
    }

    /// `a * 2^exp`, truncated to the container.
    pub fn mul_2exp(&self, exp: u32) -> Self {
        debug_assert!(self.bitlen > 0);
        let mut product = *self;
        if exp == 0 {
            return product;
        }
        debug_assert!(exp as i32 <= Self::MAX_BITS);
        let iexp = (exp as i32).min(Self::MAX_BITS);
        let bitlen_out = (self.bitlen + iexp).min(Self::MAX_BITS);
        debug_assert!(bitlen_out >= product.bitlen);
        product.set_bound(bitlen_out);
        multiple_shiftl(&mut product.limbs[..nlimbs(bitlen_out)], iexp as u32);
        product
    }

    /// `x^e`, variable time in `e`. `0^0 = 0`.
    pub fn pow(&self, mut e: u32) -> Self {
        let mut r0 = Self::one();
        let mut r1 = *self;
        while e != 0 {
            if e & 1 == 1 {
                r0 = r0.mul(&r1);
            }
            e >>= 1;
            if e != 0 {
                r1 = r1.mul(&r1);
            }
        }
        r0
    }

    /// `a * b` for a small `b`, with the bound set to `result_bound` in
    /// constant time. The caller guarantees `|prod| < 2^(result_bound-1)`.
    pub fn mul_by_int_and_set_bound(&self, b: i32, result_bound: i32) -> Self {
        debug_assert!(self.bitlen > 0);
        let xs = 0u64.wrapping_sub((b < 0) as u64);
        let xmag = ((b as i64 as u64) ^ xs).wrapping_sub(xs);
        let la = self.n();
        let a_s = self.neg_mask();
        let mut prod = self.cneg(a_s);
        let n = nlimbs(result_bound);
        debug_assert!(n <= N);
        let mut carry = 0u128;
        for i in 0..n {
            let ai = if i < la { prod.limbs[i] } else { 0 };
            let t = (ai as u128) * (xmag as u128) + carry;
            prod.limbs[i] = t as u64;
            carry = t >> LIMB_BITS;
        }
        prod.bitlen = (n as i32) * LIMB_BITS;
        let mut out = prod.cneg(a_s ^ xs);
        out.set_bound_ct(result_bound);
        out
    }

    /// `a + b` with the bound set to `result_bound` in constant time.
    #[inline]
    pub fn add_and_set_bound(&self, b: &Self, result_bound: i32) -> Self {
        let mut s = self.add(b);
        s.set_bound_ct(result_bound);
        s
    }

    /// `a + b` for a small `b`, with the bound set to `result_bound`.
    #[inline]
    pub fn add_int_and_set_bound(&self, b: i32, result_bound: i32) -> Self {
        debug_assert!(b != i32::MIN);
        let mut s = self.add(&Self::set(b as i64, 32));
        s.set_bound(result_bound);
        s
    }

    /// `a mod 2^exp`, in `[0, 2^exp)`.
    pub fn mod2exp(&self, exp: u32) -> Self {
        debug_assert!((exp as i32) < Self::MAX_BITS);
        let mut r = *self;
        r.set_bound(exp as i32 + 1);
        let nwords = nlimbs(exp as i32 + 1);
        let q = (exp / 64) as usize;
        let rem = exp % 64;
        if q < nwords {
            r.limbs[q] &= (1u64 << rem) - 1;
            for i in q + 1..nwords {
                r.limbs[i] = 0;
            }
        }
        r
    }

    /// The `[2^offset, 2^(offset+64))` window of `|a|`, independent of
    /// `offset` in timing.
    pub fn extract_u64(&self, offset: i32) -> u64 {
        debug_assert!(self.bitlen > 0);
        let abs_a = self.cneg(self.neg_mask());
        let n = abs_a.n();
        let word_off = (offset >> 6) as u64;
        let bit_off = (offset & 63) as u32;
        let mut w = [0u64; 2];
        for (j, wj) in w.iter_mut().enumerate() {
            let mut acc = 0u64;
            for (idx, &limb) in abs_a.limbs[..n].iter().enumerate() {
                let hit = ct::barrier(
                    0u64.wrapping_sub(ct::is_zero((idx as u64) ^ (word_off + j as u64)) as u64),
                );
                acc |= limb & hit;
            }
            *wj = acc;
        }
        (w[0] >> bit_off) | ((w[1] << 1) << (63 - bit_off))
    }

    /// The 63-bit magnitude window of `|a|` at `offset`, with the sign of
    /// `a`.
    pub fn extract_i64(&self, offset: i32) -> i64 {
        let u = self.extract_u64(offset) & ((1u64 << 63) - 1);
        let s: i64 = if self.is_positive() { 1 } else { -1 };
        s * (u as i64)
    }
}

/// Right shift of a limb array by `1..=63` bits.
#[inline]
pub(crate) fn shiftr(x: &mut [u64], shift: u32) {
    let n = x.len();
    for i in 0..n - 1 {
        x[i] = (x[i] >> shift) ^ (x[i + 1] << (64 - shift));
    }
    x[n - 1] >>= shift;
}

/// Right shift of a limb array by any number of bits (logical, no sign
/// extension: the caller's bound handles the sign).
pub(crate) fn multiple_shiftr(x: &mut [u64], shift: u32) {
    let mut t = shift;
    while t > 63 {
        shiftr(x, 63);
        t -= 63;
    }
    if t > 0 {
        shiftr(x, t);
    }
}

/// Left shift of a limb array by `1..=63` bits.
#[inline]
pub(crate) fn shiftl(x: &mut [u64], shift: u32) {
    let n = x.len();
    for i in (1..n).rev() {
        x[i] = (x[i] << shift) ^ (x[i - 1] >> (64 - shift));
    }
    x[0] <<= shift;
}

/// Left shift of a limb array by any number of bits.
pub(crate) fn multiple_shiftl(x: &mut [u64], shift: u32) {
    let nwords = x.len();
    let word_shift = (shift / 64) as usize;
    let bit_shift = shift % 64;
    for i in (0..nwords).rev() {
        let lo = if i >= word_shift {
            x[i - word_shift]
        } else {
            0
        };
        let hi = if i > word_shift {
            x[i - word_shift - 1]
        } else {
            0
        };
        x[i] = if bit_shift == 0 {
            lo
        } else {
            (lo << bit_shift) | (hi >> (64 - bit_shift))
        };
    }
}

/// Scratch buffer type for the division and gcd routines.
pub(crate) type Scratch = [u64; 2 * IBZ_MAX_LIMBS + 4];

/// A zeroed scratch buffer.
#[inline]
pub(crate) fn scratch() -> Scratch {
    [0u64; 2 * IBZ_MAX_LIMBS + 4]
}
