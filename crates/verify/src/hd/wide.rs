//! Integers modulo `2^192` on three `u64` limbs: the response scalars, the
//! entries of the Kani matrices and every other value the compact verifier
//! reduces modulo `2^r` (`r = 89`, `134`, `175` at levels I, III, V). All
//! arithmetic wraps modulo `2^192`; since `2^r` divides `2^192`, a wrapped
//! result masked to `r` bits is the result modulo `2^r`. Variable time;
//! everything here is public data.

#![allow(clippy::needless_range_loop)] // limb loops index several arrays at once

use core::ops::{BitAnd, BitOr, BitXor, Shl, Shr};

/// Limbs of a [`W`].
pub const W_LIMBS: usize = 3;
/// Bits of a [`W`].
pub const W_BITS: u32 = 64 * W_LIMBS as u32;

/// An unsigned integer modulo `2^192`, little-endian limbs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Hash)]
pub struct W(pub [u64; W_LIMBS]);

impl W {
    pub const ZERO: W = W([0; W_LIMBS]);
    pub const ONE: W = W([1, 0, 0]);

    #[inline]
    pub const fn from_u64(x: u64) -> W {
        W([x, 0, 0])
    }

    #[inline]
    pub const fn from_u128(x: u128) -> W {
        W([x as u64, (x >> 64) as u64, 0])
    }

    /// The low limbs of a little-endian limb slice (missing limbs are zero).
    #[inline]
    pub fn from_limbs(l: &[u64]) -> W {
        let mut r = W::ZERO;
        for (d, s) in r.0.iter_mut().zip(l.iter()) {
            *d = *s;
        }
        r
    }

    /// `2^r − 1` for `r ≤ 192`.
    #[inline]
    pub fn mask(r: u32) -> W {
        debug_assert!(r <= W_BITS);
        let mut m = W::ZERO;
        for i in 0..W_LIMBS {
            let lo = 64 * i as u32;
            if r >= lo + 64 {
                m.0[i] = u64::MAX;
            } else if r > lo {
                m.0[i] = (1u64 << (r - lo)) - 1;
            }
        }
        m
    }

    /// `2^k` for `k < 192`.
    #[inline]
    pub fn pow2(k: u32) -> W {
        let mut r = W::ZERO;
        r.0[(k / 64) as usize] = 1u64 << (k % 64);
        r
    }

    #[inline]
    pub fn low_u64(&self) -> u64 {
        self.0[0]
    }

    #[inline]
    pub fn low_u128(&self) -> u128 {
        (self.0[0] as u128) | ((self.0[1] as u128) << 64)
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|w| *w == 0)
    }

    #[inline]
    pub fn is_odd(&self) -> bool {
        self.0[0] & 1 == 1
    }

    /// Bit `i` (`false` past the width).
    #[inline]
    pub fn bit(&self, i: usize) -> bool {
        i < W_BITS as usize && (self.0[i / 64] >> (i % 64)) & 1 == 1
    }

    /// The bit length (`0` for zero).
    #[inline]
    pub fn bits(&self) -> u32 {
        for i in (0..W_LIMBS).rev() {
            if self.0[i] != 0 {
                return 64 * i as u32 + 64 - self.0[i].leading_zeros();
            }
        }
        0
    }

    /// The number of trailing zero bits (`192` for zero).
    #[inline]
    pub fn trailing_zeros(&self) -> u32 {
        for i in 0..W_LIMBS {
            if self.0[i] != 0 {
                return 64 * i as u32 + self.0[i].trailing_zeros();
            }
        }
        W_BITS
    }

    #[inline]
    pub fn wrapping_add(self, o: W) -> W {
        let mut r = W::ZERO;
        let mut carry = 0u64;
        for i in 0..W_LIMBS {
            let (s, c1) = self.0[i].overflowing_add(o.0[i]);
            let (s, c2) = s.overflowing_add(carry);
            r.0[i] = s;
            carry = (c1 | c2) as u64;
        }
        r
    }

    #[inline]
    pub fn wrapping_sub(self, o: W) -> W {
        let mut r = W::ZERO;
        let mut borrow = 0u64;
        for i in 0..W_LIMBS {
            let (d, b1) = self.0[i].overflowing_sub(o.0[i]);
            let (d, b2) = d.overflowing_sub(borrow);
            r.0[i] = d;
            borrow = (b1 | b2) as u64;
        }
        r
    }

    #[inline]
    pub fn wrapping_neg(self) -> W {
        W::ZERO.wrapping_sub(self)
    }

    /// The low 192 bits of the product.
    #[inline]
    pub fn wrapping_mul(self, o: W) -> W {
        let mut r = [0u64; W_LIMBS];
        for i in 0..W_LIMBS {
            let mut carry = 0u128;
            for j in 0..W_LIMBS - i {
                let acc = r[i + j] as u128 + (self.0[i] as u128) * (o.0[j] as u128) + carry;
                r[i + j] = acc as u64;
                carry = acc >> 64;
            }
        }
        W(r)
    }

    /// The inverse of an odd `self` modulo `2^192` (Newton from the seed 1).
    #[inline]
    pub fn inv_odd(self) -> W {
        debug_assert!(self.is_odd());
        let two = W::from_u64(2);
        let mut inv = W::ONE;
        // the precision doubles from one bit: eight rounds give 256 bits
        for _ in 0..8 {
            inv = inv.wrapping_mul(two.wrapping_sub(self.wrapping_mul(inv)));
        }
        inv
    }

    /// The little-endian bytes.
    #[inline]
    pub fn to_le_bytes(&self) -> [u8; 8 * W_LIMBS] {
        let mut out = [0u8; 8 * W_LIMBS];
        for (i, limb) in self.0.iter().enumerate() {
            out[8 * i..8 * i + 8].copy_from_slice(&limb.to_le_bytes());
        }
        out
    }

    /// From at most 24 little-endian bytes.
    #[inline]
    pub fn from_le_bytes(b: &[u8]) -> W {
        let mut r = W::ZERO;
        for (i, byte) in b.iter().enumerate().take(8 * W_LIMBS) {
            r.0[i / 8] |= (*byte as u64) << (8 * (i % 8));
        }
        r
    }
}

impl BitAnd for W {
    type Output = W;
    #[inline]
    fn bitand(self, o: W) -> W {
        W(core::array::from_fn(|i| self.0[i] & o.0[i]))
    }
}

impl BitOr for W {
    type Output = W;
    #[inline]
    fn bitor(self, o: W) -> W {
        W(core::array::from_fn(|i| self.0[i] | o.0[i]))
    }
}

impl BitXor for W {
    type Output = W;
    #[inline]
    fn bitxor(self, o: W) -> W {
        W(core::array::from_fn(|i| self.0[i] ^ o.0[i]))
    }
}

/// Hexadecimal, no leading zeros (`0` for zero).
impl core::fmt::LowerHex for W {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let top = (0..W_LIMBS).rev().find(|&i| self.0[i] != 0).unwrap_or(0);
        write!(f, "{:x}", self.0[top])?;
        for i in (0..top).rev() {
            write!(f, "{:016x}", self.0[i])?;
        }
        Ok(())
    }
}

impl Shr<u32> for W {
    type Output = W;
    #[inline]
    fn shr(self, k: u32) -> W {
        if k >= W_BITS {
            return W::ZERO;
        }
        let limbs = (k / 64) as usize;
        let bits = k % 64;
        let mut r = W::ZERO;
        for i in 0..W_LIMBS - limbs {
            let src = i + limbs;
            r.0[i] = self.0[src] >> bits;
            if bits != 0 && src + 1 < W_LIMBS {
                r.0[i] |= self.0[src + 1] << (64 - bits);
            }
        }
        r
    }
}

impl Shl<u32> for W {
    type Output = W;
    #[inline]
    fn shl(self, k: u32) -> W {
        if k >= W_BITS {
            return W::ZERO;
        }
        let limbs = (k / 64) as usize;
        let bits = k % 64;
        let mut r = W::ZERO;
        for i in (limbs..W_LIMBS).rev() {
            let src = i - limbs;
            r.0[i] = self.0[src] << bits;
            if bits != 0 && src > 0 {
                r.0[i] |= self.0[src - 1] >> (64 - bits);
            }
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::W;

    #[test]
    fn agrees_with_u128_below_128_bits() {
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s ^ (s >> 31)
        };
        for _ in 0..2000 {
            let a = ((next() as u128) << 64) | next() as u128;
            let b = ((next() as u128) << 64) | next() as u128;
            let (wa, wb) = (W::from_u128(a), W::from_u128(b));
            assert_eq!(wa.wrapping_add(wb).low_u128(), a.wrapping_add(b));
            assert_eq!(wa.wrapping_sub(wb).low_u128(), a.wrapping_sub(b));
            assert_eq!(wa.wrapping_mul(wb).low_u128(), a.wrapping_mul(b));
            assert_eq!((wa & wb).low_u128(), a & b);
            let k = (next() % 128) as u32;
            assert_eq!((wa >> k).low_u128(), a >> k);
            assert_eq!((wa << k).low_u128(), a << k);
            assert_eq!(wa.trailing_zeros().min(128), a.trailing_zeros());
            assert_eq!(wa.bits(), 128 - a.leading_zeros());
            let odd = W::from_u128(a | 1) | (W::from_u64(next()) << 128);
            assert_eq!(odd.wrapping_mul(odd.inv_odd()), W::ONE);
        }
        assert_eq!(W::mask(89), W::pow2(89).wrapping_sub(W::ONE));
        assert_eq!(W::mask(175), W::pow2(175).wrapping_sub(W::ONE));
        assert_eq!(W::mask(192).wrapping_add(W::ONE), W::ZERO);
        assert_eq!(W::ZERO.trailing_zeros(), 192);
    }
}
