//! Fixed-width unsigned integers on four `u64` limbs for the norm equation
//! `2^e − q = a_1^2 + a_2^2` of the compact format: the response degree `q`
//! (`e ≤ 255` bits), Cornacchia's descent and the square root of `−1`
//! modulo a prime by exponentiation. Variable time; everything here is
//! public data (the signature's `q`).

/// A 256-bit unsigned integer, little-endian limbs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct U4(pub [u64; 4]);

impl U4 {
    pub const ZERO: U4 = U4([0; 4]);
    pub const ONE: U4 = U4([1, 0, 0, 0]);

    pub fn from_u64(x: u64) -> U4 {
        U4([x, 0, 0, 0])
    }

    /// `2^k` for `k < 256`.
    pub fn pow2(k: u32) -> U4 {
        let mut r = U4::ZERO;
        r.0[(k / 64) as usize] = 1u64 << (k % 64);
        r
    }

    pub fn from_le_bytes(b: &[u8]) -> U4 {
        let mut r = U4::ZERO;
        for (i, byte) in b.iter().enumerate().take(32) {
            r.0[i / 8] |= (*byte as u64) << (8 * (i % 8));
        }
        r
    }

    pub fn to_le_bytes(&self) -> [u8; 32] {
        let mut out = [0u8; 32];
        for (i, limb) in self.0.iter().enumerate() {
            out[8 * i..8 * i + 8].copy_from_slice(&limb.to_le_bytes());
        }
        out
    }

    pub fn low_u128(&self) -> u128 {
        (self.0[0] as u128) | ((self.0[1] as u128) << 64)
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|w| *w == 0)
    }

    pub fn is_odd(&self) -> bool {
        self.0[0] & 1 == 1
    }

    pub fn bits(&self) -> u32 {
        for i in (0..4).rev() {
            if self.0[i] != 0 {
                return 64 * i as u32 + 64 - self.0[i].leading_zeros();
            }
        }
        0
    }

    pub fn compare(&self, other: &U4) -> core::cmp::Ordering {
        for i in (0..4).rev() {
            if self.0[i] != other.0[i] {
                return self.0[i].cmp(&other.0[i]);
            }
        }
        core::cmp::Ordering::Equal
    }

    pub fn lt(&self, other: &U4) -> bool {
        self.compare(other) == core::cmp::Ordering::Less
    }

    /// `self + other` modulo `2^256`.
    pub fn wrapping_add(&self, other: &U4) -> U4 {
        let mut r = U4::ZERO;
        let mut carry = 0u128;
        for i in 0..4 {
            let s = self.0[i] as u128 + other.0[i] as u128 + carry;
            r.0[i] = s as u64;
            carry = s >> 64;
        }
        r
    }

    /// `self − other` modulo `2^256`.
    pub fn wrapping_sub(&self, other: &U4) -> U4 {
        let mut r = U4::ZERO;
        let mut borrow = 0u64;
        for i in 0..4 {
            let (d, b1) = self.0[i].overflowing_sub(other.0[i]);
            let (d, b2) = d.overflowing_sub(borrow);
            r.0[i] = d;
            borrow = (b1 | b2) as u64;
        }
        r
    }

    /// The low 256 bits of `self * other`.
    pub fn wrapping_mul(&self, other: &U4) -> U4 {
        let w = self.mul_wide(other);
        U4([w[0], w[1], w[2], w[3]])
    }

    /// The full 512-bit product, little-endian limbs.
    pub fn mul_wide(&self, other: &U4) -> [u64; 8] {
        let mut w = [0u64; 8];
        for i in 0..4 {
            let mut carry = 0u128;
            for j in 0..4 {
                let acc = w[i + j] as u128 + (self.0[i] as u128) * (other.0[j] as u128) + carry;
                w[i + j] = acc as u64;
                carry = acc >> 64;
            }
            w[i + 4] = carry as u64;
        }
        w
    }

    pub fn shl(&self, k: u32) -> U4 {
        let mut r = U4::ZERO;
        let limbs = (k / 64) as usize;
        let bits = k % 64;
        for i in (limbs..4).rev() {
            let src = i - limbs;
            r.0[i] = self.0[src] << bits;
            if bits != 0 && src > 0 {
                r.0[i] |= self.0[src - 1] >> (64 - bits);
            }
        }
        r
    }

    pub fn shr(&self, k: u32) -> U4 {
        let mut r = U4::ZERO;
        let limbs = (k / 64) as usize;
        let bits = k % 64;
        for i in 0..4 - limbs {
            let src = i + limbs;
            r.0[i] = self.0[src] >> bits;
            if bits != 0 && src + 1 < 4 {
                r.0[i] |= self.0[src + 1] << (64 - bits);
            }
        }
        r
    }

    /// `(self / d, self % d)` for `d != 0`, by shift-and-subtract.
    pub fn div_rem(&self, d: &U4) -> (U4, U4) {
        debug_assert!(!d.is_zero());
        let mut q = U4::ZERO;
        let mut r = U4::ZERO;
        for i in (0..self.bits()).rev() {
            r = r.shl(1);
            r.0[0] |= (self.0[(i / 64) as usize] >> (i % 64)) & 1;
            if !r.lt(d) {
                r = r.wrapping_sub(d);
                q.0[(i / 64) as usize] |= 1u64 << (i % 64);
            }
        }
        (q, r)
    }

    /// `self * other mod m`, the product taken at 512 bits.
    pub fn mul_mod(&self, other: &U4, m: &U4) -> U4 {
        let w = self.mul_wide(other);
        // reduce the 512-bit value: fold limb by limb from the top
        let mut r = U4::ZERO;
        for i in (0..512).rev() {
            r = r.shl(1);
            r.0[0] |= (w[(i / 64) as usize] >> (i % 64)) & 1;
            if !r.lt(m) {
                r = r.wrapping_sub(m);
            }
        }
        r
    }

    /// `self^e mod m`.
    pub fn pow_mod(&self, e: &U4, m: &U4) -> U4 {
        let mut result = U4::ONE;
        let base = self.div_rem(m).1;
        for i in (0..e.bits()).rev() {
            result = result.mul_mod(&result, m);
            if (e.0[(i / 64) as usize] >> (i % 64)) & 1 == 1 {
                result = result.mul_mod(&base, m);
            }
        }
        result
    }

    /// `⌊√self⌋` by Newton's method.
    pub fn isqrt(&self) -> U4 {
        if self.is_zero() {
            return U4::ZERO;
        }
        let mut x = U4::ONE.shl(self.bits().div_ceil(2));
        loop {
            let (q, _) = self.div_rem(&x);
            let x_new = x.wrapping_add(&q).shr(1);
            if !x_new.lt(&x) {
                return x;
            }
            x = x_new;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::U4;

    #[test]
    fn arithmetic_small_values() {
        let a = U4::from_u64(1_000_003);
        let b = U4::from_u64(999_983);
        assert_eq!(a.wrapping_mul(&b).low_u128(), 1_000_003u128 * 999_983);
        let (q, r) = a.div_rem(&b);
        assert_eq!(q.low_u128(), 1);
        assert_eq!(r.low_u128(), 1_000_003 - 999_983);
        assert_eq!(U4::from_u64(144).isqrt().low_u128(), 12);
        assert_eq!(U4::from_u64(145).isqrt().low_u128(), 12);
        // 3^5 mod 7 = 5
        assert_eq!(
            U4::from_u64(3)
                .pow_mod(&U4::from_u64(5), &U4::from_u64(7))
                .low_u128(),
            5
        );
        let big = U4::pow2(200).wrapping_sub(&U4::ONE);
        let (q, r) = big.div_rem(&U4::pow2(100));
        assert_eq!(q, U4::pow2(100).wrapping_sub(&U4::ONE));
        assert_eq!(r, U4::pow2(100).wrapping_sub(&U4::ONE));
    }
}
