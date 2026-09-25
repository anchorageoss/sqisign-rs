//! Hexadecimal and decimal conversion (tests and vector parsing).

use super::{Ibz, IBZ_MAX_LIMBS, LIMB_BITS};

impl<const N: usize> Ibz<N> {
    /// Lowercase hexadecimal with a leading `-` for negatives, into `buf`;
    /// returns the written prefix. `buf` needs at least
    /// `2 * IBZ_MAX_LIMBS * 16 + 2` bytes.
    pub fn to_hex<'a>(&self, buf: &'a mut [u8]) -> &'a str {
        let negative = !self.is_positive();
        let absval = self.abs();
        let mut pos = 0;
        if negative {
            buf[pos] = b'-';
            pos += 1;
        }
        if absval.is_zero() {
            buf[pos] = b'0';
            pos += 1;
        } else {
            let bits = absval.bitsize();
            let nnibbles = ((bits + 3) / 4) as usize;
            for k in (0..nnibbles).rev() {
                let bitpos = k * 4;
                let limb = bitpos / 64;
                let off = bitpos % 64;
                let nibble = ((absval.limbs[limb] >> off) & 0xF) as usize;
                buf[pos] = b"0123456789abcdef"[nibble];
                pos += 1;
            }
        }
        core::str::from_utf8(&buf[..pos]).expect("ascii")
    }

    /// Parse a decimal or hexadecimal string (optional sign, no prefix).
    /// `None` on a malformed string or overflow of the container.
    pub fn from_str_radix(s: &str, base: u32) -> Option<Self> {
        if base != 10 && base != 16 {
            return None;
        }
        let bytes = s.as_bytes();
        let mut idx = 0;
        let mut negative = false;
        if idx < bytes.len() && bytes[idx] == b'-' {
            negative = true;
            idx += 1;
        } else if idx < bytes.len() && bytes[idx] == b'+' {
            idx += 1;
        }
        if idx >= bytes.len() {
            return None;
        }
        let base_ibz = Self::set(base as i64, if base == 16 { 6 } else { 5 });
        let mut v = Self::zero();
        for &c in &bytes[idx..] {
            let d = (c as char).to_digit(base)? as i64;
            let digit = Self::set(d, if base == 16 { 6 } else { 5 });
            // the bounds saturate at the container; the value is checked at the end
            v = v.mul(&base_ibz).add(&digit);
        }
        let bits = v.bitsize() + 1;
        if bits > Self::MAX_BITS || (v.bitlen == Self::MAX_BITS && !v.is_positive()) {
            return None;
        }
        v.bitlen = bits;
        if negative {
            v = v.neg();
        }
        Some(v)
    }

    /// Decimal string into `buf` (at least `IBZ_MAX_LIMBS * 20 + 2` bytes).
    pub fn to_decimal<'a>(&self, buf: &'a mut [u8]) -> &'a str {
        let negative = !self.is_positive();
        let mut absval = self.abs();
        let mut pos = 0;
        if negative {
            buf[pos] = b'-';
            pos += 1;
        }
        if absval.is_zero() {
            buf[pos] = b'0';
            pos += 1;
            return core::str::from_utf8(&buf[..pos]).expect("ascii");
        }
        let ten = Self::set(10, 5);
        let mut tmp = [0u8; IBZ_MAX_LIMBS * (LIMB_BITS as usize) + 2];
        let mut len = 0;
        while !absval.is_zero() {
            let (q, r) = absval.div(&ten);
            tmp[len] = b'0' + r.get() as u8;
            len += 1;
            absval = q;
        }
        for k in (0..len).rev() {
            buf[pos] = tmp[k];
            pos += 1;
        }
        core::str::from_utf8(&buf[..pos]).expect("ascii")
    }
}
