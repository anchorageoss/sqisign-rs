//! Small multi-precision helpers on little-endian `u64` limb slices, for the
//! scalar bookkeeping inside the ladders and the pairing discrete logarithm.
//! All slices in one call have the same length.

/// `a < b`. Variable time; used only on public scalars.
#[inline]
pub(crate) fn lt(a: &[u64], b: &[u64]) -> bool {
    for i in (0..a.len()).rev() {
        if a[i] != b[i] {
            return a[i] < b[i];
        }
    }
    false
}

/// `a == 0`.
#[inline]
pub(crate) fn is_zero(a: &[u64]) -> bool {
    a.iter().all(|&w| w == 0)
}

/// `a == 1`.
#[inline]
pub(crate) fn is_one(a: &[u64]) -> bool {
    a[0] == 1 && a[1..].iter().all(|&w| w == 0)
}

/// `a` is odd.
#[inline]
pub(crate) fn is_odd(a: &[u64]) -> bool {
    a[0] & 1 == 1
}

/// `c <- a - b mod 2^(64 n)`.
#[inline]
pub(crate) fn sub(c: &mut [u64], a: &[u64], b: &[u64]) {
    let mut borrow = 0u64;
    for i in 0..c.len() {
        let (d1, b1) = a[i].overflowing_sub(b[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        c[i] = d2;
        borrow = (b1 as u64) + (b2 as u64);
    }
}

/// `c <- a + b mod 2^(64 n)`.
#[inline]
pub(crate) fn add(c: &mut [u64], a: &[u64], b: &[u64]) {
    let mut carry = 0u64;
    for i in 0..c.len() {
        let (s1, c1) = a[i].overflowing_add(b[i]);
        let (s2, c2) = s1.overflowing_add(carry);
        c[i] = s2;
        carry = (c1 as u64) + (c2 as u64);
    }
}

/// `a <- a >> 1`, returning the bit shifted out.
#[inline]
pub(crate) fn shr1(a: &mut [u64]) -> u64 {
    let out = a[0] & 1;
    for i in 0..a.len() {
        let hi = if i + 1 < a.len() { a[i + 1] << 63 } else { 0 };
        a[i] = (a[i] >> 1) | hi;
    }
    out
}

/// `a <- a << shift` for any `shift`, dropping bits beyond the slice.
#[inline]
pub(crate) fn shl(a: &mut [u64], shift: usize) {
    let n = a.len();
    let words = shift / 64;
    let bits = shift % 64;
    if words >= n {
        a.iter_mut().for_each(|w| *w = 0);
        return;
    }
    for i in (0..n).rev() {
        let mut v = 0u64;
        if i >= words {
            v = a[i - words] << bits;
            if bits > 0 && i > words {
                v |= a[i - words - 1] >> (64 - bits);
            }
        }
        a[i] = v;
    }
}

/// Keep only the low `bits` bits of `a`.
#[inline]
pub(crate) fn mask_bits(a: &mut [u64], bits: usize) {
    for (i, w) in a.iter_mut().enumerate() {
        let lo = i * 64;
        if lo >= bits {
            *w = 0;
        } else if lo + 64 > bits {
            *w &= (1u64 << (bits - lo)) - 1;
        }
    }
}

/// Constant-time swap of `a` and `b` when `mask` is all ones.
#[inline]
pub(crate) fn cswap(a: &mut [u64], b: &mut [u64], mask: u64) {
    for i in 0..a.len() {
        let t = mask & (a[i] ^ b[i]);
        a[i] ^= t;
        b[i] ^= t;
    }
}

/// Constant-time select: `c <- if mask { b } else { a }`.
#[inline]
pub(crate) fn select(c: &mut [u64], a: &[u64], b: &[u64], mask: u64) {
    for i in 0..c.len() {
        c[i] = (a[i] & !mask) | (b[i] & mask);
    }
}

/// Bit `i` of `a`.
#[inline]
pub(crate) fn bit(a: &[u64], i: usize) -> u64 {
    (a[i / 64] >> (i % 64)) & 1
}
