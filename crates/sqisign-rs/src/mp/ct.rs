//! Constant-time fixed-point and wide-arithmetic layer: branch-free word
//! helpers, secret-amount shifts, high-product truncation, mask comparisons
//! and the Newton reciprocal the lattice reduction relies on. Constant time
//! in the limb values; bit lengths, shift bounds and loop counts are public.

use super::arith::{addc, subb};
use super::{nlimbs, Ibz, IBZ_MAX_LIMBS, LIMB_BITS};

/// Optimisation barrier: stops the compiler from turning a masked select
/// back into a branch.
#[inline(always)]
pub fn barrier(x: u64) -> u64 {
    core::hint::black_box(x)
}

/// Barrier for 32-bit values.
#[inline(always)]
pub fn barrier_u32(x: u32) -> u32 {
    core::hint::black_box(x)
}

/// Barrier for signed 32-bit values.
#[inline(always)]
pub fn barrier_i32(x: i32) -> i32 {
    core::hint::black_box(x)
}

/// 1 when `x != 0`, else 0.
#[inline(always)]
pub fn is_nonzero(x: u64) -> u32 {
    ((x | x.wrapping_neg()) >> 63) as u32
}

/// 1 when `x == 0`, else 0.
#[inline(always)]
pub fn is_zero(x: u64) -> u32 {
    1 ^ is_nonzero(x)
}

/// 1 when `x < y`, else 0.
#[inline(always)]
pub fn is_lessthan(x: u64, y: u64) -> u32 {
    ((x ^ ((x ^ y) | (x.wrapping_sub(y) ^ y))) >> 63) as u32
}

/// Bit length of a limb, zero for zero, without a value-dependent branch.
#[inline(always)]
pub fn limb_bit_length(x: u64) -> usize {
    let safe_x = x | 1;
    let len = 64 - safe_x.leading_zeros() as usize;
    let nonzero_mask = barrier(0u64.wrapping_sub((x != 0) as u64));
    len & nonzero_mask as usize
}

/// All-ones when `condition != 0`.
#[inline(always)]
pub fn mask(condition: u64) -> u64 {
    barrier(0u64.wrapping_sub((condition | condition.wrapping_neg()) >> 63))
}

/// All-ones when `b`.
#[inline(always)]
pub fn mask_bool(b: bool) -> u64 {
    mask(b as u64)
}

/// `|x|` as unsigned.
#[inline(always)]
pub fn abs64(x: i64) -> u64 {
    let ux = x as u64;
    let m = (x >> 63) as u64;
    (ux ^ m).wrapping_sub(m)
}

/// `x >> s`, zero for `s >= 64`.
#[inline(always)]
pub fn shr_sat64(x: u64, s: u32) -> u64 {
    let in_range = mask((s < 64) as u64);
    (x >> (s & 63)) & in_range
}

/// `max(e, 0)` branch-free.
#[inline(always)]
pub fn max0(e: i32) -> i32 {
    let sign_bit = (e as u32) >> 31;
    let mask_neg = 0u32.wrapping_sub(sign_bit);
    ((e as u32) & !mask_neg) as i32
}

/// Negate `x` when `mask` is all ones.
#[inline(always)]
pub fn negate_mask_i32(m: u32, x: i32) -> i32 {
    let m = barrier_u32(m);
    ((x as u32 ^ m).wrapping_sub(m)) as i32
}

/// Swap when `mask` is all ones.
#[inline(always)]
pub fn swap_u64(m: u64, a: &mut u64, b: &mut u64) {
    let t = (*a ^ *b) & barrier(m);
    *a ^= t;
    *b ^= t;
}

/// Swap when `mask` is all ones.
#[inline(always)]
pub fn swap_i32(m: u32, a: &mut i32, b: &mut i32) {
    let t = (*a ^ *b) & barrier_i32(m as i32);
    *a ^= t;
    *b ^= t;
}

/// Swap when `mask` is all ones.
#[inline(always)]
pub fn swap_i64(m: u64, a: &mut i64, b: &mut i64) {
    let t = (*a ^ *b) & barrier(m) as i64;
    *a ^= t;
    *b ^= t;
}

/// `condition ? a : b` for 32-bit values, `condition` a 0/1 flag.
#[inline(always)]
pub fn select32(condition: bool, a: i32, b: i32) -> i32 {
    let m = barrier_i32(-(condition as i32));
    (a & m) | (b & !m)
}

/// `condition ? a : b` for unsigned 32-bit values.
#[inline(always)]
pub fn select_u32(condition: bool, a: u32, b: u32) -> u32 {
    let m = barrier_u32(0u32.wrapping_sub(condition as u32));
    (a & m) | (b & !m)
}

/// `condition ? a : b` for unsigned 64-bit values.
#[inline(always)]
pub fn select_u64(condition: bool, a: u64, b: u64) -> u64 {
    let m = barrier(0u64.wrapping_sub(condition as u64));
    (a & m) | (b & !m)
}

/// Branch-free maximum.
#[inline(always)]
pub fn max(a: i32, b: i32) -> i32 {
    select32(a > b, a, b)
}

/// Branch-free minimum.
#[inline(always)]
pub fn min(a: i32, b: i32) -> i32 {
    select32(a < b, a, b)
}

/// Branch-free maximum, unsigned.
#[inline(always)]
pub fn max_u32(a: u32, b: u32) -> u32 {
    select_u32(a > b, a, b)
}

/// Branch-free minimum, unsigned.
#[inline(always)]
pub fn min_u32(a: u32, b: u32) -> u32 {
    select_u32(a < b, a, b)
}

/// Branch-free maximum, 64-bit.
#[inline(always)]
pub fn max_u64(a: u64, b: u64) -> u64 {
    select_u64(a > b, a, b)
}

/// `(|scalar|, sign)` with sign 1 when negative.
#[inline(always)]
pub fn abs_sign_i64(scalar: i64) -> (u64, u64) {
    let u = scalar as u64;
    let m = (scalar >> 63) as u64;
    ((u ^ m).wrapping_sub(m), m & 1)
}

/// `floor(sqrt(n))` for a 32-bit `n`, by a fixed 16-step masked search.
pub fn isqrt_32(n: u32) -> u64 {
    let mut res = 0u32;
    for bit in (0..16).rev() {
        let temp = res | (1u32 << bit);
        let m = barrier_u32(0u32.wrapping_sub((temp.wrapping_mul(temp) <= n) as u32));
        res = (m & temp) | (!m & res);
    }
    res as u64
}

/// High 64 bits of the signed 128-bit product.
#[inline(always)]
pub fn mul_high(a: i64, b: i64) -> i64 {
    (((a as i128) * (b as i128)) >> 64) as i64
}

/// High 64 bits of `a * a`.
#[inline(always)]
pub fn sqr_high(a: i64) -> u64 {
    (((a as i128) * (a as i128)) >> 64) as i64 as u64
}

/// `floor(x / y)` for `y != 0`, restoring division.
pub fn div_unsigned_nonzero(x: u64, y: u64) -> u64 {
    let mut q = 0u64;
    let mut r = 0u64;
    for i in (0..64).rev() {
        r = (r << 1) | ((x >> i) & 1);
        let diff = r.wrapping_sub(y);
        let borrow = ((!r & y) | ((!r | y) & diff)) >> 63;
        let m = 0u64.wrapping_sub(1 - borrow);
        r = (diff & m) | (r & !m);
        q = (q << 1) | (1 & m);
    }
    q
}

/// `floor((xhi:xlo) / y)` for `xhi < y`, `y != 0`.
fn div128_by64(xhi: u64, xlo: u64, y: u64) -> u64 {
    debug_assert!(y != 0 && xhi < y);
    let mut rlo = xhi;
    let mut q = 0u64;
    for i in (0..64).rev() {
        let rhi = rlo >> 63;
        rlo = (rlo << 1) | ((xlo >> i) & 1);
        let ge = barrier(rhi | (1 - is_lessthan(rlo, y)) as u64);
        let m = 0u64.wrapping_sub(ge);
        rlo = rlo.wrapping_sub(y & m);
        q = (q << 1) | ge;
    }
    q
}

/// Number of word-ladder stages needed for `|shift| <= max_shift`.
fn shift_word_stages(max_shift: i32) -> u32 {
    let max_ws = max_shift >> 6;
    let mut nw = 0;
    while (max_ws >> nw) != 0 {
        nw += 1;
    }
    nw
}

/// `t <<= shift` in place, `shift >= 0` secret, `max_shift` public.
fn shl_ct(t: &mut [u64], shift: i32, max_shift: i32) {
    let nwords = t.len();
    let nw = shift_word_stages(max_shift);
    let ws = shift >> 6;
    let bs = (shift & 63) as u32;
    for b in 0..nw {
        let off = 1usize << b;
        let m = 0u64.wrapping_sub(((ws >> b) & 1) as u64);
        for i in (0..nwords).rev() {
            let cand = if i >= off { t[i - off] } else { 0 };
            t[i] ^= (t[i] ^ cand) & m;
        }
    }
    for i in (0..nwords).rev() {
        let lo = t[i];
        let hi = if i >= 1 { t[i - 1] } else { 0 };
        t[i] = (lo << bs) | ((hi >> 1) >> (63 - bs));
    }
}

/// `t = shift >= 0 ? t << shift : t >> -shift` (arithmetic), shift secret
/// including its sign.
fn shift_ct(t: &mut [u64], shift: i32, max_shift: i32) {
    let nwords = t.len();
    let nw = shift_word_stages(max_shift);
    let sgn = shift >> 31;
    let mag = (shift ^ sgn) - sgn;
    let dir = sgn as i64 as u64;
    let ws = mag >> 6;
    let bs = (mag & 63) as u32;
    let mut cur = [0u64; 2 * IBZ_MAX_LIMBS + 1];
    let mut nxt = [0u64; 2 * IBZ_MAX_LIMBS + 1];
    cur[..nwords].copy_from_slice(t);
    let sign = 0u64.wrapping_sub(t[nwords - 1] >> 63);
    for b in 0..nw {
        let off = 1usize << b;
        let m = 0u64.wrapping_sub(((ws >> b) & 1) as u64);
        for i in 0..nwords {
            let cl = if i >= off { cur[i - off] } else { 0 };
            let cr = if i + off < nwords { cur[i + off] } else { sign };
            let cand = cl ^ ((cl ^ cr) & dir);
            nxt[i] = cur[i] ^ ((cur[i] ^ cand) & m);
        }
        core::mem::swap(&mut cur, &mut nxt);
    }
    for i in 0..nwords {
        let lo = cur[i];
        let hl = if i >= 1 { cur[i - 1] } else { 0 };
        let hr = if i + 1 < nwords { cur[i + 1] } else { sign };
        let vl = (lo << bs) | ((hl >> 1) >> (63 - bs));
        let vr = (lo >> bs) | ((hr << 1) << (63 - bs));
        nxt[i] = vl ^ ((vl ^ vr) & dir);
    }
    t.copy_from_slice(&nxt[..nwords]);
}

/// Arithmetic right shift of a raw two's complement buffer by a public `k`.
fn ashr_bits(out: &mut [u64], x: &[u64], nwords: usize, k: i32) {
    let sign = 0u64.wrapping_sub(x[nwords - 1] >> 63);
    let word_shift = (k / 64) as usize;
    let bit_shift = (k % 64) as u32;
    for (i, o) in out[..nwords].iter_mut().enumerate() {
        let src = i + word_shift;
        let lo = if src < nwords { x[src] } else { sign };
        let hi = if src + 1 < nwords { x[src + 1] } else { sign };
        *o = if bit_shift == 0 {
            lo
        } else {
            (lo >> bit_shift) | (hi << (64 - bit_shift))
        };
    }
}

/// Unsigned schoolbook product into a wide buffer.
fn mul_wide(out: &mut [u64], a: &[u64], na: usize, b: &[u64], nb: usize) {
    for i in 0..na {
        let mut carry = 0u128;
        for j in 0..nb {
            let sum = out[i + j] as u128 + (a[i] as u128) * (b[j] as u128) + carry;
            out[i + j] = sum as u64;
            carry = sum >> 64;
        }
        out[i + nb] = carry as u64;
    }
}

impl<const N: usize> Ibz<N> {
    /// Signed wide product `a * b` into a raw buffer of `na + nb + 1` limbs
    /// (bypassing the container bound); returns the number of limbs.
    fn wide_signed_mul(out: &mut [u64; 2 * IBZ_MAX_LIMBS + 1], a: &Self, b: &Self) -> usize {
        let sign_a = a.neg_mask();
        let sign_b = b.neg_mask();
        let sign_prod = (sign_a ^ sign_b) & a.ct_nonzero_mask() & b.ct_nonzero_mask();
        let abs_a = a.cneg(sign_a);
        let abs_b = b.cneg(sign_b);
        let na = a.n();
        let nb = b.n();
        let nwords = na + nb + 1;
        let mut prodbuf = [0u64; 2 * IBZ_MAX_LIMBS + 1];
        mul_wide(&mut prodbuf, &abs_a.limbs, na, &abs_b.limbs, nb);
        // two's complement negation, selected by the sign
        let mut negbuf = [0u64; 2 * IBZ_MAX_LIMBS + 1];
        let mut carry = 1u64;
        for i in 0..nwords {
            let (s, c) = addc(!prodbuf[i], 0, carry);
            negbuf[i] = s;
            carry = c;
        }
        for i in 0..nwords {
            out[i] = ((prodbuf[i] ^ negbuf[i]) & sign_prod) ^ prodbuf[i];
        }
        nwords
    }

    /// Narrow a raw buffer into an integer bound to `out_bitlen`. The bound
    /// is first set limb-aligned so that the constant-time decrease branch
    /// masks and sign-extends the top limb.
    fn wide_to_ibz_bound(buf: &[u64], out_bitlen: i32) -> Self {
        let nout = nlimbs(out_bitlen);
        let mut out = Self::zero();
        out.limbs[..nout].copy_from_slice(&buf[..nout]);
        out.bitlen = (nout as i32) * LIMB_BITS;
        out.set_bound_ct(out_bitlen);
        out
    }

    /// Masked shift by a secret amount: left for `shift >= 0`, arithmetic
    /// right for `shift < 0`; `|shift| <= max_shift`, truncated to
    /// `out_bitlen`.
    pub fn ct_shift(&self, shift: i32, max_shift: i32, out_bitlen: i32) -> Self {
        debug_assert!(max_shift >= 0);
        let work_bitlen = self.bitlen.max(out_bitlen);
        let nwords = nlimbs(work_bitlen);
        debug_assert!(nwords <= N);
        let mut t = *self;
        t.set_bound_ct(work_bitlen);
        shift_ct(&mut t.limbs[..nwords], shift, max_shift);
        t.bitlen = (nlimbs(out_bitlen) as i32) * LIMB_BITS;
        t.set_bound_ct(out_bitlen);
        t
    }

    /// Left shift by a secret non-negative amount, truncated to
    /// `out_bitlen`.
    pub fn ct_shl(&self, shift: i32, max_shift: i32, out_bitlen: i32) -> Self {
        debug_assert!(max_shift >= 0);
        let work_bitlen = self.bitlen.max(out_bitlen);
        let nwords = nlimbs(work_bitlen);
        debug_assert!(nwords <= N);
        let mut res = *self;
        res.set_bound_ct(work_bitlen);
        shl_ct(&mut res.limbs[..nwords], shift, max_shift);
        res.bitlen = (nlimbs(out_bitlen) as i32) * LIMB_BITS;
        res.set_bound_ct(out_bitlen);
        res
    }

    /// `floor(a * b / 2^p)` for a public `p`, bound to `out_bitlen`.
    pub fn ct_highmul_p(&self, b: &Self, p: i32, out_bitlen: i32) -> Self {
        let mut signed_buf = [0u64; 2 * IBZ_MAX_LIMBS + 1];
        let nwords = Self::wide_signed_mul(&mut signed_buf, self, b);
        let mut shifted = [0u64; 2 * IBZ_MAX_LIMBS + 1];
        ashr_bits(&mut shifted, &signed_buf, nwords, p);
        debug_assert!(nlimbs(out_bitlen) <= nwords && nlimbs(out_bitlen) <= N);
        Self::wide_to_ibz_bound(&shifted, out_bitlen)
    }

    /// `trunc(a * b * 2^shift)` for a secret `shift` bounded by `max_shift`.
    pub fn ct_highmul_s(&self, b: &Self, shift: i32, max_shift: i32, out_bitlen: i32) -> Self {
        debug_assert!(max_shift >= 0);
        let mut t = [0u64; 2 * IBZ_MAX_LIMBS + 1];
        let nwords = Self::wide_signed_mul(&mut t, self, b);
        shift_ct(&mut t[..nwords], shift, max_shift);
        debug_assert!(nlimbs(out_bitlen) <= nwords && nlimbs(out_bitlen) <= N);
        Self::wide_to_ibz_bound(&t, out_bitlen)
    }

    /// All-ones when `a < b`.
    pub fn ct_lt_mask(&self, b: &Self) -> u64 {
        let la = self.n();
        let lb = b.n();
        let n = la.max(lb);
        let sa = self.neg_mask();
        let sb = b.neg_mask();
        let mut borrow = 0u64;
        let mut top = 0u64;
        for i in 0..=n {
            let x = if i < la { self.limbs[i] } else { sa };
            let y = if i < lb { b.limbs[i] } else { sb };
            let (d, c) = subb(x, y, borrow);
            top = d;
            borrow = c;
        }
        0u64.wrapping_sub(top >> 63)
    }

    /// All-ones when `a != 0`.
    pub fn ct_nonzero_mask(&self) -> u64 {
        let mut acc = 0u64;
        for &l in &self.limbs[..self.n()] {
            acc |= l;
        }
        barrier(0u64.wrapping_sub(is_nonzero(acc) as u64))
    }

    /// Newton reciprocal: `(W, s)` with `W ~= 2^(p + s) / d`, `s =
    /// bitsize(d)` (secret), `d > 0`. `nrounds` Newton-Raphson rounds from
    /// a 64-bit seed.
    pub fn ct_fp_recip(d: &Self, p: i32, nrounds: u32) -> (Self, i32) {
        debug_assert!(d.is_positive() && !d.is_zero() && p > 63);
        let s = d.bitsize_ct();
        let shift = p - s;
        let max_shift = p.max(d.get_bound());
        let wbits = p + 4;
        let dn = d.ct_shift(shift, max_shift, wbits);
        let top64 = dn.extract_u64(p - 64);
        let q64 = div128_by64(0x7FFF_FFFF_FFFF_FFFF, 0xFFFF_FFFF_FFFF_FFFF, top64);
        let w = Self::from_digits(&[q64]);
        let mut w = w.ct_shl(p - 63, p - 63, wbits);
        let mut two_scaled = Self::set(1, 2).mul_2exp((p + 1) as u32);
        two_scaled.set_bound_ct(wbits + 2);
        for _ in 0..nrounds {
            let mut t = dn.ct_highmul_p(&w, p, wbits);
            t.set_bound_ct(wbits + 2);
            let u = two_scaled.sub(&t);
            w = w.ct_highmul_p(&u, p, wbits);
        }
        w.set_bound_ct(wbits);
        (w, s)
    }

    /// `round(x / 2^p)`, ties towards `+inf`, `p` a non-negative multiple of
    /// 64, `out_bitlen` a multiple of 64.
    pub fn ct_round_shift_limb(&self, p: i32, out_bitlen: i32) -> Self {
        debug_assert!(p >= 0 && p % 64 == 0 && out_bitlen > 0 && out_bitlen % 64 == 0);
        let p_limb = (p / 64) as usize;
        let nq = (out_bitlen / 64) as usize;
        debug_assert!(p_limb + nq <= N);
        let mut xw = *self;
        xw.set_bound_ct(((p_limb + nq) as i32) * LIMB_BITS);
        let round_bit = if p_limb > 0 {
            (xw.limbs[p_limb - 1] >> 63) & 1
        } else {
            0
        };
        let mut res = Self::zero();
        let mut carry = round_bit;
        for i in 0..nq {
            let (s, c) = addc(xw.limbs[p_limb + i], 0, carry);
            res.limbs[i] = s;
            carry = c;
        }
        res.bitlen = out_bitlen;
        res
    }
}
