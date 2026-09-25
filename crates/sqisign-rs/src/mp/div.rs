//! Variable-time Euclidean division: Knuth's schoolbook algorithm and the
//! Burnikel-Ziegler recursion of the reference.

use super::arith::{addc, scratch, subb, Scratch};
use super::{nlimbs, Ibz, LIMB_BITS};

/// Leading zeros of a limb.
#[inline]
fn clz(x: u64) -> u32 {
    x.leading_zeros()
}

/// `r = a << shift` for `shift < 64`; returns the bits shifted off the top.
fn shl_words(r: &mut [u64], a: &[u64], n: usize, shift: u32) -> u64 {
    if shift == 0 {
        r[..n].copy_from_slice(&a[..n]);
        return 0;
    }
    let mut carry = 0u64;
    for i in 0..n {
        let w = a[i];
        r[i] = (w << shift) | carry;
        carry = w >> (64 - shift);
    }
    carry
}

/// `r = a >> shift` for `shift < 64`.
fn shr_words(r: &mut [u64], a: &[u64], n: usize, shift: u32) {
    if shift == 0 {
        r[..n].copy_from_slice(&a[..n]);
        return;
    }
    let mut carry = 0u64;
    for i in (0..n).rev() {
        let w = a[i];
        r[i] = (w >> shift) | carry;
        carry = w << (64 - shift);
    }
}

/// Knuth Algorithm D: `u / v` with `v[vlen-1] != 0`, `ulen >= vlen >= 1`.
/// `q` receives `ulen - vlen + 1` words, `r` receives `vlen` words.
fn divmod_schoolbook(q: &mut [u64], r: &mut [u64], u: &[u64], ulen: usize, v: &[u64], vlen: usize) {
    debug_assert!(vlen >= 1 && v[vlen - 1] != 0 && ulen >= vlen);
    if vlen == 1 {
        let d = v[0] as u128;
        let mut rem = 0u128;
        for i in (0..ulen).rev() {
            let num = (rem << 64) | u[i] as u128;
            q[i] = (num / d) as u64;
            rem = num % d;
        }
        r[0] = rem as u64;
        return;
    }
    let shift = clz(v[vlen - 1]);
    let mut vn = scratch();
    shl_words(&mut vn, v, vlen, shift);
    let mut un = scratch();
    let topcarry = shl_words(&mut un, u, ulen, shift);
    un[ulen] = topcarry;
    let qlen = ulen - vlen + 1;
    let base: u128 = 1u128 << 64;
    for j in (0..qlen).rev() {
        let num = ((un[j + vlen] as u128) << 64) | un[j + vlen - 1] as u128;
        let mut qhat = num / vn[vlen - 1] as u128;
        let mut rhat = num % vn[vlen - 1] as u128;
        if qhat >= base {
            qhat = base - 1;
            rhat = num - qhat * vn[vlen - 1] as u128;
        }
        while rhat < base && qhat * vn[vlen - 2] as u128 > (rhat << 64) + un[j + vlen - 2] as u128 {
            qhat -= 1;
            rhat += vn[vlen - 1] as u128;
        }
        let mut borrow = 0u64;
        let mut carry = 0u128;
        for i in 0..vlen {
            let prod = qhat * vn[i] as u128 + carry;
            carry = prod >> 64;
            let prod_lo = prod as u64;
            let t = un[j + i];
            let (sub1, ba) = t.overflowing_sub(prod_lo);
            let (sub2, bb) = sub1.overflowing_sub(borrow);
            un[j + i] = sub2;
            borrow = ba as u64 + bb as u64;
        }
        {
            let t = un[j + vlen];
            let c_lo = carry as u64;
            let (sub1, ba) = t.overflowing_sub(c_lo);
            let (sub2, bb) = sub1.overflowing_sub(borrow);
            un[j + vlen] = sub2;
            borrow = ba as u64 + bb as u64;
        }
        if borrow != 0 {
            qhat -= 1;
            let mut addc_ = 0u128;
            for i in 0..vlen {
                let s = un[j + i] as u128 + vn[i] as u128 + addc_;
                un[j + i] = s as u64;
                addc_ = s >> 64;
            }
            un[j + vlen] = (un[j + vlen] as u128 + addc_) as u64;
        }
        q[j] = qhat as u64;
    }
    shr_words(r, &un, vlen, shift);
}

/// `r = a + b` on unsigned arrays, returning the final carry.
pub(crate) fn add_raw(
    r: &mut [u64],
    rlen: usize,
    a: &[u64],
    alen: usize,
    b: &[u64],
    blen: usize,
) -> u64 {
    let mut carry = 0u64;
    for i in 0..rlen {
        let av = if i < alen { a[i] } else { 0 };
        let bv = if i < blen { b[i] } else { 0 };
        let (s, c) = addc(av, bv, carry);
        r[i] = s;
        carry = c;
    }
    carry
}

/// `r = a - b` on unsigned arrays, returning 1 when the true difference is
/// negative.
pub(crate) fn sub_raw(
    r: &mut [u64],
    rlen: usize,
    a: &[u64],
    alen: usize,
    b: &[u64],
    blen: usize,
) -> u64 {
    let mut borrow = 0u64;
    for i in 0..rlen {
        let av = if i < alen { a[i] } else { 0 };
        let bv = if i < blen { b[i] } else { 0 };
        let (d, c) = subb(av, bv, borrow);
        r[i] = d;
        borrow = c;
    }
    borrow
}

/// `dst[0..n] += src[0..n] * limb`, returning the carry.
fn addmul_1_raw(dst: &mut [u64], src: &[u64], n: usize, limb: u64) -> u64 {
    let mut carry = 0u128;
    for j in 0..n {
        let t = (limb as u128) * (src[j] as u128) + dst[j] as u128 + carry;
        dst[j] = t as u64;
        carry = t >> 64;
    }
    carry as u64
}

const MUL_SWAP_RATIO_THRESHOLD: usize = 4;

/// Plain schoolbook `prod = a * b`; `prod` gets exactly `alen + blen` words.
pub(crate) fn mul_raw(prod: &mut [u64], a: &[u64], alen: usize, b: &[u64], blen: usize) {
    for p in prod[..alen + blen].iter_mut() {
        *p = 0;
    }
    let longer = alen.max(blen);
    let shorter = alen.min(blen);
    let worth_swapping = longer >= MUL_SWAP_RATIO_THRESHOLD * shorter;
    let (outer, outer_len, inner, inner_len) = if worth_swapping && alen > blen {
        (b, blen, a, alen)
    } else {
        (a, alen, b, blen)
    };
    for i in 0..outer_len {
        let oi = outer[i];
        if oi == 0 {
            continue;
        }
        let mut carry = addmul_1_raw(&mut prod[i..], inner, inner_len, oi);
        let mut k = i + inner_len;
        while carry != 0 {
            let t = prod[k] as u128 + carry as u128;
            prod[k] = t as u64;
            carry = (t >> 64) as u64;
            k += 1;
        }
    }
}

/// `r[off..] += x`.
fn add_inplace_at(r: &mut [u64], x: &[u64], xlen: usize, off: usize) {
    let mut carry = 0u128;
    let mut i = 0;
    while i < xlen {
        let s = r[off + i] as u128 + x[i] as u128 + carry;
        r[off + i] = s as u64;
        carry = s >> 64;
        i += 1;
    }
    while carry != 0 {
        let s = r[off + i] as u128 + carry;
        r[off + i] = s as u64;
        carry = s >> 64;
        i += 1;
    }
}

/// `a -= 1` on an unsigned array.
fn dec1_raw(a: &mut [u64], len: usize) {
    for i in 0..len {
        if a[i] != 0 {
            a[i] -= 1;
            return;
        }
        a[i] = u64::MAX;
    }
}

#[inline]
fn bz_should_recurse(n: usize) -> bool {
    n % 2 == 0 && n >= 2
}

/// Burnikel-Ziegler D3n2n: 3n-word `a` by 2n-word `b`, quotient `n + 1`
/// words, remainder `2n` words.
fn div_d3n2n(qhat: &mut [u64], r: &mut [u64], a: &[u64], b: &[u64], n: usize) {
    let b0 = &b[..n];
    let b1 = &b[n..2 * n];
    let a0 = &a[..n];
    let ah = &a[n..3 * n];
    let mut r1 = scratch();
    div_d2n1n(qhat, &mut r1, ah, b1, n);
    let l = 2 * n + 1;
    let mut rbuf = scratch();
    rbuf[..n].copy_from_slice(a0);
    rbuf[n..2 * n].copy_from_slice(&r1[..n]);
    let mut prod = scratch();
    mul_raw(&mut prod, qhat, n + 1, b0, n);
    let mut tmp = scratch();
    tmp[..l].copy_from_slice(&rbuf[..l]);
    let mut borrow = sub_raw(&mut rbuf, l, &tmp, l, &prod, l);
    while borrow != 0 {
        dec1_raw(qhat, n + 1);
        tmp[..l].copy_from_slice(&rbuf[..l]);
        let addcarry = add_raw(&mut rbuf, l, &tmp, l, b, 2 * n);
        borrow = if addcarry != 0 { 0 } else { 1 };
    }
    r[..2 * n].copy_from_slice(&rbuf[..2 * n]);
}

/// Burnikel-Ziegler D2n1n: 2n-word `a` by n-word `b`, quotient `n + 1`
/// words, remainder `n` words.
fn div_d2n1n(q: &mut [u64], r: &mut [u64], a: &[u64], b: &[u64], n: usize) {
    if !bz_should_recurse(n) {
        divmod_schoolbook(q, r, a, 2 * n, b, n);
        return;
    }
    let n2 = n / 2;
    let a0 = &a[..n2];
    let a1 = &a[n2..3 * n2 + n2];
    let mut q1 = scratch();
    let mut r1 = scratch();
    div_d3n2n(&mut q1, &mut r1, a1, b, n2);
    let mut aa = scratch();
    aa[..n2].copy_from_slice(a0);
    aa[n2..n2 + n].copy_from_slice(&r1[..n]);
    let mut q2 = scratch();
    div_d3n2n(&mut q2, r, &aa, b, n2);
    for x in q[..n + 1].iter_mut() {
        *x = 0;
    }
    q[..n2 + 1].copy_from_slice(&q2[..n2 + 1]);
    add_inplace_at(q, &q1, n2 + 1, n2);
}

/// Exact quotient buffer length for [`div_unsigned`].
pub(crate) fn div_qlen(ulen: usize, vlen: usize) -> usize {
    let n = vlen;
    let un_len = ulen + 1;
    let mut t = un_len / n;
    if un_len % n != 0 {
        t += 1;
    }
    t * n + 1
}

/// Burnikel-Ziegler blockwise division of `u` by `v` (`v` trimmed).
fn div_unsigned(q: &mut [u64], r: &mut [u64], u: &[u64], ulen: usize, v: &[u64], vlen: usize) {
    debug_assert!(vlen >= 1 && v[vlen - 1] != 0);
    let n = vlen;
    let shift = clz(v[n - 1]);
    let mut vn = scratch();
    shl_words(&mut vn, v, n, shift);
    let un_len = ulen + 1;
    let mut un = scratch();
    let topcarry = shl_words(&mut un, u, ulen, shift);
    un[ulen] = topcarry;
    let mut t = un_len / n;
    if un_len % n != 0 {
        t += 1;
    }
    let qbuf_len = t * n + 1;
    for x in q[..qbuf_len].iter_mut() {
        *x = 0;
    }
    let mut rem = scratch();
    for i in (0..t).rev() {
        let mut window = scratch();
        for k in 0..n {
            let j = i * n + k;
            window[k] = if j < un_len { un[j] } else { 0 };
        }
        if i == t - 1 {
            let mut diff = scratch();
            let borrow = sub_raw(&mut diff, n, &window, n, &vn, n);
            if borrow != 0 {
                rem[..n].copy_from_slice(&window[..n]);
            } else {
                rem[..n].copy_from_slice(&diff[..n]);
                q[i * n] = 1;
            }
            continue;
        }
        let rem_copy: Scratch = rem;
        window[n..2 * n].copy_from_slice(&rem_copy[..n]);
        let mut qi = scratch();
        div_d2n1n(&mut qi, &mut rem, &window, &vn, n);
        add_inplace_at(q, &qi, n + 1, i * n);
    }
    shr_words(r, &rem, n, shift);
}

/// Trimmed length of an unsigned array (at least 1).
#[inline]
pub(crate) fn trim_len(a: &[u64], mut len: usize) -> usize {
    while len > 1 && a[len - 1] == 0 {
        len -= 1;
    }
    len
}

/// Unsigned comparison of two equal-length arrays.
pub(crate) fn compare(a: &[u64], b: &[u64], nwords: usize) -> core::cmp::Ordering {
    for i in (0..nwords).rev() {
        if a[i] != b[i] {
            return a[i].cmp(&b[i]);
        }
    }
    core::cmp::Ordering::Equal
}

/// Unsigned `|a| / |b|` into `(q, r)` limb buffers with the given active
/// lengths; returns whether `|a| < |b|` (then `q = 0`, `r = |a|`).
pub(crate) fn divmod_abs(
    q: &mut Scratch,
    r: &mut Scratch,
    abs_a: &[u64],
    ulen: usize,
    abs_b: &[u64],
    vlen: usize,
    want_q: bool,
) {
    let mut a_lt_b = ulen < vlen;
    if !a_lt_b && ulen == vlen {
        a_lt_b = compare(abs_a, abs_b, ulen) == core::cmp::Ordering::Less;
    }
    if a_lt_b {
        r[..ulen].copy_from_slice(&abs_a[..ulen]);
    } else if !bz_should_recurse(vlen) {
        let mut qs = scratch();
        divmod_schoolbook(
            if want_q { q } else { &mut qs },
            r,
            abs_a,
            ulen,
            abs_b,
            vlen,
        );
    } else {
        let mut qbuf = scratch();
        let mut rbuf = scratch();
        div_unsigned(&mut qbuf, &mut rbuf, abs_a, ulen, abs_b, vlen);
        if want_q {
            let qlen = div_qlen(ulen, vlen);
            q[..qlen].copy_from_slice(&qbuf[..qlen]);
        }
        r[..vlen].copy_from_slice(&rbuf[..vlen]);
    }
}

impl<const N: usize> Ibz<N> {
    /// Euclidean division rounded towards zero: `(quotient, remainder)` with
    /// `a = quotient * b + remainder` and `|remainder| < |b|`, the remainder
    /// taking the sign of `a`. Variable time. `b` must be non-zero.
    pub fn div(&self, b: &Self) -> (Self, Self) {
        debug_assert!(self.bitlen > 0 && b.bitlen > 0);
        debug_assert!(!b.is_zero());
        let n = self.bitlen.max(b.bitlen);
        let l = nlimbs(n);
        let a_neg = self.neg_mask();
        let b_neg = b.neg_mask();
        let abs_a = self.cneg(a_neg);
        let abs_b = b.cneg(b_neg);
        let ulen = trim_len(&abs_a.limbs, l);
        let vlen = trim_len(&abs_b.limbs, l);
        let mut q = scratch();
        let mut r = scratch();
        divmod_abs(&mut q, &mut r, &abs_a.limbs, ulen, &abs_b.limbs, vlen, true);
        let mut quotient = Self::zero();
        let mut remainder = Self::zero();
        quotient.limbs[..l].copy_from_slice(&q[..l]);
        remainder.limbs[..l].copy_from_slice(&r[..l]);
        quotient.bitlen = self.bitlen;
        remainder.bitlen = b.bitlen;
        if a_neg != b_neg {
            quotient = quotient.neg();
        }
        if a_neg != 0 {
            remainder = remainder.neg();
        }
        (quotient, remainder)
    }

    /// `a mod b`, always non-negative (the sign of `b` is ignored). Variable
    /// time.
    pub fn modulo(&self, b: &Self) -> Self {
        debug_assert!(!b.is_zero());
        let n = self.bitlen.max(b.bitlen);
        let l = nlimbs(n);
        let a_neg = self.neg_mask();
        let b_neg = b.neg_mask();
        let abs_a = self.cneg(a_neg);
        let abs_b = b.cneg(b_neg);
        let ulen = trim_len(&abs_a.limbs, l);
        let vlen = trim_len(&abs_b.limbs, l);
        let mut q = scratch();
        let mut r = scratch();
        divmod_abs(
            &mut q,
            &mut r,
            &abs_a.limbs,
            ulen,
            &abs_b.limbs,
            vlen,
            false,
        );
        let mut res = Self::zero();
        res.limbs[..l].copy_from_slice(&r[..l]);
        res.bitlen = b.bitlen;
        if a_neg != 0 {
            res = res.neg();
        }
        if !res.is_positive() {
            if b_neg != 0 {
                res = res.sub(b);
            } else {
                res = res.add(b);
            }
        }
        res
    }

    /// `a mod d` for a small `d`.
    pub fn mod_u64(&self, d: u64) -> u64 {
        debug_assert!(d != 0);
        let dz = Self::from_bits(&[d], 64);
        let r = self.modulo(&dz);
        r.limbs[0]
    }

    /// Whether `b` divides `a`.
    #[inline]
    pub fn divides(&self, b: &Self) -> bool {
        self.modulo(b).is_zero()
    }
}

#[allow(dead_code)]
const _: i32 = LIMB_BITS;
