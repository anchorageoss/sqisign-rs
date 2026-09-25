//! Constant-time Lagrange-Gauss reduction of rank-2 lattices over `Z` with
//! a Lehmer window on the most significant bits, and the sum-of-two-squares
//! decomposition of a prime built on it.

use super::{
    dim2_outer_its, dim2_shift_max, dim2_val_bits, dim2_work_bits, DIM2_THRESHOLD, DIM2_W,
};
use crate::mp::{ct, Ibz};
use crate::quat::Mat2x2;

/// Partial extended gcd on 64-bit windows; returns the number of valid
/// steps and the column-convention transform.
fn partial_xgcd_msb(mut a: u64, mut b: u64, sqrta: u64, swap: bool) -> (i32, [[i32; 2]; 2]) {
    let (mut u00, mut u01, mut u10, mut u11) = (1i32, 0i32, 0i32, 1i32);
    let mut valid_mask = u64::MAX;
    let mut valid_steps = 0;
    let global_swap_mask = ct::mask_bool(swap);
    for _ in 0..DIM2_THRESHOLD + 1 {
        valid_mask &= !ct::mask_bool(a == 0 || b == 0);
        let swap_mask = ct::mask_bool(a < b) & valid_mask;
        valid_mask &= !swap_mask | global_swap_mask;
        let actual_swap = swap_mask & global_swap_mask;
        ct::swap_u64(actual_swap, &mut a, &mut b);
        ct::swap_i32(actual_swap as u32, &mut u00, &mut u10);
        ct::swap_i32(actual_swap as u32, &mut u01, &mut u11);
        valid_mask &= !ct::mask_bool(b < sqrta);
        let lza = (a | 1).leading_zeros();
        let lzb = (b | 1).leading_zeros();
        valid_mask &= !ct::mask_bool(lza > DIM2_THRESHOLD as u32 || lzb > DIM2_THRESHOLD as u32);
        let mut k = lzb.wrapping_sub(lza) & (valid_mask as u32);
        let b_shifted = b << (k & 63);
        k = k.wrapping_sub(1 & ct::mask_bool(b_shifted > a) as u32);
        k &= valid_mask as u32;
        a = a.wrapping_sub((b << (k & 63)) & valid_mask);
        u00 = u00.wrapping_sub((((u10 as u32) << (k & 31)) & (valid_mask as u32)) as i32);
        u01 = u01.wrapping_sub((((u11 as u32) << (k & 31)) & (valid_mask as u32)) as i32);
        valid_steps += (valid_mask & 1) as i32;
    }
    (valid_steps, [[u00, u01], [u10, u11]])
}

/// `x / 2^off` as a 64-bit window, for `off` of either sign.
pub(crate) fn lg_window_u64<const N: usize>(x: &Ibz<N>, off: i32) -> u64 {
    x.extract_u64(ct::max0(off)) << (ct::max0(-off) & 63)
}

/// Signed window.
pub(crate) fn lg_window_i64<const N: usize>(x: &Ibz<N>, off: i32) -> i64 {
    ((x.extract_i64(ct::max0(off)) as u64) << (ct::max0(-off) & 63)) as i64
}

/// Apply the 32-bit transform `m` to `(a, b)` with an extra secret left
/// shift on the off-diagonal term.
fn ct_apply_matrix<const N: usize>(
    a: &mut Ibz<N>,
    b: &mut Ibz<N>,
    m: &[[i32; 2]; 2],
    extra_shift: i32,
    det_bits: i32,
) {
    let wbits = dim2_work_bits(det_bits);
    let vbits = dim2_val_bits(det_bits);
    let t1 = a.mul_by_int_and_set_bound(m[1][0], wbits);
    let t2 = b.ct_shl(extra_shift, dim2_shift_max(det_bits), wbits);
    let t2 = t2.mul_by_int_and_set_bound(m[0][1], wbits);
    let na = a
        .mul_by_int_and_set_bound(m[0][0], wbits)
        .add_and_set_bound(&t2, vbits);
    let nb = b
        .mul_by_int_and_set_bound(m[1][1], wbits)
        .add_and_set_bound(&t1, vbits);
    *a = na;
    *b = nb;
}

/// `a += mu * (b << extra_shift)`.
fn ct_size_reduce<const N: usize>(
    a: &mut Ibz<N>,
    b: &Ibz<N>,
    mu: i32,
    extra_shift: i32,
    det_bits: i32,
) {
    let t1 = b.ct_shl(
        extra_shift,
        dim2_shift_max(det_bits),
        dim2_work_bits(det_bits),
    );
    let t1 = t1.mul_by_int_and_set_bound(mu, dim2_work_bits(det_bits));
    *a = a.add_and_set_bound(&t1, dim2_val_bits(det_bits));
}

/// Lehmer partial xgcd on the HNF `[[A/g, 0], [B, g]]` given as
/// `(aa, bb, cc, dd)` with `det = aa dd <= 2^det_bits`.
fn partial_xgcd_lehmer<const N: usize>(
    aa: &mut Ibz<N>,
    bb: &mut Ibz<N>,
    cc: &mut Ibz<N>,
    dd: &mut Ibz<N>,
    outer_iters: i32,
    det_bits: i32,
) {
    let mut det = aa.mul(dd);
    debug_assert!(det.bitsize_ct() <= det_bits);
    det.set_bound_ct(dim2_val_bits(det_bits));
    let ba = det.bitsize_ct();
    let shift_2k = (ba - DIM2_THRESHOLD) & !1;
    let shift_k = shift_2k / 2;
    let a_trunc = lg_window_u64(&det, shift_2k) as u32;
    let sqrta_trunc = ct::isqrt_32(a_trunc);
    for _ in 0..outer_iters {
        let diff = aa.sub(bb);
        let borrow = if diff.is_positive() { 0 } else { u64::MAX };
        Ibz::cswap(aa, bb, borrow);
        Ibz::cswap(cc, dd, borrow);
        let bits_a = aa.bitsize_ct();
        let bits_b = bb.bitsize_ct();
        let rho = bits_a - bits_b + 1;
        let bits_from_a = bits_a - DIM2_W;
        let large_diff = rho >= DIM2_THRESHOLD;
        let bits_from_b = ct::select32(
            large_diff,
            bits_b - (DIM2_W - DIM2_THRESHOLD + 1),
            bits_from_a,
        );
        let extra_shift = ct::select32(large_diff, rho - DIM2_THRESHOLD, 0);
        let pa = lg_window_u64(aa, bits_from_a);
        let pb = lg_window_u64(bb, bits_from_b);
        let diff_sh = bits_from_b - shift_k;
        let shift_right = ct::min(16, ct::max0(diff_sh));
        let mut psqrta = (sqrta_trunc >> (shift_right & 63)) as i64;
        let shift_left = ct::select32(diff_sh <= -48, 0, ct::max0(-diff_sh));
        psqrta = ((psqrta as u64) << (shift_left & 63)) as i64;
        psqrta |= ct::mask_bool(diff_sh <= -48) as i64;
        let (_, m) = partial_xgcd_msb(pa, pb, psqrta as u64, !large_diff);
        ct_apply_matrix(aa, bb, &m, extra_shift, det_bits);
        ct_apply_matrix(cc, dd, &m, extra_shift, det_bits);
        let neg_a = aa.neg_mask_pub();
        *aa = aa.cneg(neg_a);
        *cc = cc.cneg(neg_a);
        let neg_b = bb.neg_mask_pub();
        *bb = bb.cneg(neg_b);
        *dd = dd.cneg(neg_b);
        debug_assert!(aa.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(bb.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(cc.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(dd.bitsize_ct() < dim2_val_bits(det_bits));
    }
}

/// Size-reduction coefficient of `v2 = (a, c)` by `v1 = (b, d)` from
/// 64-bit windows.
fn sizered_msb(a: u64, b: u64, c: i64, d: i64) -> i32 {
    let v1_shift = ct::min(
        ((b | 1).leading_zeros() as i32) - 1,
        (d | 1).leading_zeros_signed(),
    );
    let b_shifted = (b << v1_shift) as i64;
    let d_shifted = ((d as u64) << v1_shift) as i64;
    let num = (ct::mul_high(a as i64, b_shifted) >> 1) + (ct::mul_high(c, d_shifted) >> 1);
    let denom = (ct::mul_high(b as i64, b_shifted) >> 1) + (ct::mul_high(d, d_shifted) >> 1);
    let (num_abs, num_sign) = ct::abs_sign_i64(num);
    let k_abs = ct::div_unsigned_nonzero(
        2u64.wrapping_mul(num_abs).wrapping_add(denom as u64),
        2u64.wrapping_mul(denom as u64),
    );
    let m = 0u64.wrapping_sub(num_sign & 1);
    let signed_u = (k_abs ^ m).wrapping_add(num_sign);
    -(signed_u as i64 as i32)
}

trait LeadingRedundant {
    fn leading_zeros_signed(self) -> i32;
}

impl LeadingRedundant for i64 {
    /// The reference's `__builtin_clrsbll`: redundant sign bits.
    fn leading_zeros_signed(self) -> i32 {
        let x = self;
        let y = if x < 0 { !x } else { x };
        (y.leading_zeros() as i32) - 1
    }
}

/// Post-processing size reductions after the partial xgcd.
fn sizered_lehmer<const N: usize>(
    aa: &mut Ibz<N>,
    bb: &mut Ibz<N>,
    cc: &mut Ibz<N>,
    dd: &mut Ibz<N>,
    outer_iters: i32,
    det_bits: i32,
) {
    for _ in 0..outer_iters {
        let bits_a = aa.bitsize_ct();
        let bits_b = bb.bitsize_ct();
        let bits_c = cc.bitsize_ct();
        let bits_d = dd.bitsize_ct();
        let bits_v1 = ct::max(bits_b, bits_d);
        let bits_v2 = ct::max(bits_a, bits_c);
        let bits_max = ct::max(bits_v1, bits_v2);
        let rho = bits_max - bits_v1 + 1;
        let bits_from_v2 = bits_max - DIM2_W;
        let large_diff = rho >= DIM2_THRESHOLD;
        let bits_from_v1 = ct::select32(
            large_diff,
            bits_v1 - (DIM2_W - DIM2_THRESHOLD + 1),
            bits_from_v2,
        );
        let extra_shift = ct::select32(large_diff, rho - DIM2_THRESHOLD, 0);
        let pa = lg_window_u64(aa, bits_from_v2);
        let pb = lg_window_u64(bb, bits_from_v1);
        let pc = lg_window_i64(cc, bits_from_v2);
        let pd = lg_window_i64(dd, bits_from_v1);
        let mu = sizered_msb(pa, pb, pc, pd);
        let b_copy = *bb;
        let d_copy = *dd;
        ct_size_reduce(aa, &b_copy, mu, extra_shift, det_bits);
        ct_size_reduce(cc, &d_copy, mu, extra_shift, det_bits);
        let neg_a = aa.neg_mask_pub();
        *aa = aa.cneg(neg_a);
        *cc = cc.cneg(neg_a);
        let neg_b = bb.neg_mask_pub();
        *bb = bb.cneg(neg_b);
        *dd = dd.cneg(neg_b);
        debug_assert!(aa.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(bb.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(cc.bitsize_ct() < dim2_val_bits(det_bits));
        debug_assert!(dd.bitsize_ct() < dim2_val_bits(det_bits));
    }
}

/// Constant-time Lagrange-Gauss reduction of the rank-2 basis in HNF
/// `[[A/g, 0], [B, g]]`, `det = A <= 2^det_bits` (public).
pub fn dim2_short_basis<const N: usize>(basis: &Mat2x2<N>, det_bits: i32) -> Mat2x2<N> {
    debug_assert!(dim2_work_bits(det_bits) < Ibz::<N>::MAX_BITS);
    let outer_iters = dim2_outer_its(det_bits);
    let mut r = *basis;
    for i in 0..2 {
        for j in 0..2 {
            r.0[i][j].set_bound_ct(dim2_val_bits(det_bits));
        }
    }
    let [[mut a, mut b], [mut c, mut d]] = r.0;
    partial_xgcd_lehmer(&mut a, &mut b, &mut c, &mut d, outer_iters, det_bits);
    sizered_lehmer(&mut a, &mut b, &mut c, &mut d, outer_iters, det_bits);
    sizered_lehmer(&mut b, &mut a, &mut d, &mut c, 1, det_bits);
    let bits_a = a.bitsize_ct();
    let bits_b = b.bitsize_ct();
    let bits_c = c.bitsize_ct();
    let bits_d = d.bitsize_ct();
    let bits_v1 = ct::max(bits_b, bits_d);
    let bits_v2 = ct::max(bits_a, bits_c);
    let bits_max = ct::max(bits_v1, bits_v2);
    let bits_from_vecs = ct::select32(bits_max >= DIM2_W, bits_max - DIM2_W, 0);
    let up = ct::max0(DIM2_W - 1 - bits_max) & 63;
    let pa = ((a.extract_i64(bits_from_vecs) as u64) << up) as i64;
    let pb = ((b.extract_i64(bits_from_vecs) as u64) << up) as i64;
    let pc = ((c.extract_i64(bits_from_vecs) as u64) << up) as i64;
    let pd = ((d.extract_i64(bits_from_vecs) as u64) << up) as i64;
    let norm_v1 = ct::sqr_high(pa).wrapping_add(ct::sqr_high(pc));
    let norm_v2 = ct::sqr_high(pb).wrapping_add(ct::sqr_high(pd));
    let swap_mask = ct::mask_bool(norm_v1 > norm_v2);
    Ibz::cswap(&mut a, &mut b, swap_mask);
    Ibz::cswap(&mut c, &mut d, swap_mask);
    Mat2x2([[a, b], [c, d]])
}

/// `(x, y)` with `x^2 + y^2 = p` for a prime `p` and `r^2 = -1 mod p`, in
/// constant time.
pub fn dim2_sumofsquares<const N: usize>(p: &Ibz<N>, r: &Ibz<N>) -> (Ibz<N>, Ibz<N>) {
    let det_bits = p.get_bound().max(r.get_bound());
    debug_assert!(dim2_work_bits(det_bits) < Ibz::<N>::MAX_BITS);
    let outer_iters = dim2_outer_its(det_bits);
    let mut a = p.with_bound_ct(dim2_val_bits(det_bits));
    let mut b = r.with_bound_ct(dim2_val_bits(det_bits));
    let mut c = Ibz::<N>::set(0, 2).with_bound_ct(dim2_val_bits(det_bits));
    let mut d = Ibz::<N>::set(1, 2).with_bound_ct(dim2_val_bits(det_bits));
    partial_xgcd_lehmer(&mut a, &mut b, &mut c, &mut d, outer_iters, det_bits);
    (b.abs(), d.abs())
}

impl<const N: usize> Ibz<N> {
    /// All-ones when negative (public wrapper of the internal mask).
    #[inline]
    pub fn neg_mask_pub(&self) -> u64 {
        if self.is_positive() {
            0
        } else {
            u64::MAX
        }
    }
}
