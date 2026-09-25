//! The Lehmer Lagrange-Gauss kernel on a `2 x 2` Gram matrix and the LDL
//! block it consumes.

use crate::mp::{ct, Ibz};
use crate::quat::LllConfig;

/// Positive-definite binary form `[[a, b], [b, c]]` scaled by `2^-e`.
#[derive(Clone, Copy)]
pub struct Gram2<const N: usize> {
    /// Top-left entry.
    pub a: Ibz<N>,
    /// Off-diagonal entry.
    pub b: Ibz<N>,
    /// Bottom-right entry.
    pub c: Ibz<N>,
    /// The (secret) scale: the form is `[[a, b], [b, c]] / 2^e`.
    pub e: i32,
}

/// The local `2 x 2` Gram from an LDL decomposition, truncated to a common
/// secret scale so that all three entries fit `gram_tot_bits`.
pub fn materialise_block<const N: usize>(
    dk: &Ibz<N>,
    dkp1: &Ibz<N>,
    mu: &Ibz<N>,
    p: i32,
    cfg: &LllConfig,
) -> Gram2<N> {
    let wbits = cfg.e_bits + 16;
    let b_full = mu.ct_highmul_p(dk, p, wbits);
    let mu_bfull = mu.ct_highmul_p(&b_full, p, wbits);
    let c_full = dkp1.add_and_set_bound(&mu_bfull, wbits);
    let sum_ac = dk.add_and_set_bound(&c_full, wbits + 1);
    let bs = sum_ac.bitsize_ct();
    let e = ct::max0(bs - (cfg.gram_tot_bits - 1));
    let max_e = wbits + 1;
    let out = Gram2 {
        a: dk.ct_shift(-e, max_e, cfg.gram_tot_bits),
        b: b_full.ct_shift(-e, max_e, cfg.gram_tot_bits),
        c: c_full.ct_shift(-e, max_e, cfg.gram_tot_bits),
        e,
    };
    debug_assert!(out.a.is_positive() && !out.a.is_zero());
    debug_assert!(out.c.is_positive() && !out.c.is_zero());
    out
}

const THRESHOLD_PARAM: i32 = 29;
const EXTRACT_W_PARAM: i32 = 2 * THRESHOLD_PARAM + 3;

fn lg_extra_shift_max(cfg: &LllConfig) -> i32 {
    (cfg.gram_tot_bits - THRESHOLD_PARAM + 2) / 2
}

fn fast_log2_ratio(a: u64, b: u64, clz_b: u32) -> i32 {
    let clz_a = (a | 1).leading_zeros();
    let norm_a = (a << clz_a) >> 2;
    let norm_b = (b << clz_b) >> 2;
    let three_b = 3u64.wrapping_mul(norm_b);
    let inc = ((norm_a << 1) >= three_b) as i32;
    let dec = ((norm_a << 2) < three_b) as i32;
    (clz_b as i32 - clz_a as i32) + inc - dec
}

/// MSB-based reduction of the Gram `[[a, b], [b, c]]`; column transform.
fn partial_gram_msb(mut a: u64, mut b: i64, mut c: u64, swap: bool) -> [[i32; 2]; 2] {
    let (mut u00, mut u01, mut u10, mut u11) = (1i32, 0i32, 0i32, 1i32);
    let mut valid_mask = u64::MAX;
    let global_swap_mask = ct::mask_bool(swap);
    for _ in 0..THRESHOLD_PARAM / 2 + 2 {
        valid_mask &= !ct::mask_bool(a == 0 || c == 0);
        let negate_mask = ct::mask_bool(b < 0) & valid_mask;
        let b_abs = ct::abs64(b);
        b = b_abs as i64;
        u00 = ct::negate_mask_i32(negate_mask as u32, u00);
        u10 = ct::negate_mask_i32(negate_mask as u32, u10);
        let swap_mask = ct::mask_bool(c < a) & valid_mask;
        valid_mask &= !swap_mask | global_swap_mask;
        let actual_swap = swap_mask & global_swap_mask;
        ct::swap_u64(actual_swap, &mut a, &mut c);
        ct::swap_i32(actual_swap as u32, &mut u00, &mut u01);
        ct::swap_i32(actual_swap as u32, &mut u10, &mut u11);
        valid_mask &= !ct::mask_bool(2u64.wrapping_mul(b_abs) <= a);
        let clz_a = (a | 1).leading_zeros();
        let clz_c = (c | 1).leading_zeros();
        let lim = (64 - EXTRACT_W_PARAM + THRESHOLD_PARAM) as u32;
        valid_mask &= !ct::mask_bool(clz_a >= lim || clz_c >= lim);
        let mut k = fast_log2_ratio(b_abs, a, clz_a);
        k = ct::max(0, k) & 63;
        let sub = ((b_abs << ((k + 1) & 63)) as i64).wrapping_sub((a << ((2 * k) & 63)) as i64);
        c = c.wrapping_sub((sub & valid_mask as i64) as u64);
        b = b.wrapping_sub(((a << k) & valid_mask) as i64);
        u01 = u01.wrapping_sub((((u00 as u32) << (k & 31)) & valid_mask as u32) as i32);
        u11 = u11.wrapping_sub((((u10 as u32) << (k & 31)) & valid_mask as u32) as i32);
    }
    [[u00, u01], [u10, u11]]
}

fn apply_matrix_2x2_shift<const N: usize>(
    g: &mut Gram2<N>,
    m: &[[i32; 2]; 2],
    extra_shift: i32,
    wbits: i32,
    cfg: &LllConfig,
) {
    let (m00, m01, m10, m11) = (m[0][0], m[0][1], m[1][0], m[1][1]);
    let c_a0 = m00.wrapping_mul(m00);
    let c_a1 = 2i32.wrapping_mul(m00).wrapping_mul(m10);
    let c_a2 = m10.wrapping_mul(m10);
    let c_b0 = m00.wrapping_mul(m01);
    let c_b1 = m00.wrapping_mul(m11).wrapping_add(m01.wrapping_mul(m10));
    let c_b2 = m10.wrapping_mul(m11);
    let c_c0 = m01.wrapping_mul(m01);
    let c_c1 = 2i32.wrapping_mul(m01).wrapping_mul(m11);
    let c_c2 = m11.wrapping_mul(m11);
    let max_e = lg_extra_shift_max(cfg);
    let mut na = g.a.mul_by_int_and_set_bound(c_a0, wbits);
    na = na.add_and_set_bound(&g.b.mul_by_int_and_set_bound(c_a1, wbits), wbits);
    na = na.add_and_set_bound(&g.c.mul_by_int_and_set_bound(c_a2, wbits), wbits);
    let a_shifted = g.a.ct_shl(extra_shift, max_e, wbits);
    let b_shifted = g.b.ct_shl(extra_shift, max_e, wbits);
    let mut nb = a_shifted.mul_by_int_and_set_bound(c_b0, wbits);
    nb = nb.add_and_set_bound(&g.b.mul_by_int_and_set_bound(c_b1, wbits), wbits);
    nb = nb.add_and_set_bound(&g.c.mul_by_int_and_set_bound(c_b2, wbits), wbits);
    let a_shifted2 = g.a.ct_shl(2 * extra_shift, 2 * max_e, wbits);
    let mut nc = a_shifted2.mul_by_int_and_set_bound(c_c0, wbits);
    nc = nc.add_and_set_bound(&b_shifted.mul_by_int_and_set_bound(c_c1, wbits), wbits);
    nc = nc.add_and_set_bound(&g.c.mul_by_int_and_set_bound(c_c2, wbits), wbits);
    g.a = na;
    g.b = nb;
    g.c = nc;
}

fn apply_matrix_to_u_shift<const N: usize>(
    u: &mut [[Ibz<N>; 2]; 2],
    m: &[[i32; 2]; 2],
    extra_shift: i32,
    cfg: &LllConfig,
) {
    let ub = cfg.u_bits;
    let max_e = lg_extra_shift_max(cfg);
    let nu00 = u[0][0]
        .mul_by_int_and_set_bound(m[0][0], ub)
        .add_and_set_bound(&u[0][1].mul_by_int_and_set_bound(m[1][0], ub), ub);
    let nu10 = u[1][0]
        .mul_by_int_and_set_bound(m[0][0], ub)
        .add_and_set_bound(&u[1][1].mul_by_int_and_set_bound(m[1][0], ub), ub);
    let u_shifted = u[0][0].ct_shl(extra_shift, max_e, ub);
    let nu01 = u_shifted
        .mul_by_int_and_set_bound(m[0][1], ub)
        .add_and_set_bound(&u[0][1].mul_by_int_and_set_bound(m[1][1], ub), ub);
    let u_shifted = u[1][0].ct_shl(extra_shift, max_e, ub);
    let nu11 = u_shifted
        .mul_by_int_and_set_bound(m[0][1], ub)
        .add_and_set_bound(&u[1][1].mul_by_int_and_set_bound(m[1][1], ub), ub);
    u[0][0] = nu00;
    u[1][0] = nu10;
    u[0][1] = nu01;
    u[1][1] = nu11;
}

/// Lagrange-Gauss reduce a positive-definite binary form in place and
/// return the column-convention unimodular transform. `outer_iters` is a
/// public trip count.
pub fn gram_lehmer<const N: usize>(
    g: &mut Gram2<N>,
    outer_iters: i32,
    cfg: &LllConfig,
) -> [[Ibz<N>; 2]; 2] {
    let ub = cfg.u_bits;
    let mut u = [
        [
            Ibz::<N>::set(1, 2).with_bound_ct(ub),
            Ibz::<N>::set(0, 2).with_bound_ct(ub),
        ],
        [
            Ibz::<N>::set(0, 2).with_bound_ct(ub),
            Ibz::<N>::set(1, 2).with_bound_ct(ub),
        ],
    ];
    g.a.set_bound_ct(cfg.gram_work_bits);
    g.b.set_bound_ct(cfg.gram_work_bits);
    g.c.set_bound_ct(cfg.gram_work_bits);
    for _ in 0..outer_iters {
        let diff = g.c.sub(&g.a);
        let borrow = diff.neg_mask_pub();
        Ibz::cswap(&mut g.a, &mut g.c, borrow);
        let [[mut u00, mut u01], [mut u10, mut u11]] = u;
        Ibz::cswap(&mut u00, &mut u01, borrow);
        Ibz::cswap(&mut u10, &mut u11, borrow);
        u = [[u00, u01], [u10, u11]];
        let bits_a = g.a.bitsize_ct();
        let bits_c = g.c.bitsize_ct();
        let rho = bits_c - bits_a;
        let large_diff = rho >= THRESHOLD_PARAM;
        let bits_from_c = ct::max(bits_c - EXTRACT_W_PARAM, 0);
        let extra_shift = ct::select32(large_diff, (rho - THRESHOLD_PARAM + 2) / 2, 0);
        let bits_from_a = ct::select32(large_diff, bits_from_c - 2 * extra_shift, bits_from_c);
        let bits_from_b = ct::select32(large_diff, bits_from_c - extra_shift, bits_from_c);
        let pa = g.a.extract_u64(bits_from_a);
        let pb = g.b.extract_i64(bits_from_b);
        let pc = g.c.extract_u64(bits_from_c);
        let m = partial_gram_msb(pa, pb, pc, !large_diff);
        apply_matrix_2x2_shift(g, &m, extra_shift, cfg.gram_work_bits, cfg);
        apply_matrix_to_u_shift(&mut u, &m, extra_shift, cfg);
    }
    for row in u.iter_mut() {
        for v in row.iter_mut() {
            v.set_bound_ct(ub);
        }
    }
    u
}
