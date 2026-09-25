//! Constant-time Lagrange-Gauss reduction over `Z[i]` and the reduction of
//! an `O0`-ideal basis built on it.

use super::{
    dim2_val_bits, dim2_work_bits, dim2i_herm_bits, dim2i_outer_its, dim2i_shift_max,
    DIM2I_INNER_ITS, DIM2I_THRESHOLD, DIM2I_UP_SHIFT, DIM2I_W,
};
use crate::mp::{ct, Ibz};
use crate::quat::{QuatAlg, QuatLattice};

/// A Gaussian integer with big components.
#[derive(Clone, Copy)]
struct CtComplex<const N: usize> {
    re: Ibz<N>,
    im: Ibz<N>,
}

#[derive(Clone, Copy)]
struct I32Complex {
    re: i32,
    im: i32,
}

#[derive(Clone, Copy)]
struct I64Complex {
    re: i64,
    im: i64,
}

impl<const N: usize> CtComplex<N> {
    fn zero() -> Self {
        Self {
            re: Ibz::zero(),
            im: Ibz::zero(),
        }
    }

    fn mul_si(&self, scalar: I32Complex, wbits: i32) -> Self {
        let t1 = self.re.mul_by_int_and_set_bound(scalar.re, wbits);
        let t2 = self
            .im
            .mul_by_int_and_set_bound(scalar.im.wrapping_neg(), wbits);
        let re = t1.add_and_set_bound(&t2, wbits);
        let t1 = self.re.mul_by_int_and_set_bound(scalar.im, wbits);
        let t2 = self.im.mul_by_int_and_set_bound(scalar.re, wbits);
        let im = t1.add_and_set_bound(&t2, wbits);
        Self { re, im }
    }

    fn shift_left(&self, extra_shift: i32, max_shift: i32, wbits: i32) -> Self {
        Self {
            re: self.re.ct_shl(extra_shift, max_shift, wbits),
            im: self.im.ct_shl(extra_shift, max_shift, wbits),
        }
    }

    fn add(&self, b: &Self, vbits: i32) -> Self {
        Self {
            re: self.re.add_and_set_bound(&b.re, vbits),
            im: self.im.add_and_set_bound(&b.im, vbits),
        }
    }

    fn sizeinbase_2(&self) -> u32 {
        ct::max_u32(self.re.bitsize_ct() as u32, self.im.bitsize_ct() as u32)
    }

    fn cond_swap(mask: u64, a: &mut Self, b: &mut Self) {
        Ibz::cswap(&mut a.re, &mut b.re, mask);
        Ibz::cswap(&mut a.im, &mut b.im, mask);
    }

    fn set_bound(&mut self, wbits: i32) {
        self.re.set_bound_ct(wbits);
        self.im.set_bound_ct(wbits);
    }

    fn from_rows(re: &Ibz<N>, im: &Ibz<N>, wbits: i32) -> Self {
        Self {
            re: re.with_bound_ct(wbits),
            im: im.with_bound_ct(wbits),
        }
    }
}

fn norm64_approx(a: I64Complex) -> u64 {
    ct::sqr_high(a.re).wrapping_add(ct::sqr_high(a.im))
}

fn swap_i32c(mask: u32, a: &mut I32Complex, b: &mut I32Complex) {
    ct::swap_i32(mask, &mut a.re, &mut b.re);
    ct::swap_i32(mask, &mut a.im, &mut b.im);
}

fn swap_i64c(mask: u64, a: &mut I64Complex, b: &mut I64Complex) {
    ct::swap_i64(mask, &mut a.re, &mut b.re);
    ct::swap_i64(mask, &mut a.im, &mut b.im);
}

fn mul_i32c(a: I32Complex, b: I32Complex) -> I32Complex {
    I32Complex {
        re: a
            .re
            .wrapping_mul(b.re)
            .wrapping_sub(a.im.wrapping_mul(b.im)),
        im: a
            .re
            .wrapping_mul(b.im)
            .wrapping_add(a.im.wrapping_mul(b.re)),
    }
}

fn shift_i64c(a: I64Complex, k: u32) -> I64Complex {
    I64Complex {
        re: ((a.re as u64) << (k & 63)) as i64,
        im: ((a.im as u64) << (k & 63)) as i64,
    }
}

fn shift_i32c(a: I32Complex, k: u32) -> I32Complex {
    I32Complex {
        re: ((a.re as u32) << (k & 31)) as i32,
        im: ((a.im as u32) << (k & 31)) as i32,
    }
}

fn apply_matrix_gaussian<const N: usize>(
    a: &mut CtComplex<N>,
    b: &mut CtComplex<N>,
    m: &[[I32Complex; 2]; 2],
    extra_shift: i32,
    in_bits: i32,
) {
    let wbits = dim2_work_bits(in_bits);
    let vbits = dim2_val_bits(in_bits);
    let t1 = a.mul_si(m[1][0], wbits);
    let t2 = b
        .shift_left(extra_shift, dim2i_shift_max(in_bits), wbits)
        .mul_si(m[0][1], wbits);
    let na = a.mul_si(m[0][0], wbits).add(&t2, vbits);
    let nb = b.mul_si(m[1][1], wbits).add(&t1, vbits);
    *a = na;
    *b = nb;
}

/// `ceil(log2(a) - log2(b))`, zero arguments treated as 1.
fn ceil_log2_diff_u64(a: u64, b: u64) -> i32 {
    let clz_a = (a | 1).leading_zeros() as i32;
    let clz_b = (b | 1).leading_zeros() as i32;
    let d = clz_b - clz_a;
    let sign_bit = (d as u32) >> 31;
    let mask_neg = ct::barrier_u32(0u32.wrapping_sub(sign_bit));
    let shift_b = (d as u32) & !mask_neg;
    let shift_a = ((-d) as u32) & mask_neg;
    let shifted_b = b << (shift_b & 63);
    let shifted_a = a << (shift_a & 63);
    let cond = (shifted_b.wrapping_sub(shifted_a) >> 63) as i32;
    d + cond
}

/// `round(a / b)` for `a / b` in `[-2, 2]`, `b > 0`, both below `2^62`.
fn round_ratio_62bit_ct(a: i64, b: u64) -> i32 {
    let max_abs = ct::abs64(a);
    let t1 = b >> 1;
    let t2 = b.wrapping_add(t1);
    let is_1 = (max_abs > t1) as u32;
    let is_2 = (max_abs > t2) as u32;
    let res_abs = is_1 + is_2;
    let sign_mask = (a >> 63) as u32;
    ((res_abs ^ sign_mask).wrapping_sub(sign_mask)) as i32
}

#[allow(clippy::too_many_arguments)]
fn lg_gaussian_inner(
    mut a: I64Complex,
    mut b: I64Complex,
    mut c: I64Complex,
    mut d: I64Complex,
    swap: bool,
    mut bits_from_v1: i32,
    mut bits_from_v2: i32,
) -> [[I32Complex; 2]; 2] {
    let mut u00 = I32Complex { re: 1, im: 0 };
    let mut u01 = I32Complex { re: 0, im: 0 };
    let mut u10 = I32Complex { re: 0, im: 0 };
    let mut u11 = I32Complex { re: 1, im: 0 };
    let mut valid_mask = u64::MAX;
    let global_swap_mask = ct::mask_bool(swap);
    for _ in 0..DIM2I_INNER_ITS {
        let clz_a = (ct::abs64(a.re) | ct::abs64(a.im) | 1).leading_zeros();
        let clz_b = (ct::abs64(b.re) | ct::abs64(b.im) | 1).leading_zeros();
        let clz_c = (ct::abs64(c.re) | ct::abs64(c.im) | 1).leading_zeros();
        let clz_d = (ct::abs64(d.re) | ct::abs64(d.im) | 1).leading_zeros();
        let mut base_shift_v1 = ct::max0(ct::min_u32(clz_a, clz_c) as i32 - DIM2I_UP_SHIFT);
        let mut base_shift_v2 = ct::max0(ct::min_u32(clz_b, clz_d) as i32 - DIM2I_UP_SHIFT);
        let limit = DIM2I_THRESHOLD + (64 - DIM2I_UP_SHIFT - DIM2I_W);
        let stop_v1 = ct::mask_bool(base_shift_v1 > limit) & ct::mask_bool(bits_from_v1 > 0);
        let stop_v2 = ct::mask_bool(base_shift_v2 > limit) & ct::mask_bool(bits_from_v2 > 0);
        valid_mask = ct::barrier(valid_mask & !(stop_v1 | stop_v2));
        let mut a_base = shift_i64c(a, base_shift_v1 as u32);
        let mut b_base = shift_i64c(b, base_shift_v2 as u32);
        let mut c_base = shift_i64c(c, base_shift_v1 as u32);
        let mut d_base = shift_i64c(d, base_shift_v2 as u32);
        let mut na = norm64_approx(a_base).wrapping_add(norm64_approx(c_base));
        let mut nb = norm64_approx(b_base).wrapping_add(norm64_approx(d_base));
        let base_shift_rel = base_shift_v2 - base_shift_v1;
        let v1_right = ct::select32(base_shift_rel < 0, -2 * base_shift_rel, 0) as u32;
        let v2_right = ct::select32(base_shift_rel > 0, 2 * base_shift_rel, 0) as u32;
        let lt_mask = 0u64.wrapping_sub(
            ct::shr_sat64(na, v1_right).wrapping_sub(ct::shr_sat64(nb, v2_right)) >> 63,
        );
        let swap_mask = ct::barrier(lt_mask & valid_mask);
        valid_mask = ct::barrier(valid_mask & (!swap_mask | global_swap_mask));
        let actual_swap = ct::barrier(swap_mask & global_swap_mask);
        ct::swap_u64(actual_swap, &mut na, &mut nb);
        swap_i64c(actual_swap, &mut a, &mut b);
        swap_i64c(actual_swap, &mut c, &mut d);
        swap_i32c(actual_swap as u32, &mut u00, &mut u10);
        swap_i32c(actual_swap as u32, &mut u01, &mut u11);
        swap_i64c(actual_swap, &mut a_base, &mut b_base);
        swap_i64c(actual_swap, &mut c_base, &mut d_base);
        ct::swap_i32(actual_swap as u32, &mut base_shift_v1, &mut base_shift_v2);
        ct::swap_i32(actual_swap as u32, &mut bits_from_v1, &mut bits_from_v2);
        let _ = na;
        let base_shift_rel = base_shift_v2 - base_shift_v1;
        let x = ct::mul_high(a_base.re, b_base.re)
            .wrapping_add(ct::mul_high(a_base.im, b_base.im))
            .wrapping_add(ct::mul_high(c_base.re, d_base.re))
            .wrapping_add(ct::mul_high(c_base.im, d_base.im));
        let y = ct::mul_high(a_base.im, b_base.re)
            .wrapping_sub(ct::mul_high(a_base.re, b_base.im))
            .wrapping_add(ct::mul_high(c_base.im, d_base.re))
            .wrapping_sub(ct::mul_high(c_base.re, d_base.im));
        let x_abs = ct::abs64(x);
        let y_abs = ct::abs64(y);
        let s = ceil_log2_diff_u64(ct::max_u64(x_abs, y_abs), nb) - 1;
        let k = ct::select_u32(base_shift_rel + s < 0, 0, (base_shift_rel + s) as u32);
        let b_shifted = shift_i64c(b, k);
        let d_shifted = shift_i64c(d, k);
        let shift_diff = base_shift_rel - k as i32;
        let shift_left_x = (shift_diff & !(shift_diff >> 31)) as u32;
        let shift_left_nb = ((-shift_diff) & !((-shift_diff) >> 31)) as u32;
        let x_adj = ((x as u64) << (shift_left_x & 63)) as i64;
        let y_adj = ((y as u64) << (shift_left_x & 63)) as i64;
        let nb_adj = nb << (shift_left_nb & 63);
        let xround = round_ratio_62bit_ct(x_adj, nb_adj);
        let yround = round_ratio_62bit_ct(y_adj, nb_adj);
        let zero_reduction = ct::mask_bool(xround == 0) & ct::mask_bool(yround == 0);
        valid_mask = ct::barrier(valid_mask & !zero_reduction);
        let xr = xround as i64;
        let yr = yround as i64;
        let b_final = I64Complex {
            re: xr
                .wrapping_mul(b_shifted.re)
                .wrapping_sub(yr.wrapping_mul(b_shifted.im)),
            im: xr
                .wrapping_mul(b_shifted.im)
                .wrapping_add(yr.wrapping_mul(b_shifted.re)),
        };
        let d_final = I64Complex {
            re: xr
                .wrapping_mul(d_shifted.re)
                .wrapping_sub(yr.wrapping_mul(d_shifted.im)),
            im: xr
                .wrapping_mul(d_shifted.im)
                .wrapping_add(yr.wrapping_mul(d_shifted.re)),
        };
        a.re = a.re.wrapping_sub(b_final.re & valid_mask as i64);
        a.im = a.im.wrapping_sub(b_final.im & valid_mask as i64);
        c.re = c.re.wrapping_sub(d_final.re & valid_mask as i64);
        c.im = c.im.wrapping_sub(d_final.im & valid_mask as i64);
        let mut u10_shift = shift_i32c(u10, k);
        let mut u11_shift = shift_i32c(u11, k);
        let xy = I32Complex {
            re: xround,
            im: yround,
        };
        u10_shift = mul_i32c(xy, u10_shift);
        u11_shift = mul_i32c(xy, u11_shift);
        let vm = valid_mask as u32;
        u00.re = u00.re.wrapping_sub((u10_shift.re as u32 & vm) as i32);
        u00.im = u00.im.wrapping_sub((u10_shift.im as u32 & vm) as i32);
        u01.re = u01.re.wrapping_sub((u11_shift.re as u32 & vm) as i32);
        u01.im = u01.im.wrapping_sub((u11_shift.im as u32 & vm) as i32);
    }
    [[u00, u01], [u10, u11]]
}

/// Multiply `(key, x, y, z)` by the unit of `Z[i]` that makes both parts of
/// `key` non-negative.
fn orient_canonical<const N: usize>(
    key: &mut CtComplex<N>,
    x: Option<&mut CtComplex<N>>,
    y: Option<&mut CtComplex<N>>,
    z: Option<&mut CtComplex<N>>,
) {
    let neg_re = key.re.neg_mask_pub();
    let neg_im = key.im.neg_mask_pub();
    let swap_mask = neg_re ^ neg_im;
    let mut targets: [Option<&mut CtComplex<N>>; 3] = [x, y, z];
    let apply_swap = |t: &mut CtComplex<N>| {
        Ibz::cswap(&mut t.re, &mut t.im, swap_mask);
        t.re = t.re.cneg(swap_mask);
    };
    apply_swap(key);
    for t in targets.iter_mut().flatten() {
        apply_swap(t);
    }
    let neg_mask = key.re.neg_mask_pub() | key.im.neg_mask_pub();
    let apply_neg = |t: &mut CtComplex<N>| {
        t.re = t.re.cneg(neg_mask);
        t.im = t.im.cneg(neg_mask);
    };
    apply_neg(key);
    for t in targets.iter_mut().flatten() {
        apply_neg(t);
    }
}

/// Lagrange-Gauss reduce the basis `((aa, cc), (bb, dd))` over `Z[i]`, in
/// place, carrying the optional pair `(ee, ff)` along.
#[allow(clippy::too_many_arguments)]
fn lehmer_gaussian<const N: usize>(
    aa: &mut CtComplex<N>,
    bb: &mut CtComplex<N>,
    cc: &mut CtComplex<N>,
    dd: &mut CtComplex<N>,
    mut ef: Option<(&mut CtComplex<N>, &mut CtComplex<N>)>,
    outer_iters: i32,
    in_bits: i32,
) {
    for _ in 0..outer_iters {
        let bits_a = aa.sizeinbase_2() as i32;
        let bits_b = bb.sizeinbase_2() as i32;
        let bits_c = cc.sizeinbase_2() as i32;
        let bits_d = dd.sizeinbase_2() as i32;
        let mut bits_v1 = ct::max(bits_a, bits_c);
        let mut bits_v2 = ct::max(bits_b, bits_d);
        let swap_mask = ct::mask_bool(bits_v1 < bits_v2);
        CtComplex::cond_swap(swap_mask, aa, bb);
        CtComplex::cond_swap(swap_mask, cc, dd);
        if let Some((ee, ff)) = ef.as_mut() {
            CtComplex::cond_swap(swap_mask, ee, ff);
        }
        ct::swap_i32(swap_mask as u32, &mut bits_v1, &mut bits_v2);
        let rho = bits_v1 - bits_v2 + 1;
        let bits_from_v1 = ct::select32(bits_v1 >= DIM2I_W, bits_v1 - DIM2I_W, 0);
        let large_diff = rho >= DIM2I_THRESHOLD;
        let bfv2_ideal = ct::max0(bits_v2 - (DIM2I_W - DIM2I_THRESHOLD + 1));
        let bits_from_v2 = ct::select32(large_diff, bfv2_ideal, bits_from_v1);
        let extra_shift = ct::select32(large_diff, rho - DIM2I_THRESHOLD, 0);
        let pa = I64Complex {
            re: aa.re.extract_i64(bits_from_v1),
            im: aa.im.extract_i64(bits_from_v1),
        };
        let pb = I64Complex {
            re: bb.re.extract_i64(bits_from_v2),
            im: bb.im.extract_i64(bits_from_v2),
        };
        let pc = I64Complex {
            re: cc.re.extract_i64(bits_from_v1),
            im: cc.im.extract_i64(bits_from_v1),
        };
        let pd = I64Complex {
            re: dd.re.extract_i64(bits_from_v2),
            im: dd.im.extract_i64(bits_from_v2),
        };
        let m = lg_gaussian_inner(pa, pb, pc, pd, !large_diff, bits_from_v1, bits_from_v2);
        apply_matrix_gaussian(aa, bb, &m, extra_shift, in_bits);
        apply_matrix_gaussian(cc, dd, &m, extra_shift, in_bits);
        if let Some((ee, ff)) = ef.as_mut() {
            apply_matrix_gaussian(ee, ff, &m, extra_shift, in_bits);
        }
    }
    CtComplex::cond_swap(u64::MAX, aa, bb);
    CtComplex::cond_swap(u64::MAX, cc, dd);
    if let Some((ee, ff)) = ef.as_mut() {
        CtComplex::cond_swap(u64::MAX, ee, ff);
    }
    match ef {
        Some((ee, ff)) => {
            orient_canonical(aa, Some(cc), Some(ee), None);
            orient_canonical(bb, Some(dd), Some(ff), None);
        }
        None => {
            orient_canonical(aa, Some(cc), None, None);
            orient_canonical(bb, Some(dd), None, None);
        }
    }
}

/// Exact Hermitian inner product `conj(v1) w1 + p conj(v2) w2`.
fn zi_hermitian<const N: usize>(
    v1: &CtComplex<N>,
    v2: &CtComplex<N>,
    w1: &CtComplex<N>,
    w2: &CtComplex<N>,
    p: &Ibz<N>,
    hbits: i32,
) -> CtComplex<N> {
    let b = |x: Ibz<N>| x.with_bound_ct(hbits);
    let t1 = b(v2.re.mul(&w2.re));
    let t2 = b(v2.im.mul(&w2.im));
    let t3 = t1.add_and_set_bound(&t2, hbits);
    let acc = b(t3.mul(p));
    let t1 = b(v1.re.mul(&w1.re));
    let t2 = b(v1.im.mul(&w1.im));
    let t3 = t1.add_and_set_bound(&t2, hbits);
    let re = t3.add_and_set_bound(&acc, hbits);
    let t1 = b(v2.re.mul(&w2.im));
    let t2 = b(v2.im.mul(&w2.re));
    let t3 = b(t1.sub(&t2));
    let acc = b(t3.mul(p));
    let t1 = b(v1.re.mul(&w1.im));
    let t2 = b(v1.im.mul(&w1.re));
    let t3 = b(t1.sub(&t2));
    let im = t3.add_and_set_bound(&acc, hbits);
    CtComplex { re, im }
}

fn write_zi_column_pair<const N: usize>(
    basis: &mut crate::quat::Mat4x4<N>,
    c: usize,
    z1: &CtComplex<N>,
    z2: &CtComplex<N>,
    wbits: i32,
) {
    basis.0[0][c] = z1.re;
    basis.0[1][c] = z1.im;
    basis.0[2][c] = z2.re;
    basis.0[3][c] = z2.im;
    basis.0[0][c + 1] = z1.im.neg();
    basis.0[1][c + 1] = z1.re;
    basis.0[2][c + 1] = z2.im.neg();
    basis.0[3][c + 1] = z2.re;
    for r in 0..4 {
        basis.0[r][c].set_bound_ct(wbits);
        basis.0[r][c + 1].set_bound_ct(wbits);
    }
}

/// Reduced basis of an `O0`-ideal in inert HNF (denominator 2): columns
/// `v1, i v1, v2, i v2` of a Lagrange-reduced `Z[i]`-basis, canonically
/// oriented. `bound` is a public bound on the input coefficients.
pub fn reduce_o0_ideal<const N: usize>(
    hnf: &QuatLattice<N>,
    alg: &QuatAlg<N>,
    bound: i32,
) -> QuatLattice<N> {
    let s = alg.sqrt_p.with_bound_ct(alg.sqrt_p_bits + 1);
    #[cfg(debug_assertions)]
    {
        let sq = s.mul(&s);
        let s1 = s.add_int_and_set_bound(1, alg.sqrt_p_bits + 2);
        let sq1 = s1.mul(&s1);
        debug_assert!(sq <= alg.p && alg.p < sq1);
    }
    let in_bits = bound.max(alg.sqrt_p_bits + 1);
    let wbits = dim2_val_bits(in_bits);
    let hbits = dim2i_herm_bits(in_bits, alg.p_bits);
    debug_assert!(dim2_work_bits(in_bits) < Ibz::<N>::MAX_BITS);
    debug_assert!(hbits <= Ibz::<N>::MAX_BITS);
    let hb = &hnf.basis.0;
    let mut aa = CtComplex::from_rows(&hb[0][0], &hb[1][0], wbits);
    let mut ee = CtComplex::from_rows(&hb[2][0], &hb[3][0], wbits);
    let mut bb = CtComplex::from_rows(&hb[0][2], &hb[1][2], wbits);
    let mut ff = CtComplex::from_rows(&hb[2][2], &hb[3][2], wbits);
    let mut cc = CtComplex::zero();
    let mut dd = CtComplex::zero();
    cc.re = ee.re.mul(&s).with_bound_ct(wbits);
    cc.im = ee.im.mul(&s).with_bound_ct(wbits);
    dd.re = ff.re.mul(&s).with_bound_ct(wbits);
    dd.im = ff.im.mul(&s).with_bound_ct(wbits);
    lehmer_gaussian(
        &mut aa,
        &mut bb,
        &mut cc,
        &mut dd,
        Some((&mut ee, &mut ff)),
        dim2i_outer_its(in_bits, alg.sqrt_p_bits),
        in_bits,
    );
    let mut h = zi_hermitian(&aa, &ee, &bb, &ff, &alg.p, hbits);
    orient_canonical(&mut h, Some(&mut bb), Some(&mut dd), Some(&mut ff));
    #[cfg(debug_assertions)]
    {
        let chk = zi_hermitian(&aa, &ee, &bb, &ff, &alg.p, hbits);
        debug_assert!(chk.re.is_positive() && chk.im.is_positive());
    }
    aa.set_bound(wbits);
    bb.set_bound(wbits);
    cc.set_bound(wbits);
    dd.set_bound(wbits);
    ee.set_bound(wbits);
    ff.set_bound(wbits);
    let mut reduced = QuatLattice::init();
    write_zi_column_pair(&mut reduced.basis, 0, &aa, &ee, wbits);
    write_zi_column_pair(&mut reduced.basis, 2, &bb, &ff, wbits);
    reduced.denom = hnf.denom;
    reduced
}
