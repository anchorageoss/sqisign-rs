//! Constant-time dimension-4 reduction of a general ideal in inert HNF
//! shape: a dual LDL profile, LLL tours of Lagrange-Gauss block
//! reductions, and a Minkowski finish that sorts 27 candidates with a
//! bitonic network. The output is canonical over the ideal class up to a
//! global sign, with failure probability below `2^-lambda` on honest
//! inputs.

use super::lg2::{gram_lehmer, materialise_block};
use crate::mp::{ct, Ibz};
use crate::quat::{LllConfig, Mat4x4, QuatAlg, QuatLattice};

/// LDL state: `D` raw, `L` at scale `2^P` with `L[j][j] = 2^P`, `Ainv` the
/// accumulated inverse transform.
struct FpLdl<const N: usize> {
    d: [Ibz<N>; 4],
    l: [[Ibz<N>; 4]; 4],
    ainv: [[Ibz<N>; 4]; 4],
}

struct MkResult<const N: usize> {
    t: [[i8; 4]; 4],
    ainv_total: [[Ibz<N>; 4]; 4],
    g: [[Ibz<N>; 4]; 4],
}

fn set_neg_ratio<const N: usize>(
    x: &Ibz<N>,
    r: &Ibz<N>,
    s2n: i32,
    max_shift: i32,
    cfg: &LllConfig,
) -> Ibz<N> {
    let ratio = x.ct_highmul_s(r, -s2n, max_shift, cfg.e_bits);
    ratio.neg().with_bound_ct(cfg.e_bits)
}

impl<const N: usize> FpLdl<N> {
    #[allow(clippy::too_many_arguments)]
    fn init_dual(
        two_n: &Ibz<N>,
        b: &Ibz<N>,
        c: &Ibz<N>,
        dq: &Ibz<N>,
        e: &Ibz<N>,
        p: &Ibz<N>,
        pp: i32,
        cfg: &LllConfig,
    ) -> Self {
        let mut st = Self {
            d: [Ibz::zero(); 4],
            l: [[Ibz::zero(); 4]; 4],
            ainv: [[Ibz::zero(); 4]; 4],
        };
        for i in 0..4 {
            for j in 0..4 {
                st.l[i][j] = Ibz::set(0, 2).with_bound_ct(cfg.l_bits);
            }
        }
        for i in 0..4 {
            st.l[i][i] = Ibz::set(1, 2).mul_2exp(pp as u32).with_bound_ct(cfg.l_bits);
        }
        let two_n_sq = two_n
            .mul(two_n)
            .ct_shift(-cfg.d_shift, cfg.d_shift, cfg.e_bits);
        st.d[0] = two_n_sq;
        st.d[1] = two_n_sq;
        let p_local = p.ct_shift(-cfg.d_shift, cfg.d_shift, cfg.e_bits);
        st.d[2] = p_local;
        st.d[3] = p_local;
        let (r, s2n) = Ibz::ct_fp_recip(two_n, pp, cfg.nr_rounds);
        let max_shift = two_n.get_bound();
        st.l[2][0] = set_neg_ratio(e, &r, s2n, max_shift, cfg);
        st.l[2][1] = set_neg_ratio(dq, &r, s2n, max_shift, cfg);
        st.l[3][0] = set_neg_ratio(c, &r, s2n, max_shift, cfg);
        st.l[3][1] = set_neg_ratio(b, &r, s2n, max_shift, cfg);
        for i in 0..4 {
            for j in 0..4 {
                st.ainv[i][j] = Ibz::set((i == j) as i64, 2).with_bound_ct(cfg.u_bits);
            }
        }
        st.size_reduce_from(0, pp, cfg);
        st
    }

    /// Apply the column-convention transform `u` to block `(k, k+1)`.
    fn update_block(&mut self, k: usize, u: &[[Ibz<N>; 2]; 2], pp: i32, cfg: &LllConfig) {
        let kp1 = k + 1;
        let mu_old = self.l[kp1][k].with_bound_ct(cfg.l_bits);
        let bu = cfg.block_u_bits;
        let u00 = u[0][0].with_bound_ct(bu);
        let u01 = u[0][1].with_bound_ct(bu);
        let u10 = u[1][0].with_bound_ct(bu);
        let u11 = u[1][1].with_bound_ct(bu);
        let t1 = u00.mul(&u11).with_bound_ct(2 * bu);
        let t2 = u01.mul(&u10).with_bound_ct(2 * bu);
        let det = t1.sub(&t2).with_bound_ct(32);
        debug_assert!(
            det.cmp_i32(1) == core::cmp::Ordering::Equal
                || det.cmp_i32(-1) == core::cmp::Ordering::Equal
        );
        let negate_det = det.ct_lt_mask(&Ibz::zero());
        // STEP 1: local 2x2 Gram in the orthogonalised space
        let s0 = u00.mul_2exp(pp as u32);
        let m0 = mu_old.mul(&u10);
        let eta0 = s0.add(&m0);
        debug_assert!(eta0.bitsize_ct() <= pp + 1);
        let eta0 = eta0.with_bound_ct(pp + 2);
        let s1 = u01.mul_2exp(pp as u32).with_bound_ct(cfg.l_bits);
        let m1 = mu_old.mul(&u11);
        let eta1 = s1.add(&m1).with_bound_ct(cfg.l_bits);
        debug_assert!(eta1.bitsize_ct() <= pp + 1);
        let eta0_sq = eta0.ct_highmul_p(&eta0, pp, pp + 2);
        debug_assert!(eta0_sq.bitsize_ct() <= pp + 1);
        let g00 = {
            let term1 = eta0_sq.ct_highmul_p(&self.d[k], pp, cfg.e_bits);
            let term2 = self.d[kp1].mul(&u10).with_bound_ct(cfg.e_bits);
            let term2 = term2.mul(&u10).with_bound_ct(cfg.e_bits);
            term1.add(&term2).with_bound_ct(cfg.e_bits)
        };
        let g01 = {
            let eta0_eta1 = eta0.ct_highmul_p(&eta1, pp, cfg.e_bits);
            let term1 = eta0_eta1.ct_highmul_p(&self.d[k], pp, cfg.e_bits);
            let term2 = self.d[kp1].mul(&u10).with_bound_ct(cfg.e_bits);
            let term2 = term2.mul(&u11).with_bound_ct(cfg.e_bits);
            term1.add(&term2).with_bound_ct(cfg.e_bits)
        };
        // STEP 2: diagonal norms and block coefficient
        let (inv_g00, s_g00) = Ibz::ct_fp_recip(&g00, pp, cfg.nr_rounds);
        let new_mu = g01.ct_highmul_s(&inv_g00, -s_g00, cfg.l_bits, cfg.l_bits);
        let temp_inv = self.d[k + 1].ct_highmul_s(&inv_g00, -s_g00, cfg.e_bits, cfg.e_bits);
        self.d[kp1] = temp_inv.ct_highmul_p(&self.d[k], pp, cfg.e_bits);
        self.d[k] = g00.with_bound_ct(cfg.e_bits);
        debug_assert!(self.d[k].is_positive() && !self.d[k].is_zero());
        debug_assert!(self.d[kp1].is_positive() && !self.d[kp1].is_zero());
        self.l[kp1][k] = new_mu.with_bound_ct(cfg.l_bits);
        // STEP 3: left-tail projections
        for j in 0..k {
            let a1 = self.l[k][j].mul(&u00).with_bound_ct(cfg.l_bits);
            let a2 = self.l[kp1][j].mul(&u10).with_bound_ct(cfg.l_bits);
            let new_lk = a1.add(&a2).with_bound_ct(cfg.l_bits);
            let b1 = self.l[k][j].mul(&u01).with_bound_ct(cfg.l_bits);
            let b2 = self.l[kp1][j].mul(&u11).with_bound_ct(cfg.l_bits);
            let new_lkp1 = b1.add(&b2).with_bound_ct(cfg.l_bits);
            self.l[k][j] = new_lk.with_bound_ct(cfg.l_bits);
            self.l[kp1][j] = new_lkp1.with_bound_ct(cfg.l_bits);
        }
        // STEP 4: right-tail projections
        for i in kp1 + 1..4 {
            let old_lik = self.l[i][k];
            let old_likp1 = self.l[i][kp1];
            let m_ikp1 = mu_old.ct_highmul_p(&old_likp1, pp, cfg.l_bits);
            let e0a = old_likp1.mul(&u00).with_bound_ct(cfg.l_bits);
            let e0b = m_ikp1.mul(&u10);
            debug_assert!(e0b.bitsize_ct() < cfg.l_bits);
            let e0b = e0b.with_bound_ct(cfg.l_bits);
            let eta0_l = e0a.add(&e0b).with_bound_ct(cfg.l_bits);
            let e1a = old_likp1.mul(&u01).with_bound_ct(cfg.l_bits);
            let e1b = m_ikp1.mul(&u11).with_bound_ct(cfg.l_bits);
            let eta1_l = e1a.add(&e1b).with_bound_ct(cfg.l_bits);
            let u10_lik = old_lik.mul(&u10).with_bound_ct(cfg.l_bits);
            let term_kp1 = eta0_l.sub(&u10_lik).with_bound_ct(cfg.l_bits);
            let new_likp1 = term_kp1.cneg(negate_det);
            let u11_lik = old_lik.mul(&u11).with_bound_ct(cfg.l_bits);
            let term_k = u11_lik.sub(&eta1_l).with_bound_ct(cfg.l_bits);
            let t_0 = term_k.cneg(negate_det);
            let hm = new_mu.ct_highmul_p(&new_likp1, pp, cfg.l_bits);
            let new_lik = t_0.add(&hm).with_bound_ct(cfg.l_bits);
            self.l[i][k] = new_lik.with_bound_ct(cfg.l_bits);
            self.l[i][kp1] = new_likp1.with_bound_ct(cfg.l_bits);
        }
        // Ainv column update
        for r in 0..4 {
            let x = self.ainv[r][k];
            let y = self.ainv[r][kp1];
            let p1 = u11.mul(&x).with_bound_ct(cfg.u_bits + cfg.block_u_bits);
            let p2 = u01.mul(&y);
            let newk = p1.sub(&p2).with_bound_ct(cfg.u_bits).cneg(negate_det);
            let p3 = u10.mul(&x);
            let p4 = u00.mul(&y);
            let newkp1 = p4.sub(&p3).with_bound_ct(cfg.u_bits).cneg(negate_det);
            self.ainv[r][k] = newk.with_bound_ct(cfg.u_bits);
            self.ainv[r][kp1] = newkp1.with_bound_ct(cfg.u_bits);
        }
    }

    /// Restore `|L[r][j]| <= 2^(P-1)` for rows `r >= k`, mirroring every
    /// step on `Ainv`.
    fn size_reduce_from(&mut self, k: usize, pp: i32, cfg: &LllConfig) {
        for r in k..4 {
            for j in (0..r).rev() {
                let q = self.l[r][j].ct_round_shift_limb(pp, cfg.head);
                for m in 0..j {
                    debug_assert!(self.l[j][m].bitsize_ct() <= pp);
                    self.l[j][m].set_bound_ct(cfg.p_bits);
                    let prod = self.l[j][m].mul(&q).with_bound_ct(cfg.l_bits);
                    self.l[r][m] = self.l[r][m].sub(&prod).with_bound_ct(cfg.l_bits);
                }
                let prod = q.mul_2exp(pp as u32).with_bound_ct(cfg.l_bits);
                self.l[r][j] = self.l[r][j].sub(&prod).with_bound_ct(cfg.l_bits);
                for i in 0..4 {
                    let prod = self.ainv[i][r].mul(&q).with_bound_ct(cfg.u_bits + cfg.head);
                    self.ainv[i][j] = self.ainv[i][j].add(&prod).with_bound_ct(cfg.u_bits);
                }
            }
        }
    }

    fn reduce_block(&mut self, k: usize, pp: i32, round: i32, cfg: &LllConfig) {
        let kp1 = k + 1;
        let mut g = materialise_block(&self.d[k], &self.d[kp1], &self.l[kp1][k], pp, cfg);
        let u = gram_lehmer(&mut g, cfg.lg_outer_its_for_round(round), cfg);
        self.update_block(k, &u, pp, cfg);
        self.size_reduce_from(k, pp, cfg);
    }

    fn lll_tour4(&mut self, pp: i32, cfg: &LllConfig) {
        for t in 0..cfg.tours {
            self.reduce_block(1, pp, 3 * t, cfg);
            self.reduce_block(0, pp, 3 * t + 1, cfg);
            self.reduce_block(2, pp, 3 * t + 2, cfg);
        }
    }
}

const MK_N_REAL: usize = 27;
const MK_N_TOTAL: usize = 32;

fn ct_sel8(cond: bool, a: i8, b: i8) -> i8 {
    let m = ct::barrier_i32(-(cond as i32)) as i8;
    (a & m) | (b & !m)
}

fn ct_swap_i8(mask: u64, a: &mut i8, b: &mut i8) {
    let m = ct::barrier(mask) as i8;
    let t = (*a ^ *b) & m;
    *a ^= t;
    *b ^= t;
}

fn ct_cneg_i8(mask: u64, v: i8) -> i8 {
    let m = ct::barrier(mask) as i8;
    (v.wrapping_neg() & m) | (v & !m)
}

fn mk_build_gram<const N: usize>(st: &FpLdl<N>, pp: i32, cfg: &LllConfig) -> [[Ibz<N>; 4]; 4] {
    let nb = cfg.gram_bits;
    let mut dl = [[Ibz::<N>::zero(); 4]; 4];
    for i in 0..4 {
        for k in 0..i {
            dl[k][i] = st.d[k].ct_highmul_p(&st.l[i][k], pp, nb);
        }
        dl[i][i] = st.d[i].with_bound_ct(nb);
    }
    let mut g = [[Ibz::<N>::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            let mut acc = dl[j][i].with_bound_ct(nb);
            for k in 0..j {
                let p = dl[k][i].ct_highmul_p(&st.l[j][k], pp, nb);
                acc = acc.add(&p).with_bound_ct(nb);
            }
            g[i][j] = acc;
        }
    }
    g
}

fn mk_eval_norm<const N: usize>(g: &[[Ibz<N>; 4]; 4], t: &[i8; 4], cfg: &LllConfig) -> Ibz<N> {
    let nb = cfg.gram_bits;
    let mut total = Ibz::<N>::set(0, 2).with_bound_ct(nb);
    for i in 0..4 {
        for j in 0..=i {
            let coeff = (t[i] as i32) * (t[j] as i32) * if i == j { 1 } else { 2 };
            let term = g[i][j].mul_by_int_and_set_bound(coeff, nb);
            total = total.add(&term).with_bound_ct(nb);
        }
    }
    total
}

fn mk_eval_inner<const N: usize>(
    g: &[[Ibz<N>; 4]; 4],
    t: &[i8; 4],
    u: &[i8; 4],
    cfg: &LllConfig,
) -> Ibz<N> {
    let nb = cfg.gram_bits;
    let mut total = Ibz::<N>::set(0, 2).with_bound_ct(nb);
    for i in 0..4 {
        for j in 0..=i {
            let coeff = if i == j {
                (t[i] as i32) * (u[i] as i32)
            } else {
                (t[i] as i32) * (u[j] as i32) + (t[j] as i32) * (u[i] as i32)
            };
            let term = g[i][j].mul_by_int_and_set_bound(coeff, nb);
            total = total.add(&term).with_bound_ct(nb);
        }
    }
    total
}

fn mk_t0_floor<const N: usize>(st: &FpLdl<N>, t123: &[i8; 3], pp: i32, cfg: &LllConfig) -> i32 {
    let work = cfg.l_bits;
    let t0_bits = cfg.l_bits - cfg.p_bits + 16;
    let mut shift = Ibz::<N>::set(0, 2).with_bound_ct(work);
    for j in 1..=3 {
        let term = st.l[j][0].mul_by_int_and_set_bound(t123[j - 1] as i32, work);
        shift = shift.add(&term).with_bound_ct(work);
    }
    let neg_shift = shift.neg().with_bound_ct(work);
    let t0f = neg_shift.ct_shift(-pp, work, t0_bits);
    let t0 = t0f.get();
    debug_assert!(t0 >= i8::MIN as i32 && t0 < i8::MAX as i32);
    t0
}

const MK_T_SUB: [[i8; 3]; 13] = [
    [0, 0, 1],
    [0, 1, -1],
    [0, 1, 0],
    [0, 1, 1],
    [1, -1, -1],
    [1, -1, 0],
    [1, -1, 1],
    [1, 0, -1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, -1],
    [1, 1, 0],
    [1, 1, 1],
];

fn gen_candidates<const N: usize>(
    st: &FpLdl<N>,
    pp: i32,
    cfg: &LllConfig,
) -> (
    [Ibz<N>; MK_N_TOTAL],
    [[i8; 4]; MK_N_TOTAL],
    [[Ibz<N>; 4]; 4],
) {
    let g = mk_build_gram(st, pp, cfg);
    let mut key = [Ibz::<N>::zero(); MK_N_TOTAL];
    let mut tv = [[0i8; 4]; MK_N_TOTAL];
    let mut idx = 0;
    {
        let t = [1i8, 0, 0, 0];
        key[idx] = mk_eval_norm(&g, &t, cfg);
        tv[idx] = t;
        idx += 1;
    }
    for sub in MK_T_SUB.iter() {
        let t0 = mk_t0_floor(st, sub, pp, cfg);
        for which in 0..2 {
            let t = [(t0 + which) as i8, sub[0], sub[1], sub[2]];
            key[idx] = mk_eval_norm(&g, &t, cfg);
            tv[idx] = t;
            idx += 1;
        }
    }
    debug_assert_eq!(idx, MK_N_REAL);
    let kb = cfg.gram_bits;
    let sentinel = Ibz::<N>::set(1, 2)
        .with_bound_ct(kb)
        .mul_2exp((kb - 4) as u32)
        .with_bound_ct(kb);
    for k in key[..idx].iter() {
        debug_assert!(k.bitsize_ct() < kb - 4);
    }
    for i in idx..MK_N_TOTAL {
        tv[i] = [0; 4];
        key[i] = sentinel;
    }
    (key, tv, g)
}

fn cex<const N: usize>(
    key: &mut [Ibz<N>; MK_N_TOTAL],
    tv: &mut [[i8; 4]; MK_N_TOTAL],
    i: usize,
    l: usize,
    ascending: bool,
) {
    let lt = key[l].ct_lt_mask(&key[i]);
    let swap_mask = if ascending { lt } else { !lt };
    debug_assert_eq!(key[i].get_bound(), key[l].get_bound());
    let (a, b) = if i < l {
        let (lo, hi) = key.split_at_mut(l);
        (&mut lo[i], &mut hi[0])
    } else {
        let (lo, hi) = key.split_at_mut(i);
        (&mut hi[0], &mut lo[l])
    };
    Ibz::cswap(a, b, swap_mask);
    for c in 0..4 {
        let (ti, tl) = if i < l {
            let (lo, hi) = tv.split_at_mut(l);
            (&mut lo[i][c], &mut hi[0][c])
        } else {
            let (lo, hi) = tv.split_at_mut(i);
            (&mut hi[0][c], &mut lo[l][c])
        };
        ct_swap_i8(swap_mask, ti, tl);
    }
}

fn bitonic_sort32<const N: usize>(
    key: &mut [Ibz<N>; MK_N_TOTAL],
    tv: &mut [[i8; 4]; MK_N_TOTAL],
    cfg: &LllConfig,
) {
    for k in key.iter_mut() {
        k.set_bound_ct(cfg.gram_bits);
    }
    let mut k = 2;
    while k <= MK_N_TOTAL {
        let mut j = k >> 1;
        while j > 0 {
            for i in 0..MK_N_TOTAL {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    cex(key, tv, i, l, ascending);
                }
            }
            j >>= 1;
        }
        k <<= 1;
    }
}

fn select_basis(sorted_t: &[[i8; 4]; MK_N_TOTAL]) -> [[i8; 4]; 4] {
    let mut tout = [[0i8; 4]; 4];
    let mut slot_count: i32 = 0;
    for cand in sorted_t.iter() {
        let (v0, v1, v2, v3) = (
            cand[0] as i32,
            cand[1] as i32,
            cand[2] as i32,
            cand[3] as i32,
        );
        let (m1_0, m1_1, m1_2, m1_3) = (
            tout[0][0] as i32,
            tout[0][1] as i32,
            tout[0][2] as i32,
            tout[0][3] as i32,
        );
        let t01 = m1_0 * v1 - m1_1 * v0;
        let t02 = m1_0 * v2 - m1_2 * v0;
        let t03 = m1_0 * v3 - m1_3 * v0;
        let t12 = m1_1 * v2 - m1_2 * v1;
        let t13 = m1_1 * v3 - m1_3 * v1;
        let t23 = m1_2 * v3 - m1_3 * v2;
        let test1 = (t01 | t02 | t03 | t12 | t13 | t23) != 0;
        let (m2_0, m2_1, m2_2, m2_3) = (
            tout[1][0] as i32,
            tout[1][1] as i32,
            tout[1][2] as i32,
            tout[1][3] as i32,
        );
        let m01 = m1_0 * m2_1 - m1_1 * m2_0;
        let m02 = m1_0 * m2_2 - m1_2 * m2_0;
        let m03 = m1_0 * m2_3 - m1_3 * m2_0;
        let m12 = m1_1 * m2_2 - m1_2 * m2_1;
        let m13 = m1_1 * m2_3 - m1_3 * m2_1;
        let m23 = m1_2 * m2_3 - m1_3 * m2_2;
        let s012 = v0 * m12 - v1 * m02 + v2 * m01;
        let s013 = v0 * m13 - v1 * m03 + v3 * m01;
        let s023 = v0 * m23 - v2 * m03 + v3 * m02;
        let s123 = v1 * m23 - v2 * m13 + v3 * m12;
        let test2 = (s012 | s013 | s023 | s123) != 0;
        let (m3_0, m3_1, m3_2, m3_3) = (
            tout[2][0] as i32,
            tout[2][1] as i32,
            tout[2][2] as i32,
            tout[2][3] as i32,
        );
        let c012 = m3_0 * m12 - m3_1 * m02 + m3_2 * m01;
        let c013 = m3_0 * m13 - m3_1 * m03 + m3_3 * m01;
        let c023 = m3_0 * m23 - m3_2 * m03 + m3_3 * m02;
        let c123 = m3_1 * m23 - m3_2 * m13 + m3_3 * m12;
        let det4 = v3 * c012 - v2 * c013 + v1 * c023 - v0 * c123;
        let test3 = det4 != 0;
        let c_eq0 = ct::barrier_i32((slot_count == 0) as i32) != 0;
        let c_eq1 = ct::barrier_i32((slot_count == 1) as i32) != 0;
        let c_eq2 = ct::barrier_i32((slot_count == 2) as i32) != 0;
        let c_eq3 = ct::barrier_i32((slot_count == 3) as i32) != 0;
        let cur_test = ct::select32(
            c_eq0,
            1,
            ct::select32(
                c_eq1,
                test1 as i32,
                ct::select32(c_eq2, test2 as i32, ct::select32(c_eq3, test3 as i32, 0)),
            ),
        );
        let room = ct::barrier_i32((slot_count < 4) as i32);
        let accept = room & cur_test;
        for dest in 0..4 {
            let m = ct::barrier_i32(accept & ct::barrier_i32((slot_count == dest as i32) as i32));
            for c in 0..4 {
                tout[dest][c] = ct_sel8(m != 0, cand[c], tout[dest][c]);
            }
        }
        slot_count += accept;
    }
    tout
}

fn mk_invert_transpose_int(t: &[[i8; 4]; 4]) -> [[i64; 4]; 4] {
    let mut u = [[0i64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            u[i][j] = t[i][j] as i64;
        }
    }
    let m01 = u[0][0] * u[1][1] - u[0][1] * u[1][0];
    let m02 = u[0][0] * u[1][2] - u[0][2] * u[1][0];
    let m03 = u[0][0] * u[1][3] - u[0][3] * u[1][0];
    let m12 = u[0][1] * u[1][2] - u[0][2] * u[1][1];
    let m13 = u[0][1] * u[1][3] - u[0][3] * u[1][1];
    let m23 = u[0][2] * u[1][3] - u[0][3] * u[1][2];
    let n01 = u[2][0] * u[3][1] - u[2][1] * u[3][0];
    let n02 = u[2][0] * u[3][2] - u[2][2] * u[3][0];
    let n03 = u[2][0] * u[3][3] - u[2][3] * u[3][0];
    let n12 = u[2][1] * u[3][2] - u[2][2] * u[3][1];
    let n13 = u[2][1] * u[3][3] - u[2][3] * u[3][1];
    let n23 = u[2][2] * u[3][3] - u[2][3] * u[3][2];
    let mut c = [[0i64; 4]; 4];
    c[0][0] = u[1][1] * n23 - u[1][2] * n13 + u[1][3] * n12;
    c[0][1] = -(u[1][0] * n23 - u[1][2] * n03 + u[1][3] * n02);
    c[0][2] = u[1][0] * n13 - u[1][1] * n03 + u[1][3] * n01;
    c[0][3] = -(u[1][0] * n12 - u[1][1] * n02 + u[1][2] * n01);
    c[1][0] = -(u[0][1] * n23 - u[0][2] * n13 + u[0][3] * n12);
    c[1][1] = u[0][0] * n23 - u[0][2] * n03 + u[0][3] * n02;
    c[1][2] = -(u[0][0] * n13 - u[0][1] * n03 + u[0][3] * n01);
    c[1][3] = u[0][0] * n12 - u[0][1] * n02 + u[0][2] * n01;
    c[2][0] = u[3][1] * m23 - u[3][2] * m13 + u[3][3] * m12;
    c[2][1] = -(u[3][0] * m23 - u[3][2] * m03 + u[3][3] * m02);
    c[2][2] = u[3][0] * m13 - u[3][1] * m03 + u[3][3] * m01;
    c[2][3] = -(u[3][0] * m12 - u[3][1] * m02 + u[3][2] * m01);
    c[3][0] = -(u[2][1] * m23 - u[2][2] * m13 + u[2][3] * m12);
    c[3][1] = u[2][0] * m23 - u[2][2] * m03 + u[2][3] * m02;
    c[3][2] = -(u[2][0] * m13 - u[2][1] * m03 + u[2][3] * m01);
    c[3][3] = u[2][0] * m12 - u[2][1] * m02 + u[2][2] * m01;
    let mut det_u = 0i64;
    for j in 0..4 {
        det_u += u[0][j] * c[0][j];
    }
    debug_assert!(det_u == 1 || det_u == -1);
    let mut v = [[0i64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            v[i][j] = c[i][j] * det_u;
        }
    }
    v
}

fn mk_canonicalize_t<const N: usize>(t: &mut [[i8; 4]; 4], g: &[[Ibz<N>; 4]; 4], cfg: &LllConfig) {
    for i in 1..4 {
        let ip = mk_eval_inner(g, &t[0], &t[i], cfg);
        let neg = ip.neg_mask_pub();
        for c in 0..4 {
            t[i][c] = ct_cneg_i8(neg, t[i][c]);
        }
    }
}

fn minkowski_reduce<const N: usize>(st: &FpLdl<N>, pp: i32, cfg: &LllConfig) -> MkResult<N> {
    let (mut key, mut tv, g) = gen_candidates(st, pp, cfg);
    bitonic_sort32(&mut key, &mut tv, cfg);
    let mut t = select_basis(&tv);
    mk_canonicalize_t(&mut t, &g, cfg);
    let v = mk_invert_transpose_int(&t);
    let ab = cfg.ainv_bits;
    let mut ainv_total = [[Ibz::<N>::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut acc = Ibz::<N>::set(0, 2).with_bound_ct(ab);
            for c in 0..4 {
                debug_assert!(v[j][c] >= i32::MIN as i64 && v[j][c] <= i32::MAX as i64);
                let term = st.ainv[i][c].mul_by_int_and_set_bound(v[j][c] as i32, ab);
                acc = acc.add(&term).with_bound_ct(ab);
            }
            ainv_total[i][j] = acc;
        }
    }
    MkResult { t, ainv_total, g }
}

/// Dual-Minkowski-reduce a general ideal lattice in the inert triangular
/// shape (not an `O0`-ideal). Returns `(reduced, gram_diag, ainv)`:
/// `reduced.basis = hnf.basis * ainv` with `det ainv = +-1`, and the
/// approximate dual Gram diagonal in the same reversed order, scaled by
/// `(2N)^2 p`.
pub fn dual_reduce_ideal<const N: usize>(
    hnf: &QuatLattice<N>,
    alg: &QuatAlg<N>,
) -> (QuatLattice<N>, [Ibz<N>; 4], Mat4x4<N>) {
    let cfg = &alg.lll;
    let pp = cfg.p_bits;
    let hb = &hnf.basis.0;
    let two_n = hb[0][0].with_bound_ct(cfg.hnf_bits);
    let p = alg.p.with_bound_ct(cfg.hnf_bits);
    debug_assert!(two_n.is_positive() && !two_n.is_zero());
    debug_assert!(hb[1][1] == two_n);
    debug_assert!(hb[0][0].bitsize() < cfg.hnf_bits);
    debug_assert!(hb[0][2].bitsize() < cfg.hnf_bits);
    debug_assert!(hb[0][3].bitsize() < cfg.hnf_bits);
    debug_assert!(hb[1][2].bitsize() < cfg.hnf_bits);
    debug_assert!(hb[1][3].bitsize() < cfg.hnf_bits);
    let e = hb[0][2].with_bound_ct(cfg.hnf_bits);
    let dq = hb[0][3].with_bound_ct(cfg.hnf_bits);
    let c = hb[1][2].with_bound_ct(cfg.hnf_bits);
    let b = hb[1][3].with_bound_ct(cfg.hnf_bits);
    let mut st = FpLdl::init_dual(&two_n, &b, &c, &dq, &e, &p, pp, cfg);
    st.lll_tour4(pp, cfg);
    let res = minkowski_reduce(&st, pp, cfg);
    // Aout = P * Ainv_total with the column order reversed; P = (2, 3, 0, 1)
    const SIGMA: [usize; 4] = [2, 3, 0, 1];
    let mut aout = [[Ibz::<N>::zero(); 4]; 4];
    for i in 0..4 {
        for cc in 0..4 {
            aout[i][3 - cc] = res.ainv_total[SIGMA[i]][cc].with_bound_ct(cfg.ainv_bits);
        }
    }
    let mut bout = Mat4x4::<N>::zero();
    for r in 0..4 {
        for cc in 0..4 {
            let mut acc = Ibz::<N>::set(0, 2).with_bound_ct(cfg.basis_bits);
            for m in 0..4 {
                let hm = hb[r][m].with_bound_ct(cfg.hnf_bits);
                let term = hm.mul(&aout[m][cc]);
                acc = acc.add(&term).with_bound_ct(cfg.basis_bits);
            }
            bout.0[r][cc] = acc;
        }
    }
    let mut gram_diag = [Ibz::<N>::zero(); 4];
    for i in 0..4 {
        let nrm = mk_eval_norm(&res.g, &res.t[i], cfg).with_bound_ct(cfg.gram_bits);
        gram_diag[3 - i] = nrm.mul_2exp(cfg.gram_shift as u32);
    }
    let reduced = QuatLattice {
        denom: hnf.denom,
        basis: bout,
    };
    (reduced, gram_diag, Mat4x4(aout))
}
