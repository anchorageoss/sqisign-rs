//! Qlapoty (spec Section 4.6, \[DLS26\]): two ideals equivalent to a given
//! `O0`-ideal whose norms sum to a power of two, and the endomorphism
//! `theta` the ideal-to-isogeny translation evaluates.

use super::dim2::{rounded_div, sum_two_squares};
use super::integers::cornacchia_prime;
use super::lll::lehmer::dim2_short_basis;
use super::{Mat2x2, QuatAlg, QuatAlgElem, QuatIdeal, Vec2};
use crate::mp::{DefaultDomain, Ibz, Rng};

/// Iteration counters of the Qlapoty norm equation, for measuring the
/// failure and retry behaviour (spec Section 9.4). Every counter is a total
/// over the run it was passed to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QlapotyStats {
    /// Iterations of the first loop (draws of `k`).
    pub loop_one_iterations: u64,
    /// Iterations of the second loop (values of `lambda` tried).
    pub loop_two_iterations: u64,
    /// Cornacchia targets that passed the primality test in the second loop.
    pub loop_two_primes: u64,
    /// Norm equations that ended without a solution.
    pub failures: u64,
}

/// Make the first column of a 2x2 basis the first minimum (last reduction
/// step).
fn dim2_reduce_last_step<const N: usize>(ipt: &Mat2x2<N>) -> Mat2x2<N> {
    let n0 = sum_two_squares(&ipt.0[0][0], &ipt.0[1][0]);
    let t = ipt.0[0][0]
        .mul(&ipt.0[0][1])
        .add(&ipt.0[1][0].mul(&ipt.0[1][1]));
    let tmp = rounded_div(&t, &n0);
    let mut m = Mat2x2::zero();
    m.0[0][0] = ipt.0[0][0];
    m.0[1][0] = ipt.0[1][0];
    m.0[0][1] = ipt.0[0][1].sub(&tmp.mul(&ipt.0[0][0]));
    m.0[1][1] = ipt.0[1][1].sub(&tmp.mul(&ipt.0[1][0]));
    m
}

/// The 2-dimensional lattice `{(x, y) : a x + b y = 0 mod 4N}` for the
/// generator `gen = (a + b i + j) / 2`, and optionally a target solving
/// `a x + b y = 1 mod 4N`.
fn gen_to_dim2_lattice<const N: usize>(
    gen: &QuatAlgElem<N>,
    n: &Ibz<N>,
    want_target: bool,
) -> (Mat2x2<N>, Option<Vec2<N>>) {
    debug_assert!(gen.denom == Ibz::two());
    let ab = Vec2([gen.coord.0[0], gen.coord.0[1]]);
    let n4 = n.mul_2exp(2);
    debug_assert!(ab.0[0].gcd(&ab.0[1]).gcd(&n4).is_one());
    let gcd = ab.0[0].gcd(&n4);
    let (n4div, r) = n4.div(&gcd);
    debug_assert!(r.is_zero());
    let (c, r) = ab.0[0].div(&gcd);
    debug_assert!(r.is_zero());
    let c = c.invmod(&n4div).expect("invertible");
    let c = c.mul(&ab.0[1]).neg().modulo(&n4div);
    let mut lat = Mat2x2::zero();
    lat.0[0][0] = c;
    lat.0[1][0] = gcd;
    lat.0[0][1] = n4div;
    lat.0[1][1] = Ibz::zero();
    let target = if want_target {
        let (g, u, v) = ab.0[0].xgcd(&ab.0[1]);
        let (g2, k, _l) = g.xgcd(&n4);
        debug_assert!(g2.is_one());
        Some(Vec2([k.mul(&u).modulo(&n4), k.mul(&v).modulo(&n4)]))
    } else {
        None
    };
    (lat, target)
}

/// Constant-time reduction of the 2-dimensional lattice (norm bounded by
/// `4N`).
fn dim2_ct_reduce<const N: usize>(li: &Mat2x2<N>, n: &Ibz<N>) -> Mat2x2<N> {
    let mut m = *li;
    m.0[0].swap(0, 1);
    m.0[1].swap(0, 1);
    let det_bits = n.get_bound() + 2;
    dim2_short_basis(&m, det_bits)
}

/// A short basis of the lattice for `alpha0` with the first vector `st`,
/// and the target vector.
fn get_short_basis<const N: usize>(
    st: &Vec2<N>,
    alpha0: &QuatAlgElem<N>,
    n: &Ibz<N>,
) -> (Mat2x2<N>, Vec2<N>) {
    debug_assert!(alpha0.denom == Ibz::two());
    let n4 = n.mul_2exp(2);
    let omega = *st;
    let (mut basis, target) = gen_to_dim2_lattice(alpha0, n, true);
    let target = target.expect("target");
    let c = basis.0[0][0];
    let gcd = basis.0[1][0];
    let nmod = basis.0[0][1];
    debug_assert!(gcd.mul(&nmod) == n4);
    let mut xy = Mat2x2::zero();
    let (a1, r) = omega.0[1].div(&gcd);
    debug_assert!(r.is_zero());
    xy.0[0][0] = a1;
    let mut a2 = c.mul(&omega.0[1]);
    let tmp = gcd.mul(&omega.0[0]);
    a2 = tmp.sub(&a2);
    let (a2, r) = a2.div(&n4);
    debug_assert!(r.is_zero());
    xy.0[1][0] = a2;
    let (r, x11, x01) = xy.0[0][0].xgcd(&xy.0[1][0]);
    debug_assert!(r.is_one());
    xy.0[1][1] = x11;
    xy.0[0][1] = x01.neg();
    basis = basis.mul(&xy);
    debug_assert!(basis.0[0][0] == omega.0[0] && basis.0[1][0] == omega.0[1]);
    let _ = nmod;
    let n_len = n.get_bound();
    for i in 0..2 {
        for j in 0..2 {
            basis.0[i][j].set_bound(n_len + 10);
        }
    }
    (dim2_reduce_last_step(&basis), target)
}

/// Bounds `(a, b)`, the corrected minimum `omega`, the possibly adjusted
/// generator, and `(f, big_interval)` for the first Qlapoty loop.
fn initialize_bounds<const N: usize>(
    gen: &mut QuatAlgElem<N>,
    n: &Ibz<N>,
    li: &Mat2x2<N>,
    twoe: &Ibz<N>,
    alg: &QuatAlg<N>,
) -> Option<(Vec2<N>, Vec2<N>, [u32; 2])> {
    let mut omega = Vec2([li.0[0][0], li.0[1][0]]);
    let li_norm = sum_two_squares(&omega.0[0], &omega.0[1]);
    let five = Ibz::<N>::set(5, 4);
    let mut b = Vec2::zero();
    b.0[1] = n.mul(&five).div_2exp(1);
    b.0[0] = n.div_2exp(1);
    let mut ba = twoe.mul_2exp(3);
    let tmp = n.mul_2exp(4);
    ba = ba.sub(&tmp).mul(n);
    let (ba, _) = ba.div(&alg.p);
    let r = li_norm.mul(&ba);
    let tmp = b.0[0].mul_2exp(1);
    if r <= tmp {
        return None;
    }
    let mut ab = Vec2::zero();
    let (q, r) = b.0[0].div(&li_norm);
    ab.0[0] = if r.is_zero() { q } else { q.add(&Ibz::one()) };
    let (q, _) = b.0[1].div(&li_norm);
    ab.0[1] = q;
    let tmp = ba.div_2exp(1);
    if ab.0[1] > tmp {
        ab.0[1] = tmp;
    }
    let mut cst = Ibz::<N>::set(n.bitsize() as i64, 32);
    cst = cst.mul(&cst);
    let mut f_and_big = [0u32; 2];
    f_and_big[1] = (ab.0[1].sub(&ab.0[0]) > cst) as u32;
    if f_and_big[1] == 0 {
        let mut gcd = omega.0[0].gcd(&omega.0[1]);
        f_and_big[0] = 0;
        let mut fivef = Ibz::<N>::set(1, 2);
        let cst2 = cst.mul(&cst);
        while gcd.divides(&five) && fivef < cst2 {
            let (q, r) = gcd.div(&five);
            debug_assert!(r.is_zero());
            gcd = q;
            fivef = fivef.mul(&five);
            f_and_big[0] += 1;
        }
        if f_and_big[0] != 0 {
            let (q, r) = omega.0[0].div(&fivef);
            debug_assert!(r.is_zero());
            omega.0[0] = q;
            let (q, r) = omega.0[1].div(&fivef);
            debug_assert!(r.is_zero());
            omega.0[1] = q;
            let r = omega.0[0].mul_2exp(1).sub(&omega.0[1]);
            let mut fivevec = QuatAlgElem::set(1, 2, 1, 0, 0);
            if r.divides(&five) {
                fivevec.coord.0[1] = fivevec.coord.0[1].neg();
            }
            for _ in 1..f_and_big[0] {
                let v = Vec2([fivevec.coord.0[0], fivevec.coord.0[1]]);
                let sq = v.gaussian_mul(&v);
                fivevec.coord.0[0] = sq.0[0];
                fivevec.coord.0[1] = sq.0[1];
            }
            let v = Vec2([fivevec.coord.0[0], fivevec.coord.0[1]]);
            omega = v.gaussian_mul(&omega);
            *gen = fivevec.mul(gen, alg);
        }
        ab.0[0] = fivef.mul(&ab.0[0]);
        ab.0[1] = fivef.mul(&ab.0[1]);
        f_and_big[1] = (ab.0[1].sub(&ab.0[0]) > cst) as u32;
    } else {
        f_and_big[0] = 0;
    }
    Some((ab, omega, f_and_big))
}

/// First loop: find `alpha0` and the short solution `st`.
#[allow(clippy::too_many_arguments)]
fn loop_one<const N: usize>(
    gen: &QuatAlgElem<N>,
    n: &Ibz<N>,
    absq: &Vec2<N>,
    omega: &Vec2<N>,
    f: u32,
    big_interval: bool,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
    stats: &mut QlapotyStats,
) -> Option<(QuatAlgElem<N>, Vec2<N>)> {
    let gcdomega = omega.0[0].gcd(&omega.0[1]);
    let mut onebounds = Vec2::zero();
    let tmp = absq.0[0].div_2exp(1);
    onebounds.0[0] = tmp.sqrt_floor();
    let sum = onebounds.0[0].mul(&onebounds.0[0]);
    if sum != tmp && absq.0[0].is_even() {
        onebounds.0[0] = onebounds.0[0].add(&Ibz::one());
    }
    onebounds.0[1] = absq.0[1].sqrt_floor();
    let mut k = Vec2::zero();
    let mut need_random_sampling = true;
    if onebounds.0[1] <= Ibz::set(11, 5) {
        k = Vec2([Ibz::set(4, 4), Ibz::set(3, 4)]);
        need_random_sampling = false;
    }
    if onebounds.0[1] <= Ibz::set(4, 4) {
        k = Vec2([Ibz::set(2, 3), Ibz::set(1, 3)]);
        need_random_sampling = false;
    }
    if onebounds.0[1] <= Ibz::set(2, 3) {
        k = Vec2([Ibz::set(1, 2), Ibz::set(0, 2)]);
        need_random_sampling = false;
    }
    let mut work;
    loop {
        stats.loop_one_iterations += 1;
        if need_random_sampling {
            debug_assert!(onebounds.0[0] < onebounds.0[1]);
            k.0[0] = Ibz::rand_interval(
                &onebounds.0[0],
                &onebounds.0[1],
                &mut DefaultDomain(&mut *rng),
            )?;
            let ksq = k.0[0].mul(&k.0[0]);
            let mut twobounds = Vec2::zero();
            if absq.0[0] > ksq {
                let tmp = absq.0[0].sub(&ksq);
                twobounds.0[0] = tmp.sqrt_floor();
                let sum = twobounds.0[0].mul(&twobounds.0[0]);
                if sum < tmp {
                    twobounds.0[0] = twobounds.0[0].add(&Ibz::one());
                }
            } else {
                twobounds.0[0] = Ibz::zero();
            }
            let tmp = absq.0[1].sub(&ksq);
            twobounds.0[1] = tmp.sqrt_floor();
            if twobounds.0[0] > twobounds.0[1] {
                continue;
            }
            k.0[1] = Ibz::rand_interval(
                &twobounds.0[0],
                &twobounds.0[1],
                &mut DefaultDomain(&mut *rng),
            )?;
        }
        if k.0[0].is_even() == k.0[1].is_even() {
            continue;
        }
        if !k.0[1].gcd(&k.0[0]).is_one() {
            continue;
        }
        let tmp = sum_two_squares(&k.0[0], &k.0[1]);
        if big_interval {
            if !tmp.gcd(n).is_one() {
                continue;
            }
        } else if !tmp.gcd(&gcdomega).is_one() {
            continue;
        }
        if k.0[0].is_odd() {
            k.0.swap(0, 1);
        }
        if f & 1 == 1 {
            k.0.swap(0, 1);
            k.0[0] = k.0[0].neg();
        }
        // gaussian gcd
        k.0[1] = k.0[1].neg();
        let mut gamma = k.gaussian_gcd(omega);
        k.0[1] = k.0[1].neg();
        let mut ngamma = sum_two_squares(&gamma.0[0], &gamma.0[1]);
        let mut size_k = 2 * k.0[0].get_bound();
        if size_k < 2 * k.0[1].get_bound() {
            size_k = 2 * k.0[1].get_bound();
        }
        while !ngamma.is_one() {
            gamma = gamma.gaussian_mul(&gamma);
            k = k.gaussian_mul(&gamma);
            let (q, r) = k.0[0].div(&ngamma);
            debug_assert!(r.is_zero());
            k.0[0] = q;
            let (q, r) = k.0[1].div(&ngamma);
            debug_assert!(r.is_zero());
            k.0[1] = q;
            k.0[0].set_bound(size_k);
            k.0[1].set_bound(size_k);
            k.0[1] = k.0[1].neg();
            gamma = k.gaussian_gcd(omega);
            k.0[1] = k.0[1].neg();
            ngamma = sum_two_squares(&gamma.0[0], &gamma.0[1]);
        }
        gamma = gamma.gaussian_mul(&gamma);
        k = k.gaussian_mul(&gamma);
        let (q, r) = k.0[0].div(&ngamma);
        debug_assert!(r.is_zero());
        k.0[0] = q;
        let (q, r) = k.0[1].div(&ngamma);
        debug_assert!(r.is_zero());
        k.0[1] = q;
        let kb = (n.get_bound() + alg.p.bitsize()) / 2;
        k.0[0].set_bound(kb);
        k.0[1].set_bound(kb);
        work = QuatAlgElem::set(1, 0, 0, 0, 0);
        work.coord.0[0] = k.0[0];
        work.coord.0[1] = k.0[1];
        work = work.mul(gen, alg);
        debug_assert!(work.denom == Ibz::two());
        let tmp = work.coord.0[0].gcd(&work.coord.0[1]);
        if !n.gcd(&tmp).is_one() {
            continue;
        }
        break;
    }
    let st = k.gaussian_mul(omega);
    Some((work, st))
}

/// Second loop: find `lambda`, `(s, t)` and the Cornacchia solution
/// `(sc, tc)`.
#[allow(clippy::too_many_arguments)]
fn loop_two<const N: usize>(
    alpha0: &QuatAlgElem<N>,
    r: &Ibz<N>,
    st0: &Vec2<N>,
    n: &Ibz<N>,
    twoe: &Ibz<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
    stats: &mut QlapotyStats,
) -> Option<(Ibz<N>, Ibz<N>, Ibz<N>, Ibz<N>, Ibz<N>)> {
    debug_assert!(*n > Ibz::zero());
    let delta_n = Ibz::<N>::set(1 + 2 * (r.is_even() as i64), 3).mul(n);
    let (rmat, mut target0) = get_short_basis(st0, alpha0, n);
    let (rdet, rinv) = rmat.inv_with_det_as_denom();
    let rinv = rinv.expect("invertible");
    target0.0[0] = target0.0[0].neg();
    target0.0[1] = target0.0[1].neg();
    let n_len = n.get_bound();
    let e = twoe.get_bound();
    target0.0[0].set_bound(2 * n_len + 10);
    target0.0[1].set_bound(2 * n_len + 10);
    let mut lambda = Ibz::<N>::set(1, 2).neg();
    let rtwo = r.mul_2exp(1);
    let twoemdn = twoe.sub(&delta_n);
    let quadn = n.mul_2exp(2);
    let (bound, _) = twoe.div(r);
    let bound = bound.sqrt_floor();
    let mut target = Vec2::zero();
    let mut z = Ibz::zero();
    let mut found = None;
    while lambda <= bound {
        stats.loop_two_iterations += 1;
        lambda = lambda.add(&Ibz::two());
        lambda.set_bound(50);
        if !lambda.gcd(&quadn).is_one() {
            continue;
        }
        let tmp = rtwo.mul(&lambda).mul(&lambda);
        let mut sum = twoemdn.sub(&tmp);
        sum.set_bound(e);
        let lambdainv = lambda.invmod(&quadn).expect("coprime");
        let mut tmp = sum.mul(&lambdainv).modulo(&quadn);
        tmp.set_bound(n_len + 10);
        let mut targetl = Vec2([
            target0.0[0].mul(&tmp).neg().modulo(&quadn),
            target0.0[1].mul(&tmp).neg().modulo(&quadn),
        ]);
        target = rinv.eval(&targetl);
        target.0[0] = rounded_div(&target.0[0], &rdet);
        target.0[1] = rounded_div(&target.0[1], &rdet);
        target = rmat.eval(&target);
        target = targetl.sub(&target);
        target.0[0].set_bound(2 * n_len + 10);
        target.0[1].set_bound(2 * n_len + 10);
        targetl = target;
        let _ = targetl;
        if target.0[0].is_even() == target.0[1].is_even() {
            continue;
        }
        let mut k = sum.add(&delta_n);
        let tmp = target.0[0].mul(&alpha0.coord.0[0]);
        let sum2 = target.0[1].mul(&alpha0.coord.0[1]).add(&tmp);
        let mut tmp = sum2.mul(&lambda);
        tmp.set_bound(2 * n_len + 20);
        k = k.sub(&tmp);
        let (kq, kr) = k.div(n);
        debug_assert!(kr.is_zero());
        k = kq;
        k.set_bound(e);
        debug_assert!(k.is_odd());
        z = k.add(&k);
        z = z.sub(&target.0[1].mul(&target.0[1]));
        z = z.sub(&target.0[0].mul(&target.0[0]));
        z.set_bound(2 * n_len + 10);
        if z < Ibz::zero() {
            continue;
        }
        debug_assert!((z.get() & 3) == 1);
        if (z.get() & 3) == 1 && z.probab_prime(alg.primality_num_iter, rng) {
            stats.loop_two_primes += 1;
            found = cornacchia_prime(&z);
        } else {
            found = None;
        }
        if found.is_some() {
            break;
        }
    }
    let (sc, tc) = found?;
    let mut s = target.0[0];
    let mut t = target.0[1];
    s.set_bound(n_len + 10);
    t.set_bound(n_len + 10);
    let _ = z;
    Some((lambda, s, t, sc, tc))
}

/// QlapotyNormEquation: `(mu1, mu2, theta, smallest)` with
/// `I_i = (ideal * conj(smallest) / n(ideal)) * conj(mu_i) / n(small)` two
/// ideals whose norms sum to `2^e`, and `theta = mu2 conj(mu1) / n(small)`.
/// `None` on failure (an equivalent ideal of very small norm).
pub fn normeq<const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
) -> Option<(
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
)> {
    normeq_with_stats(ideal, alg, rng, &mut QlapotyStats::default())
}

/// [`normeq`] with iteration counters added to `stats`.
pub fn normeq_with_stats<const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
    stats: &mut QlapotyStats,
) -> Option<(
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
)> {
    let r = normeq_inner(ideal, alg, rng, stats);
    if r.is_none() {
        stats.failures += 1;
    }
    r
}

fn normeq_inner<const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
    stats: &mut QlapotyStats,
) -> Option<(
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
    QuatAlgElem<N>,
)> {
    let (_gen, smallest, mut small) = ideal.shortest_equivalent(alg, rng);
    small.norm.set_bound(alg.p.bitsize() / 2 + 10);
    let mut gen = small.odd_inert_gen();
    gen.coord.0[0].set_bound(small.norm.get_bound());
    gen.coord.0[1].set_bound(small.norm.get_bound());
    gen.coord.0[2].set_bound(2);
    gen.coord.0[3].set_bound(1);
    gen.denom.set_bound(3);
    let (li, _) = gen_to_dim2_lattice(&gen, &small.norm, false);
    let li = dim2_ct_reduce(&li, &small.norm);
    let twoe = Ibz::<N>::one().mul_2exp(alg.qlapoty_used_power_of_two);
    let (ab, omega, f_and_big) = initialize_bounds(&mut gen, &small.norm, &li, &twoe, alg)?;
    #[cfg(debug_assertions)]
    {
        let sum1 = omega.0[0]
            .mul(&gen.coord.0[0])
            .add(&omega.0[1].mul(&gen.coord.0[1]));
        debug_assert!(sum1.modulo(&small.norm.mul_2exp(2)).is_zero());
    }
    let (alpha0, st0) = loop_one(
        &gen,
        &small.norm,
        &ab,
        &omega,
        f_and_big[0],
        f_and_big[1] != 0,
        alg,
        rng,
        stats,
    )?;
    let (r, t) = alpha0.norm(alg);
    debug_assert!(t.is_one());
    let (r, t) = r.div(&small.norm);
    debug_assert!(t.is_zero());
    let (lambda, s, t, mut sc, mut tc) =
        loop_two(&alpha0, &r, &st0, &small.norm, &twoe, alg, rng, stats)?;
    if s.is_even() != sc.is_even() {
        core::mem::swap(&mut sc, &mut tc);
    }
    debug_assert!(s.is_even() == sc.is_even() && t.is_even() == tc.is_even());
    let mut mu1 = QuatAlgElem::set(2, 0, 0, 0, 0);
    let mut mu2 = QuatAlgElem::set(2, 0, 0, 0, 0);
    mu1.coord.0[0] = s.add(&sc);
    mu1.coord.0[1] = t.add(&tc);
    mu2.coord.0[0] = s.sub(&sc);
    mu2.coord.0[1] = t.sub(&tc);
    let temp = alpha0.scalar_mul(&lambda);
    mu1 = mu1.scalar_mul(&small.norm);
    mu2 = mu2.scalar_mul(&small.norm);
    mu1 = mu1.add(&temp);
    mu2 = mu2.add(&temp);
    #[cfg(debug_assertions)]
    {
        let (n1, d1) = mu1.norm(alg);
        let (n2, d2) = mu2.norm(alg);
        debug_assert!(d1.is_one() && d2.is_one());
        let (q, rr) = n1.add(&n2).div(&small.norm);
        debug_assert!(rr.is_zero() && q == twoe);
    }
    let mut temp = mu1.conj();
    temp.denom = temp.denom.mul(&small.norm);
    let theta = mu2.mul(&temp, alg);
    Some((mu1, mu2, theta, smallest))
}

/// Qlapoty: `(beta1, d1, theta)` with `beta1` in the ideal,
/// `d1 = n(beta1) / n(ideal)` and `theta = conj(beta2 conj(beta1) / n)`
/// such that `n(beta1) + n(beta2) = n(ideal) 2^e`. `None` on failure.
pub fn qlapoty<const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
) -> Option<(QuatAlgElem<N>, Ibz<N>, QuatAlgElem<N>)> {
    qlapoty_with_stats(ideal, alg, rng, &mut QlapotyStats::default())
}

/// [`qlapoty`] with iteration counters added to `stats`.
pub fn qlapoty_with_stats<const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
    stats: &mut QlapotyStats,
) -> Option<(QuatAlgElem<N>, Ibz<N>, QuatAlgElem<N>)> {
    let (mu1, _mu2, mut theta, small) = normeq_with_stats(ideal, alg, rng, stats)?;
    theta.normalize();
    theta = theta.conj();
    let (n, d) = small.norm(alg);
    debug_assert!(d.is_one());
    let (mut n, d) = n.div(&ideal.norm);
    debug_assert!(d.is_zero());
    n.set_bound(alg.p.bitsize() / 2 + 1);
    let (d1, d) = mu1.norm(alg);
    debug_assert!(d.is_one());
    let (mut d1, d) = d1.div(&n);
    debug_assert!(d.is_zero());
    d1.set_bound(alg.qlapoty_used_power_of_two as i32 + 2);
    let mut beta1 = mu1.mul(&small, alg);
    beta1.denom = beta1.denom.mul(&n);
    beta1.normalize();
    for c in beta1.coord.0.iter_mut() {
        c.set_bound((ideal.norm.get_bound() + d1.get_bound()) / 2 + 3);
    }
    beta1.denom.set_bound(3);
    let tb = 2 * (beta1.coord.0[0].get_bound() + 1);
    for c in theta.coord.0.iter_mut() {
        c.set_bound(tb);
    }
    theta.denom.set_bound(4);
    Some((beta1, d1, theta))
}
