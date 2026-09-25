//! Pairings on the Kummer line by the cubical (biextension) ladder, and the
//! pairing-based two-dimensional discrete logarithm in `E[2^e]`.
//!
//! The reduced Tate pairing and the Tate-based discrete logarithm follow the
//! round-3 reference (`biextension.c`). The Weil pairing, which round 3
//! dropped, is kept because PRISM's salt-PRISM verifier checks a Weil pairing
//! identity; it matches the round-2 reference.
//!
//! Cubical addition is off by a factor of 4 from true cubical arithmetic,
//! which the final exponentiations over `Fp2` absorb.

use super::{mp, EcBasis, EcCurve, EcPoint, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};
use crate::precomp::PrimePrecomp;

/// Normalised inputs for one pairing.
struct PairingParams<L: FpBackend> {
    e: u32,
    p: EcPoint<L>,
    q: EcPoint<L>,
    pq: EcPoint<L>,
    ix_p: Fp2<L>,
    ix_q: Fp2<L>,
    a24: EcPoint<L>,
}

/// Cubical `P + Q` from cubical `P`, `Q` and `1/x(P - Q)` (the difference
/// anti-normalised as `(1 : ix)`).
#[inline]
fn cubical_add<L: FpBackend>(p: &EcPoint<L>, q: &EcPoint<L>, ix_pq: &Fp2<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let t2 = q.x.add(&q.z);
    let t3 = q.x.sub(&q.z);
    let t0 = t0.mul(&t3);
    let t1 = t1.mul(&t2);
    let t2 = t0.add(&t1);
    let t3 = t0.sub(&t1);
    EcPoint::new(ix_pq.mul(&t2.sqr()), t3.sqr())
}

/// Cubical `(P + Q, [2] Q)`; `A24` must be normalised.
#[inline]
fn cubical_dbladd<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    ix_pq: &Fp2<L>,
    a24: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let ppq_x = q.x.add(&q.z);
    let t3 = q.x.sub(&q.z);
    let t2 = ppq_x.sqr();
    let qq_z = t3.sqr();
    let t0 = t0.mul(&t3);
    let t1 = t1.mul(&ppq_x);
    let ppq_x = t0.add(&t1);
    let t3 = t0.sub(&t1);
    let ppq_z = t3.sqr();
    let ppq_x = ix_pq.mul(&ppq_x.sqr());
    let t3 = t2.sub(&qq_z);
    let qq_x = t2.mul(&qq_z);
    let t0 = t3.mul(&a24.x).add(&qq_z);
    let qq_z = t0.mul(&t3);
    (EcPoint::new(ppq_x, ppq_z), EcPoint::new(qq_x, qq_z))
}

/// `(P + [2^e] Q, [2^e] Q)` cubically.
#[inline]
fn biext_ladder_2e<L: FpBackend>(
    e: u32,
    pq: &EcPoint<L>,
    q: &EcPoint<L>,
    ix_p: &Fp2<L>,
    a24: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    let mut pnq = pq.clone();
    let mut nq = q.clone();
    for _ in 0..e {
        let (a, b) = cubical_dbladd(&pnq, &nq, ix_p, a24);
        pnq = a;
        nq = b;
    }
    (pnq, nq)
}

/// The monodromy ratio as a projective point, avoiding a division.
#[inline]
fn point_ratio<L: FpBackend>(pnq: &EcPoint<L>, nq: &EcPoint<L>, p: &EcPoint<L>) -> EcPoint<L> {
    EcPoint::new(nq.x.mul(&p.x), pnq.x.clone())
}

/// Cubical translation of `P` by the 2-torsion point `T`, constant time over
/// the three shapes `T = (A : 0)`, `(0 : B)`, `(A : B)`.
#[inline]
fn translate<L: FpBackend>(p: &mut EcPoint<L>, t: &EcPoint<L>) {
    let px = t.x.mul(&p.x).sub(&t.z.mul(&p.z));
    let pz = t.z.mul(&p.x).sub(&t.x.mul(&p.z));
    let ta_zero = t.x.ct_is_zero();
    let px = Fp2::select(&px, &p.z, ta_zero);
    let pz = Fp2::select(&pz, &p.x, ta_zero);
    let tb_zero = t.z.ct_is_zero();
    let px = Fp2::select(&px, &p.x, tb_zero);
    let pz = Fp2::select(&pz, &p.z, tb_zero);
    p.x = px;
    p.z = pz;
}

/// `g_{P,Q}^{2^e}` via the cubical arithmetic of `P + [2^e] Q`; with
/// `swap_pq` the roles of `P` and `Q` are exchanged.
#[inline]
fn monodromy_i<L: FpBackend>(params: &PairingParams<L>, swap_pq: bool) -> EcPoint<L> {
    let (p, q, ix_p) = if swap_pq {
        (&params.q, &params.p, &params.ix_q)
    } else {
        (&params.p, &params.q, &params.ix_p)
    };
    let (mut pnq, mut nq) = biext_ladder_2e(params.e - 1, &params.pq, q, ix_p, &params.a24);
    translate(&mut pnq, &nq);
    let nq_copy = nq.clone();
    translate(&mut nq, &nq_copy);
    point_ratio(&pnq, &nq, p)
}

/// Normalise `P`, `Q` to `(x : 1)` and compute `1/x(P)`, `1/x(Q)`.
#[inline]
fn cubical_normalization<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>, Fp2<L>, Fp2<L>) {
    let mut t = [p.x.clone(), p.z.clone(), q.x.clone(), q.z.clone()];
    let mut s1: [Fp2<L>; 4] = core::array::from_fn(|_| Fp2::zero());
    let mut s2: [Fp2<L>; 4] = core::array::from_fn(|_| Fp2::zero());
    Fp2::batched_inv(&mut t, &mut s1, &mut s2);
    let ix_p = p.z.mul(&t[0]);
    let ix_q = q.z.mul(&t[2]);
    let np = EcPoint::from_x(p.x.mul(&t[1]));
    let nq = EcPoint::from_x(q.x.mul(&t[3]));
    (np, nq, ix_p, ix_q)
}

fn params<L: FpBackend>(
    e: u32,
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    curve: &mut EcCurve<L>,
) -> PairingParams<L> {
    let (np, nq, ix_p, ix_q) = cubical_normalization(p, q);
    curve.normalize_a24();
    PairingParams {
        e,
        p: np,
        q: nq,
        pq: pq.clone(),
        ix_p,
        ix_q,
        a24: curve.a24.clone(),
    }
}

/// Weil pairing `e_{2^e}(P, Q)` for `P`, `Q` of order `2^e` with `pq = x(P - Q)`.
/// Neither point may be the identity.
#[inline]
pub fn weil<L: FpBackend>(
    e: u32,
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    curve: &mut EcCurve<L>,
) -> Fp2<L> {
    let params = params(e, p, q, pq, curve);
    let r0 = monodromy_i(&params, true);
    let r1 = monodromy_i(&params, false);
    let r = r0.x.mul(&r1.z).inv();
    r.mul(&r0.z).mul(&r1.x)
}

/// `a^((p + 1) / 2^f)`; the cofactor is public.
#[inline]
pub fn clear_cofac<L: FpBackend + PrimePrecomp>(a: &Fp2<L>) -> Fp2<L> {
    a.pow_vartime(&[L::COFACTOR])
}

/// Reduced Tate pairing `t_{2^e}(P, Q)` with `pq = x(P - Q)`: the unreduced
/// value raised to `(p^2 - 1) / 2^e`, split as Frobenius, inversion, cofactor
/// and `2^(f - e)` squarings.
#[inline]
pub fn reduced_tate<L: FpBackend + PrimePrecomp>(
    e: u32,
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    curve: &mut EcCurve<L>,
) -> Fp2<L> {
    let e_diff = L::TWO_ADIC_EXPONENT - e;
    let params = params(e, p, q, pq, curve);
    let r = monodromy_i(&params, true);
    // -(R.Z / R.X)^((p^2 - 1) / 2^f), with ^(p-1) as Frobenius and inverse
    let rx = r.z.mul(&r.x.frob());
    let rz = r.x.mul(&r.z.frob());
    let mut out = rx.inv().mul(&rz);
    out = clear_cofac::<L>(&out);
    for _ in 0..e_diff {
        out = out.sqr();
    }
    out
}

/// Recursive 2-power discrete logarithm: `a` with `f = g^a` given stacks of
/// `f` and `1/g`. `None` if `f` is not in the subgroup generated by `g`.
#[inline]
fn fp2_dlog_2e_rec<L: FpBackend>(
    a: &mut [u64],
    len: usize,
    pows_f: &mut [Fp2<L>],
    pows_g: &mut [Fp2<L>],
    stacklen: usize,
) -> Option<()> {
    if len == 0 {
        a.iter_mut().for_each(|w| *w = 0);
        return Some(());
    }
    if len == 1 {
        let f_is_one = pows_f[stacklen - 1].ct_is_one();
        let f_equals_g = pows_f[stacklen - 1].ct_equal(&pows_g[stacklen - 1]);
        a.iter_mut().for_each(|w| *w = 0);
        a[0] = (!f_is_one).unwrap_u8() as u64; // bit = 1 iff f != 1
        for i in 0..stacklen - 1 {
            let fg = pows_f[i].mul(&pows_g[i]);
            pows_f[i] = Fp2::select(&pows_f[i], &fg, !f_is_one);
            pows_g[i] = pows_g[i].sqr();
        }
        return if bool::from(f_is_one | f_equals_g) {
            Some(())
        } else {
            None
        };
    }
    let right = len / 2;
    let left = len - right;
    pows_f[stacklen] = pows_f[stacklen - 1].clone();
    pows_g[stacklen] = pows_g[stacklen - 1].clone();
    for _ in 0..left {
        pows_f[stacklen] = pows_f[stacklen].sqr();
        pows_g[stacklen] = pows_g[stacklen].sqr();
    }
    let n = a.len();
    let mut dlp1 = [0u64; MAX_ORDER_WORDS];
    let mut dlp2 = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e_rec(&mut dlp1[..n], right, pows_f, pows_g, stacklen + 1)?;
    fp2_dlog_2e_rec(&mut dlp2[..n], left, pows_f, pows_g, stacklen)?;
    mp::shl(&mut dlp2[..n], right);
    mp::add(a, &dlp2[..n], &dlp1[..n]);
    Some(())
}

/// `scal` with `f = g^scal` in the `2^e`-subgroup of `Fp2^*`, given `1/g`.
#[inline]
pub fn fp2_dlog_2e<L: FpBackend>(
    scal: &mut [u64],
    f: &Fp2<L>,
    g_inverse: &Fp2<L>,
    e: u32,
) -> Option<()> {
    const MAX_STACK: usize = 16;
    let mut log = 1usize;
    let mut len = e as usize;
    while len > 1 {
        len >>= 1;
        log += 1;
    }
    debug_assert!(log <= MAX_STACK);
    let mut pows_f: [Fp2<L>; MAX_STACK] = core::array::from_fn(|_| Fp2::zero());
    let mut pows_g: [Fp2<L>; MAX_STACK] = core::array::from_fn(|_| Fp2::zero());
    pows_f[0] = f.clone();
    pows_g[0] = g_inverse.clone();
    scal.iter_mut().for_each(|w| *w = 0);
    fp2_dlog_2e_rec(scal, e as usize, &mut pows_f, &mut pows_g, 1)
}

/// The biquadratic coefficients `(k00, k01, k11)` of a pair of points, whose
/// roots are `x(P +- Q)` (spec Algorithm 4.91).
#[inline]
fn compute_kappa<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    e: &EcCurve<L>,
) -> (Fp2<L>, Fp2<L>, Fp2<L>) {
    let xx = p.x.mul(&q.x);
    let xz = p.x.mul(&q.z);
    let zx = p.z.mul(&q.x);
    let zz = p.z.mul(&q.z);
    let k00 = xx.sub(&zz).sqr();
    let k11 = xz.sub(&zx).sqr();
    let k01 = xx.mul(&zz).mul(&e.a);
    let k01 = k01.add(&k01);
    let cross = xx.add(&zz).mul(&xz.add(&zx)).mul(&e.c);
    let k01 = k01.add(&cross);
    (k00, k01.add(&k01), k11)
}

/// `x(P1 - P2)` from `P1`, `P2` and their translates by a common `Q`, as the
/// shared root of the two biquadratics (spec Algorithm 4.93).
#[inline]
fn shared_difference_point<L: FpBackend>(
    p1: &EcPoint<L>,
    p2: &EcPoint<L>,
    p1q: &EcPoint<L>,
    p2q: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> EcPoint<L> {
    let (k00, k01, k11) = compute_kappa(p1, p2, curve);
    let (k00b, k01b, k11b) = compute_kappa(p1q, p2q, curve);
    let x = k01b.mul(&k00).sub(&k01.mul(&k00b));
    let z = k11b.mul(&k00).sub(&k11.mul(&k00b)).mul(&curve.c);
    EcPoint::new(x, z)
}

/// Inputs of the Tate discrete logarithm: two bases, the four cross
/// differences, inverse x-coordinates and `A24`.
struct DlogParams<L: FpBackend> {
    e: u32,
    pq: EcBasis<L>,
    rs: EcBasis<L>,
    pm_r: EcPoint<L>,
    pm_s: EcPoint<L>,
    rm_q: EcPoint<L>,
    sm_q: EcPoint<L>,
    ix_p: Fp2<L>,
    ix_q: Fp2<L>,
    a24: EcPoint<L>,
}

/// `x(P - R)`, `x(P - S)`, `x(R - Q)`, `x(S - Q)` for bases `(P, Q)` and `(R, S)`.
#[inline]
fn compute_shared_difference_points<L: FpBackend>(d: &mut DlogParams<L>, curve: &EcCurve<L>) {
    let (k00, k01, k11) = compute_kappa(&d.pq.p, &d.rs.p, curve);
    // (X : Z) solves k11 X^2 - k01 X Z + k00 Z^2 = 0
    let four_k = k00.mul(&k11);
    let four_k = four_k.add(&four_k);
    let four_k = four_k.add(&four_k);
    let mut disc = k01.sqr().sub(&four_k);
    let _ = disc.sqrt_verify();
    d.pm_r = EcPoint::new(k01.add(&disc), k11.add(&k11));
    d.rm_q = shared_difference_point(&d.rs.p, &d.pq.q, &d.pm_r, &d.pq.pmq, curve);
    d.pm_s = shared_difference_point(&d.pq.p, &d.rs.q, &d.pm_r, &d.rs.pmq, curve);
    d.sm_q = shared_difference_point(&d.rs.q, &d.pq.q, &d.pm_s, &d.pq.pmq, curve);
}

/// Normalise every point to `(x : 1)`, the curve to `(A : 1)`, and store
/// `1/x(P)`, `1/x(Q)`.
#[inline]
fn cubical_normalization_dlog<L: FpBackend>(d: &mut DlogParams<L>, curve: &mut EcCurve<L>) {
    let mut t = [
        curve.c.clone(),
        d.pq.p.x.clone(),
        d.pq.q.x.clone(),
        d.pq.p.z.clone(),
        d.pq.q.z.clone(),
        d.rs.p.z.clone(),
        d.rs.q.z.clone(),
        d.pq.pmq.z.clone(),
        d.pm_r.z.clone(),
        d.pm_s.z.clone(),
        d.rm_q.z.clone(),
        d.sm_q.z.clone(),
    ];
    let mut s1: [Fp2<L>; 12] = core::array::from_fn(|_| Fp2::zero());
    let mut s2: [Fp2<L>; 12] = core::array::from_fn(|_| Fp2::zero());
    Fp2::batched_inv(&mut t, &mut s1, &mut s2);
    curve.a = curve.a.mul(&t[0]);
    curve.c = Fp2::one();
    d.ix_p = d.pq.p.z.mul(&t[1]);
    d.ix_q = d.pq.q.z.mul(&t[2]);
    let norm = |pt: &mut EcPoint<L>, inv: &Fp2<L>| {
        pt.x = pt.x.mul(inv);
        pt.z = Fp2::one();
    };
    norm(&mut d.pq.p, &t[3]);
    norm(&mut d.pq.q, &t[4]);
    norm(&mut d.rs.p, &t[5]);
    norm(&mut d.rs.q, &t[6]);
    norm(&mut d.pq.pmq, &t[7]);
    norm(&mut d.pm_r, &t[8]);
    norm(&mut d.pm_s, &t[9]);
    norm(&mut d.rm_q, &t[10]);
    norm(&mut d.sm_q, &t[11]);
}

/// The Tate pairings `t(P, Q)`, `t(R, P)`, `t(R, Q)`, `t(S, P)`, `t(S, Q)` and
/// the four discrete logarithms.
#[inline]
fn tate_dlog_partial<L: FpBackend + PrimePrecomp>(
    d: &DlogParams<L>,
) -> Option<[[u64; MAX_ORDER_WORDS]; 4]> {
    let e_full = L::TWO_ADIC_EXPONENT;
    let e_diff = e_full - d.e;
    let n = L::ORDER_WORDS;

    let mut np = d.pq.p.clone();
    let mut nr = d.rs.p.clone();
    let mut ns = d.rs.q.clone();
    let mut npq = d.pq.pmq.clone();
    let mut pnr = d.pm_r.clone();
    let mut pns = d.pm_s.clone();
    let mut nrq = d.rm_q.clone();
    let mut nsq = d.sm_q.clone();

    for _ in 0..e_full - 1 {
        let (a, b) = cubical_dbladd(&npq, &np, &d.ix_q, &d.a24);
        npq = a;
        np = b;
    }
    for _ in 0..d.e - 1 {
        pnr = cubical_add(&pnr, &nr, &d.ix_p);
        let (a, b) = cubical_dbladd(&nrq, &nr, &d.ix_q, &d.a24);
        nrq = a;
        nr = b;
        pns = cubical_add(&pns, &ns, &d.ix_p);
        let (a, b) = cubical_dbladd(&nsq, &ns, &d.ix_q, &d.a24);
        nsq = a;
        ns = b;
    }
    translate(&mut npq, &np);
    translate(&mut pnr, &nr);
    translate(&mut nrq, &nr);
    translate(&mut pns, &ns);
    translate(&mut nsq, &ns);
    let (np_c, nr_c, ns_c) = (np.clone(), nr.clone(), ns.clone());
    translate(&mut np, &np_c);
    translate(&mut nr, &nr_c);
    translate(&mut ns, &ns_c);

    let mut w1: [Fp2<L>; 5] = core::array::from_fn(|_| Fp2::zero());
    let mut w2: [Fp2<L>; 5] = core::array::from_fn(|_| Fp2::zero());
    // t(P, Q)^(2^e_diff) = w0
    let t = point_ratio(&npq, &np, &d.pq.q);
    w1[0] = t.x;
    w2[0] = t.z;
    // t(R, P) = w0^r2
    let t = point_ratio(&pnr, &nr, &d.pq.p);
    w1[1] = t.x;
    w2[1] = t.z;
    // t(R, Q) = w0^r1 (inverted)
    let t = point_ratio(&nrq, &nr, &d.pq.q);
    w2[2] = t.x;
    w1[2] = t.z;
    // t(S, P) = w0^s2
    let t = point_ratio(&pns, &ns, &d.pq.p);
    w1[3] = t.x;
    w2[3] = t.z;
    // t(S, Q) = w0^s1 (inverted)
    let t = point_ratio(&nsq, &ns, &d.pq.q);
    w2[4] = t.x;
    w1[4] = t.z;

    // ^(p - 1) as Frobenius and inverse, projectively
    for i in 0..5 {
        let tmp = w1[i].clone();
        w1[i] = w2[i].mul(&w1[i].frob());
        w2[i] = tmp.mul(&w2[i].frob());
    }
    let mut s1: [Fp2<L>; 5] = core::array::from_fn(|_| Fp2::zero());
    let mut s2: [Fp2<L>; 5] = core::array::from_fn(|_| Fp2::zero());
    Fp2::batched_inv(&mut w2, &mut s1, &mut s2);
    for i in 0..5 {
        w1[i] = w1[i].mul(&w2[i]);
        w1[i] = clear_cofac::<L>(&w1[i]);
        for _ in 0..e_diff {
            w1[i] = w1[i].sqr();
        }
    }

    let mut out = [[0u64; MAX_ORDER_WORDS]; 4];
    fp2_dlog_2e(&mut out[1][..n], &w1[1], &w1[0], d.e)?; // r2
    fp2_dlog_2e(&mut out[0][..n], &w1[2], &w1[0], d.e)?; // r1
    fp2_dlog_2e(&mut out[3][..n], &w1[3], &w1[0], d.e)?; // s2
    fp2_dlog_2e(&mut out[2][..n], &w1[4], &w1[0], d.e)?; // s1
    Some(out)
}

/// Two-dimensional discrete logarithm by Tate pairings. `pq` is a basis of
/// the full `E[2^f]`, `rs` a basis of `E[2^e]`, `e <= f`. Returns
/// `(r1, r2, s1, s2)` with `R = [2^(f-e)] ([r1] P + [r2] Q)` and
/// `S = [2^(f-e)] ([s1] P + [s2] Q)`, each of `e` bits, **up to a common
/// sign**: x-only points determine `(R, S)` only up to `(-R, -S)` (the
/// difference `R - S` is kept), and which of the two answers comes out is
/// set by the canonical square root in the shared-difference computation
/// (the reference primes give one sign, the toy prime the other). Callers
/// build matrices from the result, for which `M` and `-M` give the same
/// Kummer points. `None` if `R` or `S` is not in the span.
#[inline]
pub fn ec_dlog_2_tate<L: FpBackend + PrimePrecomp>(
    pq: &EcBasis<L>,
    rs: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    e: u32,
) -> Option<[[u64; MAX_ORDER_WORDS]; 4]> {
    curve.normalize_a24();
    let mut d = DlogParams {
        e,
        pq: pq.clone(),
        rs: rs.clone(),
        pm_r: EcPoint::identity(),
        pm_s: EcPoint::identity(),
        rm_q: EcPoint::identity(),
        sm_q: EcPoint::identity(),
        ix_p: Fp2::zero(),
        ix_q: Fp2::zero(),
        a24: curve.a24.clone(),
    };
    compute_shared_difference_points(&mut d, curve);
    cubical_normalization_dlog(&mut d, curve);
    tate_dlog_partial(&d)
}

/// Bits of the largest 2-power subgroup a [`DlogTable`] covers: the torsion
/// `e' + 2` of every parameter set is at most 324.
pub const MAX_DLOG_BITS: usize = 324;

/// The fixed-base data of the two-adic discrete logarithm (P25): for `g`
/// of order exactly `2^t` in `F_p^2`, the powers `(g^-1)^(2^k)`, `k = 0..t`,
/// and the sixteen powers of `(g^-1)^(2^(t-4))`, the 16th roots of unity.
/// With them [`fp2_dlog_2e_fixed`] divides a found digit out of a node by
/// table products instead of the running ancestor updates of
/// [`fp2_dlog_2e`], and resolves four bits per leaf by lookup: about a
/// third of that routine's multiplications and none of its squarings of
/// the base. Built once per public key by the verifier's prepared key, or
/// once per decoding (`t` squarings and one inversion) otherwise.
pub struct DlogTable<L: FpBackend> {
    pows: [Fp2<L>; MAX_DLOG_BITS],
    roots: [Fp2<L>; 16],
    t: u32,
}

impl<L: FpBackend> DlogTable<L> {
    /// The table of `g`. `None` if `t` is outside `4..=MAX_DLOG_BITS` or if
    /// `g` does not have order exactly `2^t` (`g^(2^(t-1)) != -1`), which is
    /// how a pairing of a basis of the wrong order shows up.
    pub fn new(g: &Fp2<L>, t: u32) -> Option<Self> {
        if !(4..=MAX_DLOG_BITS as u32).contains(&t) {
            return None;
        }
        let mut pows: [Fp2<L>; MAX_DLOG_BITS] = core::array::from_fn(|_| Fp2::zero());
        pows[0] = g.inv();
        for k in 1..t as usize {
            pows[k] = pows[k - 1].sqr();
        }
        // order exactly 2^t: (g^-1)^(2^(t-1)) = -1
        if !bool::from(pows[t as usize - 1].ct_equal(&Fp2::one().neg())) {
            return None;
        }
        let zeta_inv = pows[t as usize - 4].clone();
        let mut roots: [Fp2<L>; 16] = core::array::from_fn(|_| Fp2::zero());
        roots[0] = Fp2::one();
        for j in 1..16 {
            roots[j] = roots[j - 1].mul(&zeta_inv);
        }
        Some(Self { pows, roots, t })
    }

    /// `t`.
    pub fn bits(&self) -> u32 {
        self.t
    }

    /// The digit `x < 2^len` (`len <= 4`) with `f = zeta_(2^len)^x`, where
    /// `zeta_(2^len) = g^(2^(t-len))`; `f = (zeta_16^-1)^(-x 2^(4-len))`, so
    /// the lookup index is `-x 2^(4-len) mod 16`. Variable time: the
    /// pairing values are public.
    fn leaf(&self, f: &Fp2<L>, len: usize) -> Option<u64> {
        let shift = 4 - len;
        (0..1u64 << len).find(|&x| {
            let idx = (16 - ((x << shift) & 15)) & 15;
            bool::from(f.ct_equal(&self.roots[idx as usize]))
        })
    }

    /// The bits of `x` (`len` of them, `x < 2^len`) with `f = (g^(2^k))^x`,
    /// `k + len = t`, written into `out` (little-endian limbs).
    fn rec(&self, f: &Fp2<L>, k: usize, len: usize, out: &mut [u64]) -> Option<()> {
        out.iter_mut().for_each(|w| *w = 0);
        if len <= 4 {
            out[0] = self.leaf(f, len)?;
            return Some(());
        }
        let right = len / 2;
        let left = len - right;
        // the low `right` bits from f^(2^left) = (g^(2^(k+left)))^(x_lo)
        let mut fa = f.clone();
        for _ in 0..left {
            fa = fa.sqr();
        }
        let n = out.len();
        let mut lo = [0u64; MAX_ORDER_WORDS];
        self.rec(&fa, k + left, right, &mut lo[..n])?;
        // divide x_lo out: f (g^-1)^(2^k x_lo) = (g^(2^(k+right)))^(x_hi)
        let mut fb = f.clone();
        for b in 0..right {
            if mp::bit(&lo[..n], b) == 1 {
                fb = fb.mul(&self.pows[k + b]);
            }
        }
        let mut hi = [0u64; MAX_ORDER_WORDS];
        self.rec(&fb, k + right, left, &mut hi[..n])?;
        mp::shl(&mut hi[..n], right);
        mp::add(out, &lo[..n], &hi[..n]);
        Some(())
    }
}

/// `scal` with `f = g^scal` for the `g` of `table`, `scal < 2^t`. `None` if
/// `f` is not in the subgroup generated by `g`.
pub fn fp2_dlog_2e_fixed<L: FpBackend>(
    scal: &mut [u64],
    f: &Fp2<L>,
    table: &DlogTable<L>,
) -> Option<()> {
    table.rec(f, 0, table.t as usize, scal)
}
