//! 2-power discrete logarithms on Montgomery curves over Fp2.

use super::basis::lift_basis_normalized;
use super::jacobian::jac_add;
use super::{EcBasis, EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};

type DlogResult = ([u64; 8], [u64; 8], [u64; 8], [u64; 8]);

/// Pairing computation parameters: normalized copies of P, Q, P-Q and
/// the cached inverse x-coordinates needed by the cubical ladder.
#[derive(Clone)]
struct PairingParams<L: FpBackend> {
    e: u32,
    p: EcPoint<L>,
    q: EcPoint<L>,
    pq: EcPoint<L>,
    ix_p: Fp2<L>,
    ix_q: Fp2<L>,
    a24: EcPoint<L>,
}

/// Cross-basis difference points needed for dlog: x(P-R), x(P-S), x(R-Q), x(S-Q).
#[derive(Clone)]
struct DlogDiffPoints<L: FpBackend> {
    pm_r: EcPoint<L>,
    pm_s: EcPoint<L>,
    rm_q: EcPoint<L>,
    sm_q: EcPoint<L>,
}

/// Full dlog parameter bundle for two bases {P,Q} and {R,S}.
#[derive(Clone)]
struct PairingDlogParams<L: FpBackend> {
    e: u32,
    pq: EcBasis<L>,
    rs: EcBasis<L>,
    diff: DlogDiffPoints<L>,
    ix_p: Fp2<L>,
    ix_q: Fp2<L>,
    ix_r: Fp2<L>,
    ix_s: Fp2<L>,
    a24: EcPoint<L>,
}

/// Cubical addition: given cubical reps of P, Q and `ix_pq = Z(P-Q)/X(P-Q)`,
/// compute the cubical rep of P+Q.
///
/// Cost: 3M + 2S + 3a + 3s
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
    let rz = t3.sqr();
    let t2 = t2.sqr();
    let rx = ix_pq.mul(&t2);
    EcPoint::new(rx, rz)
}

/// Cubical combined double-and-add: given cubical reps of P, Q and
/// `ix_pq = Z(P-Q)/X(P-Q)`, compute (P+Q, [2]Q).
///
/// A24 must be normalized to `((A+2C)/(4C) : 1)`.
///
/// Cost: 6M + 4S + 4a + 4s
#[inline]
fn cubical_dbladd<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    ix_pq: &Fp2<L>,
    a24: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    debug_assert!(bool::from(a24.z.ct_is_one()));

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
    let ppq_x = ppq_x.sqr();
    let ppq_x = ix_pq.mul(&ppq_x);
    let t3 = t2.sub(&qq_z);
    let qq_x = t2.mul(&qq_z);
    let t0 = t3.mul(&a24.x);
    let t0 = t0.add(&qq_z);
    let qq_z = t0.mul(&t3);
    (EcPoint::new(ppq_x, ppq_z), EcPoint::new(qq_x, qq_z))
}

/// Iterative biextension ladder: compute (P + [2ᵉ]Q, [2ᵉ]Q).
#[inline]
fn biext_ladder_2e<L: FpBackend>(
    e: u32,
    pq: &EcPoint<L>,
    q: &EcPoint<L>,
    ix_p: &Fp2<L>,
    a24: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    let mut pn_q = pq.clone();
    let mut n_q = q.clone();
    for _ in 0..e {
        let (new_pnq, new_nq) = cubical_dbladd(&pn_q, &n_q, ix_p, a24);
        pn_q = new_pnq;
        n_q = new_nq;
    }
    (pn_q, n_q)
}

/// Compute the monodromy ratio as a projective point.
#[inline]
fn point_ratio<L: FpBackend>(pn_q: &EcPoint<L>, n_q: &EcPoint<L>, p: &EcPoint<L>) -> EcPoint<L> {
    let rx = n_q.x.mul(&p.x);
    let rz = pn_q.x.clone();
    EcPoint::new(rx, rz)
}

/// Cubical translation of P by a 2-torsion point T, computed in constant time.
#[inline]
fn translate<L: FpBackend>(p: &mut EcPoint<L>, t: &EcPoint<L>) {
    // Generic case: (AX - BZ, BX - AZ)
    let ax = t.x.mul(&p.x);
    let bz = t.z.mul(&p.z);
    let px_generic = ax.sub(&bz);

    let bx = t.z.mul(&p.x);
    let az = t.x.mul(&p.z);
    let pz_generic = bx.sub(&az);

    // If T.x == 0 then result is (Z, X)
    let ta_is_zero = t.x.ct_is_zero();
    let px_new = Fp2::select(&px_generic, &p.z, ta_is_zero);
    let pz_new = Fp2::select(&pz_generic, &p.x, ta_is_zero);

    // If T.z == 0 then result is (X, Z) (identity translation)
    let tb_is_zero = t.z.ct_is_zero();
    let px_new = Fp2::select(&px_new, &p.x, tb_is_zero);
    let pz_new = Fp2::select(&pz_new, &p.z, tb_is_zero);

    p.x = px_new;
    p.z = pz_new;
}

/// Normalize P, Q for pairing computation: store `(X/Z : 1)` for each
/// point plus `Z/X` (the "inverse x-coordinate").
#[inline]
fn cubical_normalization<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>, Fp2<L>, Fp2<L>) {
    let mut t = [p.x.clone(), p.z.clone(), q.x.clone(), q.z.clone()];
    let mut s1 = [Fp2::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
    let mut s2 = [Fp2::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
    Fp2::batched_inv(&mut t, &mut s1, &mut s2);

    let ix_p = p.z.mul(&t[0]);
    let ix_q = q.z.mul(&t[2]);

    let np = EcPoint::new(p.x.mul(&t[1]), Fp2::one());
    let nq = EcPoint::new(q.x.mul(&t[3]), Fp2::one());

    (np, nq, ix_p, ix_q)
}

/// Compute the biextension monodromy via the cubical ladder.
/// When `swap_pq` is false, computes from P + [2ᵉ]Q.
/// When `swap_pq` is true, computes from Q + [2ᵉ]P.
#[inline]
fn monodromy_i<L: FpBackend>(params: &PairingParams<L>, swap_pq: bool) -> EcPoint<L> {
    let (p, q, ix_p) = if !swap_pq {
        (params.p.clone(), params.q.clone(), params.ix_p.clone())
    } else {
        (params.q.clone(), params.p.clone(), params.ix_q.clone())
    };

    let (mut pn_q, mut n_q) = biext_ladder_2e(params.e - 1, &params.pq, &q, &ix_p, &params.a24);
    translate(&mut pn_q, &n_q);
    let n_q_copy = n_q.clone();
    translate(&mut n_q, &n_q_copy);
    point_ratio(&pn_q, &n_q, &p)
}

/// Compute the Weil pairing value from normalized pairing data.
#[inline]
fn weil_n<L: FpBackend>(params: &PairingParams<L>) -> Fp2<L> {
    let r0 = monodromy_i(params, true);
    let r1 = monodromy_i(params, false);

    let r = r0.x.mul(&r1.z);
    let r = r.inv();
    let r = r.mul(&r0.z);
    r.mul(&r1.x)
}

/// Weil pairing e_{2ᵉ}(P, Q) via the biextension cubical ladder.
///
/// `pq` must be `x(P - Q)` in (X:Z) coordinates. Crashes (division by
/// zero) if either P or Q is the identity.
#[inline]
pub fn weil<L: FpBackend>(
    e: u32,
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    curve: &mut EcCurve<L>,
) -> Fp2<L> {
    let (np, nq, ix_p, ix_q) = cubical_normalization(p, q);
    curve.normalize_a24();

    let params = PairingParams {
        e,
        p: np,
        q: nq,
        pq: pq.clone(),
        ix_p,
        ix_q,
        a24: curve.a24.clone(),
    };
    weil_n(&params)
}

/// Clear the cofactor (p+1) / 2ᶠ from an 𝔽p² element by
/// exponentiation. Uses `p_cofactor_for_2f` (a small odd integer).
#[inline]
pub fn clear_cofac<L: FpBackend>(a: &Fp2<L>, cofactor: &[u64]) -> Fp2<L> {
    let mut exp = cofactor[0];
    exp >>= 1;
    let x = a.clone();
    let mut r = a.clone();
    while exp > 0 {
        r = r.sqr();
        if exp & 1 != 0 {
            r = r.mul(&x);
        }
        exp >>= 1;
    }
    r
}

/// Frobenius endomorphism on 𝔽p²: a + bi → a − bi.
/// This is conjugation since `p = 3 mod 4`.
#[inline]
pub fn fp2_frob<L: FpBackend>(x: &Fp2<L>) -> Fp2<L> {
    x.conjugate()
}

/// Reduced Tate pairing t_{2ᵉ}(P, Q) via the biextension cubical ladder.
///
/// `pq` must be `x(P - Q)` in (X:Z) coordinates.
/// Computes the unreduced pairing and applies ^((p²−1)/2ᶠ).
#[inline]
pub fn reduced_tate<L: FpBackend>(
    e: u32,
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    curve: &mut EcCurve<L>,
    torsion_even_power: u32,
    cofactor: &[u64],
) -> Fp2<L> {
    let e_diff = torsion_even_power - e;

    let (np, nq, ix_p, ix_q) = cubical_normalization(p, q);
    curve.normalize_a24();

    let params = PairingParams {
        e,
        p: np,
        q: nq,
        pq: pq.clone(),
        ix_p,
        ix_q,
        a24: curve.a24.clone(),
    };
    let r_pt = monodromy_i(&params, true);

    // Reduce: -(R.Z / R.X)^((p^2-1)/2^f)
    // Split ^(p-1) into Frobenius and ^(-1)
    let frob_rx = fp2_frob::<L>(&r_pt.x);
    let new_rx = r_pt.z.mul(&frob_rx);
    let frob_rz = fp2_frob::<L>(&r_pt.z);
    let new_rz = r_pt.x.mul(&frob_rz);
    let inv_rx = new_rx.inv();
    let r = inv_rx.mul(&new_rz);

    let mut r = clear_cofac::<L>(&r, cofactor);
    for _ in 0..e_diff {
        r = r.sqr();
    }
    r
}

/// Recursive 2-power discrete log: find `a` s.t. `f = g^a` in the
/// `2^len`-subgroup, given stacks of powers of f and g_inverse.
#[allow(clippy::needless_range_loop)]
#[inline]
fn fp2_dlog_2e_rec<L: FpBackend>(
    a: &mut [u64],
    len: usize,
    pows_f: &mut [Fp2<L>],
    pows_g: &mut [Fp2<L>],
    stacklen: usize,
) -> Option<()> {
    let nwords = a.len();
    if len == 0 {
        for w in a.iter_mut() {
            *w = 0;
        }
        return Some(());
    } else if len == 1 {
        if bool::from(pows_f[stacklen - 1].ct_is_one()) {
            for w in a.iter_mut() {
                *w = 0;
            }
            for i in 0..stacklen - 1 {
                pows_g[i] = pows_g[i].sqr();
            }
            return Some(());
        } else if bool::from(pows_f[stacklen - 1].ct_equal(&pows_g[stacklen - 1])) {
            a[0] = 1;
            for w in a[1..].iter_mut() {
                *w = 0;
            }
            for i in 0..stacklen - 1 {
                pows_f[i] = pows_f[i].mul(&pows_g[i]);
                pows_g[i] = pows_g[i].sqr();
            }
            return Some(());
        } else {
            return None;
        }
    }

    let right = (len as f64 * 0.5) as usize;
    let left = len - right;
    pows_f[stacklen] = pows_f[stacklen - 1].clone();
    pows_g[stacklen] = pows_g[stacklen - 1].clone();
    for _ in 0..left {
        pows_f[stacklen] = pows_f[stacklen].sqr();
        pows_g[stacklen] = pows_g[stacklen].sqr();
    }

    let mut dlp1 = [0u64; 8]; // max NWORDS_ORDER across all levels
    let mut dlp2 = [0u64; 8];
    let dlp1_slice = &mut dlp1[..nwords];
    let dlp2_slice = &mut dlp2[..nwords];

    fp2_dlog_2e_rec(dlp1_slice, right, pows_f, pows_g, stacklen + 1)?;
    fp2_dlog_2e_rec(dlp2_slice, left, pows_f, pows_g, stacklen)?;

    // a = dlp1 + 2^right * dlp2
    mp_shiftl_multiple(dlp2_slice, right);
    mp_add_inplace(a, dlp2_slice, dlp1_slice);

    Some(())
}

/// Compute discrete log: find `scal` such that `f = g^scal` where `g_inverse = g⁻¹`,
/// in the 2ᵉ-subgroup of 𝔽p²*.
#[inline]
fn fp2_dlog_2e<L: FpBackend>(
    scal: &mut [u64],
    f: &Fp2<L>,
    g_inverse: &Fp2<L>,
    e: u32,
) -> Option<()> {
    // Compute stack depth: ceil(log2(e)) + 1
    let mut log = 0usize;
    let mut len = e as usize;
    while len > 1 {
        len >>= 1;
        log += 1;
    }
    log += 1;

    const MAX_STACK: usize = 16;
    debug_assert!(log <= MAX_STACK);
    let mut pows_f: [Fp2<L>; MAX_STACK] = core::array::from_fn(|_| Fp2::zero());
    let mut pows_g: [Fp2<L>; MAX_STACK] = core::array::from_fn(|_| Fp2::zero());
    pows_f[0] = f.clone();
    pows_g[0] = g_inverse.clone();

    for w in scal.iter_mut() {
        *w = 0;
    }

    fp2_dlog_2e_rec(scal, e as usize, &mut pows_f, &mut pows_g, 1)
}

/// Find `scal` such that `f = g^scal` in the 2ᵉ-subgroup of `Fp2*`,
/// given `g_inverse = g⁻¹`. Returns `None` if the DLP has no solution.
#[inline]
pub fn fp2_dlog_2e_pub<L: FpBackend>(
    scal: &mut [u64],
    f: &Fp2<L>,
    g_inverse: &Fp2<L>,
    e: u32,
) -> Option<()> {
    fp2_dlog_2e(scal, f, g_inverse, e)
}

/// Left-shift a multiprecision integer by `shift` bits.
#[inline]
fn mp_shiftl_multiple(x: &mut [u64], shift: usize) {
    let mut remaining = shift;
    while remaining > 63 {
        mp_shiftl_single(x, 63);
        remaining -= 63;
    }
    if remaining > 0 {
        mp_shiftl_single(x, remaining);
    }
}

/// Left-shift by 1..63 bits.
#[inline]
fn mp_shiftl_single(x: &mut [u64], shift: usize) {
    let n = x.len();
    for i in (1..n).rev() {
        x[i] = (x[i] << shift) | (x[i - 1] >> (64 - shift));
    }
    x[0] <<= shift;
}

/// Multiprecision addition: `c = a + b`.
#[inline]
fn mp_add_inplace(c: &mut [u64], a: &[u64], b: &[u64]) {
    let n = c.len();
    let mut carry = 0u64;
    for i in 0..n {
        let (s1, c1) = a[i].overflowing_add(b[i]);
        let (s2, c2) = s1.overflowing_add(carry);
        c[i] = s2;
        carry = (c1 as u64) + (c2 as u64);
    }
}

/// Normalize both bases {P,Q} and {R,S} plus the curve coefficient
/// for dlog computation, computing inverse x-coordinates.
#[inline]
fn cubical_normalization_dlog<L: FpBackend>(
    data: &mut PairingDlogParams<L>,
    curve: &mut EcCurve<L>,
) {
    let mut t = [
        data.pq.p.x.clone(),   // 0
        data.pq.p.z.clone(),   // 1
        data.pq.q.x.clone(),   // 2
        data.pq.q.z.clone(),   // 3
        data.pq.pmq.x.clone(), // 4
        data.pq.pmq.z.clone(), // 5
        data.rs.p.x.clone(),   // 6
        data.rs.p.z.clone(),   // 7
        data.rs.q.x.clone(),   // 8
        data.rs.q.z.clone(),   // 9
        curve.c.clone(),       // 10
    ];
    let mut s1 = [
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
    ];
    let mut s2 = [
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
    ];
    Fp2::batched_inv(&mut t, &mut s1, &mut s2);

    data.ix_p = data.pq.p.z.mul(&t[0]);
    data.pq.p.x = data.pq.p.x.mul(&t[1]);
    data.pq.p.z = Fp2::one();

    data.ix_q = data.pq.q.z.mul(&t[2]);
    data.pq.q.x = data.pq.q.x.mul(&t[3]);
    data.pq.q.z = Fp2::one();

    data.pq.pmq.x = data.pq.pmq.x.mul(&t[5]);
    data.pq.pmq.z = Fp2::one();

    data.ix_r = data.rs.p.z.mul(&t[6]);
    data.rs.p.x = data.rs.p.x.mul(&t[7]);
    data.rs.p.z = Fp2::one();

    data.ix_s = data.rs.q.z.mul(&t[8]);
    data.rs.q.x = data.rs.q.x.mul(&t[9]);
    data.rs.q.z = Fp2::one();

    curve.a = curve.a.mul(&t[10]);
    curve.c = Fp2::one();
}

/// Compute the four cross-basis difference points by lifting to
/// Jacobian coordinates: x(P-R), x(P-S), x(R-Q), x(S-Q).
#[inline]
fn compute_difference_points<L: FpBackend>(data: &mut PairingDlogParams<L>, curve: &EcCurve<L>) {
    let mut pq_copy = data.pq.clone();
    let mut rs_copy = data.rs.clone();
    let (xy_p, xy_q, _) = lift_basis_normalized(&mut pq_copy, curve);
    let (xy_r, xy_s, _) = lift_basis_normalized(&mut rs_copy, curve);

    // x(P - R)
    let neg_r = xy_r.neg();
    let temp = jac_add(&neg_r, &xy_p, curve);
    data.diff.pm_r = temp.to_xz();

    // x(P - S)
    let neg_s = xy_s.neg();
    let temp = jac_add(&neg_s, &xy_p, curve);
    data.diff.pm_s = temp.to_xz();

    // x(R - Q)
    let neg_q = xy_q.neg();
    let temp = jac_add(&neg_q, &xy_r, curve);
    data.diff.rm_q = temp.to_xz();

    // x(S - Q)
    let temp = jac_add(&neg_q, &xy_s, curve);
    data.diff.sm_q = temp.to_xz();
}

/// Inline all Weil pairing computations needed for the Weil-based dlog.
#[inline]
fn weil_dlog<L: FpBackend>(data: &PairingDlogParams<L>) -> Option<DlogResult> {
    let nwords = L::NWORDS_ORDER;

    let mut n_p = data.pq.p.clone();
    let mut n_q = data.pq.q.clone();
    let mut n_r = data.rs.p.clone();
    let mut n_s = data.rs.q.clone();
    let mut n_pq = data.pq.pmq.clone();
    let mut p_nq = data.pq.pmq.clone();
    let mut n_pr = data.diff.pm_r.clone();
    let mut n_ps = data.diff.pm_s.clone();
    let mut p_nr = data.diff.pm_r.clone();
    let mut p_ns = data.diff.pm_s.clone();
    let mut n_rq = data.diff.rm_q.clone();
    let mut n_sq = data.diff.sm_q.clone();
    let mut r_nq = data.diff.rm_q.clone();
    let mut s_nq = data.diff.sm_q.clone();

    for _ in 0..data.e - 1 {
        n_pq = cubical_add(&n_pq, &n_p, &data.ix_q);
        n_pr = cubical_add(&n_pr, &n_p, &data.ix_r);
        let (new_nps, new_np) = cubical_dbladd(&n_ps, &n_p, &data.ix_s, &data.a24);
        n_ps = new_nps;
        n_p = new_np;

        p_nq = cubical_add(&p_nq, &n_q, &data.ix_p);
        r_nq = cubical_add(&r_nq, &n_q, &data.ix_r);
        let (new_snq, new_nq) = cubical_dbladd(&s_nq, &n_q, &data.ix_s, &data.a24);
        s_nq = new_snq;
        n_q = new_nq;

        p_nr = cubical_add(&p_nr, &n_r, &data.ix_p);
        let (new_nrq, new_nr) = cubical_dbladd(&n_rq, &n_r, &data.ix_q, &data.a24);
        n_rq = new_nrq;
        n_r = new_nr;

        p_ns = cubical_add(&p_ns, &n_s, &data.ix_p);
        let (new_nsq, new_ns) = cubical_dbladd(&n_sq, &n_s, &data.ix_q, &data.a24);
        n_sq = new_nsq;
        n_s = new_ns;
    }

    // Translate
    translate(&mut n_pq, &n_p);
    translate(&mut n_pr, &n_p);
    translate(&mut n_ps, &n_p);
    translate(&mut p_nq, &n_q);
    translate(&mut r_nq, &n_q);
    translate(&mut s_nq, &n_q);
    translate(&mut p_nr, &n_r);
    translate(&mut n_rq, &n_r);
    translate(&mut p_ns, &n_s);
    translate(&mut n_sq, &n_s);

    let n_p_clone = n_p.clone();
    let n_q_clone = n_q.clone();
    let n_r_clone = n_r.clone();
    let n_s_clone = n_s.clone();
    translate(&mut n_p, &n_p_clone);
    translate(&mut n_q, &n_q_clone);
    translate(&mut n_r, &n_r_clone);
    translate(&mut n_s, &n_s_clone);

    // Compute reference pairing ratios
    let mut w1 = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    let mut w2 = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];

    // e(P, Q) = w0, note: w1/w2 swapped for first element to save an inversion
    let t0 = point_ratio(&n_pq, &n_p, &data.pq.q);
    let t1 = point_ratio(&p_nq, &n_q, &data.pq.p);
    w2[0] = t0.x.mul(&t1.z);
    w1[0] = t1.x.mul(&t0.z);

    // e(P, R) = w0^r2
    let t0 = point_ratio(&n_pr, &n_p, &data.rs.p);
    let t1 = point_ratio(&p_nr, &n_r, &data.pq.p);
    w1[1] = t0.x.mul(&t1.z);
    w2[1] = t1.x.mul(&t0.z);

    // e(R, Q) = w0^r1
    let t0 = point_ratio(&n_rq, &n_r, &data.pq.q);
    let t1 = point_ratio(&r_nq, &n_q, &data.rs.p);
    w1[2] = t0.x.mul(&t1.z);
    w2[2] = t1.x.mul(&t0.z);

    // e(P, S) = w0^s2
    let t0 = point_ratio(&n_ps, &n_p, &data.rs.q);
    let t1 = point_ratio(&p_ns, &n_s, &data.pq.p);
    w1[3] = t0.x.mul(&t1.z);
    w2[3] = t1.x.mul(&t0.z);

    // e(S, Q) = w0^s1
    let t0 = point_ratio(&n_sq, &n_s, &data.pq.q);
    let t1 = point_ratio(&s_nq, &n_q, &data.rs.q);
    w1[4] = t0.x.mul(&t1.z);
    w2[4] = t1.x.mul(&t0.z);

    // Batch inversion and normalization
    let mut s1_scratch = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    let mut s2_scratch = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    Fp2::batched_inv(&mut w1, &mut s1_scratch, &mut s2_scratch);
    for i in 0..5 {
        w1[i] = w1[i].mul(&w2[i]);
    }

    let mut r1 = [0u64; 8];
    let mut r2 = [0u64; 8];
    let mut s1 = [0u64; 8];
    let mut s2 = [0u64; 8];

    fp2_dlog_2e(&mut r2[..nwords], &w1[1], &w1[0], data.e)?;
    fp2_dlog_2e(&mut r1[..nwords], &w1[2], &w1[0], data.e)?;
    fp2_dlog_2e(&mut s2[..nwords], &w1[3], &w1[0], data.e)?;
    fp2_dlog_2e(&mut s1[..nwords], &w1[4], &w1[0], data.e)?;

    Some((r1, r2, s1, s2))
}

/// Compute 2-power discrete log using the Weil pairing.
///
/// Given bases `{P, Q}` and `{R, S}` of the 2ᵉ-torsion, find
/// scalars `r1, r2, s1, s2` such that `R = [r1]P + [r2]Q` and
/// `S = [s1]P + [s2]Q`.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn ec_dlog_2_weil<L: FpBackend>(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pq: &EcBasis<L>,
    rs: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    e: u32,
) -> Option<()> {
    curve.normalize_a24();

    let mut data = PairingDlogParams {
        e,
        pq: pq.clone(),
        rs: rs.clone(),
        diff: DlogDiffPoints {
            pm_r: EcPoint::identity(),
            pm_s: EcPoint::identity(),
            rm_q: EcPoint::identity(),
            sm_q: EcPoint::identity(),
        },
        ix_p: Fp2::zero(),
        ix_q: Fp2::zero(),
        ix_r: Fp2::zero(),
        ix_s: Fp2::zero(),
        a24: curve.a24.clone(),
    };

    cubical_normalization_dlog(&mut data, curve);
    compute_difference_points(&mut data, curve);

    let (wr1, wr2, ws1, ws2) = weil_dlog(&data)?;
    let n = r1.len();
    r1.copy_from_slice(&wr1[..n]);
    r2.copy_from_slice(&wr2[..n]);
    s1.copy_from_slice(&ws1[..n]);
    s2.copy_from_slice(&ws2[..n]);
    Some(())
}

/// Inline all Tate pairing computations for partial-torsion dlog.
#[inline]
fn tate_dlog_partial<L: FpBackend>(
    data: &PairingDlogParams<L>,
    torsion_even_power: u32,
    cofactor: &[u64],
) -> Option<DlogResult> {
    let nwords = L::NWORDS_ORDER;
    let e_diff = torsion_even_power - data.e;

    let mut n_p = data.pq.p.clone();
    let mut n_q = data.pq.q.clone();
    let mut n_r = data.rs.p.clone();
    let mut n_s = data.rs.q.clone();
    let mut n_pq = data.pq.pmq.clone();
    let mut p_nr = data.diff.pm_r.clone();
    let mut p_ns = data.diff.pm_s.clone();
    let mut n_rq = data.diff.rm_q.clone();
    let mut n_sq = data.diff.sm_q.clone();

    // Full-order ladder for P, Q
    for _ in 0..torsion_even_power - 1 {
        let (new_npq, new_np) = cubical_dbladd(&n_pq, &n_p, &data.ix_q, &data.a24);
        n_pq = new_npq;
        n_p = new_np;
    }

    // Partial-order ladders for R, S
    for _ in 0..data.e - 1 {
        p_nr = cubical_add(&p_nr, &n_r, &data.ix_p);
        let (new_nrq, new_nr) = cubical_dbladd(&n_rq, &n_r, &data.ix_q, &data.a24);
        n_rq = new_nrq;
        n_r = new_nr;

        p_ns = cubical_add(&p_ns, &n_s, &data.ix_p);
        let (new_nsq, new_ns) = cubical_dbladd(&n_sq, &n_s, &data.ix_q, &data.a24);
        n_sq = new_nsq;
        n_s = new_ns;
    }

    translate(&mut n_pq, &n_p);
    translate(&mut p_nr, &n_r);
    translate(&mut n_rq, &n_r);
    translate(&mut p_ns, &n_s);
    translate(&mut n_sq, &n_s);

    let n_p_clone = n_p.clone();
    let n_q_clone = n_q.clone();
    let n_r_clone = n_r.clone();
    let n_s_clone = n_s.clone();
    translate(&mut n_p, &n_p_clone);
    translate(&mut n_q, &n_q_clone);
    translate(&mut n_r, &n_r_clone);
    translate(&mut n_s, &n_s_clone);

    // Compute Tate ratios
    let mut w1 = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    let mut w2 = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];

    // t(P, Q)
    let t0 = point_ratio(&n_pq, &n_p, &data.pq.q);
    w1[0] = t0.x.clone();
    w2[0] = t0.z.clone();

    // t(R, P) = w0^r2
    let t0 = point_ratio(&p_nr, &n_r, &data.pq.p);
    w1[1] = t0.x.clone();
    w2[1] = t0.z.clone();

    // t(R, Q) = w0^r1 , note swapped w1/w2
    let t0 = point_ratio(&n_rq, &n_r, &data.pq.q);
    w2[2] = t0.x.clone();
    w1[2] = t0.z.clone();

    // t(S, P) = w0^s2
    let t0 = point_ratio(&p_ns, &n_s, &data.pq.p);
    w1[3] = t0.x.clone();
    w2[3] = t0.z.clone();

    // t(S, Q) = w0^s1 , note swapped w1/w2
    let t0 = point_ratio(&n_sq, &n_s, &data.pq.q);
    w2[4] = t0.x.clone();
    w1[4] = t0.z.clone();

    // Batched reduction: apply ^(p-1) via Frobenius
    for i in 0..5 {
        let tmp = w1[i].clone();
        let frob = fp2_frob::<L>(&w1[i]);
        w1[i] = w2[i].mul(&frob);
        let frob = fp2_frob::<L>(&w2[i]);
        w2[i] = tmp.mul(&frob);
    }

    // Batch normalize
    let mut s1_scratch = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    let mut s2_scratch = [
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
        Fp2::<L>::zero(),
    ];
    Fp2::batched_inv(&mut w2, &mut s1_scratch, &mut s2_scratch);
    for i in 0..5 {
        w1[i] = w1[i].mul(&w2[i]);
    }

    // Clear cofactor and remaining 2^e_diff
    for item in w1.iter_mut() {
        *item = clear_cofac::<L>(item, cofactor);
        for _ in 0..e_diff {
            *item = item.sqr();
        }
    }

    let mut r1 = [0u64; 8];
    let mut r2 = [0u64; 8];
    let mut s1 = [0u64; 8];
    let mut s2 = [0u64; 8];

    fp2_dlog_2e(&mut r2[..nwords], &w1[1], &w1[0], data.e)?;
    fp2_dlog_2e(&mut r1[..nwords], &w1[2], &w1[0], data.e)?;
    fp2_dlog_2e(&mut s2[..nwords], &w1[3], &w1[0], data.e)?;
    fp2_dlog_2e(&mut s1[..nwords], &w1[4], &w1[0], data.e)?;

    Some((r1, r2, s1, s2))
}

/// Compute 2-power discrete log using the reduced Tate pairing.
///
/// `{P, Q}` must be a basis of the full `2^torsion_even_power`-torsion.
/// `{R, S}` is a basis of the 2ᵉ-torsion (where `e <= torsion_even_power`).
/// Finds scalars `r1, r2, s1, s2` such that
/// `R = [2^(f-e)]([r1]P + [r2]Q)` and `S = [2^(f-e)]([s1]P + [s2]Q)`.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn ec_dlog_2_tate<L: FpBackend>(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pq: &EcBasis<L>,
    rs: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    e: u32,
    torsion_even_power: u32,
    cofactor: &[u64],
) -> Option<()> {
    curve.normalize_a24();

    let mut data = PairingDlogParams {
        e,
        pq: pq.clone(),
        rs: rs.clone(),
        diff: DlogDiffPoints {
            pm_r: EcPoint::identity(),
            pm_s: EcPoint::identity(),
            rm_q: EcPoint::identity(),
            sm_q: EcPoint::identity(),
        },
        ix_p: Fp2::zero(),
        ix_q: Fp2::zero(),
        ix_r: Fp2::zero(),
        ix_s: Fp2::zero(),
        a24: curve.a24.clone(),
    };

    cubical_normalization_dlog(&mut data, curve);
    compute_difference_points(&mut data, curve);

    let (tr1, tr2, ts1, ts2) = tate_dlog_partial(&data, torsion_even_power, cofactor)?;
    let n = r1.len();
    r1.copy_from_slice(&tr1[..n]);
    r2.copy_from_slice(&tr2[..n]);
    s1.copy_from_slice(&ts1[..n]);
    s2.copy_from_slice(&ts2[..n]);
    Some(())
}
