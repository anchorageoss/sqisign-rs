//! x-only Montgomery arithmetic: doublings, differential addition, the
//! Montgomery ladder, the three-point ladder, the constant-time biladder
//! (`ec_biscalar_mul`), its variable-time counterpart for public data
//! (`ec_biscalar_mul_verif`), and barycentric coordinates.

use super::{mp, EcBaryCoordinates, EcBasis, EcCurve, EcPoint, MAX_ORDER_BITS, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};
use subtle::Choice;

/// `if ctl { p2 } else { p1 }`, constant time.
#[inline]
pub fn select_point<L: FpBackend>(p1: &EcPoint<L>, p2: &EcPoint<L>, ctl: Choice) -> EcPoint<L> {
    EcPoint {
        x: Fp2::select(&p1.x, &p2.x, ctl),
        z: Fp2::select(&p1.z, &p2.z, ctl),
    }
}

/// Swap `p` and `q` when `ctl` is set, constant time.
#[inline]
pub fn cswap_points<L: FpBackend>(p: &mut EcPoint<L>, q: &mut EcPoint<L>, ctl: Choice) {
    p.x.cswap(&mut q.x, ctl);
    p.z.cswap(&mut q.z, ctl);
}

impl<L: FpBackend> EcPoint<L> {
    /// The point at infinity has `Z = 0`.
    #[inline]
    pub fn is_zero(&self) -> Choice {
        self.z.ct_is_zero()
    }

    /// `X = 0` or `Z = 0`: the differential formulas break on such inputs.
    #[inline]
    pub fn has_zero_coordinate(&self) -> Choice {
        self.x.ct_is_zero() | self.z.ct_is_zero()
    }

    /// Projective equality on the Kummer line.
    #[inline]
    pub fn ct_equal(&self, other: &Self) -> Choice {
        let l_zero = self.is_zero();
        let r_zero = other.is_zero();
        let lr_equal = self.x.mul(&other.z).ct_equal(&self.z.mul(&other.x));
        (l_zero & r_zero) | (!l_zero & !r_zero & lr_equal)
    }

    /// `(X/Z : 1)`.
    #[inline]
    pub fn normalize(&mut self) {
        let z_inv = self.z.inv();
        self.x = self.x.mul(&z_inv);
        self.z = Fp2::one();
    }

    /// Whether `P` is a non-zero point of order 2: `x = 0` or
    /// `C x^2 + A x + C = 0`. Returns `Choice(0)` for the point at infinity.
    #[inline]
    pub fn is_two_torsion(&self, e: &EcCurve<L>) -> Choice {
        let t0 = self.x.add(&self.z).sqr();
        let t1 = self.x.sub(&self.z).sqr();
        let t2 = t0.sub(&t1).mul(&e.a);
        let t1 = t0.add(&t1).mul(&e.c);
        let t1 = t1.add(&t1);
        let t0 = t1.add(&t2); // 4 (C X^2 + C Z^2 + A X Z)
        !self.is_zero() & (self.x.ct_is_zero() | t0.ct_is_zero())
    }

    /// Whether `[2] P` is a non-zero point of order 2.
    #[inline]
    pub fn is_four_torsion(&self, e: &EcCurve<L>) -> Choice {
        xdbl_a24(self, &e.a24, e.is_a24_computed_and_normalized).is_two_torsion(e)
    }
}

/// `[2] P` on `E0` (`A = 0, C = 1`).
#[inline]
pub fn xdbl_e0<L: FpBackend>(p: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z).sqr();
    let t1 = p.x.sub(&p.z).sqr();
    let t2 = t0.sub(&t1);
    let t1 = t1.add(&t1);
    let x = t0.mul(&t1);
    let z = t1.add(&t2).mul(&t2);
    EcPoint { x, z }
}

/// `[2] P` from `(A : C)` given as a point `(x = A, z = C)`.
#[inline]
pub fn xdbl<L: FpBackend>(p: &EcPoint<L>, ac: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z).sqr();
    let t1 = p.x.sub(&p.z).sqr();
    let t2 = t0.sub(&t1);
    let t3 = ac.z.add(&ac.z);
    let t1 = t1.mul(&t3);
    let t1 = t1.add(&t1);
    let x = t0.mul(&t1);
    let t0 = t3.add(&ac.x).mul(&t2).add(&t1);
    let z = t0.mul(&t2);
    EcPoint { x, z }
}

/// `[2] P` from `A24 = (A + 2C : 4C)`, or `((A + 2C)/4C : 1)` when
/// `a24_normalized`.
#[inline]
pub fn xdbl_a24<L: FpBackend>(
    p: &EcPoint<L>,
    a24: &EcPoint<L>,
    a24_normalized: bool,
) -> EcPoint<L> {
    let t0 = p.x.add(&p.z).sqr();
    let t1 = p.x.sub(&p.z).sqr();
    let t2 = t0.sub(&t1);
    let t1 = if a24_normalized { t1 } else { t1.mul(&a24.z) };
    let x = t0.mul(&t1);
    let t0 = t2.mul(&a24.x).add(&t1);
    let z = t0.mul(&t2);
    EcPoint { x, z }
}

/// `P + Q` given `PQ = P - Q`.
#[inline]
pub fn xadd<L: FpBackend>(p: &EcPoint<L>, q: &EcPoint<L>, pq: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let t2 = q.x.add(&q.z);
    let t3 = q.x.sub(&q.z);
    let t0 = t0.mul(&t3);
    let t1 = t1.mul(&t2);
    let t2 = t0.add(&t1).sqr();
    let t3 = t0.sub(&t1).sqr();
    EcPoint {
        x: pq.z.mul(&t2),
        z: pq.x.mul(&t3),
    }
}

/// `([2] P, P + Q)` given `PQ = P - Q` and `A24`.
#[inline]
pub fn xdbladd<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    a24: &EcPoint<L>,
    a24_normalized: bool,
) -> (EcPoint<L>, EcPoint<L>) {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let rx = t0.sqr();
    let t2 = q.x.sub(&q.z);
    let sx = q.x.add(&q.z);
    let t0 = t0.mul(&t2);
    let rz = t1.sqr();
    let t1 = t1.mul(&sx);
    let t2 = rx.sub(&rz);
    let rz = if a24_normalized { rz } else { rz.mul(&a24.z) };
    let rx = rx.mul(&rz);
    let sx_tmp = a24.x.mul(&t2);
    let sz = t0.sub(&t1);
    let rz = rz.add(&sx_tmp);
    let sx = t0.add(&t1);
    let rz = rz.mul(&t2);
    let sz = sz.sqr().mul(&pq.x);
    let sx = sx.sqr().mul(&pq.z);
    (EcPoint { x: rx, z: rz }, EcPoint { x: sx, z: sz })
}

/// `[2] P` on `curve`, using the cached `A24` when normalised.
#[inline]
pub fn ec_dbl<L: FpBackend>(p: &EcPoint<L>, curve: &EcCurve<L>) -> EcPoint<L> {
    if curve.is_a24_computed_and_normalized {
        xdbl_a24(p, &curve.a24, true)
    } else {
        let ac = EcPoint {
            x: curve.a.clone(),
            z: curve.c.clone(),
        };
        xdbl(p, &ac)
    }
}

/// `[2^n] P`. Normalises `A24` on the curve for long chains.
#[inline]
pub fn ec_dbl_iter<L: FpBackend>(p: &EcPoint<L>, n: usize, curve: &mut EcCurve<L>) -> EcPoint<L> {
    if n == 0 {
        return p.clone();
    }
    if n > 50 {
        curve.normalize_a24();
    }
    let mut res = ec_dbl(p, curve);
    for _ in 1..n {
        res = ec_dbl(&res, curve);
    }
    res
}

/// `[2^n]` applied to each point of a basis.
#[inline]
pub fn ec_dbl_iter_basis<L: FpBackend>(
    b: &EcBasis<L>,
    n: usize,
    curve: &mut EcCurve<L>,
) -> EcBasis<L> {
    EcBasis {
        p: ec_dbl_iter(&b.p, n, curve),
        q: ec_dbl_iter(&b.q, n, curve),
        pmq: ec_dbl_iter(&b.pmq, n, curve),
    }
}

/// Whether `P` has order exactly `2^t`.
#[inline]
pub fn test_point_order_twof<L: FpBackend>(p: &EcPoint<L>, e: &EcCurve<L>, t: usize) -> bool {
    let mut curve = e.clone();
    if bool::from(p.is_zero()) {
        return false;
    }
    let test = ec_dbl_iter(p, t - 1, &mut curve);
    if bool::from(test.is_zero()) {
        return false;
    }
    bool::from(ec_dbl(&test, &curve).is_zero())
}

/// Whether all three points of a basis have order exactly `2^t`.
#[inline]
pub fn test_basis_order_twof<L: FpBackend>(b: &EcBasis<L>, e: &EcCurve<L>, t: usize) -> bool {
    test_point_order_twof(&b.p, e, t)
        && test_point_order_twof(&b.q, e, t)
        && test_point_order_twof(&b.pmq, e, t)
}

/// Montgomery ladder `[k] P` over the low `kbits` bits of the little-endian
/// limbs `k`. Constant time in `k`.
#[inline]
pub fn ec_mul<L: FpBackend>(
    p: &EcPoint<L>,
    k: &[u64],
    kbits: usize,
    curve: &mut EcCurve<L>,
) -> EcPoint<L> {
    if kbits > 50 {
        curve.normalize_a24();
    }
    let a24 = curve.ac_to_a24();
    let a24_normalized = curve.is_a24_computed_and_normalized;

    let mut r0 = EcPoint::<L>::identity();
    let mut r1 = p.clone();
    let mut prevbit = 0u8;
    for i in (0..kbits).rev() {
        let bit = mp::bit(k, i) as u8;
        let swap = bit ^ prevbit;
        prevbit = bit;
        cswap_points(&mut r0, &mut r1, Choice::from(swap));
        let (n0, n1) = xdbladd(&r0, &r1, p, &a24, a24_normalized);
        r0 = n0;
        r1 = n1;
    }
    cswap_points(&mut r0, &mut r1, Choice::from(prevbit));
    r0
}

/// Three-point ladder `P + [m] Q` given `PQ = P - Q`, over all bits of `m`.
/// Requires a normalised `A24`; `None` if it is not, or if `PQ` has a zero
/// coordinate.
#[inline]
pub fn ec_ladder3pt<L: FpBackend>(
    m: &[u64],
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    e: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    if !e.is_a24_computed_and_normalized || !bool::from(e.a24.z.ct_is_one()) {
        return None;
    }
    if bool::from(pq.has_zero_coordinate()) {
        return None;
    }
    let mut x0 = q.clone();
    let mut x1 = p.clone();
    let mut x2 = pq.clone();
    for &word in m {
        let mut t = 1u64;
        for _ in 0..64 {
            let ctl = Choice::from(((t & word) == 0) as u8);
            cswap_points(&mut x1, &mut x2, ctl);
            let (n0, n1) = xdbladd(&x0, &x1, &x2, &e.a24, true);
            x0 = n0;
            x1 = n1;
            cswap_points(&mut x1, &mut x2, ctl);
            t <<= 1;
        }
    }
    Some(x1)
}

/// `[k] P + [l] Q` for a basis `(P, Q, P - Q)`, constant time in the
/// scalars (the biladder of spec Algorithm 4.89). Scalars are taken modulo
/// `2^kbits`. `None` when the differential formulas do not apply (a zero
/// coordinate, or `kbits = 1` with points that are not 2-torsion).
#[inline]
pub fn ec_biscalar_mul<L: FpBackend>(
    scalar_p: &[u64],
    scalar_q: &[u64],
    kbits: usize,
    pq: &EcBasis<L>,
    curve: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    if bool::from(pq.pmq.z.ct_is_zero()) {
        return None;
    }
    if kbits == 1 {
        if !bool::from(
            pq.p.is_two_torsion(curve) & pq.q.is_two_torsion(curve) & pq.pmq.is_two_torsion(curve),
        ) {
            return None;
        }
        let bp = Choice::from((scalar_p[0] & 1) as u8);
        let bq = Choice::from((scalar_q[0] & 1) as u8);
        let t0 = select_point(&EcPoint::identity(), &pq.q, bq);
        let t1 = select_point(&pq.p, &pq.pmq, bq);
        return Some(select_point(&t0, &t1, bp));
    }
    let mut e = curve.clone();
    if !bool::from(e.a.ct_is_zero()) {
        e.normalize_a24();
    }
    xdblmul(&pq.p, scalar_p, &pq.q, scalar_q, &pq.pmq, kbits, &e)
}

#[inline]
fn xdblmul<L: FpBackend>(
    p: &EcPoint<L>,
    k_in: &[u64],
    q: &EcPoint<L>,
    l_in: &[u64],
    pq: &EcPoint<L>,
    kbits: usize,
    curve: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    debug_assert!(kbits > 0 && kbits <= MAX_ORDER_BITS);
    if bool::from(p.has_zero_coordinate() | q.has_zero_coordinate() | pq.has_zero_coordinate()) {
        return None;
    }
    let n = kbits.div_ceil(64).max(1);
    let mut k = [0u64; MAX_ORDER_WORDS];
    let mut l = [0u64; MAX_ORDER_WORDS];
    k[..k_in.len().min(n)].copy_from_slice(&k_in[..k_in.len().min(n)]);
    l[..l_in.len().min(n)].copy_from_slice(&l_in[..l_in.len().min(n)]);
    let k = &mut k[..n];
    let l = &mut l[..n];
    mp::mask_bits(k, kbits);
    mp::mask_bits(l, kbits);

    // sigma from the parities
    let bitk0 = k[0] & 1;
    let bitl0 = l[0] & 1;
    let maskk = 0u64.wrapping_sub(bitk0);
    let maskl = 0u64.wrapping_sub(bitl0);
    let mut sigma = [bitk0 ^ 1, bitl0 ^ 1];
    let evens = sigma[0] + sigma[1];
    let mevens = 0u64.wrapping_sub(evens & 1);
    sigma[0] &= mevens;
    sigma[1] = (sigma[1] & mevens) | (1 & !mevens);

    // even scalars become k - 1, l - 1 modulo 2^kbits
    let mut one = [0u64; MAX_ORDER_WORDS];
    one[0] = 1;
    let mut k_t = [0u64; MAX_ORDER_WORDS];
    let mut l_t = [0u64; MAX_ORDER_WORDS];
    mp::sub(&mut k_t[..n], k, &one[..n]);
    mp::sub(&mut l_t[..n], l, &one[..n]);
    mp::mask_bits(&mut k_t[..n], kbits);
    mp::mask_bits(&mut l_t[..n], kbits);
    let mut tmp = [0u64; MAX_ORDER_WORDS];
    tmp[..n].copy_from_slice(&k_t[..n]);
    mp::select(&mut k_t[..n], &tmp[..n], k, maskk);
    tmp[..n].copy_from_slice(&l_t[..n]);
    mp::select(&mut l_t[..n], &tmp[..n], l, maskl);

    // scalar recoding
    let mut r = [0u8; 2 * MAX_ORDER_BITS];
    let mut pre_sigma = 0u64;
    for i in 0..kbits {
        let mask_swap = 0u64.wrapping_sub(sigma[0] ^ pre_sigma);
        mp::cswap(&mut k_t[..n], &mut l_t[..n], mask_swap);
        let (bs1_ip1, bs2_ip1) = if i == kbits - 1 {
            (0, 0)
        } else {
            (mp::shr1(&mut k_t[..n]), mp::shr1(&mut l_t[..n]))
        };
        let bs1_i = k_t[0] & 1;
        let bs2_i = l_t[0] & 1;
        r[2 * i] = (bs1_i ^ bs1_ip1) as u8;
        r[2 * i + 1] = (bs2_i ^ bs2_ip1) as u8;
        pre_sigma = sigma[0];
        let mask_rev = 0u64.wrapping_sub(r[2 * i + 1] as u64);
        let temp = ((sigma[0] ^ sigma[1]) & mask_rev) ^ sigma[0];
        sigma[1] = ((sigma[1] ^ sigma[0]) & mask_rev) ^ sigma[1];
        sigma[0] = temp;
    }

    // point initialisation
    let choice_sigma = Choice::from(sigma[0] as u8);
    let mut r0 = EcPoint::<L>::identity();
    let mut r1 = select_point(p, q, choice_sigma);
    let mut r2 = select_point(q, p, choice_sigma);
    let mut diff1a = r1.clone();
    let mut diff1b = r2.clone();
    r2 = xadd(&r1, &r2, pq);
    if bool::from(r2.has_zero_coordinate()) {
        return None;
    }
    let mut diff2a = r2.clone();
    let mut diff2b = pq.clone();
    let a_is_zero = bool::from(curve.a.ct_is_zero());

    for i in (0..kbits).rev() {
        let h = r[2 * i] + r[2 * i + 1];
        let mut t0 = select_point(&r0, &r1, Choice::from(h & 1));
        t0 = select_point(&t0, &r2, Choice::from(h >> 1));
        t0 = if a_is_zero {
            xdbl_e0(&t0)
        } else {
            xdbl_a24(&t0, &curve.a24, true)
        };
        let ctl = Choice::from(r[2 * i + 1]);
        let t1a = select_point(&r0, &r1, ctl);
        let t2a = select_point(&r1, &r2, ctl);
        cswap_points(&mut diff1a, &mut diff1b, ctl);
        let t1 = xadd(&t1a, &t2a, &diff1a);
        let t2 = xadd(&r0, &r2, &diff2a);
        cswap_points(&mut diff2a, &mut diff2b, Choice::from(h & 1));
        r0 = t0;
        r1 = t1;
        r2 = t2;
    }

    let mut s = select_point(&r0, &r1, Choice::from((evens & 1) as u8));
    s = select_point(&s, &r2, Choice::from((bitk0 & bitl0) as u8));
    Some(s)
}

/// `[k] P + [l] Q` in variable time, for public scalars (adapted from
/// Algorithm 9 of ePrint 2017/212; the reference's `ec_biscalar_mul_verif`).
/// Scalars are taken modulo `2^kbits`. `None` when both scalars are zero or
/// a degenerate point arises.
#[inline]
pub fn ec_biscalar_mul_verif<L: FpBackend>(
    scalar_p: &[u64],
    scalar_q: &[u64],
    kbits: usize,
    pq: &EcBasis<L>,
    curve: &mut EcCurve<L>,
) -> Option<EcPoint<L>> {
    if bool::from(
        pq.p.has_zero_coordinate() | pq.q.has_zero_coordinate() | pq.pmq.has_zero_coordinate(),
    ) {
        return None;
    }
    let n = kbits.div_ceil(64).max(1);
    let mut s0 = [0u64; MAX_ORDER_WORDS];
    let mut s1 = [0u64; MAX_ORDER_WORDS];
    s0[..scalar_p.len().min(n)].copy_from_slice(&scalar_p[..scalar_p.len().min(n)]);
    s1[..scalar_q.len().min(n)].copy_from_slice(&scalar_q[..scalar_q.len().min(n)]);
    let s0 = &mut s0[..n];
    let s1 = &mut s1[..n];
    mp::mask_bits(s0, kbits);
    mp::mask_bits(s1, kbits);
    if mp::is_zero(s0) && mp::is_zero(s1) {
        return None;
    }

    let mut x0 = pq.p.clone();
    let mut x1 = pq.q.clone();
    let mut xdiff = pq.pmq.clone();
    if !bool::from(curve.a.ct_is_zero()) {
        curve.normalize_a24();
    }
    let a24 = curve.a24.clone();
    let a24n = curve.is_a24_computed_and_normalized;

    let mut t0 = [0u64; MAX_ORDER_WORDS];
    while !mp::is_zero(s0) {
        if mp::lt(s1, s0) {
            s0.swap_with_slice(s1);
            core::mem::swap(&mut x0, &mut x1);
        }
        t0[..n].copy_from_slice(s0);
        mp::shl(&mut t0[..n], 2);
        if mp::lt(s1, &t0[..n]) {
            let mut d = [0u64; MAX_ORDER_WORDS];
            mp::sub(&mut d[..n], s1, s0);
            s1.copy_from_slice(&d[..n]);
            let xtmp = x0.clone();
            x0 = xadd(&x1, &x0, &xdiff);
            xdiff = xtmp;
            if bool::from(x0.has_zero_coordinate()) {
                return None;
            }
        } else if mp::is_odd(s1) == mp::is_odd(s0) {
            let mut d = [0u64; MAX_ORDER_WORDS];
            mp::sub(&mut d[..n], s1, s0);
            s1.copy_from_slice(&d[..n]);
            mp::shr1(s1);
            x0 = xadd(&x1, &x0, &xdiff);
            x1 = xdbl_a24(&x1, &a24, a24n);
            if bool::from(x0.has_zero_coordinate() | x1.has_zero_coordinate()) {
                return None;
            }
        } else if !mp::is_odd(s1) {
            mp::shr1(s1);
            xdiff = xadd(&x1, &xdiff, &x0);
            x1 = xdbl_a24(&x1, &a24, a24n);
            if bool::from(xdiff.has_zero_coordinate() | x1.has_zero_coordinate()) {
                return None;
            }
        } else {
            mp::shr1(s0);
            xdiff = xadd(&x0, &xdiff, &x1);
            x0 = xdbl_a24(&x0, &a24, a24n);
            if bool::from(xdiff.has_zero_coordinate() | x0.has_zero_coordinate()) {
                return None;
            }
        }
    }
    while !mp::is_odd(s1) {
        mp::shr1(s1);
        x1 = xdbl_a24(&x1, &a24, a24n);
    }
    if !mp::is_one(s1) {
        x1 = ec_mul(&x1, s1, kbits, curve);
    }
    if bool::from(x1.has_zero_coordinate()) {
        return None;
    }
    Some(x1)
}

/// Barycentric coordinates of `(P, Q)`: `P + Q = (u - v : w)`, `P - Q = (u + v : w)`.
#[inline]
pub fn ec_points_to_bary_coordinates<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pmq: &EcPoint<L>,
) -> EcBaryCoordinates<L> {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let t2 = q.x.add(&q.z);
    let t3 = q.x.sub(&q.z);
    let t0 = t0.mul(&t3);
    let t1 = t1.mul(&t2);
    let lambda = t1.add(&t0);
    let mu = t1.sub(&t0);
    let t0 = pmq.z.mul(&lambda).sqr();
    let t1 = mu.mul(&pmq.x);
    let t2 = mu.mul(&pmq.z);
    let t3 = t1.sqr();
    let w = t1.mul(&t2);
    EcBaryCoordinates {
        u: t3.add(&t0),
        v: t3.sub(&t0),
        w: w.add(&w),
    }
}
