//!
//! Implements doubling (xDBL variants), differential addition (xADD),
//! combined double-and-add (xDBLADD), Montgomery ladder scalar
//! multiplication, three-point ladder, and biladder.

use super::{EcBasis, EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};
use subtle::Choice;
// Helper: constant-time select and swap for points

/// Constant-time point select: if `ctl` is clear return `p1`,
/// if `ctl` is set return `p2`.
#[inline]
pub fn select_point<L: FpBackend>(p1: &EcPoint<L>, p2: &EcPoint<L>, ctl: Choice) -> EcPoint<L> {
    EcPoint {
        x: Fp2::select(&p1.x, &p2.x, ctl),
        z: Fp2::select(&p1.z, &p2.z, ctl),
    }
}

/// Constant-time swap of two points when `ctl` is set.
#[inline]
pub fn cswap_points<L: FpBackend>(p: &mut EcPoint<L>, q: &mut EcPoint<L>, ctl: Choice) {
    p.x.cswap(&mut q.x, ctl);
    p.z.cswap(&mut q.z, ctl);
}

impl<L: FpBackend> EcPoint<L> {
    /// Test if this is the point at infinity (Z == 0).
    #[inline]
    pub fn is_zero(&self) -> Choice {
        self.z.ct_is_zero()
    }

    /// Test if either coordinate is zero.
    #[inline]
    pub fn has_zero_coordinate(&self) -> Choice {
        self.x.ct_is_zero() | self.z.ct_is_zero()
    }

    /// Projective equality: `P == Q` iff `PX * QZ == QX * PZ`.
    #[inline]
    pub fn ct_equal(&self, other: &Self) -> Choice {
        let l_zero = self.is_zero();
        let r_zero = other.is_zero();
        let t0 = self.x.mul(&other.z);
        let t1 = self.z.mul(&other.x);
        let lr_equal = t0.ct_equal(&t1);
        (l_zero & r_zero) | (!l_zero & !r_zero & lr_equal)
    }

    /// Normalize to `(X/Z : 1)` in place.
    #[inline]
    pub fn normalize(&mut self) {
        let z_inv = self.z.inv();
        self.x = self.x.mul(&z_inv);
        self.z = Fp2::one();
    }

    /// Test if P is 2-torsion but not zero on curve `E: y^2 = x^3 + (A/C)x^2 + x`.
    #[inline]
    pub fn is_two_torsion(&self, e: &EcCurve<L>) -> Choice {
        let not_zero = !self.is_zero();

        let t0 = self.x.add(&self.z);
        let t0 = t0.sqr();
        let t1 = self.x.sub(&self.z);
        let t1 = t1.sqr();
        let t2 = t0.sub(&t1);
        let t1 = t0.add(&t1);
        let t2 = t2.mul(&e.a);
        let t1 = t1.mul(&e.c);
        let t1 = t1.add(&t1);
        let t0 = t1.add(&t2); // 4(CX^2 + CZ^2 + AXZ)

        let x_is_zero = self.x.ct_is_zero();
        let tmp_is_zero = t0.ct_is_zero();

        not_zero & (x_is_zero | tmp_is_zero)
    }

    /// Test if P is 4-torsion (i.e. \[2\]P is 2-torsion but not zero).
    #[inline]
    pub fn is_four_torsion(&self, e: &EcCurve<L>) -> Choice {
        let test = xdbl_a24(self, &e.a24, e.is_a24_computed_and_normalized);
        test.is_two_torsion(e)
    }
}

impl<L: FpBackend> EcBasis<L> {
    /// Check if basis points (P, Q) form a full 4-torsion basis.
    #[inline]
    pub fn is_four_torsion(&self, e: &EcCurve<L>) -> Choice {
        let p2 = xdbl_a24(&self.p, &e.a24, e.is_a24_computed_and_normalized);
        let q2 = xdbl_a24(&self.q, &e.a24, e.is_a24_computed_and_normalized);
        p2.is_two_torsion(e) & q2.is_two_torsion(e) & !p2.ct_equal(&q2)
    }
}

/// Check if an (X:Z) point has order exactly 2ᵗ.
#[inline]
pub fn test_point_order_twof<L: FpBackend>(p: &EcPoint<L>, e: &EcCurve<L>, t: usize) -> Choice {
    let mut curve = e.clone();
    let mut test = p.clone();

    if bool::from(test.is_zero()) {
        return Choice::from(0u8);
    }
    test = ec_dbl_iter(&test, t - 1, &mut curve);
    if bool::from(test.is_zero()) {
        return Choice::from(0u8);
    }
    test = ec_dbl(&test, &curve);
    test.is_zero()
}

/// Check if all three basis points (P, Q, P-Q) have order exactly 2ᵗ.
#[inline]
pub fn test_basis_order_twof<L: FpBackend>(b: &EcBasis<L>, e: &EcCurve<L>, t: usize) -> Choice {
    let check_p = test_point_order_twof(&b.p, e, t);
    let check_q = test_point_order_twof(&b.q, e, t);
    let check_pmq = test_point_order_twof(&b.pmq, e, t);
    check_p & check_q & check_pmq
}

/// Check if a Jacobian point has order exactly 2ᵗ.
#[inline]
pub fn test_jac_order_twof<L: FpBackend>(
    p: &super::JacPoint<L>,
    e: &EcCurve<L>,
    t: usize,
) -> Choice {
    let mut test = p.clone();
    if bool::from(test.z.ct_is_zero()) {
        return Choice::from(0u8);
    }
    for _ in 0..t - 1 {
        test = super::jacobian::jac_dbl(&test, e);
    }
    if bool::from(test.z.ct_is_zero()) {
        return Choice::from(0u8);
    }
    test = super::jacobian::jac_dbl(&test, e);
    test.z.ct_is_zero()
}

/// Doubling on the special curve E0 with `(A:C) = (0:1)`.
#[inline]
pub fn xdbl_e0<L: FpBackend>(p: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t0 = t0.sqr();
    let t1 = p.x.sub(&p.z);
    let t1 = t1.sqr();
    let t2 = t0.sub(&t1);
    let t1 = t1.add(&t1);
    let qx = t0.mul(&t1);
    let qz_base = t1.add(&t2);
    let qz = qz_base.mul(&t2);
    EcPoint { x: qx, z: qz }
}

/// Doubling using `(A:C)` directly (computing `(A+2C : 4C)` on-the-fly).
/// The `ac` parameter is `(A : C)` packed as `(x=A, z=C)`.
#[inline]
pub fn xdbl<L: FpBackend>(p: &EcPoint<L>, ac: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t0 = t0.sqr();
    let t1 = p.x.sub(&p.z);
    let t1 = t1.sqr();
    let t2 = t0.sub(&t1);
    let t3 = ac.z.add(&ac.z);
    let t1 = t1.mul(&t3);
    let t1 = t1.add(&t1);
    let qx = t0.mul(&t1);
    let t0 = t3.add(&ac.x);
    let t0 = t0.mul(&t2);
    let t0 = t0.add(&t1);
    let qz = t0.mul(&t2);
    EcPoint { x: qx, z: qz }
}

/// Doubling using the precomputed A24 = `(A+2C : 4C)` or `((A+2C)/(4C) : 1)`
/// if `a24_normalized` is true.
#[inline]
pub fn xdbl_a24<L: FpBackend>(
    p: &EcPoint<L>,
    a24: &EcPoint<L>,
    a24_normalized: bool,
) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t0 = t0.sqr();
    let t1 = p.x.sub(&p.z);
    let t1 = t1.sqr();
    let t2 = t0.sub(&t1);
    let t1 = if !a24_normalized { t1.mul(&a24.z) } else { t1 };
    let qx = t0.mul(&t1);
    let t0 = t2.mul(&a24.x);
    let t0 = t0.add(&t1);
    let qz = t0.mul(&t2);
    EcPoint { x: qx, z: qz }
}

/// Differential addition: `R = P + Q` given `PQ = P - Q`.
#[inline]
pub fn xadd<L: FpBackend>(p: &EcPoint<L>, q: &EcPoint<L>, pq: &EcPoint<L>) -> EcPoint<L> {
    let t0 = p.x.add(&p.z);
    let t1 = p.x.sub(&p.z);
    let t2 = q.x.add(&q.z);
    let t3 = q.x.sub(&q.z);
    let t0 = t0.mul(&t3);
    let t1 = t1.mul(&t2);
    let t2 = t0.add(&t1);
    let t3 = t0.sub(&t1);
    let t2 = t2.sqr();
    let t3 = t3.sqr();
    let rx = pq.z.mul(&t2);
    let rz = pq.x.mul(&t3);
    EcPoint { x: rx, z: rz }
}

/// Simultaneous doubling and differential addition.
/// Returns `(2*P, P+Q)` given `PQ = P-Q`.
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
    let rz = if !a24_normalized { rz.mul(&a24.z) } else { rz };
    let rx = rx.mul(&rz);
    let sx_tmp = a24.x.mul(&t2);
    let sz = t0.sub(&t1);
    let rz = rz.add(&sx_tmp);
    let sx = t0.add(&t1);
    let rz = rz.mul(&t2);
    let sz = sz.sqr();
    let sx = sx.sqr();
    let sz = sz.mul(&pq.x);
    let sx = sx.mul(&pq.z);
    (EcPoint { x: rx, z: rz }, EcPoint { x: sx, z: sz })
}

/// Double a point on the given curve.
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

/// Iterated doubling: compute \[2ⁿ\]P.
#[inline]
pub fn ec_dbl_iter<L: FpBackend>(p: &EcPoint<L>, n: usize, curve: &mut EcCurve<L>) -> EcPoint<L> {
    if n == 0 {
        return p.clone();
    }

    if n > 50 {
        curve.normalize_a24();
    }

    if curve.is_a24_computed_and_normalized {
        let mut res = xdbl_a24(p, &curve.a24, true);
        for _ in 0..n - 1 {
            res = xdbl_a24(&res, &curve.a24, true);
        }
        res
    } else {
        let ac = EcPoint {
            x: curve.a.clone(),
            z: curve.c.clone(),
        };
        let mut res = xdbl(p, &ac);
        for _ in 0..n - 1 {
            res = xdbl(&res, &ac);
        }
        res
    }
}

/// Iterated doubling of a full basis {P, Q, P-Q}.
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

/// Montgomery ladder: compute `[k]P`.
///
/// `k` is a little-endian array of 64-bit limbs; `kbits` is the
/// bit-length of the scalar.
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
    xmul(p, k, kbits, curve)
}

#[inline]
fn xmul<L: FpBackend>(p: &EcPoint<L>, k: &[u64], kbits: usize, curve: &EcCurve<L>) -> EcPoint<L> {
    let a24 = if !curve.is_a24_computed_and_normalized {
        curve.ac_to_a24()
    } else {
        curve.a24.clone()
    };
    let a24_normalized = curve.is_a24_computed_and_normalized;

    let mut r0 = EcPoint::<L>::identity();
    let mut r1 = p.clone();
    let mut prevbit: u32 = 0;

    for i in (0..kbits).rev() {
        let bit = ((k[i >> 6] >> (i & 63)) & 1) as u32;
        let swap = bit ^ prevbit;
        prevbit = bit;
        cswap_points(&mut r0, &mut r1, Choice::from(swap as u8));
        let (new_r0, new_r1) = xdbladd(&r0, &r1, p, &a24, a24_normalized);
        r0 = new_r0;
        r1 = new_r1;
    }
    cswap_points(&mut r0, &mut r1, Choice::from(prevbit as u8));
    r0
}

/// Three-point Montgomery ladder: compute `P + [m]*Q` given `PQ = P - Q`.
///
/// Requires that `E.a24` is normalized (i.e. `((A+2C)/(4C) : 1)`).
/// Returns `None` if preconditions are not met.
#[inline]
pub fn ec_ladder3pt<L: FpBackend>(
    m: &[u64],
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq: &EcPoint<L>,
    e: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    if !e.is_a24_computed_and_normalized {
        return None;
    }
    if !bool::from(e.a24.z.ct_is_one()) {
        return None;
    }
    if bool::from(pq.has_zero_coordinate()) {
        return None;
    }

    let mut x0 = q.clone();
    let mut x1 = p.clone();
    let mut x2 = pq.clone();

    for &mi in m {
        let mut t: u64 = 1;
        for _ in 0..64 {
            let ctl = Choice::from(((t & mi) == 0) as u8);
            cswap_points(&mut x1, &mut x2, ctl);
            let (new_x0, new_x1) = xdbladd(&x0, &x1, &x2, &e.a24, true);
            x0 = new_x0;
            x1 = new_x1;
            cswap_points(&mut x1, &mut x2, ctl);
            t <<= 1;
        }
    }
    Some(x1)
}

/// Multiprecision subtraction: `c = a - b` (mod 2^(64*n)).
#[inline]
fn mp_sub(a: &[u64], b: &[u64], c: &mut [u64]) {
    let n = a.len();
    let mut borrow: u64 = 0;
    for i in 0..n {
        let (diff, b1) = a[i].overflowing_sub(b[i]);
        let (diff2, b2) = diff.overflowing_sub(borrow);
        c[i] = diff2;
        borrow = (b1 as u64) + (b2 as u64);
    }
}

/// Constant-time select: if mask == 0 then `c = a`, else `c = b`.
#[inline]
fn select_ct(c: &mut [u64], a: &[u64], b: &[u64], mask: u64) {
    for i in 0..c.len() {
        c[i] = (a[i] & !mask) | (b[i] & mask);
    }
}

/// Constant-time swap: if option != 0 then swap a and b.
#[inline]
fn swap_ct(a: &mut [u64], b: &mut [u64], option: u64) {
    for i in 0..a.len() {
        let t = option & (a[i] ^ b[i]);
        a[i] ^= t;
        b[i] ^= t;
    }
}

/// Right-shift a multiprecision integer by `shift` bits (shift < 64).
/// Returns the bits shifted out from the bottom.
#[inline]
fn mp_shiftr(x: &mut [u64], shift: u32) -> u64 {
    let n = x.len();
    let mut carry: u64 = 0;
    for i in (0..n).rev() {
        let new_carry = x[i] << (64 - shift);
        x[i] = (x[i] >> shift) | carry;
        carry = new_carry;
    }
    carry >> (64 - shift)
}

/// Biladder: compute `[k]*P + [l]*Q` given `PQ = P - Q`.
///
/// Returns `None` if formulas are not valid (e.g. any input has a zero
/// coordinate, or kbits == 1 degenerate case fails).
#[inline]
pub fn ec_biscalar_mul<L: FpBackend>(
    scalar_p: &[u64],
    scalar_q: &[u64],
    kbits: usize,
    pq_basis: &EcBasis<L>,
    curve: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    if bool::from(pq_basis.pmq.z.ct_is_zero()) {
        return None;
    }

    if kbits == 1 {
        // 2-torsion case
        if !bool::from(pq_basis.p.is_two_torsion(curve))
            || !bool::from(pq_basis.q.is_two_torsion(curve))
            || !bool::from(pq_basis.pmq.is_two_torsion(curve))
        {
            return None;
        }
        let bp = scalar_p[0] & 1;
        let bq = scalar_q[0] & 1;
        let res = match (bp, bq) {
            (0, 0) => EcPoint::identity(),
            (1, 0) => pq_basis.p.clone(),
            (0, 1) => pq_basis.q.clone(),
            (1, 1) => pq_basis.pmq.clone(),
            // SAFETY: bp and bq are produced by `& 1`, so they can only be 0 or 1
            _ => unreachable!("bp and bq are produced by `& 1`, so they can only be 0 or 1"),
        };
        return Some(res);
    }

    let mut e = curve.clone();
    if !bool::from(e.a.ct_is_zero()) {
        e.normalize_a24();
    }
    xdblmul(
        &pq_basis.p,
        scalar_p,
        &pq_basis.q,
        scalar_q,
        &pq_basis.pmq,
        kbits,
        &e,
    )
}

#[inline]
fn xdblmul<L: FpBackend>(
    p: &EcPoint<L>,
    k: &[u64],
    q: &EcPoint<L>,
    l: &[u64],
    pq: &EcPoint<L>,
    kbits: usize,
    curve: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    if bool::from(p.has_zero_coordinate())
        || bool::from(q.has_zero_coordinate())
        || bool::from(pq.has_zero_coordinate())
    {
        return None;
    }

    let nwords = k.len();

    // Derive sigma according to parity
    let bitk0 = k[0] & 1;
    let bitl0 = l[0] & 1;
    let maskk: u64 = 0u64.wrapping_sub(bitk0);
    let maskl: u64 = 0u64.wrapping_sub(bitl0);
    let mut sigma = [bitk0 ^ 1, bitl0 ^ 1];
    let evens = sigma[0] + sigma[1];
    let mevens: u64 = 0u64.wrapping_sub(evens & 1);

    // If both even or both odd, pick sigma = (0, 1)
    sigma[0] &= mevens;
    sigma[1] = (sigma[1] & mevens) | (1 & !mevens);

    // Convert even scalars to odd
    let mut one = [0u64; 16]; // max nwords
    one[0] = 1;
    let mut k_t = [0u64; 16];
    let mut l_t = [0u64; 16];
    let mut tmp = [0u64; 16];
    mp_sub(&k[..nwords], &one[..nwords], &mut k_t[..nwords]);
    mp_sub(&l[..nwords], &one[..nwords], &mut l_t[..nwords]);
    // select: if mask is all-ones (scalar was odd), keep original k/l
    tmp[..nwords].copy_from_slice(&k_t[..nwords]);
    select_ct(&mut k_t[..nwords], &tmp[..nwords], &k[..nwords], maskk);
    tmp[..nwords].copy_from_slice(&l_t[..nwords]);
    select_ct(&mut l_t[..nwords], &tmp[..nwords], &l[..nwords], maskl);

    // Scalar recoding
    let mut r = [0u64; 1024]; // 2 * BITS max
    let mut pre_sigma: u64 = 0;

    for i in 0..kbits {
        let mask_swap = 0u64.wrapping_sub(sigma[0] ^ pre_sigma);
        swap_ct(&mut k_t[..nwords], &mut l_t[..nwords], mask_swap);

        let (bs1_ip1, bs2_ip1) = if i == kbits - 1 {
            (0u64, 0u64)
        } else {
            let a = mp_shiftr(&mut k_t[..nwords], 1);
            let b = mp_shiftr(&mut l_t[..nwords], 1);
            (a, b)
        };
        let bs1_i = k_t[0] & 1;
        let bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        pre_sigma = sigma[0];
        let mask_rev = 0u64.wrapping_sub(r[2 * i + 1]);
        let mut temp = 0u64;
        select_ct(
            core::slice::from_mut(&mut temp),
            core::slice::from_ref(&sigma[0]),
            core::slice::from_ref(&sigma[1]),
            mask_rev,
        );
        let mut s1_new = 0u64;
        select_ct(
            core::slice::from_mut(&mut s1_new),
            core::slice::from_ref(&sigma[1]),
            core::slice::from_ref(&sigma[0]),
            mask_rev,
        );
        sigma[0] = temp;
        sigma[1] = s1_new;
    }

    // Point initialization
    let mut pts = [
        EcPoint::<L>::identity(),
        EcPoint::<L>::identity(),
        EcPoint::<L>::identity(),
    ];
    let choice_sigma = Choice::from(sigma[0] as u8);
    pts[1] = select_point(p, q, choice_sigma);
    pts[2] = select_point(q, p, choice_sigma);

    let mut diff1a = pts[1].clone();
    let mut diff1b = pts[2].clone();

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    pts[2] = xadd(&pts[1], &pts[2], pq);
    if bool::from(pts[2].has_zero_coordinate()) {
        return None;
    }

    let mut diff2a = pts[2].clone();
    let mut diff2b = pq.clone();

    let a_is_zero = curve.a.ct_is_zero();

    // Main loop
    for i in (0..kbits).rev() {
        let h = r[2 * i] + r[2 * i + 1];
        let choice1 = Choice::from((h & 1) as u8);
        let mut t0 = select_point(&pts[0], &pts[1], choice1);
        let choice2 = Choice::from((h >> 1) as u8);
        t0 = select_point(&t0, &pts[2], choice2);

        t0 = if bool::from(a_is_zero) {
            xdbl_e0(&t0)
        } else {
            xdbl_a24(&t0, &curve.a24, true)
        };

        let choice_r1 = Choice::from(r[2 * i + 1] as u8);
        let t1_a = select_point(&pts[0], &pts[1], choice_r1);
        let t2_a = select_point(&pts[1], &pts[2], choice_r1);

        cswap_points(&mut diff1a, &mut diff1b, choice_r1);
        let t1 = xadd(&t1_a, &t2_a, &diff1a);
        let t2 = xadd(&pts[0], &pts[2], &diff2a);

        let choice_h = Choice::from((h & 1) as u8);
        cswap_points(&mut diff2a, &mut diff2b, choice_h);

        pts[0] = t0;
        pts[1] = t1;
        pts[2] = t2;
    }

    // Output R[evens]
    let choice_evens = Choice::from((evens & 1) as u8);
    let mut s = select_point(&pts[0], &pts[1], choice_evens);
    let choice_both = Choice::from((bitk0 & bitl0) as u8);
    s = select_point(&s, &pts[2], choice_both);

    Some(s)
}
