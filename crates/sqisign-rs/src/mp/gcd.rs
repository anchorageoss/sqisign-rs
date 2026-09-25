//! Variable-time gcd, extended gcd, modular inverse, CRT and the Jacobi
//! symbol, with Jebelean's double-digit Lehmer inner loop.

use super::arith::{scratch, Scratch};
use super::div::{add_raw, compare, divmod_abs, mul_raw, sub_raw, trim_len};
use super::{nlimbs, Ibz};
use core::cmp::Ordering;

/// Lehmer's inner loop on the top two limbs; returns the number of steps
/// and the cofactor matrix.
fn lehmer_inner(x: u128, y: u128) -> (i32, [i64; 4]) {
    let (mut a, mut b, mut c, mut d): (i128, i128, i128, i128) = (1, 0, 0, 1);
    let (mut xx, mut yy) = (x, y);
    const BOUND: i128 = 1i128 << 62;
    const Q_BOUND: u128 = 1u128 << 64;
    const Q_SMALL: u32 = 32;
    let mut steps = 0;
    loop {
        if yy == 0 {
            break;
        }
        let (q, rem) = {
            let mut r = xx;
            let mut cnt = 0u32;
            while cnt < Q_SMALL && r >= yy {
                r -= yy;
                cnt += 1;
            }
            if r < yy {
                (cnt as u128, r)
            } else {
                let extra_q = r / yy;
                (Q_SMALL as u128 + extra_q, r - extra_q * yy)
            }
        };
        if q >= Q_BOUND {
            break;
        }
        let w = q as i128;
        let nc = a - w * c;
        let nd = b - w * d;
        if nc > BOUND || nc < -BOUND || nd > BOUND || nd < -BOUND {
            break;
        }
        let (na, nb) = (c, d);
        let (nx, ny) = (yy, rem);
        let v_diff = (nd - d).abs();
        let abs_nd = nd.abs();
        if !(ny >= abs_nd as u128 && nx - ny >= v_diff as u128) {
            break;
        }
        a = na;
        b = nb;
        c = nc;
        d = nd;
        xx = nx;
        yy = ny;
        steps += 1;
        if b == 0 {
            break;
        }
    }
    (steps, [a as i64, b as i64, c as i64, d as i64])
}

/// `result = a * u + b * v` for signed single-limb `a, b` and non-negative
/// arrays `u, v`; the result is known to be non-negative.
fn lehmer_combine(
    result: &mut Scratch,
    reslen: usize,
    a: i64,
    u: &[u64],
    len1: usize,
    b: i64,
    v: &[u64],
    len2: usize,
) {
    let neg1 = a < 0;
    let neg2 = b < 0;
    debug_assert!(!(neg1 && neg2));
    let abs1 = a.unsigned_abs();
    let abs2 = b.unsigned_abs();
    let mut term1 = scratch();
    let mut term2 = scratch();
    let mut tmp1 = scratch();
    let mut tmp2 = scratch();
    mul_raw(&mut tmp1, u, len1, &[abs1], 1);
    mul_raw(&mut tmp2, v, len2, &[abs2], 1);
    let c1 = (len1 + 1).min(reslen);
    let c2 = (len2 + 1).min(reslen);
    term1[..c1].copy_from_slice(&tmp1[..c1]);
    term2[..c2].copy_from_slice(&tmp2[..c2]);
    if neg1 == neg2 {
        let carry = add_raw(result, reslen, &term1, reslen, &term2, reslen);
        debug_assert_eq!(carry, 0);
    } else {
        let (pos, neg) = if neg1 {
            (&term2, &term1)
        } else {
            (&term1, &term2)
        };
        let borrow = sub_raw(result, reslen, pos, reslen, neg, reslen);
        debug_assert_eq!(borrow, 0);
    }
}

/// Shared state of the Lehmer loops: the two magnitudes with their active
/// lengths.
struct LehmerState {
    u: Scratch,
    v: Scratch,
    l: usize,
    ulen: usize,
    vlen: usize,
}

impl LehmerState {
    fn new(abs_a: &[u64], abs_b: &[u64], l: usize) -> (Self, bool) {
        let mut u = scratch();
        let mut v = scratch();
        let a_first = compare(abs_a, abs_b, l) != Ordering::Less;
        if a_first {
            u[..l].copy_from_slice(&abs_a[..l]);
            v[..l].copy_from_slice(&abs_b[..l]);
        } else {
            u[..l].copy_from_slice(&abs_b[..l]);
            v[..l].copy_from_slice(&abs_a[..l]);
        }
        let ulen = trim_len(&u, l);
        let vlen = trim_len(&v, l);
        (
            Self {
                u,
                v,
                l,
                ulen,
                vlen,
            },
            a_first,
        )
    }

    fn done(&self) -> bool {
        self.v[self.vlen - 1] == 0
    }

    /// Either a Lehmer step (returns `Some(matrix)`) or a full division step
    /// (returns `None` and the quotient in `q` with its length).
    fn step(&mut self, q: &mut Scratch) -> Option<[i64; 4]> {
        let (steps, m) = if self.ulen - self.vlen >= 2 {
            (0, [0; 4])
        } else {
            let has_lo = self.ulen >= 2;
            let x = ((self.u[self.ulen - 1] as u128) << 64)
                | if has_lo {
                    self.u[self.ulen - 2] as u128
                } else {
                    0
                };
            let y = ((self.v[self.ulen - 1] as u128) << 64)
                | if has_lo {
                    self.v[self.ulen - 2] as u128
                } else {
                    0
                };
            lehmer_inner(x, y)
        };
        if steps == 0 {
            let mut rbuf = scratch();
            divmod_abs(q, &mut rbuf, &self.u, self.ulen, &self.v, self.vlen, true);
            let v_copy: Scratch = self.v;
            self.u = scratch();
            self.u[..self.vlen].copy_from_slice(&v_copy[..self.vlen]);
            self.v = scratch();
            self.v[..self.vlen].copy_from_slice(&rbuf[..self.vlen]);
            None
        } else {
            let reslen = self.ulen.max(self.vlen) + 2;
            let mut new_u = scratch();
            let mut new_v = scratch();
            lehmer_combine(
                &mut new_u, reslen, m[0], &self.u, self.ulen, m[1], &self.v, self.vlen,
            );
            lehmer_combine(
                &mut new_v, reslen, m[2], &self.u, self.ulen, m[3], &self.v, self.vlen,
            );
            let ccount = reslen.min(self.l);
            self.u = scratch();
            self.v = scratch();
            self.u[..ccount].copy_from_slice(&new_u[..ccount]);
            self.v[..ccount].copy_from_slice(&new_v[..ccount]);
            Some(m)
        }
    }

    fn retrim(&mut self) {
        self.ulen = trim_len(&self.u, self.l);
        self.vlen = trim_len(&self.v, self.l);
        debug_assert!(compare(&self.u, &self.v, self.l) != Ordering::Less);
    }
}

impl<const N: usize> Ibz<N> {
    /// Signed single-limb coefficient times an integer (variable time).
    fn mul_small_signed(coef: i64, x: &Self) -> Self {
        if coef == 0 || x.is_zero() {
            return Self::zero();
        }
        let x_neg = x.neg_mask();
        let result_neg = ((coef < 0) as u64 ^ x_neg) & 1;
        let abs_coef = coef.unsigned_abs();
        let abs_x = x.cneg(x_neg);
        let xlen = abs_x.n();
        let mut result = scratch();
        mul_raw(&mut result, &abs_x.limbs, xlen, &[abs_coef], 1);
        let mut result_bitlen = abs_x.bitlen + (64 - abs_coef.leading_zeros() as i32) - 1;
        result_bitlen = result_bitlen.min(Self::MAX_BITS - 1);
        let mut prod = Self::from_bits(&result[..xlen + 1], result_bitlen);
        if result_neg != 0 {
            prod = prod.neg();
        }
        prod
    }

    /// Greatest common divisor (non-negative; `gcd(0, 0) = 0`). Variable
    /// time.
    pub fn gcd(&self, b: &Self) -> Self {
        if self.is_zero() && b.is_zero() {
            return Self::zero();
        }
        let n = self.bitlen.max(b.bitlen);
        let l = nlimbs(n);
        let abs_a = self.abs();
        let abs_b = b.abs();
        let (mut st, _) = LehmerState::new(&abs_a.limbs, &abs_b.limbs, l);
        while !st.done() {
            let mut q = scratch();
            let _ = st.step(&mut q);
            st.retrim();
        }
        let mut g = Self::zero();
        g.bitlen = n;
        g.limbs[..l].copy_from_slice(&st.u[..l]);
        g
    }

    /// Extended gcd: `(g, u, v)` with `u a + v b = g`. `a` and `b` must be
    /// at least 64 bits below the container bound. Variable time.
    pub fn xgcd(&self, b: &Self) -> (Self, Self, Self) {
        let n = self.bitlen.max(b.bitlen);
        let l = nlimbs(n);
        let a_neg = self.neg_mask();
        let b_neg = b.neg_mask();
        let abs_a = self.cneg(a_neg);
        let abs_b = b.cneg(b_neg);
        let (mut st, a_first) = LehmerState::new(&abs_a.limbs, &abs_b.limbs, l);
        let (mut cu_a, mut cu_b, mut cv_a, mut cv_b) = if a_first {
            (
                Self::set(1, 2),
                Self::set(0, 2),
                Self::set(0, 2),
                Self::set(1, 2),
            )
        } else {
            (
                Self::set(0, 2),
                Self::set(1, 2),
                Self::set(1, 2),
                Self::set(0, 2),
            )
        };
        while !st.done() {
            let mut qbuf = scratch();
            let (ulen, vlen) = (st.ulen, st.vlen);
            match st.step(&mut qbuf) {
                None => {
                    let qlen = trim_len(&qbuf, super::div_qlen_pub(ulen, vlen));
                    let qz = Self::from_digits(&qbuf[..qlen]);
                    let new_cv_a = cu_a.sub(&qz.mul(&cv_a));
                    let new_cv_b = cu_b.sub(&qz.mul(&cv_b));
                    cu_a = cv_a;
                    cu_b = cv_b;
                    cv_a = new_cv_a;
                    cv_b = new_cv_b;
                }
                Some(m) => {
                    let new_cu_a = Self::mul_small_signed(m[0], &cu_a)
                        .add(&Self::mul_small_signed(m[1], &cv_a));
                    let new_cv_a = Self::mul_small_signed(m[2], &cu_a)
                        .add(&Self::mul_small_signed(m[3], &cv_a));
                    let new_cu_b = Self::mul_small_signed(m[0], &cu_b)
                        .add(&Self::mul_small_signed(m[1], &cv_b));
                    let new_cv_b = Self::mul_small_signed(m[2], &cu_b)
                        .add(&Self::mul_small_signed(m[3], &cv_b));
                    cu_a = new_cu_a;
                    cv_a = new_cv_a;
                    cu_b = new_cu_b;
                    cv_b = new_cv_b;
                }
            }
            st.retrim();
        }
        let mut g = Self::zero();
        g.bitlen = n;
        g.limbs[..l].copy_from_slice(&st.u[..l]);
        cu_a.bitlen = cu_a.bitsize() + 1;
        cu_b.bitlen = cu_b.bitsize() + 1;
        (g, cu_a.cneg(a_neg), cu_b.cneg(b_neg))
    }

    /// Modular inverse in `[0, m)`, `None` if it does not exist. `m` must be
    /// positive. Variable time.
    pub fn invmod(&self, m: &Self) -> Option<Self> {
        debug_assert!(m.is_positive() && !m.is_zero());
        let n = self.bitlen.max(m.bitlen);
        let l = nlimbs(n);
        let a_neg = self.neg_mask();
        let abs_a = self.cneg(a_neg);
        let (mut st, a_first) = LehmerState::new(&abs_a.limbs, &m.limbs, l);
        let (mut cu, mut cv) = if a_first {
            (Self::set(1, 2), Self::set(0, 2))
        } else {
            (Self::set(0, 2), Self::set(1, 2))
        };
        while !st.done() {
            let mut qbuf = scratch();
            let (ulen, vlen) = (st.ulen, st.vlen);
            match st.step(&mut qbuf) {
                None => {
                    let qlen = trim_len(&qbuf, super::div_qlen_pub(ulen, vlen));
                    let qz = Self::from_digits(&qbuf[..qlen]);
                    let new_cv = cu.sub(&qz.mul(&cv));
                    cu = cv;
                    cv = new_cv;
                }
                Some(mm) => {
                    let new_cu =
                        Self::mul_small_signed(mm[0], &cu).add(&Self::mul_small_signed(mm[1], &cv));
                    let new_cv =
                        Self::mul_small_signed(mm[2], &cu).add(&Self::mul_small_signed(mm[3], &cv));
                    cu = new_cu;
                    cv = new_cv;
                }
            }
            st.retrim();
        }
        let coeff = cu.cneg(a_neg);
        let inv = coeff.modulo(m);
        if st.ulen == 1 && st.u[0] == 1 {
            Some(inv)
        } else {
            None
        }
    }

    /// Chinese remainder: `x` in `[0, m1 m2)` with `x = a1 mod m1`,
    /// `x = a2 mod m2`, for coprime positive moduli. Variable time.
    pub fn crt(a1: &Self, a2: &Self, m1: &Self, m2: &Self) -> Self {
        let u = m1.invmod(m2).expect("coprime moduli");
        let diff_a = a2.sub(a1);
        debug_assert!(diff_a.bitlen + u.bitlen <= Self::MAX_BITS);
        let r = diff_a.mul(&u).modulo(m2);
        let x0 = a1.add(&m1.mul(&r));
        debug_assert!(m1.bitlen + m2.bitlen <= Self::MAX_BITS);
        x0.modulo(&m1.mul(m2))
    }

    /// Jacobi symbol `(a / p)` for odd `p`: `-1`, `0` or `1`. Variable
    /// time.
    pub fn legendre(&self, p: &Self) -> i32 {
        if self.is_zero() {
            return 0;
        }
        let mut u = *self;
        let mut v = *p;
        let mut sym = true;
        while !u.is_one() {
            debug_assert!(!u.is_zero() && !v.is_zero());
            let u4 = u.limbs[0] % 4;
            if u4 == 0 {
                u = u.div_2exp(2);
                continue;
            }
            if u4 == 2 {
                u = u.div_2exp(1);
                if v.limbs[0] % 8 == 3 || v.limbs[0] % 8 == 5 {
                    sym = !sym;
                }
                continue;
            }
            if u >= v {
                u = u.sub(&v);
                continue;
            }
            if u4 == 3 && v.limbs[0] % 4 == 3 {
                sym = !sym;
            }
            core::mem::swap(&mut u, &mut v);
        }
        if sym {
            1
        } else {
            -1
        }
    }

    /// Inverse modulo `2^e` by Newton-Hensel lifting; `None` for even `a`.
    pub(crate) fn invmod_2e(&self, e: i32) -> Option<Self> {
        debug_assert!(e >= 0 && e < Self::MAX_BITS);
        if e == 0 {
            return Some(Self::zero());
        }
        let exists = self.is_odd();
        let mut x = Self::one();
        let mut prec = 1;
        while prec < e {
            let next = (2 * prec).min(e);
            let t = Self::two().sub(&self.mul(&x));
            x = x.mul(&t).mod2exp(next as u32);
            prec = next;
        }
        let inv = x.mod2exp(e as u32);
        if exists {
            Some(inv)
        } else {
            None
        }
    }

    /// Inverse of `[[r1, r2], [s1, s2]]` modulo `2^e`, in place. Returns
    /// whether the matrix was invertible (odd determinant).
    pub fn invmat(r1: &mut Self, r2: &mut Self, s1: &mut Self, s2: &mut Self, e: i32) -> bool {
        let mut a = *r1;
        let mut b = *r2;
        let mut c = *s1;
        let mut d = *s2;
        for v in [&mut a, &mut b, &mut c, &mut d] {
            if v.bitlen > e + 1 {
                v.bitlen = e + 1;
            }
        }
        let det = a.mul(&d).sub(&b.mul(&c));
        let exists = det.is_odd();
        let det_inv = match det.invmod_2e(e) {
            Some(v) => v,
            None => det.mod2exp(e as u32),
        };
        // when the inverse does not exist the outputs are meaningless, as in the reference
        let det_inv = if exists { det_inv } else { Self::zero() };
        *r1 = d.mul(&det_inv).mod2exp(e as u32);
        *r2 = b.neg().mul(&det_inv).mod2exp(e as u32);
        *s1 = c.neg().mul(&det_inv).mod2exp(e as u32);
        *s2 = a.mul(&det_inv).mod2exp(e as u32);
        exists
    }
}
