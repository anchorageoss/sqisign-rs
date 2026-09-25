//! Two-dimensional integer vectors and matrices, Gaussian integers as
//! `Vec2`, and the modular `2 x 2` operations the protocol uses.

use super::{Mat2x2, Vec2};
use crate::mp::Ibz;

impl<const N: usize> Mat2x2<N> {
    /// Determinant of `[[a11, a12], [a21, a22]]`.
    pub fn det_from(a11: &Ibz<N>, a12: &Ibz<N>, a21: &Ibz<N>, a22: &Ibz<N>) -> Ibz<N> {
        a11.mul(a22).sub(&a12.mul(a21))
    }

    /// `mat * vec`.
    pub fn eval(&self, v: &Vec2<N>) -> Vec2<N> {
        let m = &self.0;
        Vec2([
            m[0][0].mul(&v.0[0]).add(&m[0][1].mul(&v.0[1])),
            m[1][0].mul(&v.0[0]).add(&m[1][1].mul(&v.0[1])),
        ])
    }

    /// `scalar * mat`.
    pub fn scalar_mul(&self, scalar: &Ibz<N>) -> Self {
        let mut out = *self;
        for r in 0..2 {
            for c in 0..2 {
                out.0[r][c] = scalar.mul(&self.0[r][c]);
            }
        }
        out
    }

    /// `a * b`.
    pub fn mul(&self, b: &Self) -> Self {
        let mut sums = Self::zero();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    sums.0[i][j] = sums.0[i][j].add(&self.0[i][k].mul(&b.0[k][j]));
                }
            }
        }
        sums
    }

    /// `(det, adjugate)`: `adjugate / det` is the inverse when `det != 0`.
    /// Returns `None` for a singular matrix (the determinant is then zero).
    pub fn inv_with_det_as_denom(&self) -> (Ibz<N>, Option<Self>) {
        let det = self.0[0][0]
            .mul(&self.0[1][1])
            .sub(&self.0[0][1].mul(&self.0[1][0]));
        if det.is_zero() {
            return (det, None);
        }
        let inv = Self([
            [self.0[1][1], self.0[0][1].neg()],
            [self.0[1][0].neg(), self.0[0][0]],
        ]);
        (det, Some(inv))
    }

    /// `a * b mod m`.
    pub fn mul_mod(&self, b: &Self, m: &Ibz<N>) -> Self {
        let mut sums = Self::zero();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    sums.0[i][j] = sums.0[i][j].add(&self.0[i][k].mul(&b.0[k][j])).modulo(m);
                }
            }
        }
        sums
    }

    /// Inverse modulo `m`; the zero matrix and `false` when the determinant
    /// is not invertible.
    pub fn inv_mod(&self, m: &Ibz<N>) -> (Self, bool) {
        let mut det = self.0[0][0].mul(&self.0[1][1]).modulo(m);
        det = det.sub(&self.0[0][1].mul(&self.0[1][0])).modulo(m);
        let (det, ok) = match det.invmod(m) {
            Some(d) => (d, true),
            None => (Ibz::zero(), false),
        };
        let det = det.mul(&Ibz::set(ok as i64, 2));
        let mut inv = Self([
            [self.0[1][1], self.0[0][1].neg()],
            [self.0[1][0].neg(), self.0[0][0]],
        ]);
        for r in 0..2 {
            for c in 0..2 {
                inv.0[r][c] = inv.0[r][c].mul(&det).modulo(m);
            }
        }
        (inv, ok)
    }

    /// NormalizeMatrix: for a matrix with entries in `[0, 2^e)`, negate
    /// modulo `2^e` so that the first non-zero entry is below `2^(e-1)`.
    pub fn normalize(&mut self, e: u32) {
        let mut tmp = Ibz::one().mul_2exp(e - 1);
        for i in 0..2 {
            for j in 0..2 {
                debug_assert!(self.0[i][j].is_positive());
                if self.0[i][j].is_zero() {
                    continue;
                }
                let r = self.0[i][j].cmp(&tmp);
                if r == core::cmp::Ordering::Less {
                    return;
                }
                if r == core::cmp::Ordering::Greater {
                    tmp = tmp.mul_2exp(1);
                    for k in 0..2 {
                        for l in 0..2 {
                            self.0[k][l] = tmp.sub(&self.0[k][l]);
                            debug_assert!(self.0[k][l].is_positive());
                        }
                    }
                    return;
                }
            }
        }
    }
}

impl<const N: usize> Vec2<N> {
    /// `a + b`.
    pub fn add(&self, b: &Self) -> Self {
        Self([self.0[0].add(&b.0[0]), self.0[1].add(&b.0[1])])
    }

    /// `a - b`.
    pub fn sub(&self, b: &Self) -> Self {
        Self([self.0[0].sub(&b.0[0]), self.0[1].sub(&b.0[1])])
    }

    /// Whether both coordinates are zero.
    pub fn is_zero(&self) -> bool {
        self.0[0].is_zero() && self.0[1].is_zero()
    }

    /// Whether the Gaussian integer is a unit.
    pub fn gaussian_is_unit(&self) -> bool {
        sum_two_squares(&self.0[0], &self.0[1]).is_one()
    }

    /// Product of Gaussian integers `a * b` in `Z[i]`.
    pub fn gaussian_mul(&self, b: &Self) -> Self {
        let r = self.0[0].add(&self.0[1]).mul(&b.0[0]);
        let s = self.0[0].sub(&self.0[1]).mul(&b.0[1]);
        let t = b.0[0].add(&b.0[1]).mul(&self.0[1]);
        Self([r.sub(&t), t.add(&s)])
    }

    /// Euclidean division of Gaussian integers: `(q, r)` with
    /// `a = q b + r`, `N(r) < N(b)`. `b` must be non-zero.
    pub fn gaussian_euclidean_division(&self, b: &Self) -> (Self, Self) {
        let n = sum_two_squares(&b.0[0], &b.0[1]);
        let mut tmp_q = *b;
        tmp_q.0[1] = tmp_q.0[1].neg();
        tmp_q = tmp_q.gaussian_mul(self);
        tmp_q.0[1] = rounded_div(&tmp_q.0[1], &n);
        tmp_q.0[0] = rounded_div(&tmp_q.0[0], &n);
        let prod = tmp_q.gaussian_mul(b);
        let mut tmp_r = self.sub(&prod);
        let mut max_a = self.0[1].get_bound();
        if self.0[0].get_bound() > max_a {
            max_a = self.0[0].get_bound();
        }
        tmp_r.0[0].set_bound(b.0[0].get_bound());
        tmp_r.0[1].set_bound(b.0[1].get_bound());
        tmp_q.0[0].set_bound(max_a);
        tmp_q.0[1].set_bound(max_a);
        #[cfg(debug_assertions)]
        {
            let nr = sum_two_squares(&tmp_r.0[0], &tmp_r.0[1]);
            debug_assert!(nr < n);
            let chk = tmp_q.gaussian_mul(b).add(&tmp_r);
            debug_assert!(chk.0[0] == self.0[0] && chk.0[1] == self.0[1]);
        }
        (tmp_q, tmp_r)
    }

    /// A gcd of two Gaussian integers, not both zero.
    pub fn gaussian_gcd(&self, b: &Self) -> Self {
        let (mut q, mut r) = if self.is_zero() {
            (*b, *self)
        } else {
            (*self, *b)
        };
        while !r.is_zero() {
            let (_, rem) = q.gaussian_euclidean_division(&r);
            q = rem;
            core::mem::swap(&mut q, &mut r);
        }
        q
    }
}

/// `a^2 + b^2`.
pub fn sum_two_squares<const N: usize>(a: &Ibz<N>, b: &Ibz<N>) -> Ibz<N> {
    a.mul(a).add(&b.mul(b))
}

/// `round(a / b)`, ties away from zero.
pub fn rounded_div<const N: usize>(a: &Ibz<N>, b: &Ibz<N>) -> Ibz<N> {
    let abs_b = b.abs();
    let sa = if a.is_positive() { 1i64 } else { -1 };
    let sb = if b.is_positive() { 1i64 } else { -1 };
    let sign_q = Ibz::<N>::set(sa * sb, 2);
    let (q, r) = a.div(b);
    let r2 = r.abs();
    let r2 = r2.add(&r2);
    let sign_neg = sign_q < Ibz::zero();
    let adj = Ibz::<N>::set((if sign_neg { -1 } else { 1 }) * (r2 > abs_b) as i64, 2);
    q.add(&adj)
}
