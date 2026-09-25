//! Four-dimensional integer vectors and matrices.

use super::{Mat2x2, Mat4x4, Vec4};
use crate::mp::Ibz;

impl<const N: usize> Vec4<N> {
    /// `a + b`.
    pub fn add(&self, b: &Self) -> Self {
        Self([
            self.0[0].add(&b.0[0]),
            self.0[1].add(&b.0[1]),
            self.0[2].add(&b.0[2]),
            self.0[3].add(&b.0[3]),
        ])
    }

    /// `a - b`.
    pub fn sub(&self, b: &Self) -> Self {
        Self([
            self.0[0].sub(&b.0[0]),
            self.0[1].sub(&b.0[1]),
            self.0[2].sub(&b.0[2]),
            self.0[3].sub(&b.0[3]),
        ])
    }

    /// Whether all coordinates are zero.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|v| v.is_zero())
    }

    /// gcd of the coordinates.
    pub fn content(&self) -> Ibz<N> {
        let mut c = self.0[0].gcd(&self.0[1]);
        c = self.0[2].gcd(&c);
        self.0[3].gcd(&c)
    }

    /// `scalar * vec`.
    pub fn scalar_mul(&self, scalar: &Ibz<N>) -> Self {
        Self([
            self.0[0].mul(scalar),
            self.0[1].mul(scalar),
            self.0[2].mul(scalar),
            self.0[3].mul(scalar),
        ])
    }

    /// Divide every coordinate by `scalar` (rounded towards zero); returns
    /// whether the division was exact everywhere.
    pub fn scalar_div_inplace(&mut self, scalar: &Ibz<N>) -> bool {
        let mut exact = true;
        for v in self.0.iter_mut() {
            let (q, r) = v.div(scalar);
            *v = q;
            exact = exact && r.is_zero();
        }
        exact
    }

    /// `(vec / scalar, exact)`.
    pub fn scalar_div(&self, scalar: &Ibz<N>) -> (Self, bool) {
        let mut q = *self;
        let exact = q.scalar_div_inplace(scalar);
        (q, exact)
    }
}

impl<const N: usize> Mat4x4<N> {
    /// `a * b`.
    pub fn mul(&self, b: &Self) -> Self {
        let mut m = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    m.0[i][j] = m.0[i][j].add(&self.0[i][k].mul(&b.0[k][j]));
                }
            }
        }
        m
    }

    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut m = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                m.0[i][j] = self.0[j][i];
            }
        }
        m
    }

    /// Equality of values.
    pub fn equals(&self, b: &Self) -> bool {
        (0..4).all(|i| (0..4).all(|j| self.0[i][j] == b.0[i][j]))
    }

    /// `scalar * mat`.
    pub fn scalar_mul(&self, scalar: &Ibz<N>) -> Self {
        let mut m = *self;
        for i in 0..4 {
            for j in 0..4 {
                m.0[i][j] = self.0[i][j].mul(scalar);
            }
        }
        m
    }

    /// gcd of all entries.
    pub fn gcd(&self) -> Ibz<N> {
        let mut d = self.0[0][0];
        for i in 0..4 {
            for j in 0..4 {
                d = d.gcd(&self.0[i][j]);
            }
        }
        d
    }

    /// `(mat / scalar, exact)`.
    pub fn scalar_div(&self, scalar: &Ibz<N>) -> (Self, bool) {
        let mut m = *self;
        let mut exact = true;
        for i in 0..4 {
            for j in 0..4 {
                let (q, r) = self.0[i][j].div(scalar);
                m.0[i][j] = q;
                exact = exact && r.is_zero();
            }
        }
        (m, exact)
    }

    /// `mat * vec`.
    pub fn eval(&self, v: &Vec4<N>) -> Vec4<N> {
        let mut s = Vec4::zero();
        for i in 0..4 {
            for j in 0..4 {
                s.0[i] = s.0[i].add(&self.0[i][j].mul(&v.0[j]));
            }
        }
        s
    }

    /// `vec * mat` (row vector on the left).
    pub fn eval_t(&self, v: &Vec4<N>) -> Vec4<N> {
        let mut s = Vec4::zero();
        for i in 0..4 {
            for j in 0..4 {
                s.0[i] = s.0[i].add(&self.0[j][i].mul(&v.0[j]));
            }
        }
        s
    }

    /// Whether the matrix is upper triangular.
    pub fn is_triangular(&self) -> bool {
        self.0[1][0].is_zero()
            && self.0[2][0].is_zero()
            && self.0[2][1].is_zero()
            && self.0[3][0].is_zero()
            && self.0[3][1].is_zero()
            && self.0[3][2].is_zero()
    }

    /// Whether the matrix is in Hermite normal form.
    pub fn is_hnf(&self) -> bool {
        let mut res = true;
        let zero = Ibz::<N>::zero();
        for i in 0..4 {
            for j in 0..i {
                res = res && self.0[i][j].is_zero();
            }
            let mut found = false;
            let mut ind = 0;
            for j in i..4 {
                if found {
                    res = res && self.0[i][j] >= zero;
                    res = res && self.0[i][ind] > self.0[i][j];
                } else if !self.0[i][j].is_zero() {
                    found = true;
                    ind = j;
                    res = res && self.0[i][j] > zero;
                }
            }
        }
        let linestart: i32 = -1;
        for j in 0..4 {
            let mut i = 0;
            while i < 4 && self.0[i][j].is_zero() {
                i += 1;
            }
            if i != 4 {
                res = res && linestart < i as i32;
            }
        }
        res
    }

    /// Determinant and the transposed adjugate (`adj / det` is the inverse).
    /// The adjugate is zero when the determinant is.
    pub fn inv_with_det_as_denom(&self) -> (Ibz<N>, Self) {
        let m = &self.0;
        let mut s = [Ibz::<N>::zero(); 6];
        let mut c = [Ibz::<N>::zero(); 6];
        for i in 0..3 {
            s[i] = Mat2x2::det_from(&m[0][0], &m[0][i + 1], &m[1][0], &m[1][i + 1]);
            c[i] = Mat2x2::det_from(&m[2][0], &m[2][i + 1], &m[3][0], &m[3][i + 1]);
        }
        for i in 0..2 {
            s[3 + i] = Mat2x2::det_from(&m[0][1], &m[0][2 + i], &m[1][1], &m[1][2 + i]);
            c[3 + i] = Mat2x2::det_from(&m[2][1], &m[2][2 + i], &m[3][1], &m[3][2 + i]);
        }
        s[5] = Mat2x2::det_from(&m[0][2], &m[0][3], &m[1][2], &m[1][3]);
        c[5] = Mat2x2::det_from(&m[2][2], &m[2][3], &m[3][2], &m[3][3]);
        let mut det = Ibz::<N>::zero();
        for i in 0..6 {
            let prod = s[i].mul(&c[5 - i]);
            if i != 1 && i != 4 {
                det = det.add(&prod);
            } else {
                det = det.sub(&prod);
            }
        }
        let pmp = |a1: &Ibz<N>, a2: &Ibz<N>, b1: &Ibz<N>, b2: &Ibz<N>, c1: &Ibz<N>, c2: &Ibz<N>| {
            a1.mul(a2).sub(&b1.mul(b2)).add(&c1.mul(c2))
        };
        let mpm = |a1: &Ibz<N>, a2: &Ibz<N>, b1: &Ibz<N>, b2: &Ibz<N>, c1: &Ibz<N>, c2: &Ibz<N>| {
            b1.mul(b2).sub(&a1.mul(a2)).sub(&c1.mul(c2))
        };
        let mut work = Self::zero();
        for j in 0..4usize {
            let j0 = (j == 0) as usize;
            let j1 = (j == 1) as usize;
            let j2 = (j == 2) as usize;
            let j3 = (j == 3) as usize;
            let jg1 = (j > 1) as usize;
            for k in 0..2usize {
                let row = 1 - k;
                let args = (
                    &m[row][j0],
                    &c[6 - j - j0],
                    &m[row][2 - jg1],
                    &c[4 - j - j1],
                    &m[row][3 - j3],
                    &c[3 - j - j1 - j2],
                );
                work.0[j][k] = if (k + j + 1) % 2 == 1 {
                    pmp(args.0, args.1, args.2, args.3, args.4, args.5)
                } else {
                    mpm(args.0, args.1, args.2, args.3, args.4, args.5)
                };
            }
            for k in 2..4usize {
                let row = 3 - (k == 3) as usize;
                let args = (
                    &m[row][j0],
                    &s[6 - j - j0],
                    &m[row][2 - jg1],
                    &s[4 - j - j1],
                    &m[row][3 - j3],
                    &s[3 - j - j1 - j2],
                );
                work.0[j][k] = if (k + j + 1) % 2 == 1 {
                    pmp(args.0, args.1, args.2, args.3, args.4, args.5)
                } else {
                    mpm(args.0, args.1, args.2, args.3, args.4, args.5)
                };
            }
        }
        let res = Ibz::<N>::set(!det.is_zero() as i64, 2);
        (det, work.scalar_mul(&res))
    }
}

/// Quadratic form evaluation `coord^T qf coord`.
pub fn qf_eval<const N: usize>(qf: &Mat4x4<N>, coord: &Vec4<N>) -> Ibz<N> {
    let sum = qf.eval(coord);
    let mut acc = sum.0[0].mul(&coord.0[0]);
    for i in 1..4 {
        acc = acc.add(&sum.0[i].mul(&coord.0[i]));
    }
    acc
}
