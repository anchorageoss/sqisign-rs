//! Lattices of the quaternion algebra: membership, inclusion, Gram
//! matrices and sampling from a ball.

use super::dim4::qf_eval;
use super::lll::applications::bound_parallelogram;
use super::{Mat4x4, QuatAlg, QuatAlgElem, QuatLattice, Vec4};
use crate::mp::{Ibz, Rng};

impl<const N: usize> QuatLattice<N> {
    /// Divide basis and denominator by their gcd; the denominator becomes
    /// positive.
    pub fn reduce_denom(&self) -> Self {
        let mut gcd = self.basis.gcd();
        gcd = gcd.gcd(&self.denom);
        let (basis, _) = self.basis.scalar_div(&gcd);
        let (denom, _) = self.denom.div(&gcd);
        Self {
            denom: denom.abs(),
            basis,
        }
    }

    /// Whether `elem` lies in the lattice (which must have a full-rank
    /// upper-triangular basis), and its coordinates in the basis.
    pub fn contains(&self, elem: &QuatAlgElem<N>) -> (bool, Vec4<N>) {
        debug_assert!(self.basis.is_triangular());
        debug_assert!((0..4).all(|i| !self.basis.0[i][i].is_zero()));
        let mut res = true;
        let mut elem_vec = elem.coord.scalar_mul(&self.denom);
        let prod = self.basis.scalar_mul(&elem.denom);
        let mut coords = Vec4::zero();
        for i in (0..4).rev() {
            let (q, r) = elem_vec.0[i].div(&prod.0[i][i]);
            coords.0[i] = q;
            res = res && r.is_zero();
            let col = Vec4([prod.0[0][i], prod.0[1][i], prod.0[2][i], prod.0[3][i]]);
            elem_vec = elem_vec.sub(&col.scalar_mul(&coords.0[i]));
        }
        (res, coords)
    }

    /// Whether `self` is included in `overlat` (upper-triangular basis).
    pub fn inclusion(&self, overlat: &Self) -> bool {
        let mut res = true;
        for i in 0..4 {
            let elem = QuatAlgElem {
                denom: self.denom,
                coord: Vec4([
                    self.basis.0[0][i],
                    self.basis.0[1][i],
                    self.basis.0[2][i],
                    self.basis.0[3][i],
                ]),
            };
            res = res && overlat.contains(&elem).0;
        }
        res
    }

    /// Equality of lattices with triangular bases.
    pub fn equals(&self, b: &Self) -> bool {
        self.inclusion(b) && b.inclusion(self)
    }

    /// Gram matrix of the trace form `denom^2 Tr(a conj(b))` on the basis.
    pub fn gram(&self, alg: &QuatAlg<N>) -> Mat4x4<N> {
        let mut g = Mat4x4::zero();
        for i in 0..4 {
            for j in 0..=i {
                let mut acc = Ibz::zero();
                for k in 0..4 {
                    let mut tmp = self.basis.0[k][i].mul(&self.basis.0[k][j]);
                    if k >= 2 {
                        tmp = tmp.mul(&alg.p);
                    }
                    acc = acc.add(&tmp);
                }
                g.0[i][j] = acc.mul(&Ibz::two());
            }
        }
        for i in 0..4 {
            for j in i + 1..4 {
                g.0[i][j] = g.0[j][i];
            }
        }
        g
    }

    /// Sample a uniform non-zero element of norm at most `radius` (relative
    /// to the ideal norm) from a product lattice as output by
    /// [`super::QuatIdeal::mul_o0`]. Returns the element and the norm
    /// divided by the ideal norm, or `None` if the ball is trivial or the
    /// randomness failed.
    pub fn sample_from_ball(
        &self,
        radius: &Ibz<N>,
        alg: &QuatAlg<N>,
        rng: &mut impl Rng,
    ) -> Option<(QuatAlgElem<N>, Ibz<N>)> {
        debug_assert!(*radius > Ibz::zero());
        let lat = *self;
        let g = super::ideal::product_gram_matrix(&lat, alg);
        let rad = radius.mul_2exp(1);
        let (box_, u) = bound_parallelogram(&lat, &rad, alg)?;
        let mut x = Vec4::zero();
        let mut tmp;
        loop {
            for i in 0..4 {
                if box_.0[i].is_zero() {
                    x.0[i] = Ibz::zero();
                } else {
                    let two_b = box_.0[i].mul_2exp(1);
                    let v = Ibz::rand_interval(&Ibz::zero(), &two_b, rng)?;
                    let mut v = v.sub(&box_.0[i]);
                    v.set_bound(box_.0[i].get_bound());
                    x.0[i] = v;
                }
            }
            let y = u.eval_t(&x);
            tmp = qf_eval(&g, &y);
            let rejected = (tmp.get() % 4 == 0) || tmp > rad;
            if !rejected {
                x = y;
                break;
            }
        }
        debug_assert!(tmp.is_even());
        let ideal_norm = tmp.div_2exp(1);
        let mut res = QuatAlgElem {
            denom: lat.denom,
            coord: lat.basis.eval(&x),
        };
        res.normalize();
        res.denom.set_bound(3);
        for i in 0..4 {
            res.coord.0[i].set_bound(
                (radius.get_bound() + lat.basis.0[0][0].get_bound()) / 2
                    + res.denom.get_bound()
                    + 3,
            );
        }
        let mut ideal_norm = ideal_norm;
        ideal_norm.set_bound(radius.get_bound());
        Some((res, ideal_norm))
    }
}
