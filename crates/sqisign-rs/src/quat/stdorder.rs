//! Functions specific to the special maximal order `O0`.

use super::integers::cornacchia_prime;
use super::{QuatAlg, QuatAlgElem, Vec4};
use crate::mp::{Ibz, Rng};

impl<const N: usize> QuatAlgElem<N> {
    /// Coordinates of an element of `O0` in the basis
    /// `(1, i, (i+j)/2, (1+ij)/2)`.
    pub fn to_o0_basis(&self) -> Vec4<N> {
        let mut v = Vec4::zero();
        v.0[2] = self.coord.0[2].mul_2exp(1);
        v.0[3] = self.coord.0[3].mul_2exp(1);
        v.0[0] = self.coord.0[0].sub(&self.coord.0[3]);
        v.0[1] = self.coord.0[1].sub(&self.coord.0[2]);
        for c in v.0.iter_mut() {
            debug_assert!(c.divides(&self.denom));
            let (q, _) = c.div(&self.denom);
            *c = q;
        }
        v
    }

    /// `x mod (m O0)` for `x` in `O0`, with the coordinate bounds set to
    /// `bound(m) + 2`.
    pub fn mod_o0(&self, m: &Ibz<N>, alg: &QuatAlg<N>) -> Self {
        let mut coords = self.to_o0_basis();
        for c in coords.0.iter_mut() {
            *c = c.modulo(m);
        }
        let mut red = Self {
            denom: alg.o0.denom,
            coord: alg.o0.basis.eval(&coords),
        };
        let len_n = m.get_bound();
        for c in red.coord.0.iter_mut() {
            c.set_bound(len_n + 2);
        }
        red
    }
}

/// RepresentInteger: a primitive element of `O0` of odd norm `n_gamma`.
/// `None` on failure (likely for small `n_gamma`) or randomness failure.
pub fn represent_integer<const N: usize>(
    n_gamma: &Ibz<N>,
    alg: &QuatAlg<N>,
    rng: &mut impl Rng,
) -> Option<QuatAlgElem<N>> {
    debug_assert!(n_gamma.is_odd());
    let bit_bound = alg.p.bitsize();
    let p = &alg.p;
    // search space size roughly n_gamma / p
    let tmp = p.mul(p).sqrt_floor();
    let (mut counter, _) = n_gamma.div(&tmp);
    counter.set_bound(bit_bound);
    // first bound sqrt(n_gamma / p - 1)
    let (mut sq_bound, _) = n_gamma.div(p);
    sq_bound = sq_bound.sub(&Ibz::one());
    sq_bound.set_bound(bit_bound);
    let bound = sq_bound.sqrt_floor();
    let mut coeffs = Vec4::<N>::zero();
    while counter != Ibz::zero() {
        counter = counter.sub(&Ibz::one());
        counter.set_bound(bit_bound);
        coeffs.0[2] = Ibz::rand_interval(&Ibz::one(), &bound, rng)?;
        let mut cornacchia_target = coeffs.0[2].mul(&coeffs.0[2]);
        let mut tmp = cornacchia_target.mul(p);
        tmp = n_gamma.sub(&tmp);
        let (mut t2, _) = tmp.div(p);
        t2.set_bound(bit_bound);
        let t2 = t2.sqrt_floor();
        if t2 == Ibz::zero() {
            continue;
        }
        coeffs.0[3] = Ibz::rand_interval(&Ibz::one(), &t2, rng)?;
        let tmp = coeffs.0[3].mul(&coeffs.0[3]);
        cornacchia_target = cornacchia_target.add(&tmp);
        cornacchia_target = cornacchia_target.mul(p);
        cornacchia_target = n_gamma.sub(&cornacchia_target);
        cornacchia_target.set_bound(n_gamma.get_bound());
        debug_assert!(cornacchia_target > Ibz::zero());
        let found = if cornacchia_target.probab_prime(alg.primality_num_iter, rng) {
            cornacchia_prime(&cornacchia_target)
        } else {
            None
        };
        if let Some((x, y)) = found {
            coeffs.0[0] = x;
            coeffs.0[1] = y;
            for c in coeffs.0[..2].iter_mut() {
                c.set_bound(cornacchia_target.get_bound() / 2 + 2);
            }
            let gamma = QuatAlgElem {
                denom: Ibz::set(1, 2),
                coord: coeffs,
            };
            #[cfg(debug_assertions)]
            {
                let (n, d) = gamma.norm(alg);
                debug_assert!(d.is_one() && n == *n_gamma);
                debug_assert!(alg.o0.contains(&gamma).0);
            }
            return Some(gamma);
        }
    }
    None
}
