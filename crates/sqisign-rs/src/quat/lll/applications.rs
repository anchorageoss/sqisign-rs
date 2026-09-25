//! Lattice-reduction entry points the ideal layer uses.

use super::dim4::dual_reduce_ideal;
use super::gaussian::reduce_o0_ideal;
use crate::mp::Ibz;
use crate::quat::{Mat4x4, QuatAlg, QuatIdeal, QuatLattice, Vec4};

/// Reduced basis of an `O0`-ideal (constant time): columns `v1, i v1, v2,
/// i v2` of a Lagrange-reduced `Z[i]`-basis, the same for every equivalent
/// input up to one global `Z[i]` unit. Denominator 2.
pub fn reduce_basis<const N: usize>(ideal: &QuatIdeal<N>, alg: &QuatAlg<N>) -> QuatLattice<N> {
    let ipt = ideal.to_lattice();
    let bound = ideal.norm.get_bound() + 2;
    let mut reduced = reduce_o0_ideal(&ipt, alg, bound);
    debug_assert!(reduced.denom == Ibz::two());
    reduced.denom = Ibz::set(2, 3);
    let pb = alg.p.bitsize();
    let nb = ideal.norm.get_bound();
    for i in 0..4 {
        for j in 0..4 {
            let b = 1 + ((((i < 2) as i32) * pb + nb) / 2);
            debug_assert!(reduced.basis.0[i][j].bitsize() <= b);
            reduced.basis.0[i][j].set_bound_ct(b);
        }
    }
    for j in 0..2 {
        let b = 1 + ((pb / 2 + nb) / 2) + 1;
        debug_assert!(reduced.basis.0[j][0].bitsize() <= b);
        reduced.basis.0[j][0].set_bound_ct(b);
    }
    for j in 2..4 {
        let mut b = 1 + (-pb / 2 + nb) + 1;
        if b < 1 {
            debug_assert!(reduced.basis.0[j][0].is_zero());
            b = 1;
        }
        debug_assert!(reduced.basis.0[j][0].bitsize() <= b);
        reduced.basis.0[j][0].set_bound_ct(b);
    }
    reduced
}

/// A parallelogram containing the ball of the given radius in the product
/// lattice `lat` (shape of [`QuatIdeal::mul_o0`]): `(box, U)` such that the
/// vectors are `(x1 x2 x3 x4) U` with `|x_i| <= box_i`. `None` when only
/// the origin lies in it.
pub fn bound_parallelogram<const N: usize>(
    lat: &QuatLattice<N>,
    radius: &Ibz<N>,
    alg: &QuatAlg<N>,
) -> Option<(Vec4<N>, Mat4x4<N>)> {
    let (_reduced, gram_diag, ainv) = dual_reduce_ideal(lat, alg);
    let mut u = Mat4x4::<N>::zero();
    for i in 0..4 {
        for j in 0..4 {
            u.0[i][j] = ainv.0[j][i];
        }
    }
    let den = lat.basis.0[0][0].mul(&alg.p);
    let mut box_ = Vec4::<N>::zero();
    let mut trivial = true;
    for i in 0..4 {
        let num = gram_diag[i].mul(radius);
        let (q, _) = num.div(&den);
        box_.0[i] = q.sqrt_floor();
        trivial = trivial && box_.0[i].is_zero();
    }
    for i in 0..4 {
        for j in 0..4 {
            let b = u.0[i][j].bitsize() + 1;
            u.0[i][j].set_bound(b);
        }
    }
    for v in box_.0.iter_mut() {
        let b = v.bitsize() + 1;
        v.set_bound(b);
    }
    if trivial {
        None
    } else {
        Some((box_, u))
    }
}
