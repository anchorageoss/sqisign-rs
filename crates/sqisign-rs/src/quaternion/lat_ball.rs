//!
//! Samples short vectors from a lattice within a ball of specified radius.
//! Uses LLL reduction of the dual lattice to compute a bounding
//! parallelogram, then rejection-samples from that parallelogram.

use super::algebra::quat_alg_normalize;
use super::dim4::{
    ibz_mat_4x4_eval, ibz_mat_4x4_eval_t, ibz_mat_4x4_identity, ibz_mat_4x4_inv_with_det_as_denom,
    ibz_mat_4x4_scalar_mul, quat_qf_eval,
};
use super::intbig::{ibz_div, ibz_rand_interval, ibz_sqrt_floor, Ibz};
use super::lattice::quat_lattice_gram;
use super::lll::quat_lll_core;
use super::types::{IbzMat4x4, IbzVec4, QuatAlg, QuatAlgElem, QuatLattice};
use num_bigint::BigInt;
use num_traits::Zero;
use rand::Rng;

/// Compute a bounding parallelogram for a ball of given radius in a lattice.
///
/// Given a Gram matrix `G` and a radius, computes a bounding box `box_out`
/// and a transformation matrix `U` such that any lattice point within the
/// ball is contained in the parallelogram defined by `box_out` after
/// applying `U`.
///
/// Returns `(box_out, U, non_trivial)` where `non_trivial` is true
/// if at least one dimension is non-zero.
pub fn quat_lattice_bound_parallelogram(g: &IbzMat4x4, radius: &Ibz) -> (IbzVec4, IbzMat4x4, bool) {
    let (mut dual_g, denom, ok) = ibz_mat_4x4_inv_with_det_as_denom(g);
    debug_assert!(ok);

    let mut u = ibz_mat_4x4_identity();
    quat_lll_core(&mut dual_g, &mut u);

    let mut box_out = IbzVec4::default();
    let mut trivial = true;
    for i in 0..4 {
        let prod = &dual_g[i][i] * radius;
        let (q, _) = ibz_div(&prod, &denom);
        box_out[i] = ibz_sqrt_floor(&q);
        if !box_out[i].is_zero() {
            trivial = false;
        }
    }

    // Compute transpose transformation matrix: U = inv(U)^T * denom
    let (inv_u, det_u, _) = ibz_mat_4x4_inv_with_det_as_denom(&u);
    u = ibz_mat_4x4_scalar_mul(&det_u, &inv_u);
    // U is unitary so det(U) = ±1 and det_u * inv gives the adjugate directly.

    (box_out, u, !trivial)
}

/// Sample a lattice vector of norm at most `radius`.
///
/// Uses rejection sampling from a bounding parallelogram computed via
/// LLL reduction of the dual lattice.
///
/// Returns `Some(element)` on success, `None` if the bounding box only
/// contains the origin.
pub fn quat_lattice_sample_from_ball(
    lattice: &QuatLattice,
    alg: &QuatAlg,
    radius: &Ibz,
    rng: &mut impl Rng,
) -> Option<QuatAlgElem> {
    assert!(radius > &Ibz::zero());
    let g = quat_lattice_gram(lattice, alg);

    // Correct ball radius by the denominator and factor of 2
    let mut rad = radius * &lattice.denom;
    rad = &rad * &lattice.denom;
    rad = &rad * &BigInt::from(2);

    let (box_out, u, ok) = quat_lattice_bound_parallelogram(&g, &rad);
    if !ok {
        return None;
    }

    let mut x = IbzVec4::default();

    loop {
        // Sample vector in bounding box
        for i in 0..4 {
            if box_out[i].is_zero() {
                x[i] = Ibz::zero();
            } else {
                let tmp = &box_out[i] + &box_out[i];
                let raw = ibz_rand_interval(rng, &Ibz::zero(), &tmp);
                x[i] = &raw - &box_out[i];
            }
        }

        // Map to parallelogram: x = x^T * U
        x = ibz_mat_4x4_eval_t(&x, &u);

        // Evaluate quadratic form
        let tmp = quat_qf_eval(&g, &x);
        if !tmp.is_zero() && tmp <= rad {
            break;
        }
    }

    // Evaluate linear combination: res = lattice.basis * x
    let coord = ibz_mat_4x4_eval(&lattice.basis, &x);
    let mut res = QuatAlgElem {
        coord,
        denom: lattice.denom.clone(),
    };
    quat_alg_normalize(&mut res);

    Some(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::algebra::{quat_alg_init_set_ui, quat_alg_norm};
    use crate::quaternion::lattice::quat_lattice_contains;
    use num_bigint::BigInt;

    #[test]
    fn test_sample_from_ball_identity() {
        let mut rng = rand::thread_rng();
        let alg = quat_alg_init_set_ui(11);
        let lattice = QuatLattice {
            basis: ibz_mat_4x4_identity(),
            denom: BigInt::from(1),
        };

        for r in 1..=20 {
            let radius = BigInt::from(r);
            let res = quat_lattice_sample_from_ball(&lattice, &alg, &radius, &mut rng);
            assert!(res.is_some(), "should find a vector for radius {}", r);
            let res = res.unwrap();

            // Check it's in the lattice
            assert!(quat_lattice_contains(&lattice, &res).is_some());

            // Check norm <= radius
            let (norm_n, norm_d) = quat_alg_norm(&res, &alg);
            let bound = &norm_d * &radius;
            assert!(norm_n <= bound, "norm should be <= radius");
        }
    }

    #[test]
    fn test_sample_from_ball_with_denom() {
        let mut rng = rand::thread_rng();
        let alg = quat_alg_init_set_ui(11);
        let lattice = QuatLattice {
            basis: ibz_mat_4x4_identity(),
            denom: BigInt::from(13),
        };

        for r in 1..=20 {
            let radius = BigInt::from(r);
            let res = quat_lattice_sample_from_ball(&lattice, &alg, &radius, &mut rng);
            assert!(res.is_some());
            let res = res.unwrap();
            assert!(quat_lattice_contains(&lattice, &res).is_some());
            let (norm_n, norm_d) = quat_alg_norm(&res, &alg);
            let bound = &norm_d * &radius;
            assert!(norm_n <= bound);
        }
    }
}
