//!
//! Provides `RepresentInteger`, the core algorithm that finds quaternion
//! elements of prescribed norm in a given order. Also provides setup for
//! the standard maximal order O₀ and functions for creating ideals of
//! specified norm.

use super::algebra::{
    quat_alg_add, quat_alg_elem_is_zero, quat_alg_make_primitive, quat_alg_mul, quat_alg_norm,
    quat_alg_scalar,
};
use super::dim4::ibz_mat_4x4_eval;
use super::ideal::quat_lideal_create;
use super::intbig::{
    ibz_div, ibz_divides, ibz_gcd, ibz_get, ibz_probab_prime, ibz_rand_interval, ibz_sqrt_floor,
    ibz_sqrt_mod_p, Ibz,
};
use super::integers::ibz_cornacchia_prime;
use super::types::{
    IbzVec4, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal, QuatPExtremalMaximalOrder,
    QuatRepresentIntegerParams,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::Rng;
use zeroize::Zeroize;

/// Initialize the standard maximal order O₀ lattice.
///
/// The order has basis `{1, i, (i+j)/2, (1+ij)/2}` with denominator 2.
pub fn quat_lattice_o0_set() -> QuatLattice {
    let mut o0 = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            o0.basis[i][j] = BigInt::from(0);
        }
    }
    o0.denom = BigInt::from(2);
    o0.basis[0][0] = BigInt::from(2);
    o0.basis[1][1] = BigInt::from(2);
    o0.basis[2][2] = BigInt::from(1);
    o0.basis[1][2] = BigInt::from(1);
    o0.basis[3][3] = BigInt::from(1);
    o0.basis[0][3] = BigInt::from(1);
    o0
}

/// Initialize the standard p-extremal maximal order O₀.
///
/// Sets `z = i`, `t = j`, `q = 1`, and the order lattice to the standard O₀.
pub fn quat_lattice_o0_set_extremal() -> QuatPExtremalMaximalOrder {
    let mut o0 = QuatPExtremalMaximalOrder::default();
    o0.z.coord[1] = BigInt::from(1);
    o0.z.denom = BigInt::from(1);
    o0.t.coord[2] = BigInt::from(1);
    o0.t.denom = BigInt::from(1);
    o0.q = 1;
    o0.order = quat_lattice_o0_set();
    o0
}

/// Create an element of a p-extremal maximal order from lattice coordinates.
///
/// Given coefficients `[x, y, z, t]` and order generators, computes
/// `x + y·order.z + z·order.t + t·order.t·order.z`.
pub fn quat_order_elem_create(
    order: &QuatPExtremalMaximalOrder,
    coeffs: &IbzVec4,
    alg: &QuatAlg,
) -> QuatAlgElem {
    let one = Ibz::one();

    // elem = x
    let mut elem = quat_alg_scalar(&coeffs[0], &one);

    // quat_temp = y * order.z
    let mut quat_temp = quat_alg_scalar(&coeffs[1], &one);
    quat_temp = quat_alg_mul(&order.z, &quat_temp, alg);
    elem = quat_alg_add(&elem, &quat_temp);

    // quat_temp = z * order.t
    quat_temp = quat_alg_scalar(&coeffs[2], &one);
    quat_temp = quat_alg_mul(&order.t, &quat_temp, alg);
    elem = quat_alg_add(&elem, &quat_temp);

    // quat_temp = t * order.t * order.z
    quat_temp = quat_alg_scalar(&coeffs[3], &one);
    quat_temp = quat_alg_mul(&order.t, &quat_temp, alg);
    quat_temp = quat_alg_mul(&quat_temp, &order.z, alg);
    elem = quat_alg_add(&elem, &quat_temp);

    elem
}

/// Represent an integer as the norm of a quaternion element in a given order.
///
/// Finds `gamma` in the order such that `nrd(gamma) = n_gamma`. Uses
/// random sampling for the j and ij coordinates, then Cornacchia's algorithm
/// for the 1 and i coordinates.
///
/// Returns `true` if successful, `false` otherwise. On success, `gamma`
/// contains the quaternion element with the requested norm.
pub fn quat_represent_integer(
    gamma: &mut QuatAlgElem,
    n_gamma: &Ibz,
    non_diag: bool,
    params: &QuatRepresentIntegerParams,
    rng: &mut impl Rng,
) -> bool {
    if n_gamma.is_even() {
        return false;
    }

    let standard_order = params.order.q == 1;
    let q = BigInt::from(params.order.q);

    // Adjust norm (multiply by 4 for standard or non_diag orders)
    let adjusted_n_gamma = if non_diag || standard_order {
        n_gamma * &BigInt::from(4)
    } else {
        n_gamma.clone()
    };

    // Compute bound = sqrt(adjusted_n_gamma / (p * q))
    let (sq_bound, _) = ibz_div(&adjusted_n_gamma, &params.algebra.p);
    let sq_bound = &sq_bound - &q;
    let bound = ibz_sqrt_floor(&sq_bound);

    // Compute search space size ~ n_gamma / (p^2 * q)
    let mut temp = &q * &params.algebra.p;
    temp = &temp * &params.algebra.p;
    temp = ibz_sqrt_floor(&temp);
    let (counter_init, _) = ibz_div(&adjusted_n_gamma, &temp);

    let mut counter = counter_init.clone();
    let mut found = false;
    let mut coeffs = IbzVec4::default();

    while !found && counter > Ibz::zero() {
        counter = &counter - &Ibz::one();

        // Sample j-coordinate
        coeffs[2] = ibz_rand_interval(rng, &Ibz::one(), &bound);

        // Compute second bound
        let mut cornacchia_target = &coeffs[2] * &coeffs[2];
        temp = &cornacchia_target * &params.algebra.p;
        temp = &adjusted_n_gamma - &temp;
        let sq_bound2 = &q * &params.algebra.p;
        let (temp2, _) = ibz_div(&temp, &sq_bound2);
        temp = ibz_sqrt_floor(&temp2);

        if temp.is_zero() {
            continue;
        }

        // Sample ij-coordinate
        coeffs[3] = ibz_rand_interval(rng, &Ibz::one(), &temp);

        // Compute Cornacchia target: adjusted_n_gamma - p*(z^2 + q*t^2)
        temp = &coeffs[3] * &coeffs[3];
        temp = &q * &temp;
        cornacchia_target = &cornacchia_target + &temp;
        cornacchia_target = &cornacchia_target * &params.algebra.p;
        cornacchia_target = &adjusted_n_gamma - &cornacchia_target;
        debug_assert!(cornacchia_target > Ibz::zero());

        // Apply Cornacchia: find x, y such that x^2 + q*y^2 = cornacchia_target
        let is_prime = ibz_probab_prime(&cornacchia_target, params.primality_test_iterations) > 0;

        if is_prime {
            let corn_result = ibz_cornacchia_prime(&q, &cornacchia_target);
            if let Some((x, y)) = corn_result {
                coeffs[0] = x;
                coeffs[1] = y;
                found = true;
            }
        }

        if found && non_diag && standard_order {
            if coeffs[0].is_odd() != coeffs[3].is_odd() {
                coeffs.0.swap(0, 1);
            }
            let c0 = ibz_get(&coeffs[0]);
            let c3 = ibz_get(&coeffs[3]);
            let c1 = ibz_get(&coeffs[1]);
            let c2 = ibz_get(&coeffs[2]);
            found = found && ((c0 - c3).rem_euclid(4) == 2) && ((c1 - c2).rem_euclid(4) == 2);
        }

        if found {
            // Create gamma from coefficients
            *gamma = quat_order_elem_create(params.order, &coeffs, params.algebra);

            let (mut prim_coeffs, mut content) =
                quat_alg_make_primitive(gamma, &params.order.order);

            if non_diag || standard_order {
                found = content == BigInt::from(2);
            } else {
                found = content.is_one();
            }

            if found {
                // Set gamma from the primitive coefficients in the order basis
                let new_coord = ibz_mat_4x4_eval(&params.order.order.basis, &prim_coeffs);
                gamma.coord = new_coord;
                gamma.denom = params.order.order.denom.clone();
            }

            prim_coeffs.zeroize();
            super::intbig::ibz_zeroize(&mut content);
        }
    }

    // Scrub the accepted secret coordinates; `gamma` (the output) is the
    // caller's to zeroize. Per-iteration Cornacchia/bound temporaries are left
    // to the zeroing allocator (Tier 3) on std.
    coeffs.zeroize();
    found
}

/// Create a random left ideal of O₀ with prescribed norm.
///
/// If `is_prime` is true, uses a trace-zero element approach.
/// Otherwise, uses `quat_represent_integer` with a prime cofactor.
pub fn quat_sampling_random_ideal_o0_given_norm(
    norm: &Ibz,
    is_prime: bool,
    params: &QuatRepresentIntegerParams,
    prime_cofactor: Option<&Ibz>,
    rng: &mut impl Rng,
) -> Option<QuatLeftIdeal> {
    let mut gen = QuatAlgElem::default();
    let mut found = false;

    if is_prime {
        while !found {
            gen.coord[0] = Ibz::zero();
            let n_minus_1 = norm - &Ibz::one();
            for i in 1..4 {
                gen.coord[i] = ibz_rand_interval(rng, &Ibz::zero(), &n_minus_1);
            }

            let (n_temp, norm_d) = quat_alg_norm(&gen, params.algebra);
            debug_assert!(norm_d.is_one());

            let disc = -&n_temp;
            let disc_mod = super::intbig::ibz_mod(&disc, norm);

            let sqrt_result = ibz_sqrt_mod_p(&disc_mod, norm);
            if let Some(sqrt) = sqrt_result {
                gen.coord[0] = sqrt;
                found = !quat_alg_elem_is_zero(&gen);
            }
        }
        // is_prime path found
    } else {
        let cofactor =
            prime_cofactor.expect("invariant: prime_cofactor required for non-prime norm");
        assert!(!norm.is_zero());
        let n_temp = cofactor * norm;
        found = quat_represent_integer(&mut gen, &n_temp, false, params, rng);
        debug_assert!(!found || !quat_alg_elem_is_zero(&gen));
    }

    // Rerandomize the class
    found = false;
    let mut gen_rerand = QuatAlgElem::default();
    while !found {
        for i in 0..4 {
            gen_rerand.coord[i] = ibz_rand_interval(rng, &Ibz::one(), norm);
        }
        let (n_temp, norm_d) = quat_alg_norm(&gen_rerand, params.algebra);
        debug_assert!(norm_d.is_one());
        let disc = ibz_gcd(&n_temp, norm);
        found = disc.is_one() && !quat_alg_elem_is_zero(&gen_rerand);
    }
    // rerandomize done

    gen = quat_alg_mul(&gen, &gen_rerand, params.algebra);
    let lideal = quat_lideal_create(&gen, norm, &params.order.order, params.algebra)
        .expect("invariant: norm-equation ideal has perfect-square index");
    debug_assert_eq!(&lideal.norm, norm);

    // Scrub the secret connecting/rerandomization elements: they determine the
    // returned ideal (which the caller zeroizes) but are not themselves part of
    // it. On no_std there is no zeroing allocator to catch them on free.
    gen.zeroize();
    gen_rerand.zeroize();

    Some(lideal)
}

/// Change coordinates from the standard quaternion basis to the O₀ basis.
///
/// Given an element `coord/denom` in the basis `{1, i, j, ij}`, computes
/// its coordinates in the O₀ basis `{1, i, (i+j)/2, (1+ij)/2}`.
pub fn quat_change_to_o0_basis(el: &QuatAlgElem) -> IbzVec4 {
    let mut vec = IbzVec4::default();

    // v[2] = 2 * el.coord[2]
    vec[2] = &el.coord[2] + &el.coord[2];
    // v[3] = 2 * el.coord[3]
    vec[3] = &el.coord[3] + &el.coord[3];
    // v[0] = el.coord[0] - el.coord[3]
    vec[0] = &el.coord[0] - &el.coord[3];
    // v[1] = el.coord[1] - el.coord[2]
    vec[1] = &el.coord[1] - &el.coord[2];

    debug_assert!(ibz_divides(&vec[0], &el.denom));
    debug_assert!(ibz_divides(&vec[1], &el.denom));
    debug_assert!(ibz_divides(&vec[2], &el.denom));
    debug_assert!(ibz_divides(&vec[3], &el.denom));

    let (q0, _) = ibz_div(&vec[0], &el.denom);
    vec[0] = q0;
    let (q1, _) = ibz_div(&vec[1], &el.denom);
    vec[1] = q1;
    let (q2, _) = ibz_div(&vec[2], &el.denom);
    vec[2] = q2;
    let (q3, _) = ibz_div(&vec[3], &el.denom);
    vec[3] = q3;

    vec
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::algebra::{quat_alg_elem_equal, quat_alg_init_set_ui};
    use num_bigint::BigInt;

    #[test]
    fn test_o0_set() {
        let o0 = quat_lattice_o0_set();
        assert_eq!(o0.denom, BigInt::from(2));
        assert_eq!(o0.basis[0][0], BigInt::from(2));
        assert_eq!(o0.basis[1][1], BigInt::from(2));
        assert_eq!(o0.basis[2][2], BigInt::from(1));
        assert_eq!(o0.basis[1][2], BigInt::from(1));
        assert_eq!(o0.basis[3][3], BigInt::from(1));
        assert_eq!(o0.basis[0][3], BigInt::from(1));
    }

    #[test]
    fn test_o0_set_extremal() {
        let o0 = quat_lattice_o0_set_extremal();
        assert_eq!(o0.q, 1);
        assert_eq!(o0.z.coord[1], BigInt::from(1));
        assert_eq!(o0.t.coord[2], BigInt::from(1));
    }

    #[test]
    fn test_order_elem_create() {
        let alg = quat_alg_init_set_ui(103);
        let o0 = quat_lattice_o0_set_extremal();
        let coeffs = IbzVec4([
            BigInt::from(1),
            BigInt::from(7),
            BigInt::from(2),
            BigInt::from(-2),
        ]);
        let elem = quat_order_elem_create(&o0, &coeffs, &alg);
        assert_eq!(elem.denom, BigInt::from(1));
        // x + y·i + z·j + t·ji = 1 + 7i + 2j + (-2)ji
        // ji = -ij, so (-2)ji = 2ij → coord = [1, 7, 2, 2]
        assert_eq!(elem.coord[0], BigInt::from(1));
        assert_eq!(elem.coord[1], BigInt::from(7));
        assert_eq!(elem.coord[2], BigInt::from(2));
        assert_eq!(elem.coord[3], BigInt::from(2));
    }

    #[test]
    fn test_change_to_o0_basis() {
        // elem = (2, 7, 1, -4)/2
        let elem = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(2),
                BigInt::from(7),
                BigInt::from(1),
                BigInt::from(-4),
            ]),
            denom: BigInt::from(2),
        };

        let out = quat_change_to_o0_basis(&elem);
        // Verify: O0.basis * out / O0.denom should give back elem
        let o0 = quat_lattice_o0_set();
        let recon_coord = ibz_mat_4x4_eval(&o0.basis, &out);
        let recon = QuatAlgElem {
            coord: recon_coord,
            denom: o0.denom.clone(),
        };
        assert!(quat_alg_elem_equal(&recon, &elem));
    }
}
