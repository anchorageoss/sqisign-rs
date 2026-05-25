//!
//! Implements the L² algorithm of Nguyen-Stehlé for reducing 4-dimensional
//! lattice bases. Uses DPE (double-precision extended) floating point for
//! approximate Gram-Schmidt orthogonalization and exact integer arithmetic
//! for basis updates.
//!
//! Also provides LLL verification routines using exact rational arithmetic,
//! and higher-level applications: ideal basis reduction, ideal multiplication
//! with reduction, and finding prime-norm equivalent ideals.

use super::algebra::quat_alg_conj;
use super::dim4::{ibz_mat_4x4_copy, ibz_mat_4x4_eval, ibz_mat_4x4_scalar_mul, quat_qf_eval};
use super::dpe::{
    dpe_cmp, dpe_cmp_d, dpe_div, dpe_get_z, dpe_mul, dpe_round, dpe_set, dpe_set_d, dpe_set_z,
    dpe_sub, Dpe,
};
use super::ideal::{quat_lideal_class_gram, quat_lideal_mul, quat_lideal_norm};
use super::intbig::{ibz_div, ibz_probab_prime, ibz_rand_interval_minm_m, Ibz};
use super::lattice::{quat_lattice_gram, quat_lattice_mul};
use super::rational::{
    ibq_abs, ibq_add, ibq_cmp, ibq_inv, ibq_mul, ibq_set, ibq_sub, ibq_vec_4_copy_ibz, Ibq,
    IbqMat4x4, IbqVec4,
};
use super::types::{IbzMat4x4, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use rand::Rng;

/// Floating-point LLL parameter `δ̄`. Must satisfy `1/4 < δ̄ ≤ 1`;
/// closer to 1 gives stronger reduction at higher cost.
const DELTABAR: f64 = 0.995;
/// Floating-point LLL parameter `η̄` (size-reduction threshold).
/// Must satisfy `1/2 ≤ η̄ < sqrt(δ̄)`.
const ETABAR: f64 = 0.505;
/// Exact rational LLL `δ = 99/100`, used for the integer Lovász condition.
const DELTA_NUM: i32 = 99;
const DELTA_DENOM: i32 = 100;
/// Exact rational `ε = 1/100`, the gap between `δ̄` and `δ` that absorbs
/// floating-point rounding errors in the Gram-Schmidt approximation.
const EPSILON_NUM: i32 = 1;
const EPSILON_DENOM: i32 = 100;

/// Access entry of symmetric matrix: `G[max(i,j)][min(i,j)]`.
fn sym(g: &IbzMat4x4, i: usize, j: usize) -> &Ibz {
    if i < j {
        &g[j][i]
    } else {
        &g[i][j]
    }
}

fn sym_mut(g: &mut IbzMat4x4, i: usize, j: usize) -> &mut Ibz {
    if i < j {
        &mut g[j][i]
    } else {
        &mut g[i][j]
    }
}

/// In-place L² lattice reduction.
///
/// Given a 4×4 Gram matrix `G` (lower-triangular, symmetric) and a
/// corresponding basis matrix (columns are basis vectors), reduces the
/// basis in-place using the L² algorithm and updates `G` accordingly.
pub fn quat_lll_core(g: &mut IbzMat4x4, basis: &mut IbzMat4x4) {
    // Floating-point variables for Gram-Schmidt and Lovász conditions
    let mut r = [[Dpe::default(); 4]; 4];
    let mut u = [[Dpe::default(); 4]; 4];
    let mut lovasz = [Dpe::default(); 4];

    let delta_bar = dpe_set_d(DELTABAR);

    let mut xf: Dpe;
    let mut tmp_f: Dpe;
    let mut x: Ibz;
    let mut tmp_i: Ibz;

    // Main L² loop
    r[0][0] = dpe_set_z(&g[0][0]);
    let mut kappa: i32 = 1;

    while kappa < 4 {
        let ku = kappa as usize;
        let mut done = false;

        while !done {
            // Recompute the κ-th row of the Cholesky factorization
            for j in 0..=ku {
                r[ku][j] = dpe_set_z(&g[ku][j]);
                for k in 0..j {
                    tmp_f = dpe_mul(&r[ku][k], &u[j][k]);
                    r[ku][j] = dpe_sub(&r[ku][j], &tmp_f);
                }
                if j < ku {
                    u[ku][j] = dpe_div(&r[ku][j], &r[j][j]);
                }
            }

            done = true;
            // Size reduce
            for i in (0..ku).rev() {
                if dpe_cmp_d(&u[ku][i], ETABAR) > 0 || dpe_cmp_d(&u[ku][i], -ETABAR) < 0 {
                    done = false;
                    xf = dpe_set(&u[ku][i]);
                    xf = dpe_round(&xf);
                    x = dpe_get_z(&xf);

                    // Update basis: b_κ ← b_κ - X·b_i
                    for j in 0..4 {
                        tmp_i = &x * &basis[j][i];
                        basis[j][ku] = &basis[j][ku] - &tmp_i;
                    }

                    // Update lower half of the Gram matrix
                    tmp_i = &x * sym(g, ku, i);
                    g[ku][ku] = &g[ku][ku] - &tmp_i;
                    for j in 0..4 {
                        let sij = sym(g, i, j).clone();
                        tmp_i = &x * &sij;
                        let skj = sym_mut(g, ku, j);
                        *skj = &*skj - &tmp_i;
                    }

                    // Update u[kappa][j]
                    #[allow(clippy::needless_range_loop)]
                    for j in 0..i {
                        let xf_copy = xf;
                        tmp_f = dpe_mul(&xf_copy, &u[i][j]);
                        u[ku][j] = dpe_sub(&u[ku][j], &tmp_f);
                    }
                }
            }
        }

        // Check Lovász conditions
        lovasz[0] = dpe_set_z(&g[ku][ku]);
        for i in 1..ku {
            tmp_f = dpe_mul(&u[ku][i - 1], &r[ku][i - 1]);
            lovasz[i] = dpe_sub(&lovasz[i - 1], &tmp_f);
        }

        let mut swap = kappa;
        while swap > 0 {
            let su = swap as usize;
            tmp_f = dpe_mul(&delta_bar, &r[su - 1][su - 1]);
            if dpe_cmp(&tmp_f, &lovasz[su - 1]) < 0 {
                break;
            }
            swap -= 1;
        }

        if kappa != swap {
            let su = swap as usize;
            for j in (su + 1..=ku).rev() {
                for i in 0..4 {
                    // Swap basis[i][j] and basis[i][j-1]
                    let tmp = basis[i][j].clone();
                    basis[i][j] = basis[i][j - 1].clone();
                    basis[i][j - 1] = tmp;

                    if i == j - 1 {
                        // Swap g[i][i] and g[j][j]
                        let tmp = g[i][i].clone();
                        g[i][i] = g[j][j].clone();
                        g[j][j] = tmp;
                    } else if i != j {
                        // Swap sym(G, i, j) and sym(G, i, j-1)
                        let (i1, j1) = if i < j { (j, i) } else { (i, j) };
                        let (i2, j2) = if i < j - 1 { (j - 1, i) } else { (i, j - 1) };
                        let tmp = g[i1][j1].clone();
                        g[i1][j1] = g[i2][j2].clone();
                        g[i2][j2] = tmp;
                    }
                }
            }
            // Copy row u[κ] and r[κ] in swap position
            for i in 0..su {
                u[su][i] = dpe_set(&u[ku][i]);
                r[su][i] = dpe_set(&r[ku][i]);
            }
            r[su][su] = lovasz[su];
            kappa = swap;
        }

        kappa += 1;
    }

    // Fill in the upper half of the Gram matrix
    for i in 0..4 {
        for j in (i + 1)..4 {
            g[i][j] = g[j][i].clone();
        }
    }
}

/// LLL-reduce a lattice basis.
///
/// Computes the Gram matrix of the lattice with respect to the quaternion
/// norm, then applies in-place L² reduction to the basis.
pub fn quat_lattice_lll(lattice: &QuatLattice, alg: &QuatAlg) -> IbzMat4x4 {
    let mut g = quat_lattice_gram(lattice, alg);
    let mut red = ibz_mat_4x4_copy(&lattice.basis);
    quat_lll_core(&mut g, &mut red);
    red
}

/// Set the rational LLL parameters delta and eta from the L² constants.
pub fn quat_lll_set_ibq_parameters() -> (Ibq, Ibq) {
    // Denominators are non-zero constants (2, EPSILON_DENOM, DELTA_DENOM).
    let half =
        ibq_set(&Ibz::one(), &BigInt::from(2)).expect("invariant: denominator 2 is non-zero");
    let epsilon = ibq_set(&BigInt::from(EPSILON_NUM), &BigInt::from(EPSILON_DENOM))
        .expect("invariant: EPSILON_DENOM is non-zero");
    let eta = ibq_add(&half, &epsilon);
    let delta = ibq_set(&BigInt::from(DELTA_NUM), &BigInt::from(DELTA_DENOM))
        .expect("invariant: DELTA_DENOM is non-zero");
    (delta, eta)
}

/// Bilinear form for the quaternion norm.
///
/// Computes `vec0[0]*vec1[0] + vec0[1]*vec1[1] + q*(vec0[2]*vec1[2] + vec0[3]*vec1[3])`.
pub fn quat_lll_bilinear(vec0: &IbqVec4, vec1: &IbqVec4, q: &Ibz) -> Ibq {
    let one = Ibz::one();
    // Denominator is 1, which is non-zero.
    let norm_q = ibq_set(q, &one).expect("invariant: denominator 1 is non-zero");

    let mut sum = ibq_mul(&vec0[0], &vec1[0]);
    let prod = ibq_mul(&vec0[1], &vec1[1]);
    sum = ibq_add(&sum, &prod);

    let mut prod = ibq_mul(&vec0[2], &vec1[2]);
    prod = ibq_mul(&prod, &norm_q);
    sum = ibq_add(&sum, &prod);

    let mut prod = ibq_mul(&vec0[3], &vec1[3]);
    prod = ibq_mul(&prod, &norm_q);
    ibq_add(&sum, &prod)
}

/// Gram-Schmidt orthogonalization (transposed) using exact rational arithmetic.
pub fn quat_lll_gram_schmidt_transposed_with_ibq(mat: &IbzMat4x4, q: &Ibz) -> IbqMat4x4 {
    let mut work = IbqMat4x4::default();

    // Transpose: work[i] = column i of mat, stored as rational vector
    for i in 0..4 {
        work[i] = ibq_vec_4_copy_ibz(&mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
    }

    for i in 0..4 {
        let norm = quat_lll_bilinear(&work[i], &work[i], q);
        // Gram-Schmidt norm of a nonzero vector is always positive.
        debug_assert!(!norm.num.is_zero());
        let inv_norm = ibq_inv(&norm).expect("invariant: Gram-Schmidt norm is non-zero");

        for j in (i + 1)..4 {
            let vec = ibq_vec_4_copy_ibz(&mat[0][j], &mat[1][j], &mat[2][j], &mat[3][j]);
            let b = quat_lll_bilinear(&work[i], &vec, q);
            let coeff = ibq_mul(&inv_norm, &b);

            for k in 0..4 {
                let prod = ibq_mul(&coeff, &work[i][k]);
                work[j][k] = ibq_sub(&work[j][k], &prod);
            }
        }
    }

    work
}

/// Verify that a basis is LLL-reduced with given parameters delta and eta.
pub fn quat_lll_verify(mat: &IbzMat4x4, delta: &Ibq, eta: &Ibq, alg: &QuatAlg) -> bool {
    let ot = quat_lll_gram_schmidt_transposed_with_ibq(mat, &alg.p);
    let mut res = true;

    // Check small bilinear products/norms (size-reducedness)
    for i in 0..4 {
        for j in 0..i {
            let tmp_vec = ibq_vec_4_copy_ibz(&mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
            let b = quat_lll_bilinear(&ot[j], &tmp_vec, &alg.p);
            let norm = quat_lll_bilinear(&ot[j], &ot[j], &alg.p);
            // Gram-Schmidt norm of a nonzero vector is always positive.
            debug_assert!(!norm.num.is_zero());
            let inv = ibq_inv(&norm).expect("invariant: Gram-Schmidt norm is non-zero");
            let mu = ibq_mul(&b, &inv);
            let mu_abs = ibq_abs(&mu);
            res = res && (ibq_cmp(&mu_abs, eta) <= 0);
        }
    }

    // Check Lovász conditions
    for i in 1..4 {
        let tmp_vec = ibq_vec_4_copy_ibz(&mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
        let b = quat_lll_bilinear(&ot[i - 1], &tmp_vec, &alg.p);
        let norm = quat_lll_bilinear(&ot[i - 1], &ot[i - 1], &alg.p);
        // Gram-Schmidt norm of a nonzero vector is always positive.
        debug_assert!(!norm.num.is_zero());
        let inv = ibq_inv(&norm).expect("invariant: Gram-Schmidt norm is non-zero");
        let mu = ibq_mul(&b, &inv);
        // mu^2
        let mu_sq = ibq_mul(&mu, &mu);
        // delta - mu^2
        let delta_minus_mu_sq = ibq_sub(delta, &mu_sq);
        // norm(b_i*)
        let norm_i = quat_lll_bilinear(&ot[i], &ot[i], &alg.p);
        // (delta - mu^2) * norm(b_{i-1}*)
        let threshold = ibq_mul(&norm, &delta_minus_mu_sq);
        res = res && (ibq_cmp(&norm_i, &threshold) >= 0);
    }

    res
}

/// Reduce the basis of a left ideal using LLL, returning the reduced basis
/// and a corrected Gram matrix.
///
/// The Gram matrix is scaled so that `gram[i][i]` represents the reduced
/// norm of the i-th basis vector divided by the ideal norm.
pub fn quat_lideal_reduce_basis(lideal: &QuatLeftIdeal, alg: &QuatAlg) -> (IbzMat4x4, IbzMat4x4) {
    let gram_corrector = &lideal.lattice.denom * &lideal.lattice.denom;
    let mut gram = quat_lideal_class_gram(lideal, alg);
    let mut reduced = ibz_mat_4x4_copy(&lideal.lattice.basis);
    quat_lll_core(&mut gram, &mut reduced);

    // Scale gram by denom^2
    gram = ibz_mat_4x4_scalar_mul(&gram_corrector, &gram);
    // Halve diagonal and zero upper triangle
    for i in 0..4 {
        gram[i][i] = &gram[i][i] >> 1usize;
        for j in (i + 1)..4 {
            gram[i][j] = Ibz::zero();
        }
    }
    (reduced, gram)
}

/// Multiply two left ideals and LLL-reduce the result.
pub fn quat_lideal_lideal_mul_reduced(
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    alg: &QuatAlg,
) -> (QuatLeftIdeal, IbzMat4x4) {
    let mut prod = QuatLeftIdeal {
        lattice: quat_lattice_mul(&lideal1.lattice, &lideal2.lattice, alg),
        parent_order: lideal1.parent_order.clone(),
        ..Default::default()
    };
    quat_lideal_norm(&mut prod)
        .expect("invariant: product of valid ideals has perfect-square index");

    let (red, gram) = quat_lideal_reduce_basis(&prod, alg);
    prod.lattice.basis = red;

    (prod, gram)
}

/// Find an equivalent ideal with prime norm by random sampling.
///
/// Given a left ideal, finds a random element in the ideal whose conjugate
/// produces an equivalent ideal of prime norm.
pub fn quat_lideal_prime_norm_reduced_equivalent(
    lideal: &mut QuatLeftIdeal,
    alg: &QuatAlg,
    primality_num_iter: u32,
    equiv_bound_coeff: i32,
    rng: &mut impl Rng,
) -> bool {
    use crate::quaternion::lattice::quat_lattice_contains;
    let (red, gram) = quat_lideal_reduce_basis(lideal, alg);

    let adjusted_norm = &lideal.lattice.denom * &lideal.lattice.denom;

    let mut equiv_num_iter = 2 * equiv_bound_coeff as i64 + 1;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;

    let mut found = false;
    let mut ctr = 0i64;

    while !found && ctr < equiv_num_iter {
        ctr += 1;

        let mut new_alpha = QuatAlgElem::default();
        new_alpha.coord[0] = ibz_rand_interval_minm_m(rng, equiv_bound_coeff);
        new_alpha.coord[1] = ibz_rand_interval_minm_m(rng, equiv_bound_coeff);
        new_alpha.coord[2] = ibz_rand_interval_minm_m(rng, equiv_bound_coeff);
        new_alpha.coord[3] = ibz_rand_interval_minm_m(rng, equiv_bound_coeff);

        let mut tmp = quat_qf_eval(&gram, &new_alpha.coord);
        let (quotient, remainder) = ibz_div(&tmp, &adjusted_norm);
        debug_assert!(remainder.is_zero());
        tmp = quotient;

        if ibz_probab_prime(&tmp, primality_num_iter) > 0 {
            new_alpha.coord = ibz_mat_4x4_eval(&red, &new_alpha.coord);
            new_alpha.denom = lideal.lattice.denom.clone();
            debug_assert!(quat_lattice_contains(&lideal.lattice, &new_alpha).is_some());

            new_alpha = quat_alg_conj(&new_alpha);
            new_alpha.denom = &new_alpha.denom * &lideal.norm;
            *lideal = quat_lideal_mul(lideal, &new_alpha, alg);
            debug_assert!(ibz_probab_prime(&lideal.norm, primality_num_iter) > 0);

            found = true;
        }
    }

    found
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::algebra::quat_alg_init_set_ui;
    use crate::quaternion::lattice::{quat_lattice_equal, quat_lattice_hnf_inplace};

    #[test]
    fn test_lll_set_ibq_parameters() {
        let (delta, eta) = quat_lll_set_ibq_parameters();
        // delta = 99/100
        assert_eq!(delta.num, BigInt::from(99));
        assert_eq!(delta.den, BigInt::from(100));
        // eta = 1/2 + 1/100 = 51/100
        let expected_eta = ibq_set(&BigInt::from(51), &BigInt::from(100)).unwrap();
        assert_eq!(ibq_cmp(&eta, &expected_eta), 0);
    }

    #[test]
    fn test_lll_bilinear() {
        let v0 = ibq_vec_4_copy_ibz(
            &BigInt::from(1),
            &BigInt::from(2),
            &BigInt::from(3),
            &BigInt::from(4),
        );
        // Invert each element: 1/1, 1/2, 1/3, 1/4
        let mut v0_inv = IbqVec4::default();
        for i in 0..4 {
            v0_inv[i] = ibq_inv(&v0[i]).unwrap();
        }
        let v1 = ibq_vec_4_copy_ibz(
            &BigInt::from(9),
            &BigInt::from(-8),
            &BigInt::from(7),
            &BigInt::from(-6),
        );
        let q = BigInt::from(3);
        let b = quat_lll_bilinear(&v0_inv, &v1, &q);
        // Expected: 15/2
        let expected = ibq_set(&BigInt::from(15), &BigInt::from(2)).unwrap();
        assert_eq!(ibq_cmp(&b, &expected), 0);
    }

    #[test]
    fn test_lll_verify_reduced_basis() {
        let alg = quat_alg_init_set_ui(3);
        let (delta, eta) = quat_lll_set_ibq_parameters();

        // Reduced basis from C test
        let mut mat = IbzMat4x4::default();
        mat[0][0] = BigInt::from(0);
        mat[0][1] = BigInt::from(2);
        mat[0][2] = BigInt::from(3);
        mat[0][3] = BigInt::from(-14);
        mat[1][0] = BigInt::from(2);
        mat[1][1] = BigInt::from(-1);
        mat[1][2] = BigInt::from(-4);
        mat[1][3] = BigInt::from(-8);
        mat[2][0] = BigInt::from(1);
        mat[2][1] = BigInt::from(-2);
        mat[2][2] = BigInt::from(1);
        mat[2][3] = BigInt::from(0);
        mat[3][0] = BigInt::from(1);
        mat[3][1] = BigInt::from(1);
        mat[3][2] = BigInt::from(0);
        mat[3][3] = BigInt::from(7);

        assert!(quat_lll_verify(&mat, &delta, &eta, &alg));
    }

    #[test]
    fn test_lattice_lll() {
        let alg = quat_alg_init_set_ui(103);
        let (delta, eta) = quat_lll_set_ibq_parameters();

        // Set lattice from C test
        let mut lat = QuatLattice {
            denom: BigInt::from(60),
            ..QuatLattice::default()
        };
        for i in 0..4 {
            for j in 0..4 {
                lat.basis[i][j] = BigInt::from(0);
            }
        }
        lat.basis[0][0] = BigInt::from(3);
        lat.basis[1][0] = BigInt::from(7);
        lat.basis[0][1] = BigInt::from(1);
        lat.basis[3][1] = BigInt::from(-6);
        lat.basis[1][2] = BigInt::from(12);
        lat.basis[2][2] = BigInt::from(5);
        lat.basis[0][3] = BigInt::from(-19);
        lat.basis[3][3] = BigInt::from(3);

        quat_lattice_hnf_inplace(&mut lat);

        let red = quat_lattice_lll(&lat, &alg);

        // Verify reduced
        assert!(quat_lll_verify(&red, &delta, &eta, &alg));

        // Verify same lattice
        let mut test = QuatLattice {
            denom: lat.denom.clone(),
            basis: red,
        };
        quat_lattice_hnf_inplace(&mut test);
        assert!(quat_lattice_equal(&test, &lat));
    }
}
