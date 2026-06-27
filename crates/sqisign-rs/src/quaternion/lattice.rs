//!
//! Operations on rank-4 lattices in `B_{p,∞}`, represented by an integer
//! basis matrix divided by a common denominator. Includes HNF reduction,
//! lattice addition (sum), intersection via duality, containment testing,
//! index computation, and Gram matrix computation.

use super::algebra::quat_alg_coord_mul;
use super::dim4::{
    ibz_mat_4x4_copy, ibz_mat_4x4_equal, ibz_mat_4x4_eval, ibz_mat_4x4_gcd,
    ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_scalar_div, ibz_mat_4x4_scalar_mul,
    ibz_mat_4x4_transpose, ibz_vec_4_copy_ibz, ibz_vec_4_scalar_div, ibz_vec_4_scalar_mul,
};
use super::hnf::ibz_mat_4xn_hnf_mod_core;
use super::intbig::{ibz_div, ibz_gcd, Ibz};
use super::types::{IbzMat4x4, IbzVec4, QuatAlg, QuatAlgElem, QuatLattice};
use alloc::vec::Vec;
use num_traits::{Signed, Zero};

/// Reduce the denominator of a lattice by dividing out common factors
/// between the basis entries and the denominator.
pub fn quat_lattice_reduce_denom(lat: &QuatLattice) -> QuatLattice {
    let mut gcd = ibz_mat_4x4_gcd(&lat.basis);
    gcd = ibz_gcd(&gcd, &lat.denom);
    let (basis, _) = ibz_mat_4x4_scalar_div(&gcd, &lat.basis);
    let (denom, _) = ibz_div(&lat.denom, &gcd);
    QuatLattice {
        basis,
        denom: denom.abs(),
    }
}

/// Test equality of two lattices. Reduces both to canonical HNF form
/// before comparing.
pub fn quat_lattice_equal(lat1: &QuatLattice, lat2: &QuatLattice) -> bool {
    let mut a = quat_lattice_reduce_denom(lat1);
    let mut b = quat_lattice_reduce_denom(lat2);
    a.denom = a.denom.abs();
    b.denom = b.denom.abs();
    quat_lattice_hnf_inplace(&mut a);
    quat_lattice_hnf_inplace(&mut b);
    a.denom == b.denom && ibz_mat_4x4_equal(&a.basis, &b.basis)
}

/// Test whether `sublat` is a sublattice of `overlat`.
pub fn quat_lattice_inclusion(sublat: &QuatLattice, overlat: &QuatLattice) -> bool {
    let sum = quat_lattice_add(overlat, sublat);
    quat_lattice_equal(&sum, overlat)
}

/// Conjugation of a lattice (negate rows 1,2,3 of the basis matrix).
///
/// The result is NOT in HNF.
pub fn quat_lattice_conjugate_without_hnf(lat: &QuatLattice) -> QuatLattice {
    let mut basis = ibz_mat_4x4_copy(&lat.basis);
    for row in 1..4 {
        for col in 0..4 {
            basis[row][col] = -&basis[row][col];
        }
    }
    QuatLattice {
        basis,
        denom: lat.denom.clone(),
    }
}

/// Dual lattice (not in HNF).
///
/// Stores `dual_basis = denom * adjugate(basis)^T` with
/// `dual_denom = det(basis)`, representing `denom * basis^{-T}`.
pub fn quat_lattice_dual_without_hnf(lat: &QuatLattice) -> QuatLattice {
    let (inv, det, _ok) = ibz_mat_4x4_inv_with_det_as_denom(&lat.basis);
    let inv_t = ibz_mat_4x4_transpose(&inv);
    let dual_basis = ibz_mat_4x4_scalar_mul(&lat.denom, &inv_t);
    QuatLattice {
        basis: dual_basis,
        denom: det,
    }
}

/// Lattice sum (addition): the smallest lattice containing both `lat1` and `lat2`.
///
/// Computes HNF of the combined generators modulo the GCD of their determinants.
pub fn quat_lattice_add(lat1: &QuatLattice, lat2: &QuatLattice) -> QuatLattice {
    let mut generators = Vec::with_capacity(8);

    // Scale lat2 columns by lat1.denom
    let tmp1 = ibz_mat_4x4_scalar_mul(&lat1.denom, &lat2.basis);
    for j in 0..4 {
        generators.push(IbzVec4([
            tmp1[0][j].clone(),
            tmp1[1][j].clone(),
            tmp1[2][j].clone(),
            tmp1[3][j].clone(),
        ]));
    }
    let (_, det1, _) = ibz_mat_4x4_inv_with_det_as_denom(&tmp1);

    // Scale lat1 columns by lat2.denom
    let tmp2 = ibz_mat_4x4_scalar_mul(&lat2.denom, &lat1.basis);
    for j in 0..4 {
        generators.push(IbzVec4([
            tmp2[0][j].clone(),
            tmp2[1][j].clone(),
            tmp2[2][j].clone(),
            tmp2[3][j].clone(),
        ]));
    }
    let (_, det2, _) = ibz_mat_4x4_inv_with_det_as_denom(&tmp2);

    debug_assert!(!det1.is_zero());
    debug_assert!(!det2.is_zero());
    let detprod = ibz_gcd(&det1, &det2);

    let basis = ibz_mat_4xn_hnf_mod_core(&generators, &detprod);
    let denom = &lat1.denom * &lat2.denom;
    let res = QuatLattice { basis, denom };
    quat_lattice_reduce_denom(&res)
}

/// Lattice intersection via duality: `lat1 ∩ lat2 = (lat1* + lat2*)*`.
pub fn quat_lattice_intersect(lat1: &QuatLattice, lat2: &QuatLattice) -> QuatLattice {
    let dual1 = quat_lattice_dual_without_hnf(lat1);
    let dual2 = quat_lattice_dual_without_hnf(lat2);
    let dual_sum = quat_lattice_add(&dual1, &dual2);
    let mut res = quat_lattice_dual_without_hnf(&dual_sum);
    quat_lattice_hnf_inplace(&mut res);
    res
}

/// Multiply each column of the basis matrix by an algebra element (as coordinates).
///
/// Result is NOT in HNF. This computes `basis * elem` in the quaternion algebra
/// at the coordinate level.
pub fn quat_lattice_mat_alg_coord_mul_without_hnf(
    lat_basis: &IbzMat4x4,
    coord: &IbzVec4,
    alg: &QuatAlg,
) -> IbzMat4x4 {
    let mut prod = IbzMat4x4::default();
    for col in 0..4 {
        let a = ibz_vec_4_copy_ibz(
            &lat_basis[0][col],
            &lat_basis[1][col],
            &lat_basis[2][col],
            &lat_basis[3][col],
        );
        let p = quat_alg_coord_mul(&a, coord, alg);
        for row in 0..4 {
            prod[row][col] = p[row].clone();
        }
    }
    prod
}

/// Right-multiply a lattice by an algebra element, then HNF-reduce.
pub fn quat_lattice_alg_elem_mul(
    lat: &QuatLattice,
    elem: &QuatAlgElem,
    alg: &QuatAlg,
) -> QuatLattice {
    let basis = quat_lattice_mat_alg_coord_mul_without_hnf(&lat.basis, &elem.coord, alg);
    let denom = &lat.denom * &elem.denom;
    let mut result = QuatLattice { basis, denom };
    quat_lattice_hnf_inplace(&mut result);
    result
}

/// Lattice multiplication in the quaternion algebra.
///
/// Computes the product lattice generated by all pairwise products of basis
/// elements from `lat1` and `lat2`.
pub fn quat_lattice_mul(lat1: &QuatLattice, lat2: &QuatLattice, alg: &QuatAlg) -> QuatLattice {
    let mut generators = Vec::with_capacity(16);
    let mut detmat = IbzMat4x4::default();

    for k in 0..4 {
        let elem1 = ibz_vec_4_copy_ibz(
            &lat1.basis[0][k],
            &lat1.basis[1][k],
            &lat1.basis[2][k],
            &lat1.basis[3][k],
        );
        for i in 0..4 {
            let elem2 = ibz_vec_4_copy_ibz(
                &lat2.basis[0][i],
                &lat2.basis[1][i],
                &lat2.basis[2][i],
                &lat2.basis[3][i],
            );
            let elem_res = quat_alg_coord_mul(&elem1, &elem2, alg);
            if k == 0 {
                for j in 0..4 {
                    detmat[i][j] = elem_res[j].clone();
                }
            }
            generators.push(elem_res);
        }
    }

    let (_, det, _) = ibz_mat_4x4_inv_with_det_as_denom(&detmat);
    let det_abs = det.abs();
    let basis = ibz_mat_4xn_hnf_mod_core(&generators, &det_abs);
    let denom = &lat1.denom * &lat2.denom;
    let res = QuatLattice { basis, denom };
    quat_lattice_reduce_denom(&res)
}

/// Test whether algebra element `x` is contained in lattice `lat`.
///
/// Returns `Some(coord)` with the coordinate vector if contained, `None` otherwise.
/// The lattice must be full-rank.
pub fn quat_lattice_contains(lat: &QuatLattice, x: &QuatAlgElem) -> Option<IbzVec4> {
    let (inv, det, _) = ibz_mat_4x4_inv_with_det_as_denom(&lat.basis);
    debug_assert!(!det.is_zero());

    let work_coord = ibz_mat_4x4_eval(&inv, &x.coord);
    let work_coord = ibz_vec_4_scalar_mul(&lat.denom, &work_coord);
    let divisor = &x.denom * &det;
    let (result, ok) = ibz_vec_4_scalar_div(&divisor, &work_coord);

    if ok {
        Some(result)
    } else {
        None
    }
}

/// Compute the index `[overlat : sublat]` (ratio of lattice volumes).
pub fn quat_lattice_index(sublat: &QuatLattice, overlat: &QuatLattice) -> Ibz {
    let (_, det_sub, _) = ibz_mat_4x4_inv_with_det_as_denom(&sublat.basis);
    let over_denom4 = {
        let d2 = &overlat.denom * &overlat.denom;
        &d2 * &d2
    };
    let mut index = &det_sub * &over_denom4;

    let sub_denom4 = {
        let d2 = &sublat.denom * &sublat.denom;
        &d2 * &d2
    };
    let (_, det_over, _) = ibz_mat_4x4_inv_with_det_as_denom(&overlat.basis);
    let divisor = &sub_denom4 * &det_over;

    let (q, r) = ibz_div(&index, &divisor);
    debug_assert!(r.is_zero());
    index = q;

    index.abs()
}

/// Put a lattice into Hermite Normal Form in place.
pub fn quat_lattice_hnf_inplace(lat: &mut QuatLattice) {
    let (_, modulus, _) = ibz_mat_4x4_inv_with_det_as_denom(&lat.basis);
    let modulus = modulus.abs();

    let mut generators = Vec::with_capacity(4);
    for j in 0..4 {
        generators.push(IbzVec4([
            lat.basis[0][j].clone(),
            lat.basis[1][j].clone(),
            lat.basis[2][j].clone(),
            lat.basis[3][j].clone(),
        ]));
    }

    lat.basis = ibz_mat_4xn_hnf_mod_core(&generators, &modulus);
    let reduced = quat_lattice_reduce_denom(lat);
    lat.basis = reduced.basis;
    lat.denom = reduced.denom;
}

/// Put a lattice into Hermite Normal Form (returns new lattice).
pub fn quat_lattice_hnf(lat: &QuatLattice) -> QuatLattice {
    let mut result = lat.clone();
    quat_lattice_hnf_inplace(&mut result);
    result
}

/// Compute the Gram matrix of the reduced norm form on a lattice.
///
/// `G[i][j] = 2 * sum_k (basis[k][i] * basis[k][j] * (k >= 2 ? p : 1))`
///
/// The factor of 2 comes from using the reduced norm (not the bilinear form).
pub fn quat_lattice_gram(lattice: &QuatLattice, alg: &QuatAlg) -> IbzMat4x4 {
    let mut g = IbzMat4x4::default();

    for i in 0..4 {
        for j in 0..=i {
            let mut val = Ibz::zero();
            for k in 0..4 {
                let mut tmp = &lattice.basis[k][i] * &lattice.basis[k][j];
                if k >= 2 {
                    tmp = &tmp * &alg.p;
                }
                val = &val + &tmp;
            }
            g[i][j] = &val * &Ibz::from(2);
        }
    }
    // Symmetrize
    for i in 0..4 {
        for j in (i + 1)..4 {
            g[i][j] = g[j][i].clone();
        }
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn make_order_p7() -> (QuatAlg, QuatLattice) {
        let alg = QuatAlg::new(&BigInt::from(7));
        // Standard maximal order O0 for p=7: basis is identity/1
        // Actually for B_{7,∞}, the standard order has basis (1, i, (i+j)/2, (1+ij)/2)
        // But for simplicity, let's use the order Z[1,i,j,ij] with denom=1
        // That order has discriminant 4p^2, not maximal. For testing lattice ops
        // we can use any full-rank lattice.
        let mut lat = QuatLattice::default();
        for i in 0..4 {
            lat.basis[i][i] = BigInt::from(1);
        }
        lat.denom = BigInt::from(1);
        (alg, lat)
    }

    #[test]
    fn test_reduce_denom() {
        let mut lat = QuatLattice::default();
        for i in 0..4 {
            lat.basis[i][i] = BigInt::from(6);
        }
        lat.denom = BigInt::from(3);
        let reduced = quat_lattice_reduce_denom(&lat);
        assert_eq!(reduced.denom, BigInt::from(1));
        for i in 0..4 {
            assert_eq!(reduced.basis[i][i], BigInt::from(2));
        }
    }

    #[test]
    fn test_lattice_equal() {
        let (_, lat) = make_order_p7();
        assert!(quat_lattice_equal(&lat, &lat));

        // Scale by 2/2, same lattice
        let mut scaled = QuatLattice::default();
        for i in 0..4 {
            scaled.basis[i][i] = BigInt::from(2);
        }
        scaled.denom = BigInt::from(2);
        assert!(quat_lattice_equal(&lat, &scaled));
    }

    #[test]
    fn test_lattice_contains() {
        let (_, lat) = make_order_p7();
        // Element (1,0,0,0)/1 should be contained
        let x = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(1),
        };
        let coord = quat_lattice_contains(&lat, &x);
        assert!(coord.is_some());
        let c = coord.unwrap();
        assert_eq!(c[0], BigInt::from(1));

        // Element (1,0,0,0)/2 should NOT be contained in Z^4
        let y = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(2),
        };
        assert!(quat_lattice_contains(&lat, &y).is_none());
    }

    #[test]
    fn test_lattice_index() {
        let (_, lat) = make_order_p7();
        // sublat = 2 * lat (index should be 2^4 = 16)
        let mut sublat = QuatLattice::default();
        for i in 0..4 {
            sublat.basis[i][i] = BigInt::from(2);
        }
        sublat.denom = BigInt::from(1);
        let idx = quat_lattice_index(&sublat, &lat);
        assert_eq!(idx, BigInt::from(16));
    }

    #[test]
    fn test_lattice_add() {
        let (_, lat) = make_order_p7();
        // lat + lat = lat
        let sum = quat_lattice_add(&lat, &lat);
        assert!(quat_lattice_equal(&sum, &lat));
    }

    #[test]
    fn test_conjugate() {
        let (_, lat) = make_order_p7();
        let conj = quat_lattice_conjugate_without_hnf(&lat);
        // Row 0 unchanged, rows 1-3 negated
        assert_eq!(conj.basis[0][0], BigInt::from(1));
        // For identity matrix, conjugate negates rows 1,2,3
        // So basis[1][1] = -1 etc.
        assert_eq!(conj.basis[1][1], BigInt::from(-1));
    }
}
