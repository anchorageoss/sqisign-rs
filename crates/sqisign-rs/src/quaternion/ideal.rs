//!
//! Provides creation, multiplication, addition, intersection, norm computation,
//! generator search, conjugation, right order/transporter, and maximality
//! testing for left ideals of maximal orders in `B_{p,∞}`.

use super::algebra::{quat_alg_elem_is_zero, quat_alg_norm};
use super::dim4::{
    ibz_mat_4x4_eval, ibz_mat_4x4_identity, ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_mul,
    ibz_mat_4x4_scalar_mul, ibz_mat_4x4_transpose, ibz_vec_4_content, ibz_vec_4_set,
};
use super::intbig::{ibz_div, ibz_divides, ibz_gcd, ibz_sqrt, Ibz};
use super::lattice::{
    quat_lattice_add, quat_lattice_alg_elem_mul, quat_lattice_conjugate_without_hnf,
    quat_lattice_equal, quat_lattice_gram, quat_lattice_index, quat_lattice_intersect,
    quat_lattice_mul, quat_lattice_reduce_denom,
};
use super::types::{IbzMat4x4, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal};
use num_traits::{One, Zero};

/// Compute and set the norm of a left ideal from its lattice and parent order.
///
/// Returns `None` if the lattice index is not a perfect square (i.e. the
/// lattice does not represent a valid ideal).
pub fn quat_lideal_norm(lideal: &mut QuatLeftIdeal) -> Option<()> {
    let index = quat_lattice_index(&lideal.lattice, &lideal.parent_order);
    lideal.norm = ibz_sqrt(&index)?;
    Some(())
}

/// Verify that the stored norm of a left ideal is correct.
pub fn quat_lideal_norm_verify(lideal: &QuatLeftIdeal) -> bool {
    let index = quat_lattice_index(&lideal.lattice, &lideal.parent_order);
    if let Some(sqrt_idx) = ibz_sqrt(&index) {
        sqrt_idx == lideal.norm
    } else {
        false
    }
}

pub fn quat_lideal_copy(src: &QuatLeftIdeal) -> QuatLeftIdeal {
    src.clone()
}

/// Create a principal left ideal `O * x` where `O` is a maximal order.
///
/// The element `x` must belong to the order. Returns the ideal with
/// norm set to `nrd(x)`.
pub fn quat_lideal_create_principal(
    x: &QuatAlgElem,
    order: &QuatLattice,
    alg: &QuatAlg,
) -> QuatLeftIdeal {
    debug_assert!(quat_order_is_maximal(order, alg));

    let lattice = quat_lattice_alg_elem_mul(order, x, alg);
    let lattice = quat_lattice_reduce_denom(&lattice);

    let (norm_n, norm_d) = quat_alg_norm(x, alg);
    debug_assert!(norm_d.is_one());

    QuatLeftIdeal {
        lattice,
        norm: norm_n,
        parent_order: order.clone(),
    }
}

/// Create a left ideal `O*x + O*N` where `O` is a maximal order.
///
/// The element `x` must be nonzero and belong to the order.
pub fn quat_lideal_create(
    x: &QuatAlgElem,
    n: &Ibz,
    order: &QuatLattice,
    alg: &QuatAlg,
) -> Option<QuatLeftIdeal> {
    debug_assert!(quat_order_is_maximal(order, alg));
    debug_assert!(!quat_alg_elem_is_zero(x));

    let mut lideal = quat_lideal_create_principal(x, order, alg);

    // O*N lattice
    let on = QuatLattice {
        basis: ibz_mat_4x4_scalar_mul(n, &order.basis),
        denom: order.denom.clone(),
    };

    lideal.lattice = quat_lattice_add(&lideal.lattice, &on);
    lideal.parent_order = order.clone();
    quat_lideal_norm(&mut lideal)?;

    Some(lideal)
}

/// Find a generator for the left ideal.
///
/// Searches for a primitive element `gen` such that the ideal equals
/// `O*gen + O*N(ideal)`, where `gcd(nrd(gen)/N(ideal), N(ideal)) = 1`.
///
/// Returns `Some(gen)` if found (always succeeds for valid ideals).
pub fn quat_lideal_generator(lideal: &QuatLeftIdeal, alg: &QuatAlg) -> Option<QuatAlgElem> {
    let mut int_norm = 0i32;
    loop {
        int_norm += 1;
        for a in -int_norm..=int_norm {
            for b in (-int_norm + a.abs())..=(int_norm - a.abs()) {
                for c in (-int_norm + a.abs() + b.abs())..=(int_norm - a.abs() - b.abs()) {
                    let d = int_norm - a.abs() - b.abs() - c.abs();
                    let vec = ibz_vec_4_set(a, b, c, d);
                    let content = ibz_vec_4_content(&vec);
                    if content.is_one() {
                        let coord = ibz_mat_4x4_eval(&lideal.lattice.basis, &vec);
                        let gen = QuatAlgElem {
                            coord,
                            denom: lideal.lattice.denom.clone(),
                        };
                        let (norm_int, norm_denom) = quat_alg_norm(&gen, alg);
                        debug_assert!(norm_denom.is_one());
                        let (q, r) = ibz_div(&norm_int, &lideal.norm);
                        debug_assert!(r.is_zero());
                        let gcd = ibz_gcd(&lideal.norm, &q);
                        if gcd.is_one() {
                            return Some(gen);
                        }
                    }
                }
            }
        }
    }
}

/// Multiply an ideal on the right by an algebra element: `I * alpha`.
pub fn quat_lideal_mul(
    lideal: &QuatLeftIdeal,
    alpha: &QuatAlgElem,
    alg: &QuatAlg,
) -> QuatLeftIdeal {
    let lattice = quat_lattice_alg_elem_mul(&lideal.lattice, alpha, alg);
    let (norm_alpha, norm_d) = quat_alg_norm(alpha, alg);
    let new_norm = &lideal.norm * &norm_alpha;
    debug_assert!(ibz_divides(&new_norm, &norm_d));
    let (final_norm, _) = ibz_div(&new_norm, &norm_d);

    QuatLeftIdeal {
        lattice,
        norm: final_norm,
        parent_order: lideal.parent_order.clone(),
    }
}

/// Intersection of two left ideals: `I1 ∩ I2`.
///
/// Both ideals must share the same parent order.
pub fn quat_lideal_inter(i1: &QuatLeftIdeal, i2: &QuatLeftIdeal) -> QuatLeftIdeal {
    let lattice = quat_lattice_intersect(&i1.lattice, &i2.lattice);
    let mut result = QuatLeftIdeal {
        lattice,
        norm: Ibz::zero(),
        parent_order: i1.parent_order.clone(),
    };
    quat_lideal_norm(&mut result)
        .expect("invariant: intersection of valid ideals has perfect-square index");
    result
}

/// Test equality of two left ideals.
pub fn quat_lideal_equals(i1: &QuatLeftIdeal, i2: &QuatLeftIdeal) -> bool {
    quat_lattice_equal(&i1.parent_order, &i2.parent_order)
        && i1.norm == i2.norm
        && quat_lattice_equal(&i1.lattice, &i2.lattice)
}

/// Inverse lattice `\overline{I} / N(I)` (not in HNF).
///
/// The inverse of a left ideal `I` of norm `N` is the lattice
/// `{conjugate(x) : x ∈ I} / N`.
pub fn quat_lideal_inverse_lattice_without_hnf(lideal: &QuatLeftIdeal) -> QuatLattice {
    let mut inv = quat_lattice_conjugate_without_hnf(&lideal.lattice);
    inv.denom = &inv.denom * &lideal.norm;
    inv
}

/// Right transporter: `{x ∈ B : I1 * x ⊆ I2}`.
///
/// Both ideals must share the same parent order.
pub fn quat_lideal_right_transporter(
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    alg: &QuatAlg,
) -> QuatLattice {
    let inv = quat_lideal_inverse_lattice_without_hnf(lideal1);
    quat_lattice_mul(&inv, &lideal2.lattice, alg)
}

/// Right order of a left ideal: `O_R(I) = {x ∈ B : I * x ⊆ I}`.
pub fn quat_lideal_right_order(lideal: &QuatLeftIdeal, alg: &QuatAlg) -> QuatLattice {
    quat_lideal_right_transporter(lideal, lideal, alg)
}

/// Gram matrix of the reduced norm form on an ideal, divided by
/// `norm * denom²`.
pub fn quat_lideal_class_gram(lideal: &QuatLeftIdeal, alg: &QuatAlg) -> IbzMat4x4 {
    let mut g = quat_lattice_gram(&lideal.lattice, alg);

    let divisor = {
        let d2 = &lideal.lattice.denom * &lideal.lattice.denom;
        &d2 * &lideal.norm
    };

    for i in 0..4 {
        for j in 0..=i {
            let (q, r) = ibz_div(&g[i][j], &divisor);
            debug_assert!(r.is_zero());
            g[i][j] = q;
        }
    }
    for i in 0..4 {
        for j in (i + 1)..4 {
            g[i][j] = g[j][i].clone();
        }
    }
    g
}

/// Conjugate of a left ideal (not in HNF).
///
/// Returns `(conjugate_ideal, new_parent_order)` where the new parent order
/// is the right order of the original ideal.
pub fn quat_lideal_conjugate_without_hnf(
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) -> (QuatLeftIdeal, QuatLattice) {
    let new_parent_order = quat_lideal_right_order(lideal, alg);
    let conj_lattice = quat_lattice_conjugate_without_hnf(&lideal.lattice);
    let conj = QuatLeftIdeal {
        lattice: conj_lattice,
        norm: lideal.norm.clone(),
        parent_order: new_parent_order.clone(),
    };
    (conj, new_parent_order)
}

/// Compute the discriminant of an order in the quaternion algebra.
///
/// Returns `(ok, disc)` where `ok` is true if the discriminant could be computed.
pub fn quat_order_discriminant(order: &QuatLattice, alg: &QuatAlg) -> (bool, Ibz) {
    let transposed = ibz_mat_4x4_transpose(&order.basis);

    let mut norm_mat = ibz_mat_4x4_identity();
    norm_mat[2][2] = alg.p.clone();
    norm_mat[3][3] = alg.p.clone();
    let norm_mat = ibz_mat_4x4_scalar_mul(&Ibz::from(2), &norm_mat);

    let prod = ibz_mat_4x4_mul(&transposed, &norm_mat);
    let prod = ibz_mat_4x4_mul(&prod, &order.basis);

    let (_, det, _) = ibz_mat_4x4_inv_with_det_as_denom(&prod);

    let d2 = &order.denom * &order.denom;
    let d4 = &d2 * &d2;
    let d8 = &d4 * &d4;

    let (sqr, remainder) = ibz_div(&det, &d8);
    let ok1 = remainder.is_zero();
    let ok2 = if let Some(disc) = ibz_sqrt(&sqr) {
        return (ok1, disc);
    } else {
        false
    };

    (ok1 && ok2, Ibz::zero())
}

/// Test whether a lattice represents a maximal order in the quaternion algebra.
///
/// An order is maximal iff its discriminant equals `p`.
pub fn quat_order_is_maximal(order: &QuatLattice, alg: &QuatAlg) -> bool {
    let (ok, disc) = quat_order_discriminant(order, alg);
    ok && disc == alg.p
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::types::IbzVec4;
    use num_bigint::BigInt;

    // Build the standard maximal order O0 for p=7:
    // O0 has basis (1, i, (i+j)/2, (1+ij)/2) with denom=2
    fn make_maximal_order_p7() -> (QuatAlg, QuatLattice) {
        let alg = QuatAlg::new(&BigInt::from(7));
        // Basis columns (in basis 1,i,j,ij):
        // col0 = (1,0,0,0)*2 = (2,0,0,0)  → represents 1 with denom 2
        // col1 = (0,1,0,0)*2 = (0,2,0,0)  → represents i with denom 2
        // col2 = (0,1,1,0)                 → represents (i+j)/2 with denom 2
        // col3 = (1,0,0,1)                 → represents (1+ij)/2 with denom 2
        let lat = QuatLattice {
            basis: IbzMat4x4([
                [
                    BigInt::from(2),
                    BigInt::from(0),
                    BigInt::from(0),
                    BigInt::from(1),
                ],
                [
                    BigInt::from(0),
                    BigInt::from(2),
                    BigInt::from(1),
                    BigInt::from(0),
                ],
                [
                    BigInt::from(0),
                    BigInt::from(0),
                    BigInt::from(1),
                    BigInt::from(0),
                ],
                [
                    BigInt::from(0),
                    BigInt::from(0),
                    BigInt::from(0),
                    BigInt::from(1),
                ],
            ]),
            denom: BigInt::from(2),
        };
        (alg, lat)
    }

    #[test]
    fn test_order_is_maximal() {
        let (alg, order) = make_maximal_order_p7();
        assert!(quat_order_is_maximal(&order, &alg));
    }

    #[test]
    fn test_order_discriminant() {
        let (alg, order) = make_maximal_order_p7();
        let (ok, disc) = quat_order_discriminant(&order, &alg);
        assert!(ok);
        assert_eq!(disc, BigInt::from(7));
    }

    #[test]
    fn test_lideal_create_principal() {
        let (alg, order) = make_maximal_order_p7();
        // x = 1+i (norm = 2)
        let x = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(1),
        };
        let lideal = quat_lideal_create_principal(&x, &order, &alg);
        assert_eq!(lideal.norm, BigInt::from(2));
    }

    #[test]
    fn test_lideal_create() {
        let (alg, order) = make_maximal_order_p7();
        let x = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(1),
        };
        let lideal = quat_lideal_create(&x, &BigInt::from(2), &order, &alg).unwrap();
        // norm should divide 2
        assert!(ibz_divides(&BigInt::from(2), &lideal.norm));
    }

    #[test]
    fn test_lideal_equals() {
        let (alg, order) = make_maximal_order_p7();
        let x = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(1),
        };
        let i1 = quat_lideal_create(&x, &BigInt::from(2), &order, &alg).unwrap();
        let i2 = quat_lideal_create(&x, &BigInt::from(2), &order, &alg).unwrap();
        assert!(quat_lideal_equals(&i1, &i2));
    }

    #[test]
    fn test_lideal_norm_verify() {
        let (alg, order) = make_maximal_order_p7();
        let x = QuatAlgElem {
            coord: IbzVec4([
                BigInt::from(1),
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(0),
            ]),
            denom: BigInt::from(1),
        };
        let lideal = quat_lideal_create_principal(&x, &order, &alg);
        assert!(quat_lideal_norm_verify(&lideal));
    }
}
