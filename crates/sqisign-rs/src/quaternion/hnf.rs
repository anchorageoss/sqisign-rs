//!
//! Implements Cohen's Algorithm 2.4.8 for computing the HNF of a lattice
//! given as column generators modulo a known determinant multiple. Also
//! provides helper functions used internally by the HNF reduction:
//! centered modular reduction, extended GCD with u ≠ 0 guarantee, and
//! modular vector operations.

use super::intbig::{ibz_div, ibz_div_floor, ibz_mod, ibz_xgcd, Ibz};
use super::types::{IbzMat4x4, IbzVec4};
use alloc::vec::Vec;
use num_traits::{One, Signed, Zero};

/// Non-zero modular reduction: `x mod m`, but if the result is zero, returns `m`.
pub fn ibz_mod_not_zero(x: &Ibz, modulus: &Ibz) -> Ibz {
    let m = ibz_mod(x, modulus);
    if m.is_zero() {
        modulus.clone()
    } else {
        m
    }
}

/// Centered modular reduction: result is in `(-m/2, m/2]`.
///
/// Ties go positive (i.e., the result favors positive values).
pub fn ibz_centered_mod(a: &Ibz, modulus: &Ibz) -> Ibz {
    assert!(modulus > &Ibz::zero());
    let (d, _) = ibz_div_floor(modulus, &Ibz::from(2));
    let tmp = ibz_mod_not_zero(a, modulus);
    if tmp > d {
        &tmp - modulus
    } else {
        tmp
    }
}

/// Conditional assignment: if `c` is true, returns `x`; otherwise returns `y`.
pub fn ibz_conditional_assign(x: &Ibz, y: &Ibz, c: bool) -> Ibz {
    if c {
        x.clone()
    } else {
        y.clone()
    }
}

/// Extended GCD with the guarantee that `u ≠ 0` and `x * u > 0` when `x ≠ 0`.
///
/// Returns `(d, u, v)` where `d = gcd(x, y)` and `d = x*u + y*v`.
/// Additional postconditions: `0 < x*u ≤ |y * (x/d)|` when both `x` and `y`
/// are nonzero.
pub fn ibz_xgcd_with_u_not_0(x: &Ibz, y: &Ibz) -> (Ibz, Ibz, Ibz) {
    if x.is_zero() && y.is_zero() {
        return (Ibz::one(), Ibz::one(), Ibz::zero());
    }

    let (mut d, mut u, mut v) = ibz_xgcd(x, y);

    if d.is_negative() {
        d = -d;
        u = -u;
        v = -v;
    }

    // Ensure u != 0
    if u.is_zero() {
        if !x.is_zero() {
            if y.is_zero() {
                let q = ibz_div(x, &Ibz::one()).0;
                v = &v - &q;
            } else {
                let (q, _) = ibz_div(x, y);
                v = &v - &q;
            }
        }
        u = Ibz::one();
    }

    // Ensure x * u > 0 when x != 0
    if !x.is_zero() {
        let r = x * y;
        let neg = r < Ibz::zero();
        let mut prod = x * &u;
        while prod <= Ibz::zero() {
            let (q_y, _) = ibz_div(y, &d);
            let (q_x, _) = ibz_div(x, &d);
            if neg {
                u = &u + &(-&q_y);
                v = &v - &(-&q_x);
            } else {
                u = &u + &q_y;
                v = &v - &q_x;
            }
            prod = x * &u;
        }
    }

    (d, u, v)
}

/// Test whether a 4×4 integer matrix is in Hermite Normal Form.
///
/// Checks: upper-triangular, each pivot (first nonzero entry per row)
/// is positive and strictly larger than all entries to its right.
pub fn ibz_mat_4x4_is_hnf(mat: &IbzMat4x4) -> bool {
    let zero = Ibz::zero();
    let mut res = true;

    for i in 0..4 {
        // Upper triangular check
        for j in 0..i {
            res = res && mat[i][j].is_zero();
        }
        // Find first nonzero element in row
        let mut found = false;
        let mut ind = 0;
        for j in i..4 {
            if found {
                res = res && (mat[i][j] >= zero);
                res = res && (mat[i][ind] > mat[i][j]);
            } else if !mat[i][j].is_zero() {
                found = true;
                ind = j;
                res = res && (mat[i][j] > zero);
            }
        }
    }

    // Check that first nonzero element index per column is strictly increasing
    let linestart: i32 = -1;
    for j in 0..4 {
        let mut i = 0;
        while i < 4 && mat[i][j].is_zero() {
            i += 1;
        }
        if i != 4 {
            res = res && (linestart < i as i32);
        }
    }
    res
}

/// Linear combination with centered modular reduction:
/// `lc[i] = centered_mod(coeff_a * vec_a[i] + coeff_b * vec_b[i], mod)`.
pub fn ibz_vec_4_linear_combination_mod(
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
    modulus: &Ibz,
) -> IbzVec4 {
    let mut lc = IbzVec4::default();
    for i in 0..4 {
        let s = coeff_a * &vec_a[i] + coeff_b * &vec_b[i];
        lc[i] = ibz_centered_mod(&s, modulus);
    }
    lc
}

/// Copy a vector with centered modular reduction on each component.
pub fn ibz_vec_4_copy_mod(vec: &IbzVec4, modulus: &Ibz) -> IbzVec4 {
    let mut res = IbzVec4::default();
    for i in 0..4 {
        res[i] = ibz_centered_mod(&vec[i], modulus);
    }
    res
}

/// Scalar multiplication of a vector with modular reduction:
/// `prod[i] = (scalar * vec[i]) mod m`.
pub fn ibz_vec_4_scalar_mul_mod(scalar: &Ibz, vec: &IbzVec4, modulus: &Ibz) -> IbzVec4 {
    let mut prod = IbzVec4::default();
    for i in 0..4 {
        prod[i] = ibz_mod(&(&vec[i] * scalar), modulus);
    }
    prod
}

/// Hermite Normal Form of a lattice given by column generators modulo a
/// known determinant multiple.
///
/// Implements Cohen's Algorithm 2.4.8. The `generators` slice must contain
/// more than 3 vectors (each representing a column of the generator matrix,
/// stored as rows of `generators`). The `modulus` must be a positive multiple
/// of the lattice volume.
pub fn ibz_mat_4xn_hnf_mod_core(generators: &[IbzVec4], modulus: &Ibz) -> IbzMat4x4 {
    let n = generators.len();
    assert!(n > 3);
    assert!(modulus > &Ibz::zero());

    let mut a: Vec<IbzVec4> = generators.to_vec();
    let mut w = [
        IbzVec4::default(),
        IbzVec4::default(),
        IbzVec4::default(),
        IbzVec4::default(),
    ];
    let mut m = modulus.clone();

    let mut i: i32 = 3;
    let mut k = (n - 1) as i32;
    let mut j = k;

    while i != -1 {
        let iu = i as usize;
        while j != 0 {
            j -= 1;
            let ju = j as usize;
            let ku = k as usize;
            if !a[ju][iu].is_zero() {
                let (d, u, v) = ibz_xgcd_with_u_not_0(&a[ku][iu], &a[ju][iu]);
                let c = {
                    let mut tmp = IbzVec4::default();
                    for idx in 0..4 {
                        tmp[idx] = &u * &a[ku][idx] + &v * &a[ju][idx];
                    }
                    tmp
                };
                let (coeff_1, _) = ibz_div(&a[ku][iu], &d);
                let coeff_2 = {
                    let (q, _) = ibz_div(&a[ju][iu], &d);
                    -q
                };
                a[ju] = ibz_vec_4_linear_combination_mod(&coeff_1, &a[ju], &coeff_2, &a[ku], &m);
                a[ku] = ibz_vec_4_copy_mod(&c, &m);
            }
        }
        let ku = k as usize;
        let (d, u, _v) = ibz_xgcd_with_u_not_0(&a[ku][iu], &m);
        w[iu] = ibz_vec_4_scalar_mul_mod(&u, &a[ku], &m);
        if w[iu][iu].is_zero() {
            w[iu][iu] = m.clone();
        }
        for h in (iu + 1)..4 {
            let (q, _) = ibz_div_floor(&w[h][iu], &w[iu][iu]);
            let neg_q = -q;
            let w_i_clone = w[iu].clone();
            for idx in 0..4 {
                w[h][idx] = &w[h][idx] + &neg_q * &w_i_clone[idx];
            }
        }
        let (new_m, r) = ibz_div(&m, &d);
        assert!(r.is_zero());
        m = new_m;

        if i != 0 {
            k -= 1;
            i -= 1;
            j = k;
            let ku = k as usize;
            let iu = i as usize;
            if a[ku][iu].is_zero() {
                a[ku][iu] = m.clone();
            }
        } else {
            k -= 1;
            i -= 1;
            j = k;
        }
    }

    // Transpose w into hnf: hnf[i][j] = w[j][i]
    let mut hnf = IbzMat4x4::default();
    for j in 0..4 {
        for i in 0..4 {
            hnf[i][j] = w[j][i].clone();
        }
    }
    hnf
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn test_mod_not_zero() {
        assert_eq!(
            ibz_mod_not_zero(&BigInt::from(0), &BigInt::from(5)),
            BigInt::from(5)
        );
        assert_eq!(
            ibz_mod_not_zero(&BigInt::from(3), &BigInt::from(5)),
            BigInt::from(3)
        );
        assert_eq!(
            ibz_mod_not_zero(&BigInt::from(10), &BigInt::from(5)),
            BigInt::from(5)
        );
        assert_eq!(
            ibz_mod_not_zero(&BigInt::from(7), &BigInt::from(5)),
            BigInt::from(2)
        );
    }

    #[test]
    fn test_centered_mod() {
        // mod=7: d=floor(7/2)=3, mod_not_zero maps 0→7, then 7>3→7-7=0
        assert_eq!(
            ibz_centered_mod(&BigInt::from(0), &BigInt::from(7)),
            BigInt::from(0)
        );
        assert_eq!(
            ibz_centered_mod(&BigInt::from(1), &BigInt::from(7)),
            BigInt::from(1)
        );
        assert_eq!(
            ibz_centered_mod(&BigInt::from(4), &BigInt::from(7)),
            BigInt::from(-3)
        );
        assert_eq!(
            ibz_centered_mod(&BigInt::from(3), &BigInt::from(7)),
            BigInt::from(3)
        );
        assert_eq!(
            ibz_centered_mod(&BigInt::from(-1), &BigInt::from(7)),
            BigInt::from(-1)
        );
    }

    #[test]
    fn test_xgcd_with_u_not_0() {
        let (d, u, v) = ibz_xgcd_with_u_not_0(&BigInt::from(12), &BigInt::from(8));
        assert_eq!(d, BigInt::from(4));
        assert_eq!(
            &BigInt::from(12) * &u + &BigInt::from(8) * &v,
            BigInt::from(4)
        );
        assert!(!u.is_zero());
        assert!(&BigInt::from(12) * &u > BigInt::zero());

        // Both zero
        let (d, u, v) = ibz_xgcd_with_u_not_0(&BigInt::from(0), &BigInt::from(0));
        assert_eq!(d, BigInt::from(1));
        assert_eq!(u, BigInt::from(1));
        assert_eq!(v, BigInt::from(0));

        // One zero
        let (d, u, v) = ibz_xgcd_with_u_not_0(&BigInt::from(0), &BigInt::from(5));
        assert_eq!(d, BigInt::from(5));
        assert_eq!(
            &BigInt::from(0) * &u + &BigInt::from(5) * &v,
            BigInt::from(5)
        );

        let (d, u, v) = ibz_xgcd_with_u_not_0(&BigInt::from(5), &BigInt::from(0));
        assert_eq!(d, BigInt::from(5));
        assert_eq!(
            &BigInt::from(5) * &u + &BigInt::from(0) * &v,
            BigInt::from(5)
        );
        assert!(!u.is_zero());
    }

    #[test]
    fn test_hnf_identity() {
        use crate::quaternion::dim4::ibz_mat_4x4_identity;
        // HNF of identity columns (with extra generators)
        let id = ibz_mat_4x4_identity();
        let mut gens = Vec::new();
        // Extract columns as vectors
        for j in 0..4 {
            let mut v = IbzVec4::default();
            for i in 0..4 {
                v[i] = id[i][j].clone();
            }
            gens.push(v);
        }
        // Add an extra generator (column sum) to have > 3 generators
        // Use 5 generators: the 4 identity columns plus one extra
        let mut extra = IbzVec4::default();
        extra[0] = BigInt::from(1);
        extra[1] = BigInt::from(1);
        extra[2] = BigInt::from(1);
        extra[3] = BigInt::from(1);
        gens.push(extra);

        let hnf = ibz_mat_4xn_hnf_mod_core(&gens, &BigInt::from(1));
        assert!(ibz_mat_4x4_is_hnf(&hnf));
        // Should be identity
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert_eq!(hnf[i][j], BigInt::from(1), "hnf[{i}][{j}]");
                } else {
                    assert_eq!(hnf[i][j], BigInt::from(0), "hnf[{i}][{j}]");
                }
            }
        }
    }
}
