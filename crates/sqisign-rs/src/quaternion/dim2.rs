use super::intbig::{ibz_invmod, ibz_mod, Ibz};
use super::types::{IbzMat2x2, IbzVec2};
use num_bigint::BigInt;

pub fn ibz_vec_2_set(a0: i32, a1: i32) -> IbzVec2 {
    IbzVec2([BigInt::from(a0), BigInt::from(a1)])
}

pub fn ibz_mat_2x2_set(a00: i32, a01: i32, a10: i32, a11: i32) -> IbzMat2x2 {
    IbzMat2x2([
        [BigInt::from(a00), BigInt::from(a01)],
        [BigInt::from(a10), BigInt::from(a11)],
    ])
}

pub fn ibz_mat_2x2_copy(src: &IbzMat2x2) -> IbzMat2x2 {
    src.clone()
}

/// Element-wise addition of two 2x2 matrices.
pub fn ibz_mat_2x2_add(a: &IbzMat2x2, b: &IbzMat2x2) -> IbzMat2x2 {
    IbzMat2x2([
        [&a[0][0] + &b[0][0], &a[0][1] + &b[0][1]],
        [&a[1][0] + &b[1][0], &a[1][1] + &b[1][1]],
    ])
}

/// 2x2 determinant from four individual elements: `a11*a22 - a12*a21`.
pub fn ibz_mat_2x2_det_from_ibz(a11: &Ibz, a12: &Ibz, a21: &Ibz, a22: &Ibz) -> Ibz {
    &(a11 * a22) - &(a12 * a21)
}

/// Matrix-vector product `mat * vec`.
pub fn ibz_mat_2x2_eval(mat: &IbzMat2x2, vec: &IbzVec2) -> IbzVec2 {
    IbzVec2([
        &(&mat[0][0] * &vec[0]) + &(&mat[0][1] * &vec[1]),
        &(&mat[1][0] * &vec[0]) + &(&mat[1][1] * &vec[1]),
    ])
}

/// Modular matrix multiplication: `(mat_a * mat_b) mod m`.
///
/// Uses a triple loop with reduction after each inner product accumulation
/// to keep intermediate values bounded.
pub fn ibz_2x2_mul_mod(mat_a: &IbzMat2x2, mat_b: &IbzMat2x2, m: &Ibz) -> IbzMat2x2 {
    let zero = BigInt::from(0);
    let mut res = IbzMat2x2([[zero.clone(), zero.clone()], [zero.clone(), zero]]);
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let prod = &mat_a[i][k] * &mat_b[k][j];
                res[i][j] = ibz_mod(&(&res[i][j] + &prod), m);
            }
        }
    }
    res
}

/// Modular matrix inverse via adjugate and determinant inversion.
///
/// Returns `(inverse_matrix, true)` if the determinant is invertible mod `m`,
/// or `(zero_matrix, false)` otherwise.
///
/// When the determinant is not invertible, returns the zero matrix (all
/// entries zeroed).
pub fn ibz_mat_2x2_inv_mod(mat: &IbzMat2x2, m: &Ibz) -> (IbzMat2x2, bool) {
    // Compute det = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) mod m
    let mut det = ibz_mod(&(&mat[0][0] * &mat[1][1]), m);
    let prod = &mat[0][1] * &mat[1][0];
    det = ibz_mod(&(&det - &prod), m);

    // Try to invert the determinant mod m
    let (res, det_inv) = match ibz_invmod(&det, m) {
        Some(inv) => (true, inv),
        None => (false, BigInt::from(0)),
    };

    // Build the adjugate: swap diagonal, negate off-diagonal
    let mut inv = IbzMat2x2([
        [mat[1][1].clone(), -&mat[0][1]],
        [-&mat[1][0], mat[0][0].clone()],
    ]);

    // Multiply every entry by det_inv (which is 0 if non-invertible) and reduce
    for i in 0..2 {
        for j in 0..2 {
            inv[i][j] = ibz_mod(&(&inv[i][j] * &det_inv), m);
        }
    }

    (inv, res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::Zero;

    #[test]
    fn test_vec_2_set() {
        let vec = ibz_vec_2_set(2, 5);
        assert_eq!(vec[0], BigInt::from(2));
        assert_eq!(vec[1], BigInt::from(5));
    }

    #[test]
    fn test_mat_2x2_set() {
        let mat = ibz_mat_2x2_set(2, 7, -1, 5);
        assert_eq!(mat[0][0], BigInt::from(2));
        assert_eq!(mat[0][1], BigInt::from(7));
        assert_eq!(mat[1][0], BigInt::from(-1));
        assert_eq!(mat[1][1], BigInt::from(5));
    }

    #[test]
    fn test_mat_2x2_copy() {
        let mat = ibz_mat_2x2_set(1, -1, 2, 4);
        let copy = ibz_mat_2x2_copy(&mat);
        assert_eq!(copy[0][0], BigInt::from(1));
        assert_eq!(copy[0][1], BigInt::from(-1));
        assert_eq!(copy[1][0], BigInt::from(2));
        assert_eq!(copy[1][1], BigInt::from(4));
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(mat[i][j], copy[i][j]);
            }
        }
    }

    #[test]
    fn test_mat_2x2_det_from_ibz() {
        // Identity: det = 1
        let det = ibz_mat_2x2_det_from_ibz(
            &BigInt::from(1),
            &BigInt::from(0),
            &BigInt::from(0),
            &BigInt::from(1),
        );
        assert_eq!(det, BigInt::from(1));

        // det([[2,3],[1,-2]]) = -4-3 = -7
        let det = ibz_mat_2x2_det_from_ibz(
            &BigInt::from(2),
            &BigInt::from(3),
            &BigInt::from(1),
            &BigInt::from(-2),
        );
        assert_eq!(det, BigInt::from(-7));

        // det([[0,3],[-1,0]]) = 0-(-3) = 3
        let det = ibz_mat_2x2_det_from_ibz(
            &BigInt::from(0),
            &BigInt::from(3),
            &BigInt::from(-1),
            &BigInt::from(0),
        );
        assert_eq!(det, BigInt::from(3));

        // det([[a,a],[a,a]]) = 0 for any a
        let a = BigInt::from(2);
        let det = ibz_mat_2x2_det_from_ibz(&a, &a, &a, &a);
        assert!(det.is_zero());
    }

    #[test]
    fn test_mat_2x2_eval() {
        // [[1,-1],[2,4]] * [1,-1] = [2, -2]
        let mat = ibz_mat_2x2_set(1, -1, 2, 4);
        let vec = ibz_vec_2_set(1, -1);
        let res = ibz_mat_2x2_eval(&mat, &vec);
        assert_eq!(res[0], BigInt::from(2));
        assert_eq!(res[1], BigInt::from(-2));

        // [[2,-2],[1,3]] * [2,4] = [-4, 14]
        let mat = ibz_mat_2x2_set(2, -2, 1, 3);
        let vec = ibz_vec_2_set(2, 4);
        let res = ibz_mat_2x2_eval(&mat, &vec);
        assert_eq!(res[0], BigInt::from(-4));
        assert_eq!(res[1], BigInt::from(14));
    }

    #[test]
    fn test_2x2_mul_mod() {
        let m = BigInt::from(7);
        let a = ibz_mat_2x2_set(2, -2, 1, 3);
        let b = ibz_mat_2x2_set(5, 3, 4, 1);

        // a * b mod 7
        let prod = ibz_2x2_mul_mod(&a, &b, &m);
        let cmp = ibz_mat_2x2_set(2, 4, 3, 6);
        assert_eq!(prod, cmp);

        // b * a mod 7 (non-commutative)
        let prod = ibz_2x2_mul_mod(&b, &a, &m);
        let cmp = ibz_mat_2x2_set(6, 6, 2, 2);
        assert_eq!(prod, cmp);

        // Squaring: [[2,7],[1,-2]]^2 mod 12 = [[11,0],[0,11]]
        let m = BigInt::from(12);
        let a = ibz_mat_2x2_set(2, 7, 1, -2);
        let sq = ibz_2x2_mul_mod(&a, &a, &m);
        let cmp = ibz_mat_2x2_set(11, 0, 0, 11);
        assert_eq!(sq, cmp);
    }

    #[test]
    fn test_mat_2x2_inv_mod() {
        // Invertible case: m=7, a=[[2,-3],[1,3]]
        let m = BigInt::from(7);
        let a = ibz_mat_2x2_set(2, -3, 1, 3);
        let (inv, ok) = ibz_mat_2x2_inv_mod(&a, &m);
        assert!(ok);
        // Verify: inv * a mod m = identity
        let prod = ibz_2x2_mul_mod(&inv, &a, &m);
        let id = ibz_mat_2x2_set(1, 0, 0, 1);
        assert_eq!(prod, id);

        // Invertible case: m=12, a=[[2,7],[1,-2]]
        let m = BigInt::from(12);
        let a = ibz_mat_2x2_set(2, 7, 1, -2);
        let (inv, ok) = ibz_mat_2x2_inv_mod(&a, &m);
        assert!(ok);
        let prod = ibz_2x2_mul_mod(&a, &inv, &m);
        let id = ibz_mat_2x2_set(1, 0, 0, 1);
        assert_eq!(prod, id);

        // Non-invertible: det=0 mod 25 for [[2,-2],[-1,1]]
        let m = BigInt::from(25);
        let a = ibz_mat_2x2_set(2, -2, -1, 1);
        let (_, ok) = ibz_mat_2x2_inv_mod(&a, &m);
        assert!(!ok);

        // Non-invertible: det=2*(-2)-3*1=-7≡0 mod 7
        let m = BigInt::from(7);
        let a = ibz_mat_2x2_set(2, 3, 1, -2);
        let (_, ok) = ibz_mat_2x2_inv_mod(&a, &m);
        assert!(!ok);

        // Non-invertible: det=2*(-2)-1*1=-5≡0 mod 5, but m=25 and -5 is not invertible mod 25
        let m = BigInt::from(25);
        let a = ibz_mat_2x2_set(2, 1, 1, -2);
        let (_, ok) = ibz_mat_2x2_inv_mod(&a, &m);
        assert!(!ok);
    }
}
