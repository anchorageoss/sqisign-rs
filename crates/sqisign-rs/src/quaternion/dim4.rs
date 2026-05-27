use super::dim2::ibz_mat_2x2_det_from_ibz;
use super::intbig::{ibz_div, ibz_gcd, Ibz};
use super::types::{IbzMat4x4, IbzVec4};
use num_bigint::BigInt;
use num_traits::{One, Zero};

pub fn ibz_vec_4_set(c0: i32, c1: i32, c2: i32, c3: i32) -> IbzVec4 {
    IbzVec4([
        BigInt::from(c0),
        BigInt::from(c1),
        BigInt::from(c2),
        BigInt::from(c3),
    ])
}

pub fn ibz_vec_4_copy(src: &IbzVec4) -> IbzVec4 {
    src.clone()
}

/// Construct a 4-vector from four big integer references.
pub fn ibz_vec_4_copy_ibz(c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) -> IbzVec4 {
    IbzVec4([c0.clone(), c1.clone(), c2.clone(), c3.clone()])
}

/// GCD of all four components (the "content" of the vector).
pub fn ibz_vec_4_content(v: &IbzVec4) -> Ibz {
    let d = ibz_gcd(&v[0], &v[1]);
    let d = ibz_gcd(&d, &v[2]);
    ibz_gcd(&d, &v[3])
}

pub fn ibz_vec_4_negate(vec: &IbzVec4) -> IbzVec4 {
    IbzVec4([-&vec[0], -&vec[1], -&vec[2], -&vec[3]])
}

pub fn ibz_vec_4_add(a: &IbzVec4, b: &IbzVec4) -> IbzVec4 {
    IbzVec4([&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2], &a[3] + &b[3]])
}

pub fn ibz_vec_4_sub(a: &IbzVec4, b: &IbzVec4) -> IbzVec4 {
    IbzVec4([&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2], &a[3] - &b[3]])
}

pub fn ibz_vec_4_is_zero(x: &IbzVec4) -> bool {
    x[0].is_zero() && x[1].is_zero() && x[2].is_zero() && x[3].is_zero()
}

/// Linear combination: `coeff_a * vec_a + coeff_b * vec_b`.
pub fn ibz_vec_4_linear_combination(
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
) -> IbzVec4 {
    IbzVec4([
        &(coeff_a * &vec_a[0]) + &(coeff_b * &vec_b[0]),
        &(coeff_a * &vec_a[1]) + &(coeff_b * &vec_b[1]),
        &(coeff_a * &vec_a[2]) + &(coeff_b * &vec_b[2]),
        &(coeff_a * &vec_a[3]) + &(coeff_b * &vec_b[3]),
    ])
}

pub fn ibz_vec_4_scalar_mul(scalar: &Ibz, vec: &IbzVec4) -> IbzVec4 {
    IbzVec4([
        scalar * &vec[0],
        scalar * &vec[1],
        scalar * &vec[2],
        scalar * &vec[3],
    ])
}

/// Component-wise exact division by a scalar.
///
/// Returns `(quotient, true)` if every division is exact (zero remainder),
/// or `(quotient_with_truncation, false)` otherwise.
pub fn ibz_vec_4_scalar_div(scalar: &Ibz, vec: &IbzVec4) -> (IbzVec4, bool) {
    let mut res = true;
    let mut quot = IbzVec4::default();
    for i in 0..4 {
        let (q, r) = ibz_div(&vec[i], scalar);
        quot[i] = q;
        res = res && r.is_zero();
    }
    (quot, res)
}

pub fn ibz_mat_4x4_mul(a: &IbzMat4x4, b: &IbzMat4x4) -> IbzMat4x4 {
    let mut mat = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            let mut acc = Ibz::zero();
            for k in 0..4 {
                acc += &a[i][k] * &b[k][j];
            }
            mat[i][j] = acc;
        }
    }
    mat
}

pub fn ibz_mat_4x4_copy(src: &IbzMat4x4) -> IbzMat4x4 {
    src.clone()
}

pub fn ibz_mat_4x4_negate(mat: &IbzMat4x4) -> IbzMat4x4 {
    let mut res = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            res[i][j] = -&mat[i][j];
        }
    }
    res
}

pub fn ibz_mat_4x4_transpose(mat: &IbzMat4x4) -> IbzMat4x4 {
    let mut res = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            res[i][j] = mat[j][i].clone();
        }
    }
    res
}

pub fn ibz_mat_4x4_zero() -> IbzMat4x4 {
    IbzMat4x4::default()
}

/// 4×4 identity matrix.
pub fn ibz_mat_4x4_identity() -> IbzMat4x4 {
    let mut mat = IbzMat4x4::default();
    for i in 0..4 {
        mat[i][i] = Ibz::one();
    }
    mat
}

/// Test whether the matrix is the identity.
///
/// Diagonal entries must be 1; off-diagonal entries must be 0.
pub fn ibz_mat_4x4_is_identity(mat: &IbzMat4x4) -> bool {
    for i in 0..4 {
        for j in 0..4 {
            if mat[i][j].is_one() != (i == j) {
                return false;
            }
        }
    }
    true
}

pub fn ibz_mat_4x4_equal(a: &IbzMat4x4, b: &IbzMat4x4) -> bool {
    for i in 0..4 {
        for j in 0..4 {
            if a[i][j] != b[i][j] {
                return false;
            }
        }
    }
    true
}

/// Scalar-matrix multiplication.
pub fn ibz_mat_4x4_scalar_mul(scalar: &Ibz, mat: &IbzMat4x4) -> IbzMat4x4 {
    let mut res = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            res[i][j] = scalar * &mat[i][j];
        }
    }
    res
}

/// GCD of all 16 entries.
pub fn ibz_mat_4x4_gcd(mat: &IbzMat4x4) -> Ibz {
    let mut d = mat[0][0].clone();
    for i in 0..4 {
        for j in 0..4 {
            d = ibz_gcd(&d, &mat[i][j]);
        }
    }
    d
}

/// Element-wise exact division by a scalar.
///
/// Returns `(quotient, true)` if every division is exact,
/// or `(quotient_with_truncation, false)` otherwise.
pub fn ibz_mat_4x4_scalar_div(scalar: &Ibz, mat: &IbzMat4x4) -> (IbzMat4x4, bool) {
    let mut res = true;
    let mut quot = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            let (q, r) = ibz_div(&mat[i][j], scalar);
            quot[i][j] = q;
            res = res && r.is_zero();
        }
    }
    (quot, res)
}

/// Cofactor helper: `a1*a2 - b1*b2 + c1*c2`.
pub fn ibz_inv_dim4_make_coeff_pmp(
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) -> Ibz {
    &(a1 * a2) - &(b1 * b2) + &(c1 * c2)
}

/// Cofactor helper: `-a1*a2 + b1*b2 - c1*c2`.
pub fn ibz_inv_dim4_make_coeff_mpm(
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) -> Ibz {
    &(b1 * b2) - &(a1 * a2) - &(c1 * c2)
}

/// 4×4 matrix inverse via Laplace expansion.
///
/// Returns `(adjugate, determinant, invertible)` where `adjugate * mat == det * I`.
/// If the determinant is zero, `adjugate` is the zero matrix and `invertible` is false.
pub fn ibz_mat_4x4_inv_with_det_as_denom(mat: &IbzMat4x4) -> (IbzMat4x4, Ibz, bool) {
    // Compute 2x2 minors from rows 0,1 (s) and rows 2,3 (c).
    let mut s = [
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
    ];
    let mut c = [
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
        Ibz::zero(),
    ];

    // s[i] = det of columns (0, i+1) from rows 0,1 for i=0,1,2
    // c[i] = det of columns (0, i+1) from rows 2,3 for i=0,1,2
    for i in 0..3 {
        s[i] = ibz_mat_2x2_det_from_ibz(&mat[0][0], &mat[0][i + 1], &mat[1][0], &mat[1][i + 1]);
        c[i] = ibz_mat_2x2_det_from_ibz(&mat[2][0], &mat[2][i + 1], &mat[3][0], &mat[3][i + 1]);
    }

    // s[3+i] = det of columns (1, 2+i) from rows 0,1 for i=0,1
    // c[3+i] = det of columns (1, 2+i) from rows 2,3 for i=0,1
    for i in 0..2 {
        s[3 + i] = ibz_mat_2x2_det_from_ibz(&mat[0][1], &mat[0][2 + i], &mat[1][1], &mat[1][2 + i]);
        c[3 + i] = ibz_mat_2x2_det_from_ibz(&mat[2][1], &mat[2][2 + i], &mat[3][1], &mat[3][2 + i]);
    }

    // s[5] = det of columns (2,3) from rows 0,1
    // c[5] = det of columns (2,3) from rows 2,3
    s[5] = ibz_mat_2x2_det_from_ibz(&mat[0][2], &mat[0][3], &mat[1][2], &mat[1][3]);
    c[5] = ibz_mat_2x2_det_from_ibz(&mat[2][2], &mat[2][3], &mat[3][2], &mat[3][3]);

    // Compute determinant: sum of s[i]*c[5-i] with signs +,-,+,+,-,+
    let mut work_det = Ibz::zero();
    for i in 0..6 {
        let prod = &s[i] * &c[5 - i];
        if i != 1 && i != 4 {
            work_det += &prod;
        } else {
            work_det -= &prod;
        }
    }

    // Compute the transposed adjugate matrix.
    let mut work = IbzMat4x4::default();

    // Helper to convert bool to usize for index arithmetic (C boolean -> 0 or 1).
    let b = |cond: bool| -> usize { usize::from(cond) };

    for j in 0..4usize {
        // k = 0, 1: use c[] array and rows 0,1 of mat (row index = 1-k)
        for k in 0..2usize {
            let row = 1 - k;
            let col_a = b(j == 0);
            let col_b = 2 - b(j > 1);
            let col_c = 3 - b(j == 3);
            let idx_a = 6 - j - b(j == 0);
            let idx_b = 4 - j - b(j == 1);
            let idx_c = 3 - j - b(j == 1) - b(j == 2);

            if (k + j + 1) % 2 == 1 {
                work[j][k] = ibz_inv_dim4_make_coeff_pmp(
                    &mat[row][col_a],
                    &c[idx_a],
                    &mat[row][col_b],
                    &c[idx_b],
                    &mat[row][col_c],
                    &c[idx_c],
                );
            } else {
                work[j][k] = ibz_inv_dim4_make_coeff_mpm(
                    &mat[row][col_a],
                    &c[idx_a],
                    &mat[row][col_b],
                    &c[idx_b],
                    &mat[row][col_c],
                    &c[idx_c],
                );
            }
        }

        // k = 2, 3: use s[] array and rows 2,3 of mat (row index = 3-(k==3))
        for k in 2..4usize {
            let row = 3 - b(k == 3);
            let col_a = b(j == 0);
            let col_b = 2 - b(j > 1);
            let col_c = 3 - b(j == 3);
            let idx_a = 6 - j - b(j == 0);
            let idx_b = 4 - j - b(j == 1);
            let idx_c = 3 - j - b(j == 1) - b(j == 2);

            if (k + j + 1) % 2 == 1 {
                work[j][k] = ibz_inv_dim4_make_coeff_pmp(
                    &mat[row][col_a],
                    &s[idx_a],
                    &mat[row][col_b],
                    &s[idx_b],
                    &mat[row][col_c],
                    &s[idx_c],
                );
            } else {
                work[j][k] = ibz_inv_dim4_make_coeff_mpm(
                    &mat[row][col_a],
                    &s[idx_a],
                    &mat[row][col_b],
                    &s[idx_b],
                    &mat[row][col_c],
                    &s[idx_c],
                );
            }
        }
    }

    let invertible = !work_det.is_zero();

    // If det is zero, return the zero matrix for adjugate.
    let adj = if invertible {
        work
    } else {
        IbzMat4x4::default()
    };

    (adj, work_det, invertible)
}

/// Matrix-vector product `mat * vec`.
pub fn ibz_mat_4x4_eval(mat: &IbzMat4x4, vec: &IbzVec4) -> IbzVec4 {
    let mut res = IbzVec4::default();
    for i in 0..4 {
        let mut acc = Ibz::zero();
        for j in 0..4 {
            acc += &mat[i][j] * &vec[j];
        }
        res[i] = acc;
    }
    res
}

/// Transposed matrix-vector product `vec^T * mat` (equivalently `mat^T * vec`).
pub fn ibz_mat_4x4_eval_t(vec: &IbzVec4, mat: &IbzMat4x4) -> IbzVec4 {
    let mut res = IbzVec4::default();
    for i in 0..4 {
        let mut acc = Ibz::zero();
        for j in 0..4 {
            acc += &mat[j][i] * &vec[j];
        }
        res[i] = acc;
    }
    res
}

/// Quadratic form evaluation: `coord^T * qf * coord`.
pub fn quat_qf_eval(qf: &IbzMat4x4, coord: &IbzVec4) -> Ibz {
    let sum = ibz_mat_4x4_eval(qf, coord);
    let mut result = Ibz::zero();
    for i in 0..4 {
        result += &sum[i] * &coord[i];
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::Zero;

    #[test]
    fn test_vec_4_set() {
        let a = ibz_vec_4_set(1, 2, 3, 4);
        assert_eq!(a[0], BigInt::from(1));
        assert_eq!(a[1], BigInt::from(2));
        assert_eq!(a[2], BigInt::from(3));
        assert_eq!(a[3], BigInt::from(4));
    }

    #[test]
    fn test_vec_4_copy() {
        let a = ibz_vec_4_set(1, 2, 3, 4);
        let b = ibz_vec_4_copy(&a);
        for i in 0..4 {
            assert_eq!(a[i], b[i]);
        }
    }

    #[test]
    fn test_vec_4_negate() {
        let a = ibz_vec_4_set(1, 2, 3, 4);
        let b = ibz_vec_4_negate(&a);
        assert_eq!(b[0], BigInt::from(-1));
        assert_eq!(b[1], BigInt::from(-2));
        assert_eq!(b[2], BigInt::from(-3));
        assert_eq!(b[3], BigInt::from(-4));
    }

    #[test]
    fn test_vec_4_copy_ibz() {
        let (a, b, c, d) = (
            BigInt::from(1),
            BigInt::from(2),
            BigInt::from(3),
            BigInt::from(4),
        );
        let coord = ibz_vec_4_copy_ibz(&a, &b, &c, &d);
        assert_eq!(coord[0], a);
        assert_eq!(coord[1], b);
        assert_eq!(coord[2], c);
        assert_eq!(coord[3], d);
    }

    #[test]
    fn test_vec_4_add() {
        let a = ibz_vec_4_copy_ibz(
            &BigInt::from(1),
            &BigInt::from(-2),
            &BigInt::from(7),
            &BigInt::from(199),
        );
        let b = ibz_vec_4_copy_ibz(
            &BigInt::from(-6),
            &BigInt::from(2),
            &BigInt::from(67),
            &BigInt::from(-22),
        );
        let c = ibz_vec_4_add(&a, &b);
        assert_eq!(c[0], BigInt::from(-5));
        assert_eq!(c[1], BigInt::from(0));
        assert_eq!(c[2], BigInt::from(74));
        assert_eq!(c[3], BigInt::from(177));

        let a = ibz_vec_4_copy_ibz(
            &BigInt::from(-122),
            &BigInt::from(0),
            &BigInt::from(-7),
            &BigInt::from(1889),
        );
        let b = ibz_vec_4_copy_ibz(
            &BigInt::from(-6),
            &BigInt::from(2),
            &BigInt::from(67),
            &BigInt::from(-1889),
        );
        let c = ibz_vec_4_add(&a, &b);
        assert_eq!(c[0], BigInt::from(-128));
        assert_eq!(c[1], BigInt::from(2));
        assert_eq!(c[2], BigInt::from(60));
        assert_eq!(c[3], BigInt::from(0));

        // Self-add: a + a = 2*a
        let a = ibz_vec_4_copy_ibz(
            &BigInt::from(-1),
            &BigInt::from(2),
            &BigInt::from(-7),
            &BigInt::from(19),
        );
        let c = ibz_vec_4_add(&a, &a);
        assert_eq!(c[0], BigInt::from(-2));
        assert_eq!(c[1], BigInt::from(4));
        assert_eq!(c[2], BigInt::from(-14));
        assert_eq!(c[3], BigInt::from(38));
    }

    #[test]
    fn test_vec_4_sub() {
        let a = ibz_vec_4_copy_ibz(
            &BigInt::from(1),
            &BigInt::from(-2),
            &BigInt::from(7),
            &BigInt::from(199),
        );
        let b = ibz_vec_4_copy_ibz(
            &BigInt::from(-6),
            &BigInt::from(2),
            &BigInt::from(67),
            &BigInt::from(-22),
        );
        let c = ibz_vec_4_sub(&a, &b);
        assert_eq!(c[0], BigInt::from(7));
        assert_eq!(c[1], BigInt::from(-4));
        assert_eq!(c[2], BigInt::from(-60));
        assert_eq!(c[3], BigInt::from(221));

        // Self-sub: a - a = 0
        let a = ibz_vec_4_copy_ibz(
            &BigInt::from(-1),
            &BigInt::from(2),
            &BigInt::from(-7),
            &BigInt::from(19),
        );
        let c = ibz_vec_4_sub(&a, &a);
        assert!(ibz_vec_4_is_zero(&c));
    }

    #[test]
    fn test_vec_4_is_zero() {
        assert!(ibz_vec_4_is_zero(&ibz_vec_4_set(0, 0, 0, 0)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(20, 0, 0, 0)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(0, -1, 0, 0)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(0, 0, 2, 0)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(0, 0, 0, 1)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(1, 1, 1, 1)));
        assert!(!ibz_vec_4_is_zero(&ibz_vec_4_set(-1, 1, 1, -1)));
    }

    #[test]
    fn test_vec_4_linear_combination() {
        let a = ibz_vec_4_set(1, 2, 3, 4);
        let b = ibz_vec_4_set(-2, 1, 3, -3);
        let ca = BigInt::from(2);
        let cb = BigInt::from(-1);
        // 2*[1,2,3,4] + (-1)*[-2,1,3,-3] = [2+2, 4-1, 6-3, 8+3] = [4, 3, 3, 11]
        let lc = ibz_vec_4_linear_combination(&ca, &a, &cb, &b);
        assert_eq!(lc[0], BigInt::from(4));
        assert_eq!(lc[1], BigInt::from(3));
        assert_eq!(lc[2], BigInt::from(3));
        assert_eq!(lc[3], BigInt::from(11));

        // 2*a + (-1)*a = a
        let lc = ibz_vec_4_linear_combination(&ca, &a, &cb, &a);
        assert_eq!(lc[0], BigInt::from(1));
        assert_eq!(lc[1], BigInt::from(2));
        assert_eq!(lc[2], BigInt::from(3));
        assert_eq!(lc[3], BigInt::from(4));
    }

    #[test]
    fn test_vec_4_scalar_mul() {
        let s = BigInt::from(5);
        let vec = ibz_vec_4_set(0, 1, 2, 3);
        let prod = ibz_vec_4_scalar_mul(&s, &vec);
        for i in 0..4 {
            assert_eq!(prod[i], BigInt::from(i as i32 * 5));
        }
    }

    #[test]
    fn test_vec_4_scalar_div() {
        let s = BigInt::from(5);
        let vec = ibz_vec_4_set(0, 5, 10, 15);
        let (quot, exact) = ibz_vec_4_scalar_div(&s, &vec);
        assert!(exact);
        for i in 0..4 {
            assert_eq!(quot[i], BigInt::from(i as i32));
        }

        // Non-exact division
        let vec = ibz_vec_4_set(1, 5, 10, 15);
        let (quot, exact) = ibz_vec_4_scalar_div(&s, &vec);
        assert!(!exact);
        assert_eq!(quot[0], BigInt::from(0)); // 1/5 truncated = 0
    }

    #[test]
    fn test_mat_4x4_mul() {
        // Set up matrices
        let mut a = IbzMat4x4::default();
        a[0][0] = BigInt::from(1);
        a[0][1] = BigInt::from(2);
        a[0][2] = BigInt::from(1);
        a[1][1] = BigInt::from(1);
        a[1][2] = BigInt::from(3);
        a[2][2] = BigInt::from(1);
        a[2][3] = BigInt::from(4);
        a[3][3] = BigInt::from(1);

        let mut b = IbzMat4x4::default();
        b[0][0] = BigInt::from(-1);
        b[1][0] = BigInt::from(1);
        b[1][1] = BigInt::from(-2);
        b[2][0] = BigInt::from(1);
        b[2][2] = BigInt::from(3);
        b[3][0] = BigInt::from(1);
        b[3][1] = BigInt::from(5);
        b[3][2] = BigInt::from(-1);
        b[3][3] = BigInt::from(2);

        let prod = ibz_mat_4x4_mul(&a, &b);

        // Verify the product manually:
        // Row 0: [1*-1+2*1+1*1, 1*0+2*-2+1*0, 1*0+2*0+1*3, 0] = [2, -4, 3, 0]
        assert_eq!(prod[0][0], BigInt::from(2));
        assert_eq!(prod[0][1], BigInt::from(-4));
        assert_eq!(prod[0][2], BigInt::from(3));
        assert_eq!(prod[0][3], BigInt::from(0));
        // Row 1: [0+1+3, 0-2+0, 0+0+9, 0] = [4, -2, 9, 0]
        assert_eq!(prod[1][0], BigInt::from(4));
        assert_eq!(prod[1][1], BigInt::from(-2));
        assert_eq!(prod[1][2], BigInt::from(9));
        assert_eq!(prod[1][3], BigInt::from(0));
        // Row 2: [0+0+1+4, 0+0+0+20, 0+0+3-4, 0+0+0+8] = [5, 20, -1, 8]
        assert_eq!(prod[2][0], BigInt::from(5));
        assert_eq!(prod[2][1], BigInt::from(20));
        assert_eq!(prod[2][2], BigInt::from(-1));
        assert_eq!(prod[2][3], BigInt::from(8));
        // Row 3: [1, 5, -1, 2]
        assert_eq!(prod[3][0], BigInt::from(1));
        assert_eq!(prod[3][1], BigInt::from(5));
        assert_eq!(prod[3][2], BigInt::from(-1));
        assert_eq!(prod[3][3], BigInt::from(2));
    }

    #[test]
    fn test_mat_4x4_copy() {
        let mut mat = IbzMat4x4::default();
        mat[0][0] = BigInt::from(1);
        mat[0][1] = BigInt::from(2);
        mat[0][2] = BigInt::from(-7);
        mat[0][3] = BigInt::from(77);
        mat[2][0] = BigInt::from(13);
        mat[1][1] = BigInt::from(20);
        mat[3][2] = BigInt::from(-77);
        mat[3][3] = BigInt::from(7);
        let copy = ibz_mat_4x4_copy(&mat);
        assert!(ibz_mat_4x4_equal(&copy, &mat));
    }

    #[test]
    fn test_mat_4x4_negate() {
        let mut mat = IbzMat4x4::default();
        mat[0][0] = BigInt::from(1);
        mat[0][1] = BigInt::from(2);
        mat[0][2] = BigInt::from(-7);
        mat[0][3] = BigInt::from(77);
        mat[2][0] = BigInt::from(13);
        mat[1][1] = BigInt::from(20);
        mat[3][2] = BigInt::from(-77);
        mat[3][3] = BigInt::from(7);

        let mut cmp = IbzMat4x4::default();
        cmp[0][0] = BigInt::from(-1);
        cmp[0][1] = BigInt::from(-2);
        cmp[0][2] = BigInt::from(7);
        cmp[0][3] = BigInt::from(-77);
        cmp[2][0] = BigInt::from(-13);
        cmp[1][1] = BigInt::from(-20);
        cmp[3][2] = BigInt::from(77);
        cmp[3][3] = BigInt::from(-7);

        let neg = ibz_mat_4x4_negate(&mat);
        assert!(ibz_mat_4x4_equal(&neg, &cmp));
    }

    #[test]
    fn test_mat_4x4_transpose() {
        let mut mat = IbzMat4x4::default();
        mat[0][0] = BigInt::from(1);
        mat[0][1] = BigInt::from(2);
        mat[0][2] = BigInt::from(-7);
        mat[0][3] = BigInt::from(77);
        mat[2][0] = BigInt::from(13);
        mat[1][1] = BigInt::from(20);
        mat[3][2] = BigInt::from(-77);
        mat[3][3] = BigInt::from(7);

        let mut cmp = IbzMat4x4::default();
        cmp[0][0] = BigInt::from(1);
        cmp[1][0] = BigInt::from(2);
        cmp[2][0] = BigInt::from(-7);
        cmp[3][0] = BigInt::from(77);
        cmp[0][2] = BigInt::from(13);
        cmp[1][1] = BigInt::from(20);
        cmp[2][3] = BigInt::from(-77);
        cmp[3][3] = BigInt::from(7);

        let t = ibz_mat_4x4_transpose(&mat);
        assert!(ibz_mat_4x4_equal(&t, &cmp));
    }

    #[test]
    fn test_mat_4x4_zero() {
        let z = ibz_mat_4x4_zero();
        for i in 0..4 {
            for j in 0..4 {
                assert!(z[i][j].is_zero());
            }
        }
    }

    #[test]
    fn test_mat_4x4_identity() {
        let id = ibz_mat_4x4_identity();
        assert!(ibz_mat_4x4_is_identity(&id));

        // Modify and verify it's no longer identity
        let mut m = id.clone();
        m[0][1] = BigInt::from(1);
        assert!(!ibz_mat_4x4_is_identity(&m));

        let mut m = ibz_mat_4x4_identity();
        m[3][3] = BigInt::from(0);
        assert!(!ibz_mat_4x4_is_identity(&m));
    }

    #[test]
    fn test_mat_4x4_equal() {
        let mut a = IbzMat4x4::default();
        let mut b = IbzMat4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                a[i][j] = BigInt::from((i + j) as i32);
                b[i][j] = BigInt::from((i + j) as i32);
            }
        }
        assert!(ibz_mat_4x4_equal(&a, &b));
        b[2][2] = BigInt::from(2); // was 4, now 2
        assert!(!ibz_mat_4x4_equal(&a, &b));
    }

    #[test]
    fn test_mat_4x4_scalar_mul() {
        let s = BigInt::from(5);
        let mut mat = IbzMat4x4::default();
        let mut cmp = IbzMat4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from((i + j) as i32);
                cmp[i][j] = BigInt::from((i + j) as i32 * 5);
            }
        }
        let prod = ibz_mat_4x4_scalar_mul(&s, &mat);
        assert!(ibz_mat_4x4_equal(&prod, &cmp));
    }

    #[test]
    fn test_mat_4x4_gcd() {
        // GCD of matrix with d=2
        let mut mat = IbzMat4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from(2 * i as i32 * j as i32);
            }
        }
        let g = ibz_mat_4x4_gcd(&mat);
        assert_eq!(g, BigInt::from(2));

        // GCD with d=21
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from(21 * i as i32 * j as i32);
            }
        }
        let g = ibz_mat_4x4_gcd(&mat);
        assert_eq!(g, BigInt::from(21));
    }

    #[test]
    fn test_mat_4x4_scalar_div() {
        let s = BigInt::from(5);
        let mut mat = IbzMat4x4::default();
        let mut cmp = IbzMat4x4::default();
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from((i + j) as i32 * 5);
                cmp[i][j] = BigInt::from((i + j) as i32);
            }
        }
        let (quot, exact) = ibz_mat_4x4_scalar_div(&s, &mat);
        assert!(exact);
        assert!(ibz_mat_4x4_equal(&quot, &cmp));
    }

    #[test]
    fn test_inv_dim4_make_coeff_pmp() {
        // 0*3 - (-1)*0 + (-1)*0 = 0
        let c = ibz_inv_dim4_make_coeff_pmp(
            &BigInt::from(0),
            &BigInt::from(3),
            &BigInt::from(-1),
            &BigInt::from(0),
            &BigInt::from(-1),
            &BigInt::from(0),
        );
        assert_eq!(c, BigInt::from(0));

        // 2*3 - (-1)*1 + (-4)*2 = 6+1-8 = -1
        let c = ibz_inv_dim4_make_coeff_pmp(
            &BigInt::from(2),
            &BigInt::from(3),
            &BigInt::from(-1),
            &BigInt::from(1),
            &BigInt::from(-4),
            &BigInt::from(2),
        );
        assert_eq!(c, BigInt::from(-1));

        // a*a - a*a + a*a = a^2 for a=2 → 4
        let a = BigInt::from(2);
        let c = ibz_inv_dim4_make_coeff_pmp(&a, &a, &a, &a, &a, &a);
        assert_eq!(c, BigInt::from(4));
    }

    #[test]
    fn test_inv_dim4_make_coeff_mpm() {
        // -0*3 + (-1)*0 - (-1)*0 = 0
        let c = ibz_inv_dim4_make_coeff_mpm(
            &BigInt::from(0),
            &BigInt::from(3),
            &BigInt::from(-1),
            &BigInt::from(0),
            &BigInt::from(-1),
            &BigInt::from(0),
        );
        assert_eq!(c, BigInt::from(0));

        // -2*3 + (-1)*1 - (-4)*2 = -6-1+8 = 1
        let c = ibz_inv_dim4_make_coeff_mpm(
            &BigInt::from(2),
            &BigInt::from(3),
            &BigInt::from(-1),
            &BigInt::from(1),
            &BigInt::from(-4),
            &BigInt::from(2),
        );
        assert_eq!(c, BigInt::from(1));

        // -a*a + a*a - a*a = -a^2 for a=2 → -4
        let a = BigInt::from(2);
        let c = ibz_inv_dim4_make_coeff_mpm(&a, &a, &a, &a, &a, &a);
        assert_eq!(c, BigInt::from(-4));
    }

    fn validate_inv_with_det(mat: &IbzMat4x4, inv: &IbzMat4x4, det: &Ibz) -> bool {
        let det_id = ibz_mat_4x4_scalar_mul(det, &ibz_mat_4x4_identity());
        let prod1 = ibz_mat_4x4_mul(inv, mat);
        let prod2 = ibz_mat_4x4_mul(mat, inv);
        ibz_mat_4x4_equal(&prod1, &det_id) && ibz_mat_4x4_equal(&prod2, &det_id)
    }

    #[test]
    fn test_mat_4x4_inv_with_det_as_denom() {
        // Zero matrix: not invertible
        let mat = ibz_mat_4x4_zero();
        let (_, det, ok) = ibz_mat_4x4_inv_with_det_as_denom(&mat);
        assert!(!ok);
        assert!(det.is_zero());

        // Identity: det=1, inv=I
        let mat = ibz_mat_4x4_identity();
        let (inv, det, ok) = ibz_mat_4x4_inv_with_det_as_denom(&mat);
        assert!(ok);
        assert!(det.is_one());
        assert!(validate_inv_with_det(&mat, &inv, &det));

        // Upper triangular matrix
        let mut mat = IbzMat4x4::default();
        mat[0][0] = BigInt::from(2);
        mat[0][1] = BigInt::from(-17);
        mat[0][2] = BigInt::from(3);
        mat[0][3] = BigInt::from(5);
        mat[1][1] = BigInt::from(-2);
        mat[1][2] = BigInt::from(3);
        mat[1][3] = BigInt::from(2);
        mat[2][2] = BigInt::from(-3);
        mat[2][3] = BigInt::from(0);
        mat[3][3] = BigInt::from(1);
        let (inv, det, ok) = ibz_mat_4x4_inv_with_det_as_denom(&mat);
        assert!(ok);
        assert_eq!(det, BigInt::from(12));
        assert!(validate_inv_with_det(&mat, &inv, &det));

        // Full matrix with lower triangular additions
        mat[3][0] = BigInt::from(1);
        mat[3][1] = BigInt::from(8);
        mat[3][2] = BigInt::from(-9);
        mat[2][0] = BigInt::from(3);
        mat[2][1] = BigInt::from(0);
        mat[1][0] = BigInt::from(4);
        let (inv, det, ok) = ibz_mat_4x4_inv_with_det_as_denom(&mat);
        assert!(ok);
        assert_eq!(det, BigInt::from(-1503));
        assert!(validate_inv_with_det(&mat, &inv, &det));
    }

    #[test]
    fn test_mat_4x4_eval() {
        // mat[i][j] = i*j, vec[i] = i
        let mut mat = IbzMat4x4::default();
        let mut vec = IbzVec4::default();
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from(i as i32 * j as i32);
            }
            vec[i] = BigInt::from(i as i32);
        }
        // Expected: row i = sum_j (i*j)*j = i*sum_j j^2 = i*14
        let res = ibz_mat_4x4_eval(&mat, &vec);
        assert_eq!(res[0], BigInt::from(0));
        assert_eq!(res[1], BigInt::from(14));
        assert_eq!(res[2], BigInt::from(28));
        assert_eq!(res[3], BigInt::from(42));

        // eval_t: mat^T * vec should give same result when called with transposed mat
        let mat_t = ibz_mat_4x4_transpose(&mat);
        let res_t = ibz_mat_4x4_eval_t(&vec, &mat_t);
        assert_eq!(res_t[0], BigInt::from(0));
        assert_eq!(res_t[1], BigInt::from(14));
        assert_eq!(res_t[2], BigInt::from(28));
        assert_eq!(res_t[3], BigInt::from(42));

        // Second test: mat[i][j] = i*(j-1)+1, vec[i] = i*i-2
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = BigInt::from(i as i32 * (j as i32 - 1) + 1);
            }
            vec[i] = BigInt::from(i as i32 * i as i32 - 2);
        }
        let res = ibz_mat_4x4_eval(&mat, &vec);
        assert_eq!(res[0], BigInt::from(6));
        assert_eq!(res[1], BigInt::from(24));
        assert_eq!(res[2], BigInt::from(42));
        assert_eq!(res[3], BigInt::from(60));
    }

    #[test]
    fn test_qf_eval() {
        // qf[i][j] = i*j, coord[i] = i
        // Expected: coord^T * qf * coord
        // qf*coord: row i = sum_j (i*j)*j = i*14
        // Then coord^T * result = sum_i i*(i*14) = 14*(0+1+4+9) = 14*14 = 196
        let mut qf = IbzMat4x4::default();
        let mut coord = IbzVec4::default();
        for i in 0..4 {
            for j in 0..4 {
                qf[i][j] = BigInt::from(i as i32 * j as i32);
            }
            coord[i] = BigInt::from(i as i32);
        }
        let res = quat_qf_eval(&qf, &coord);
        assert_eq!(res, BigInt::from(196));

        // qf[i][j] = (i+1)*(j+1)-4, coord[i] = (i-1)*2-2
        for i in 0..4 {
            for j in 0..4 {
                qf[i][j] = BigInt::from((i as i32 + 1) * (j as i32 + 1) - 4);
            }
            coord[i] = BigInt::from((i as i32 - 1) * 2 - 2);
        }
        let res = quat_qf_eval(&qf, &coord);
        assert_eq!(res, BigInt::from(-64));
    }

    #[test]
    fn test_vec_4_content() {
        // All zeros
        let x = ibz_vec_4_set(0, 0, 0, 0);
        assert_eq!(ibz_vec_4_content(&x), BigInt::from(0));

        // [5, 25, 125, 30] → GCD = 5
        let x = ibz_vec_4_set(5, 25, 125, 30);
        assert_eq!(ibz_vec_4_content(&x), BigInt::from(5));

        // [5, 2, 125, 30] → GCD = 1
        let x = ibz_vec_4_set(5, 2, 125, 30);
        assert_eq!(ibz_vec_4_content(&x), BigInt::from(1));

        // [5, -2, 125, 0] → GCD = 1
        let x = ibz_vec_4_set(5, -2, 125, 0);
        assert_eq!(ibz_vec_4_content(&x), BigInt::from(1));

        // [0, -2, 0, 0] → GCD = 2
        let x = ibz_vec_4_set(0, -2, 0, 0);
        assert_eq!(ibz_vec_4_content(&x), BigInt::from(2));
    }
}
