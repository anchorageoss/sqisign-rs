use crate::fp::{Fp2, FpBackend};
use subtle::Choice;

use super::{BasisChangeMatrix, PrecompBasisChangeMatrix, ThetaPoint};

/// Select a basis change matrix in constant time.
///
/// If `option` is clear, result = `m1` (regular matrix).
/// If `option` is set, result = expanded `m2` (precomputed matrix).
#[inline]
pub fn select_base_change_matrix<L: FpBackend>(
    m1: &BasisChangeMatrix<L>,
    m2: &PrecompBasisChangeMatrix,
    fp2_constants: &[Fp2<L>; 5],
    option: Choice,
) -> BasisChangeMatrix<L> {
    let mut out = BasisChangeMatrix::default();
    for i in 0..4 {
        for j in 0..4 {
            out.m[i][j] = Fp2::select(&m1.m[i][j], &fp2_constants[m2.m[i][j] as usize], option);
        }
    }
    out
}

/// Expand a precomputed basis change matrix (u8 indices) into full 𝔽p² elements.
#[inline]
pub fn set_base_change_matrix_from_precomp<L: FpBackend>(
    m: &PrecompBasisChangeMatrix,
    fp2_constants: &[Fp2<L>; 5],
) -> BasisChangeMatrix<L> {
    let mut out = BasisChangeMatrix::default();
    for i in 0..4 {
        for j in 0..4 {
            out.m[i][j] = fp2_constants[m.m[i][j] as usize].clone();
        }
    }
    out
}

/// Extract a coordinate from a theta point by index (mod 4).
///
/// 0 = x, 1 = y, 2 = z, 3 = t
#[inline]
pub fn choose_index_theta_point<L: FpBackend>(ind: usize, p: &ThetaPoint<L>) -> Fp2<L> {
    match ind % 4 {
        0 => p.x.clone(),
        1 => p.y.clone(),
        2 => p.z.clone(),
        3 => p.t.clone(),
        // SAFETY: x % 4 can only produce 0, 1, 2, or 3
        _ => unreachable!("x % 4 can only produce 0, 1, 2, or 3"),
    }
}

/// Apply a 4×4 basis change matrix to a theta point.
///
/// When `pt_not_zero` is true, all 4 columns of M are used (general case).
/// When false, the t-coordinate column (column 3) is skipped as an
/// optimization when P.t is known to be zero.
#[inline]
pub fn apply_isomorphism_general<L: FpBackend>(
    m: &BasisChangeMatrix<L>,
    p: &ThetaPoint<L>,
    pt_not_zero: bool,
) -> ThetaPoint<L> {
    let out_x = p.x.mul(&m.m[0][0]);
    let out_x = out_x.add(&p.y.mul(&m.m[0][1]));
    let out_x = out_x.add(&p.z.mul(&m.m[0][2]));

    let out_y = p.x.mul(&m.m[1][0]);
    let out_y = out_y.add(&p.y.mul(&m.m[1][1]));
    let out_y = out_y.add(&p.z.mul(&m.m[1][2]));

    let out_z = p.x.mul(&m.m[2][0]);
    let out_z = out_z.add(&p.y.mul(&m.m[2][1]));
    let out_z = out_z.add(&p.z.mul(&m.m[2][2]));

    let out_t = p.x.mul(&m.m[3][0]);
    let out_t = out_t.add(&p.y.mul(&m.m[3][1]));
    let out_t = out_t.add(&p.z.mul(&m.m[3][2]));

    let out_x = if pt_not_zero {
        out_x.add(&p.t.mul(&m.m[0][3]))
    } else {
        out_x
    };
    let out_y = if pt_not_zero {
        out_y.add(&p.t.mul(&m.m[1][3]))
    } else {
        out_y
    };
    let out_z = if pt_not_zero {
        out_z.add(&p.t.mul(&m.m[2][3]))
    } else {
        out_z
    };
    let out_t = if pt_not_zero {
        out_t.add(&p.t.mul(&m.m[3][3]))
    } else {
        out_t
    };

    ThetaPoint {
        x: out_x,
        y: out_y,
        z: out_z,
        t: out_t,
    }
}

/// Apply a 4×4 basis change matrix to a theta point (full, all 4 columns).
#[inline]
pub fn apply_isomorphism<L: FpBackend>(
    m: &BasisChangeMatrix<L>,
    p: &ThetaPoint<L>,
) -> ThetaPoint<L> {
    apply_isomorphism_general(m, p, true)
}

/// Multiply two 4×4 basis change matrices: res = m1 * m2.
#[inline]
pub fn base_change_matrix_multiplication<L: FpBackend>(
    m1: &BasisChangeMatrix<L>,
    m2: &BasisChangeMatrix<L>,
) -> BasisChangeMatrix<L> {
    let mut res = BasisChangeMatrix::default();
    for i in 0..4 {
        for j in 0..4 {
            let mut sum = Fp2::<L>::zero();
            for k in 0..4 {
                let prod = m1.m[i][k].mul(&m2.m[k][j]);
                sum = sum.add(&prod);
            }
            res.m[i][j] = sum;
        }
    }
    res
}
