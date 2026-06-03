//!
//! Operations on elements of the quaternion algebra `B_{p,∞}` with basis
//! `{1, i, j, ij}` where `i² = -1` and `j² = -p`. Elements are represented
//! as `coord / denom` where `coord` is a 4-vector of integer numerators.

use super::dim4::{
    ibz_vec_4_add, ibz_vec_4_content, ibz_vec_4_copy_ibz, ibz_vec_4_is_zero, ibz_vec_4_scalar_div,
    ibz_vec_4_scalar_mul, ibz_vec_4_set, ibz_vec_4_sub,
};
use super::intbig::{ibz_div, ibz_gcd, Ibz};
use super::types::{IbzVec4, QuatAlg, QuatAlgElem};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

/// Initialize a quaternion algebra from a small unsigned prime.
pub fn quat_alg_init_set_ui(p: u32) -> QuatAlg {
    QuatAlg::new(&BigInt::from(p))
}

/// Coordinate multiplication in the quaternion algebra.
///
/// Computes the coordinate vector of the product `a * b` in the basis
/// `{1, i, j, ij}` of `B_{p,∞}`, without considering denominators.
pub fn quat_alg_coord_mul(a: &IbzVec4, b: &IbzVec4, alg: &QuatAlg) -> IbzVec4 {
    let mut sum = IbzVec4::default();

    // 1-coordinate: a0*b0 - a1*b1 - p*(a2*b2 + a3*b3)
    sum[0] = -(&a[2] * &b[2]) - &a[3] * &b[3];
    sum[0] = &sum[0] * &alg.p;
    sum[0] = &sum[0] + &a[0] * &b[0];
    sum[0] = &sum[0] - &a[1] * &b[1];

    // i-coordinate: a0*b1 + a1*b0 + p*(a2*b3 - a3*b2)
    sum[1] = &a[2] * &b[3] - &a[3] * &b[2];
    sum[1] = &sum[1] * &alg.p;
    sum[1] = &sum[1] + &a[0] * &b[1];
    sum[1] = &sum[1] + &a[1] * &b[0];

    // j-coordinate: a0*b2 + a2*b0 - a1*b3 + a3*b1
    sum[2] = &a[0] * &b[2] + &a[2] * &b[0];
    sum[2] = &sum[2] - &a[1] * &b[3];
    sum[2] = &sum[2] + &a[3] * &b[1];

    // ij-coordinate: a0*b3 + a3*b0 - a2*b1 + a1*b2
    sum[3] = &a[0] * &b[3] + &a[3] * &b[0];
    sum[3] = &sum[3] - &a[2] * &b[1];
    sum[3] = &sum[3] + &a[1] * &b[2];

    sum
}

/// Put two quaternion algebra elements on a common denominator.
///
/// Returns `(res_a, res_b)` where both have the same denominator `lcm(a.denom, b.denom)`.
pub fn quat_alg_equal_denom(a: &QuatAlgElem, b: &QuatAlgElem) -> (QuatAlgElem, QuatAlgElem) {
    let gcd = ibz_gcd(&a.denom, &b.denom);
    let (red_a_denom, _) = ibz_div(&a.denom, &gcd);
    let (red_b_denom, _) = ibz_div(&b.denom, &gcd);

    let mut res_a = QuatAlgElem::default();
    let mut res_b = QuatAlgElem::default();

    for i in 0..4 {
        res_a.coord[i] = &a.coord[i] * &red_b_denom;
        res_b.coord[i] = &b.coord[i] * &red_a_denom;
    }

    let common_reduced = &red_a_denom * &red_b_denom;
    res_a.denom = &common_reduced * &gcd;
    res_b.denom = &common_reduced * &gcd;

    (res_a, res_b)
}

/// Addition of two quaternion algebra elements.
pub fn quat_alg_add(a: &QuatAlgElem, b: &QuatAlgElem) -> QuatAlgElem {
    let (ra, rb) = quat_alg_equal_denom(a, b);
    QuatAlgElem {
        denom: ra.denom,
        coord: ibz_vec_4_add(&ra.coord, &rb.coord),
    }
}

/// Subtraction of two quaternion algebra elements.
pub fn quat_alg_sub(a: &QuatAlgElem, b: &QuatAlgElem) -> QuatAlgElem {
    let (ra, rb) = quat_alg_equal_denom(a, b);
    QuatAlgElem {
        denom: ra.denom,
        coord: ibz_vec_4_sub(&ra.coord, &rb.coord),
    }
}

/// Multiplication of two quaternion algebra elements.
pub fn quat_alg_mul(a: &QuatAlgElem, b: &QuatAlgElem, alg: &QuatAlg) -> QuatAlgElem {
    QuatAlgElem {
        denom: &a.denom * &b.denom,
        coord: quat_alg_coord_mul(&a.coord, &b.coord, alg),
    }
}

/// Reduced norm of a quaternion algebra element.
///
/// Returns `(numerator, denominator)` of `nrd(a)`. The result is always
/// non-negative.
pub fn quat_alg_norm(a: &QuatAlgElem, alg: &QuatAlg) -> (Ibz, Ibz) {
    let conj = quat_alg_conj(a);
    let norm_elem = quat_alg_mul(a, &conj, alg);
    let g = ibz_gcd(&norm_elem.coord[0], &norm_elem.denom);
    let (num, _) = ibz_div(&norm_elem.coord[0], &g);
    let (den, _) = ibz_div(&norm_elem.denom, &g);
    (num.abs(), den.abs())
}

/// Create a scalar quaternion element `numerator / denominator`.
pub fn quat_alg_scalar(numerator: &Ibz, denominator: &Ibz) -> QuatAlgElem {
    QuatAlgElem {
        denom: denominator.clone(),
        coord: ibz_vec_4_copy_ibz(numerator, &Ibz::zero(), &Ibz::zero(), &Ibz::zero()),
    }
}

/// Standard involution (conjugation): negates the `i`, `j`, `ij` components.
pub fn quat_alg_conj(x: &QuatAlgElem) -> QuatAlgElem {
    QuatAlgElem {
        denom: x.denom.clone(),
        coord: IbzVec4([x.coord[0].clone(), -&x.coord[1], -&x.coord[2], -&x.coord[3]]),
    }
}

/// Factor an element of an order into its primitive and content parts.
///
/// Given `x ∈ order`, returns `(primitive_coord, content)` where `content`
/// is the GCD of `x`'s coordinates in the order basis, and `primitive_coord`
/// is the coordinate vector divided by `content`.
pub fn quat_alg_make_primitive(
    x: &QuatAlgElem,
    order: &super::types::QuatLattice,
) -> (IbzVec4, Ibz) {
    let coord = super::lattice::quat_lattice_contains(order, x)
        .expect("invariant: element must be contained in the order");
    let content = ibz_vec_4_content(&coord);
    let (primitive, _) = ibz_vec_4_scalar_div(&content, &coord);
    (primitive, content)
}

/// Normalize a quaternion element so that `gcd(denom, content(coord)) = 1`
/// and `denom > 0`.
pub fn quat_alg_normalize(x: &mut QuatAlgElem) {
    let mut gcd = ibz_vec_4_content(&x.coord);
    gcd = ibz_gcd(&gcd, &x.denom);
    let (new_denom, _) = ibz_div(&x.denom, &gcd);
    x.denom = new_denom;
    let (new_coord, _) = ibz_vec_4_scalar_div(&gcd, &x.coord);
    x.coord = new_coord;
    // Ensure denom > 0
    if x.denom < Ibz::zero() {
        x.denom = -&x.denom;
        x.coord = IbzVec4([-&x.coord[0], -&x.coord[1], -&x.coord[2], -&x.coord[3]]);
    }
}

/// Compares as rational values (normalizes denominators).
pub fn quat_alg_elem_equal(a: &QuatAlgElem, b: &QuatAlgElem) -> bool {
    let diff = quat_alg_sub(a, b);
    quat_alg_elem_is_zero(&diff)
}

pub fn quat_alg_elem_is_zero(x: &QuatAlgElem) -> bool {
    ibz_vec_4_is_zero(&x.coord)
}

pub fn quat_alg_elem_set(
    denom: i32,
    coord0: i32,
    coord1: i32,
    coord2: i32,
    coord3: i32,
) -> QuatAlgElem {
    QuatAlgElem {
        denom: BigInt::from(denom),
        coord: ibz_vec_4_set(coord0, coord1, coord2, coord3),
    }
}

pub fn quat_alg_elem_copy(src: &QuatAlgElem) -> QuatAlgElem {
    src.clone()
}

/// Construct a quaternion element from big integer references.
pub fn quat_alg_elem_copy_ibz(
    denom: &Ibz,
    coord0: &Ibz,
    coord1: &Ibz,
    coord2: &Ibz,
    coord3: &Ibz,
) -> QuatAlgElem {
    QuatAlgElem {
        denom: denom.clone(),
        coord: ibz_vec_4_copy_ibz(coord0, coord1, coord2, coord3),
    }
}

/// Multiply a quaternion element by a scalar (multiply all coordinates).
pub fn quat_alg_elem_mul_by_scalar(scalar: &Ibz, elem: &QuatAlgElem) -> QuatAlgElem {
    QuatAlgElem {
        denom: elem.denom.clone(),
        coord: ibz_vec_4_scalar_mul(scalar, &elem.coord),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_alg() -> QuatAlg {
        QuatAlg::new(&BigInt::from(7))
    }

    #[test]
    fn test_coord_mul() {
        let alg = make_alg();
        // 1 * 1 = 1
        let one = ibz_vec_4_set(1, 0, 0, 0);
        let prod = quat_alg_coord_mul(&one, &one, &alg);
        assert_eq!(prod[0], BigInt::from(1));
        assert!(prod[1].is_zero() && prod[2].is_zero() && prod[3].is_zero());

        // i * i = -1
        let i_vec = ibz_vec_4_set(0, 1, 0, 0);
        let prod = quat_alg_coord_mul(&i_vec, &i_vec, &alg);
        assert_eq!(prod[0], BigInt::from(-1));
        assert!(prod[1].is_zero() && prod[2].is_zero() && prod[3].is_zero());

        // j * j = -p = -7
        let j_vec = ibz_vec_4_set(0, 0, 1, 0);
        let prod = quat_alg_coord_mul(&j_vec, &j_vec, &alg);
        assert_eq!(prod[0], BigInt::from(-7));
        assert!(prod[1].is_zero() && prod[2].is_zero() && prod[3].is_zero());

        // i * j = ij
        let prod = quat_alg_coord_mul(&i_vec, &j_vec, &alg);
        assert!(prod[0].is_zero() && prod[1].is_zero() && prod[2].is_zero());
        assert_eq!(prod[3], BigInt::from(1));

        // j * i = -ij
        let prod = quat_alg_coord_mul(&j_vec, &i_vec, &alg);
        assert!(prod[0].is_zero() && prod[1].is_zero() && prod[2].is_zero());
        assert_eq!(prod[3], BigInt::from(-1));
    }

    #[test]
    fn test_norm() {
        let alg = make_alg();
        // nrd(1) = 1
        let one = quat_alg_elem_set(1, 1, 0, 0, 0);
        let (n, d) = quat_alg_norm(&one, &alg);
        assert_eq!(n, BigInt::from(1));
        assert_eq!(d, BigInt::from(1));

        // nrd(i) = 1
        let i_elem = quat_alg_elem_set(1, 0, 1, 0, 0);
        let (n, d) = quat_alg_norm(&i_elem, &alg);
        assert_eq!(n, BigInt::from(1));
        assert_eq!(d, BigInt::from(1));

        // nrd(j) = p = 7
        let j_elem = quat_alg_elem_set(1, 0, 0, 1, 0);
        let (n, d) = quat_alg_norm(&j_elem, &alg);
        assert_eq!(n, BigInt::from(7));
        assert_eq!(d, BigInt::from(1));

        // nrd(1+i) = 2
        let elem = quat_alg_elem_set(1, 1, 1, 0, 0);
        let (n, d) = quat_alg_norm(&elem, &alg);
        assert_eq!(n, BigInt::from(2));
        assert_eq!(d, BigInt::from(1));
    }

    #[test]
    fn test_conj() {
        let x = quat_alg_elem_set(2, 1, 3, 5, 7);
        let c = quat_alg_conj(&x);
        assert_eq!(c.denom, BigInt::from(2));
        assert_eq!(c.coord[0], BigInt::from(1));
        assert_eq!(c.coord[1], BigInt::from(-3));
        assert_eq!(c.coord[2], BigInt::from(-5));
        assert_eq!(c.coord[3], BigInt::from(-7));
    }

    #[test]
    fn test_add_sub() {
        let a = quat_alg_elem_set(1, 1, 2, 3, 4);
        let b = quat_alg_elem_set(1, 5, 6, 7, 8);
        let sum = quat_alg_add(&a, &b);
        let diff = quat_alg_sub(&a, &b);

        // Normalize to check values
        let mut sum_n = sum.clone();
        quat_alg_normalize(&mut sum_n);
        assert_eq!(sum_n.coord[0], BigInt::from(6));
        assert_eq!(sum_n.coord[1], BigInt::from(8));

        let mut diff_n = diff.clone();
        quat_alg_normalize(&mut diff_n);
        assert_eq!(diff_n.coord[0], BigInt::from(-4));
        assert_eq!(diff_n.coord[1], BigInt::from(-4));
    }

    #[test]
    fn test_mul() {
        let alg = make_alg();
        // (1+i)(1-i) = 1 - i + i - i^2 = 1 + 1 = 2
        let a = quat_alg_elem_set(1, 1, 1, 0, 0);
        let b = quat_alg_elem_set(1, 1, -1, 0, 0);
        let prod = quat_alg_mul(&a, &b, &alg);
        let mut prod_n = prod;
        quat_alg_normalize(&mut prod_n);
        assert_eq!(prod_n.coord[0], BigInt::from(2));
        assert!(prod_n.coord[1].is_zero());
        assert!(prod_n.coord[2].is_zero());
        assert!(prod_n.coord[3].is_zero());
        assert_eq!(prod_n.denom, BigInt::from(1));
    }

    #[test]
    fn test_normalize() {
        // (6, [12, 18, 24, 30]) → normalize → gcd = gcd(6, gcd(12,18,24,30)) = gcd(6,6) = 6
        // → (1, [2, 3, 4, 5])
        let mut x = QuatAlgElem {
            denom: BigInt::from(6),
            coord: IbzVec4([
                BigInt::from(12),
                BigInt::from(18),
                BigInt::from(24),
                BigInt::from(30),
            ]),
        };
        quat_alg_normalize(&mut x);
        assert_eq!(x.denom, BigInt::from(1));
        assert_eq!(x.coord[0], BigInt::from(2));
        assert_eq!(x.coord[1], BigInt::from(3));
        assert_eq!(x.coord[2], BigInt::from(4));
        assert_eq!(x.coord[3], BigInt::from(5));

        // Negative denom gets flipped
        let mut y = QuatAlgElem {
            denom: BigInt::from(-3),
            coord: IbzVec4([
                BigInt::from(6),
                BigInt::from(9),
                BigInt::from(12),
                BigInt::from(15),
            ]),
        };
        quat_alg_normalize(&mut y);
        assert!(y.denom > Ibz::zero());
    }

    #[test]
    fn test_elem_equal() {
        let a = quat_alg_elem_set(1, 1, 2, 3, 4);
        let b = quat_alg_elem_set(1, 1, 2, 3, 4);
        assert!(quat_alg_elem_equal(&a, &b));

        // Same value different denom: (2, [2,4,6,8]) == (1, [1,2,3,4])
        let c = quat_alg_elem_set(2, 2, 4, 6, 8);
        assert!(quat_alg_elem_equal(&a, &c));

        let d = quat_alg_elem_set(1, 1, 2, 3, 5);
        assert!(!quat_alg_elem_equal(&a, &d));
    }

    #[test]
    fn test_scalar() {
        let s = quat_alg_scalar(&BigInt::from(5), &BigInt::from(3));
        assert_eq!(s.coord[0], BigInt::from(5));
        assert!(s.coord[1].is_zero());
        assert!(s.coord[2].is_zero());
        assert!(s.coord[3].is_zero());
        assert_eq!(s.denom, BigInt::from(3));
    }
}
