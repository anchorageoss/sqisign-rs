//! Cross-validation test: compare Rust dim2/dim4/integers operations
//! byte-for-byte against the reference output.
//!
//! Expected values come from running `tools/c-validate/dim_cval` which
//! links against the reference quaternion library.

use num_bigint::BigInt;
use num_traits::Zero;
use sqisign_rs::quaternion::dim2::*;
use sqisign_rs::quaternion::dim4::*;
use sqisign_rs::quaternion::intbig::{ibz_to_str, Ibz};
use sqisign_rs::quaternion::integers::*;
use sqisign_rs::quaternion::types::*;

fn to_hex(x: &Ibz) -> String {
    ibz_to_str(x, 16)
}

fn assert_hex(label: &str, got: &Ibz, expected: &str) {
    let h = to_hex(got);
    assert_eq!(h, expected, "{label}: got {h} expected {expected}");
}

// -----------------------------------------------------------------------
// Section 1: dim2
// -----------------------------------------------------------------------

#[test]
fn crossval_dim2_vec2_set() {
    let v = ibz_vec_2_set(42, -17);
    assert_hex("vec2_set_0", &v[0], "2a");
    assert_hex("vec2_set_1", &v[1], "-11");
}

#[test]
fn crossval_dim2_mat2_set_copy() {
    let m = ibz_mat_2x2_set(3, -7, 11, 5);
    let m2 = ibz_mat_2x2_copy(&m);
    assert_hex("mat2_set_00", &m2[0][0], "3");
    assert_hex("mat2_set_01", &m2[0][1], "-7");
    assert_hex("mat2_set_10", &m2[1][0], "b");
    assert_hex("mat2_set_11", &m2[1][1], "5");
}

#[test]
fn crossval_dim2_mat2_add() {
    let a = ibz_mat_2x2_set(10, 20, 30, 40);
    let b = ibz_mat_2x2_set(5, -3, 7, -9);
    let sum = ibz_mat_2x2_add(&a, &b);
    assert_hex("mat2_add_00", &sum[0][0], "f");
    assert_hex("mat2_add_01", &sum[0][1], "11");
    assert_hex("mat2_add_10", &sum[1][0], "25");
    assert_hex("mat2_add_11", &sum[1][1], "1f");
}

#[test]
fn crossval_dim2_det() {
    let det = ibz_mat_2x2_det_from_ibz(
        &BigInt::from(3),
        &BigInt::from(7),
        &BigInt::from(-2),
        &BigInt::from(5),
    );
    assert_hex("det2x2", &det, "1d");
}

#[test]
fn crossval_dim2_eval() {
    let m = ibz_mat_2x2_set(3, -7, 11, 5);
    let v = ibz_vec_2_set(4, -2);
    let res = ibz_mat_2x2_eval(&m, &v);
    assert_hex("mat2_eval_0", &res[0], "1a");
    assert_hex("mat2_eval_1", &res[1], "22");
}

#[test]
fn crossval_dim2_mulmod() {
    let a = ibz_mat_2x2_set(5, 3, -2, 7);
    let b = ibz_mat_2x2_set(1, -4, 6, 2);
    let m = BigInt::from(13);
    let prod = ibz_2x2_mul_mod(&a, &b, &m);
    assert_hex("mat2_mulmod_00", &prod[0][0], "a");
    assert_hex("mat2_mulmod_01", &prod[0][1], "c");
    assert_hex("mat2_mulmod_10", &prod[1][0], "1");
    assert_hex("mat2_mulmod_11", &prod[1][1], "9");
}

#[test]
fn crossval_dim2_inv_mod() {
    let m = ibz_mat_2x2_set(5, 3, -2, 7);
    let modulus = BigInt::from(13);
    let (inv, ok) = ibz_mat_2x2_inv_mod(&m, &modulus);
    assert!(ok, "mat2_inv_ok");
    assert_hex("mat2_inv_00", &inv[0][0], "a");
    assert_hex("mat2_inv_01", &inv[0][1], "5");
    assert_hex("mat2_inv_10", &inv[1][0], "1");
    assert_hex("mat2_inv_11", &inv[1][1], "9");
}

#[test]
fn crossval_dim2_inv_mod_noninvertible() {
    let m = ibz_mat_2x2_set(2, 3, 1, -2);
    let modulus = BigInt::from(7);
    let (_inv, ok) = ibz_mat_2x2_inv_mod(&m, &modulus);
    assert!(!ok, "mat2_inv_noninv_ok should be false");
}

// -----------------------------------------------------------------------
// Section 2: dim4
// -----------------------------------------------------------------------

#[test]
fn crossval_dim4_vec4_set() {
    let v = ibz_vec_4_set(100, -200, 300, -400);
    assert_hex("vec4_set_0", &v[0], "64");
    assert_hex("vec4_set_1", &v[1], "-c8");
    assert_hex("vec4_set_2", &v[2], "12c");
    assert_hex("vec4_set_3", &v[3], "-190");
}

#[test]
fn crossval_dim4_vec4_negate() {
    let v = ibz_vec_4_set(1, -2, 3, -4);
    let neg = ibz_vec_4_negate(&v);
    assert_hex("vec4_neg_0", &neg[0], "-1");
    assert_hex("vec4_neg_1", &neg[1], "2");
    assert_hex("vec4_neg_2", &neg[2], "-3");
    assert_hex("vec4_neg_3", &neg[3], "4");
}

#[test]
fn crossval_dim4_vec4_add() {
    let a = ibz_vec_4_set(10, 20, 30, 40);
    let b = ibz_vec_4_set(-5, 15, -25, 35);
    let sum = ibz_vec_4_add(&a, &b);
    assert_hex("vec4_add_0", &sum[0], "5");
    assert_hex("vec4_add_1", &sum[1], "23");
    assert_hex("vec4_add_2", &sum[2], "5");
    assert_hex("vec4_add_3", &sum[3], "4b");
}

#[test]
fn crossval_dim4_vec4_sub() {
    let a = ibz_vec_4_set(10, 20, 30, 40);
    let b = ibz_vec_4_set(-5, 15, -25, 35);
    let diff = ibz_vec_4_sub(&a, &b);
    assert_hex("vec4_sub_0", &diff[0], "f");
    assert_hex("vec4_sub_1", &diff[1], "5");
    assert_hex("vec4_sub_2", &diff[2], "37");
    assert_hex("vec4_sub_3", &diff[3], "5");
}

#[test]
fn crossval_dim4_vec4_scalar_mul() {
    let v = ibz_vec_4_set(3, -7, 11, -13);
    let scalar = BigInt::from(5);
    let prod = ibz_vec_4_scalar_mul(&scalar, &v);
    assert_hex("vec4_smul_0", &prod[0], "f");
    assert_hex("vec4_smul_1", &prod[1], "-23");
    assert_hex("vec4_smul_2", &prod[2], "37");
    assert_hex("vec4_smul_3", &prod[3], "-41");
}

#[test]
fn crossval_dim4_vec4_scalar_div() {
    let v = ibz_vec_4_set(15, -35, 55, -65);
    let scalar = BigInt::from(5);
    let (quot, ok) = ibz_vec_4_scalar_div(&scalar, &v);
    assert!(ok, "vec4_sdiv_ok");
    assert_hex("vec4_sdiv_0", &quot[0], "3");
    assert_hex("vec4_sdiv_1", &quot[1], "-7");
    assert_hex("vec4_sdiv_2", &quot[2], "b");
    assert_hex("vec4_sdiv_3", &quot[3], "-d");
}

#[test]
fn crossval_dim4_vec4_content() {
    let v = ibz_vec_4_set(12, -18, 24, -30);
    let content = ibz_vec_4_content(&v);
    assert_hex("vec4_content", &content, "6");
}

#[test]
fn crossval_dim4_vec4_linear_combination() {
    let a = ibz_vec_4_set(1, 2, 3, 4);
    let b = ibz_vec_4_set(5, 6, 7, 8);
    let ca = BigInt::from(3);
    let cb = BigInt::from(-2);
    let lc = ibz_vec_4_linear_combination(&ca, &a, &cb, &b);
    assert_hex("vec4_lincomb_0", &lc[0], "-7");
    assert_hex("vec4_lincomb_1", &lc[1], "-6");
    assert_hex("vec4_lincomb_2", &lc[2], "-5");
    assert_hex("vec4_lincomb_3", &lc[3], "-4");
}

#[test]
fn crossval_dim4_vec4_is_zero() {
    let v_zero = ibz_vec_4_set(0, 0, 0, 0);
    assert!(ibz_vec_4_is_zero(&v_zero), "vec4_iszero_yes");

    let v_nonzero = ibz_vec_4_set(0, 0, 1, 0);
    assert!(!ibz_vec_4_is_zero(&v_nonzero), "vec4_iszero_no");
}

#[test]
fn crossval_dim4_mat4_identity() {
    let id = ibz_mat_4x4_identity();
    assert!(ibz_mat_4x4_is_identity(&id), "mat4_isid_yes");

    let z = ibz_mat_4x4_zero();
    assert!(!ibz_mat_4x4_is_identity(&z), "mat4_isid_no");
}

fn make_seq_mat() -> IbzMat4x4 {
    let zero = BigInt::zero();
    let mut m = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = BigInt::from((i * 4 + j + 1) as i32);
        }
    }
    m
}

#[test]
fn crossval_dim4_mat4_mul() {
    let a = make_seq_mat();

    let zero = BigInt::zero();
    let mut b = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    b[0][0] = BigInt::from(2);
    b[0][2] = BigInt::from(1);
    b[1][1] = BigInt::from(3);
    b[1][3] = BigInt::from(1);
    b[2][0] = BigInt::from(1);
    b[2][2] = BigInt::from(2);
    b[3][1] = BigInt::from(1);
    b[3][3] = BigInt::from(3);

    let prod = ibz_mat_4x4_mul(&a, &b);

    let expected: [[&str; 4]; 4] = [
        ["5", "a", "7", "e"],
        ["11", "1a", "13", "1e"],
        ["1d", "2a", "1f", "2e"],
        ["29", "3a", "2b", "3e"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_mul_{i}{j}"), &prod[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_mat4_transpose() {
    let m = make_seq_mat();
    let mt = ibz_mat_4x4_transpose(&m);

    let expected: [[&str; 4]; 4] = [
        ["1", "5", "9", "d"],
        ["2", "6", "a", "e"],
        ["3", "7", "b", "f"],
        ["4", "8", "c", "10"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_trans_{i}{j}"), &mt[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_mat4_scalar_mul() {
    let m = make_seq_mat();
    let scalar = BigInt::from(-3);
    let prod = ibz_mat_4x4_scalar_mul(&scalar, &m);

    let expected: [[&str; 4]; 4] = [
        ["-3", "-6", "-9", "-c"],
        ["-f", "-12", "-15", "-18"],
        ["-1b", "-1e", "-21", "-24"],
        ["-27", "-2a", "-2d", "-30"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_smul_{i}{j}"), &prod[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_mat4_scalar_div() {
    let zero = BigInt::zero();
    let mut m = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = BigInt::from(((i * 4 + j + 1) * 6) as i32);
        }
    }
    let scalar = BigInt::from(3);
    let (quot, ok) = ibz_mat_4x4_scalar_div(&scalar, &m);
    assert!(ok, "mat4_sdiv_ok");

    let expected: [[&str; 4]; 4] = [
        ["2", "4", "6", "8"],
        ["a", "c", "e", "10"],
        ["12", "14", "16", "18"],
        ["1a", "1c", "1e", "20"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_sdiv_{i}{j}"), &quot[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_mat4_negate() {
    let m = make_seq_mat();
    let neg = ibz_mat_4x4_negate(&m);

    let expected: [[&str; 4]; 4] = [
        ["-1", "-2", "-3", "-4"],
        ["-5", "-6", "-7", "-8"],
        ["-9", "-a", "-b", "-c"],
        ["-d", "-e", "-f", "-10"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_neg_{i}{j}"), &neg[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_mat4_gcd() {
    let zero = BigInt::zero();
    let mut m = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = BigInt::from(((i * 4 + j + 1) * 6) as i32);
        }
    }
    let g = ibz_mat_4x4_gcd(&m);
    assert_hex("mat4_gcd", &g, "6");
}

#[test]
fn crossval_dim4_mat4_eval() {
    let m = make_seq_mat();
    let v = ibz_vec_4_set(1, -1, 2, -2);
    let res = ibz_mat_4x4_eval(&m, &v);
    assert_hex("mat4_eval_0", &res[0], "-3");
    assert_hex("mat4_eval_1", &res[1], "-3");
    assert_hex("mat4_eval_2", &res[2], "-3");
    assert_hex("mat4_eval_3", &res[3], "-3");
}

#[test]
fn crossval_dim4_mat4_eval_t() {
    let m = make_seq_mat();
    let v = ibz_vec_4_set(1, -1, 2, -2);
    let res = ibz_mat_4x4_eval_t(&v, &m);
    assert_hex("mat4_evalt_0", &res[0], "-c");
    assert_hex("mat4_evalt_1", &res[1], "-c");
    assert_hex("mat4_evalt_2", &res[2], "-c");
    assert_hex("mat4_evalt_3", &res[3], "-c");
}

#[test]
fn crossval_dim4_mat4_equal() {
    let a = make_seq_mat();
    let mut b = make_seq_mat();
    assert!(ibz_mat_4x4_equal(&a, &b), "mat4_equal_yes");
    b[2][3] = BigInt::from(999);
    assert!(!ibz_mat_4x4_equal(&a, &b), "mat4_equal_no");
}

#[test]
fn crossval_dim4_mat4_inv() {
    let zero = BigInt::zero();
    let mut m = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    // Upper triangular: [[2,1,3,0],[0,4,0,1],[0,0,3,2],[0,0,0,2]]
    m[0][0] = BigInt::from(2);
    m[0][1] = BigInt::from(1);
    m[0][2] = BigInt::from(3);
    m[1][1] = BigInt::from(4);
    m[1][3] = BigInt::from(1);
    m[2][2] = BigInt::from(3);
    m[2][3] = BigInt::from(2);
    m[3][3] = BigInt::from(2);

    let (inv, det, ok) = ibz_mat_4x4_inv_with_det_as_denom(&m);
    assert!(ok, "mat4_inv_ok");
    assert_hex("mat4_inv_det", &det, "30");

    let expected: [[&str; 4]; 4] = [
        ["18", "-6", "-18", "1b"],
        ["0", "c", "0", "-6"],
        ["0", "0", "10", "-10"],
        ["0", "0", "0", "18"],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_hex(&format!("mat4_inv_{i}{j}"), &inv[i][j], expected[i][j]);
        }
    }
}

#[test]
fn crossval_dim4_qf_eval() {
    let zero = BigInt::zero();
    let mut qf = IbzMat4x4([
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero.clone()],
        [zero.clone(), zero.clone(), zero.clone(), zero],
    ]);
    qf[0][0] = BigInt::from(1);
    qf[1][1] = BigInt::from(1);
    qf[2][2] = BigInt::from(3);
    qf[3][3] = BigInt::from(3);

    let coord = ibz_vec_4_set(2, 3, 1, -1);
    let result = quat_qf_eval(&qf, &coord);
    assert_hex("qf_eval", &result, "13");
}

// -----------------------------------------------------------------------
// Section 3: integers (Cornacchia)
// -----------------------------------------------------------------------

#[test]
fn crossval_cornacchia_1_5() {
    let n = BigInt::from(1);
    let p = BigInt::from(5);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_1_5 should succeed");
    assert_hex("corn_1_5_x", &x, "2");
    assert_hex("corn_1_5_y", &y, "1");
}

#[test]
fn crossval_cornacchia_1_2() {
    let n = BigInt::from(1);
    let p = BigInt::from(2);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_1_2 should succeed");
    assert_hex("corn_1_2_x", &x, "1");
    assert_hex("corn_1_2_y", &y, "1");
}

#[test]
fn crossval_cornacchia_1_41() {
    let n = BigInt::from(1);
    let p = BigInt::from(41);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_1_41 should succeed");
    assert_hex("corn_1_41_x", &x, "5");
    assert_hex("corn_1_41_y", &y, "4");
}

#[test]
fn crossval_cornacchia_2_3() {
    let n = BigInt::from(2);
    let p = BigInt::from(3);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_2_3 should succeed");
    assert_hex("corn_2_3_x", &x, "1");
    assert_hex("corn_2_3_y", &y, "1");
}

#[test]
fn crossval_cornacchia_3_7() {
    let n = BigInt::from(3);
    let p = BigInt::from(7);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_3_7 should succeed");
    assert_hex("corn_3_7_x", &x, "2");
    assert_hex("corn_3_7_y", &y, "1");
}

#[test]
fn crossval_cornacchia_1_7_no_solution() {
    let n = BigInt::from(1);
    let p = BigInt::from(7);
    assert!(
        ibz_cornacchia_prime(&n, &p).is_none(),
        "corn_1_7 should have no solution"
    );
}

#[test]
fn crossval_cornacchia_1_104729() {
    let n = BigInt::from(1);
    let p = BigInt::from(104729);
    let (x, y) = ibz_cornacchia_prime(&n, &p).expect("corn_1_104729 should succeed");
    assert_hex("corn_1_104729_x", &x, "143");
    assert_hex("corn_1_104729_y", &y, "14");
    let check = &x * &x + &y * &y;
    assert_eq!(
        check, p,
        "corn_1_104729 verification: x^2+y^2 should equal p"
    );
}
