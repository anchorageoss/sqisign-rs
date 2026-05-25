//! Cross-validation test: compare Rust LLL/normeq/lat_ball operations
//! byte-for-byte against the reference output.
//!
//! Expected values come from running `tools/c-validate/lll_cval` which
//! links against the reference quaternion library.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use sqisign_rs::quaternion::algebra::*;
use sqisign_rs::quaternion::dim4::*;
use sqisign_rs::quaternion::intbig::{ibz_to_str, Ibz};
use sqisign_rs::quaternion::lat_ball::*;
use sqisign_rs::quaternion::lll::*;
use sqisign_rs::quaternion::normeq::*;
use sqisign_rs::quaternion::rational::*;
use sqisign_rs::quaternion::types::*;

fn to_hex(x: &Ibz) -> String {
    ibz_to_str(x, 16)
}

fn assert_hex(label: &str, got: &Ibz, expected: &str) {
    let h = to_hex(got);
    assert_eq!(h, expected, "{label}: got {h} expected {expected}");
}

fn make_o0() -> QuatLattice {
    let mut lat = QuatLattice {
        denom: BigInt::from(2),
        ..QuatLattice::default()
    };
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = BigInt::zero();
        }
    }
    lat.basis[0][0] = BigInt::from(2);
    lat.basis[1][1] = BigInt::from(2);
    lat.basis[2][2] = BigInt::from(1);
    lat.basis[1][2] = BigInt::from(1);
    lat.basis[3][3] = BigInt::from(1);
    lat.basis[0][3] = BigInt::from(1);
    lat
}

// -----------------------------------------------------------------------
// Section 1: quat_lattice_lll
// -----------------------------------------------------------------------

#[test]
fn crossval_lattice_lll_p103() {
    let alg = quat_alg_init_set_ui(103);
    let lat = make_o0();
    let red = quat_lattice_lll(&lat, &alg);

    assert_hex("red_00", &red[0][0], "2");
    assert_hex("red_01", &red[0][1], "0");
    assert_hex("red_02", &red[0][2], "0");
    assert_hex("red_03", &red[0][3], "1");
    assert_hex("red_10", &red[1][0], "0");
    assert_hex("red_11", &red[1][1], "2");
    assert_hex("red_12", &red[1][2], "1");
    assert_hex("red_13", &red[1][3], "0");
    assert_hex("red_20", &red[2][0], "0");
    assert_hex("red_21", &red[2][1], "0");
    assert_hex("red_22", &red[2][2], "1");
    assert_hex("red_23", &red[2][3], "0");
    assert_hex("red_30", &red[3][0], "0");
    assert_hex("red_31", &red[3][1], "0");
    assert_hex("red_32", &red[3][2], "0");
    assert_hex("red_33", &red[3][3], "1");
}

#[test]
fn crossval_lattice_lll_p11() {
    let alg = quat_alg_init_set_ui(11);
    let lat = make_o0();
    let red = quat_lattice_lll(&lat, &alg);

    assert_hex("red_00", &red[0][0], "2");
    assert_hex("red_01", &red[0][1], "0");
    assert_hex("red_02", &red[0][2], "0");
    assert_hex("red_03", &red[0][3], "1");
    assert_hex("red_10", &red[1][0], "0");
    assert_hex("red_11", &red[1][1], "2");
    assert_hex("red_12", &red[1][2], "1");
    assert_hex("red_13", &red[1][3], "0");
    assert_hex("red_20", &red[2][0], "0");
    assert_hex("red_21", &red[2][1], "0");
    assert_hex("red_22", &red[2][2], "1");
    assert_hex("red_23", &red[2][3], "0");
    assert_hex("red_30", &red[3][0], "0");
    assert_hex("red_31", &red[3][1], "0");
    assert_hex("red_32", &red[3][2], "0");
    assert_hex("red_33", &red[3][3], "1");
}

// -----------------------------------------------------------------------
// Section 2: quat_lll_verify
// -----------------------------------------------------------------------

#[test]
fn crossval_lll_verify() {
    let alg = quat_alg_init_set_ui(103);
    let lat = make_o0();
    let red = quat_lattice_lll(&lat, &alg);

    let delta = ibq_set(&BigInt::from(99), &BigInt::from(100)).unwrap();
    let eta = ibq_set(&BigInt::from(51), &BigInt::from(100)).unwrap();

    let verified = quat_lll_verify(&red, &delta, &eta, &alg);
    assert!(verified, "LLL-reduced basis should pass verification");
}

// -----------------------------------------------------------------------
// Section 3: quat_lideal_reduce_basis
// -----------------------------------------------------------------------

#[test]
fn crossval_reduce_basis() {
    let alg = quat_alg_init_set_ui(103);
    let parent = make_o0();

    let lideal = QuatLeftIdeal {
        lattice: make_o0(),
        parent_order: parent,
        norm: BigInt::one(),
    };

    let (red, gram) = quat_lideal_reduce_basis(&lideal, &alg);

    assert_hex("red_00", &red[0][0], "2");
    assert_hex("red_01", &red[0][1], "0");
    assert_hex("red_02", &red[0][2], "0");
    assert_hex("red_03", &red[0][3], "1");
    assert_hex("red_10", &red[1][0], "0");
    assert_hex("red_11", &red[1][1], "2");
    assert_hex("red_12", &red[1][2], "1");
    assert_hex("red_13", &red[1][3], "0");
    assert_hex("red_20", &red[2][0], "0");
    assert_hex("red_21", &red[2][1], "0");
    assert_hex("red_22", &red[2][2], "1");
    assert_hex("red_23", &red[2][3], "0");
    assert_hex("red_30", &red[3][0], "0");
    assert_hex("red_31", &red[3][1], "0");
    assert_hex("red_32", &red[3][2], "0");
    assert_hex("red_33", &red[3][3], "1");

    assert_hex("gram_00", &gram[0][0], "4");
    assert_hex("gram_01", &gram[0][1], "0");
    assert_hex("gram_02", &gram[0][2], "0");
    assert_hex("gram_03", &gram[0][3], "0");
    assert_hex("gram_10", &gram[1][0], "0");
    assert_hex("gram_11", &gram[1][1], "4");
    assert_hex("gram_12", &gram[1][2], "0");
    assert_hex("gram_13", &gram[1][3], "0");
    assert_hex("gram_20", &gram[2][0], "0");
    assert_hex("gram_21", &gram[2][1], "4");
    assert_hex("gram_22", &gram[2][2], "68");
    assert_hex("gram_23", &gram[2][3], "0");
    assert_hex("gram_30", &gram[3][0], "4");
    assert_hex("gram_31", &gram[3][1], "0");
    assert_hex("gram_32", &gram[3][2], "0");
    assert_hex("gram_33", &gram[3][3], "68");
}

// -----------------------------------------------------------------------
// Section 4: quat_lattice_O0_set
// -----------------------------------------------------------------------

#[test]
fn crossval_lattice_o0_set() {
    let o0 = quat_lattice_o0_set();

    assert_hex("o0_denom", &o0.denom, "2");
    assert_hex("o0_00", &o0.basis[0][0], "2");
    assert_hex("o0_01", &o0.basis[0][1], "0");
    assert_hex("o0_02", &o0.basis[0][2], "0");
    assert_hex("o0_03", &o0.basis[0][3], "1");
    assert_hex("o0_10", &o0.basis[1][0], "0");
    assert_hex("o0_11", &o0.basis[1][1], "2");
    assert_hex("o0_12", &o0.basis[1][2], "1");
    assert_hex("o0_13", &o0.basis[1][3], "0");
    assert_hex("o0_20", &o0.basis[2][0], "0");
    assert_hex("o0_21", &o0.basis[2][1], "0");
    assert_hex("o0_22", &o0.basis[2][2], "1");
    assert_hex("o0_23", &o0.basis[2][3], "0");
    assert_hex("o0_30", &o0.basis[3][0], "0");
    assert_hex("o0_31", &o0.basis[3][1], "0");
    assert_hex("o0_32", &o0.basis[3][2], "0");
    assert_hex("o0_33", &o0.basis[3][3], "1");
}

// -----------------------------------------------------------------------
// Section 5: quat_change_to_O0_basis
// -----------------------------------------------------------------------

#[test]
fn crossval_change_to_o0_basis() {
    let el = QuatAlgElem {
        coord: ibz_vec_4_set(2, 7, 1, -4),
        denom: BigInt::from(2),
    };

    let vec = quat_change_to_o0_basis(&el);

    assert_hex("vec_0", &vec[0], "3");
    assert_hex("vec_1", &vec[1], "3");
    assert_hex("vec_2", &vec[2], "1");
    assert_hex("vec_3", &vec[3], "-4");
}

// -----------------------------------------------------------------------
// Section 6: quat_lattice_bound_parallelogram
// -----------------------------------------------------------------------

#[test]
fn crossval_bound_parallelogram() {
    let mut g = IbzMat4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            g[i][j] = if i == j {
                BigInt::from(4)
            } else {
                BigInt::zero()
            };
        }
    }

    let radius = BigInt::from(100);
    let (box_out, u, ok) = quat_lattice_bound_parallelogram(&g, &radius);

    assert!(ok, "bound_para should return true (non-trivial)");

    assert_hex("box_0", &box_out[0], "5");
    assert_hex("box_1", &box_out[1], "5");
    assert_hex("box_2", &box_out[2], "5");
    assert_hex("box_3", &box_out[3], "5");

    assert_hex("U_00", &u[0][0], "1");
    assert_hex("U_01", &u[0][1], "0");
    assert_hex("U_02", &u[0][2], "0");
    assert_hex("U_03", &u[0][3], "0");
    assert_hex("U_10", &u[1][0], "0");
    assert_hex("U_11", &u[1][1], "1");
    assert_hex("U_12", &u[1][2], "0");
    assert_hex("U_13", &u[1][3], "0");
    assert_hex("U_20", &u[2][0], "0");
    assert_hex("U_21", &u[2][1], "0");
    assert_hex("U_22", &u[2][2], "1");
    assert_hex("U_23", &u[2][3], "0");
    assert_hex("U_30", &u[3][0], "0");
    assert_hex("U_31", &u[3][1], "0");
    assert_hex("U_32", &u[3][2], "0");
    assert_hex("U_33", &u[3][3], "1");
}
