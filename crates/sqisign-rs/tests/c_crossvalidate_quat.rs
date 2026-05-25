//! Cross-validation test: compare Rust algebra/hnf/lattice/ideal operations
//! byte-for-byte against the reference output.
//!
//! Expected values come from running `tools/c-validate/quat_cval` which
//! links against the reference quaternion library.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use sqisign_rs::quaternion::algebra::*;
use sqisign_rs::quaternion::dim4::*;
use sqisign_rs::quaternion::hnf::*;
use sqisign_rs::quaternion::ideal::*;
use sqisign_rs::quaternion::intbig::{ibz_to_str, Ibz};
use sqisign_rs::quaternion::lattice::*;
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
// Section 1: algebra
// -----------------------------------------------------------------------

#[test]
fn crossval_coord_mul_p7() {
    let alg = quat_alg_init_set_ui(7);
    let a = ibz_vec_4_set(152, 57, 190, 28);
    let b = ibz_vec_4_set(165, 35, 231, 770);
    let c = quat_alg_coord_mul(&a, &b, &alg);
    assert_hex("coord_mul_p7_0", &c[0], "-6a379");
    assert_hex("coord_mul_p7_1", &c[1], "f290d");
    assert_hex("coord_mul_p7_2", &c[2], "5c00");
    assert_hex("coord_mul_p7_3", &c[3], "1f4b1");
}

#[test]
fn crossval_coord_mul_p11() {
    let mut alg = quat_alg_init_set_ui(7);
    alg.p = BigInt::from(11);
    let a = ibz_vec_4_set(152, 57, 190, 28);
    let b = ibz_vec_4_set(165, 35, 231, 770);
    let c = quat_alg_coord_mul(&a, &b, &alg);
    assert_hex("coord_mul_p11_0", &c[0], "-aa221");
    assert_hex("coord_mul_p11_1", &c[1], "17b1ed");
    assert_hex("coord_mul_p11_2", &c[2], "5c00");
    assert_hex("coord_mul_p11_3", &c[3], "1f4b1");
}

#[test]
fn crossval_mul_p7() {
    let alg = quat_alg_init_set_ui(7);
    let a = quat_alg_elem_copy_ibz(
        &BigInt::from(76),
        &BigInt::from(152),
        &BigInt::from(57),
        &BigInt::from(190),
        &BigInt::from(28),
    );
    let b = quat_alg_elem_copy_ibz(
        &BigInt::from(385),
        &BigInt::from(165),
        &BigInt::from(35),
        &BigInt::from(231),
        &BigInt::from(770),
    );
    let c = quat_alg_mul(&a, &b, &alg);
    assert_hex("mul_p7_denom", &c.denom, "724c");
    assert_hex("mul_p7_coord0", &c.coord[0], "-6a379");
    assert_hex("mul_p7_coord1", &c.coord[1], "f290d");
    assert_hex("mul_p7_coord2", &c.coord[2], "5c00");
    assert_hex("mul_p7_coord3", &c.coord[3], "1f4b1");
}

#[test]
fn crossval_add() {
    let a = quat_alg_elem_copy_ibz(
        &BigInt::from(9),
        &BigInt::from(-12),
        &BigInt::from(0),
        &BigInt::from(-7),
        &BigInt::from(19),
    );
    let b1 = quat_alg_elem_copy_ibz(
        &BigInt::from(3),
        &BigInt::from(-6),
        &BigInt::from(2),
        &BigInt::from(7),
        &BigInt::from(-19),
    );
    let c1 = quat_alg_add(&a, &b1);
    assert_hex("add_1_denom", &c1.denom, "9");
    assert_hex("add_1_coord0", &c1.coord[0], "-1e");
    assert_hex("add_1_coord1", &c1.coord[1], "6");
    assert_hex("add_1_coord2", &c1.coord[2], "e");
    assert_hex("add_1_coord3", &c1.coord[3], "-26");

    let b2 = quat_alg_elem_copy_ibz(
        &BigInt::from(6),
        &BigInt::from(-6),
        &BigInt::from(2),
        &BigInt::from(7),
        &BigInt::from(-19),
    );
    let c2 = quat_alg_add(&a, &b2);
    assert_hex("add_2_denom", &c2.denom, "12");
    assert_hex("add_2_coord0", &c2.coord[0], "-2a");
    assert_hex("add_2_coord1", &c2.coord[1], "6");
    assert_hex("add_2_coord2", &c2.coord[2], "7");
    assert_hex("add_2_coord3", &c2.coord[3], "-13");
}

#[test]
fn crossval_sub() {
    let a = quat_alg_elem_copy_ibz(
        &BigInt::from(9),
        &BigInt::from(-12),
        &BigInt::from(0),
        &BigInt::from(-7),
        &BigInt::from(19),
    );
    let b = quat_alg_elem_copy_ibz(
        &BigInt::from(3),
        &BigInt::from(-6),
        &BigInt::from(2),
        &BigInt::from(7),
        &BigInt::from(-19),
    );
    let c = quat_alg_sub(&a, &b);
    assert_hex("sub_1_denom", &c.denom, "9");
    assert_hex("sub_1_coord0", &c.coord[0], "6");
    assert_hex("sub_1_coord1", &c.coord[1], "-6");
    assert_hex("sub_1_coord2", &c.coord[2], "-1c");
    assert_hex("sub_1_coord3", &c.coord[3], "4c");
}

#[test]
fn crossval_norm() {
    let alg_11 = quat_alg_init_set_ui(11);
    let a1 = quat_alg_elem_copy_ibz(
        &BigInt::from(2),
        &BigInt::from(1),
        &BigInt::from(5),
        &BigInt::from(7),
        &BigInt::from(2),
    );
    let (num1, den1) = quat_alg_norm(&a1, &alg_11);
    assert_hex("norm1_num", &num1, "261");
    assert_hex("norm1_denom", &den1, "4");

    let a2 = quat_alg_elem_copy_ibz(
        &BigInt::from(76),
        &BigInt::from(152),
        &BigInt::from(57),
        &BigInt::from(190),
        &BigInt::from(28),
    );
    let (num2, den2) = quat_alg_norm(&a2, &alg_11);
    assert_hex("norm2_num", &num2, "697cd");
    assert_hex("norm2_denom", &den2, "1690");

    let alg_7 = quat_alg_init_set_ui(7);
    let (num3, den3) = quat_alg_norm(&a2, &alg_7);
    assert_hex("norm3_num", &num3, "4577d");
    assert_hex("norm3_denom", &den3, "1690");
}

#[test]
fn crossval_conj() {
    let a = quat_alg_elem_copy_ibz(
        &BigInt::from(25),
        &BigInt::from(-125),
        &BigInt::from(2),
        &BigInt::from(0),
        &BigInt::from(-30),
    );
    let c = quat_alg_conj(&a);
    assert_hex("conj_1_denom", &c.denom, "19");
    assert_hex("conj_1_coord0", &c.coord[0], "-7d");
    assert_hex("conj_1_coord1", &c.coord[1], "-2");
    assert_hex("conj_1_coord2", &c.coord[2], "0");
    assert_hex("conj_1_coord3", &c.coord[3], "1e");
}

#[test]
fn crossval_normalize() {
    let mut x1 = quat_alg_elem_copy_ibz(
        &BigInt::from(48),
        &BigInt::from(-36),
        &BigInt::from(18),
        &BigInt::from(0),
        &BigInt::from(-300),
    );
    quat_alg_normalize(&mut x1);
    assert_hex("normalize_1_denom", &x1.denom, "8");
    assert_hex("normalize_1_coord0", &x1.coord[0], "-6");
    assert_hex("normalize_1_coord1", &x1.coord[1], "3");
    assert_hex("normalize_1_coord2", &x1.coord[2], "0");
    assert_hex("normalize_1_coord3", &x1.coord[3], "-32");

    let mut x2 = quat_alg_elem_copy_ibz(
        &BigInt::from(-6),
        &BigInt::from(-36),
        &BigInt::from(18),
        &BigInt::from(0),
        &BigInt::from(-300),
    );
    quat_alg_normalize(&mut x2);
    assert_hex("normalize_2_denom", &x2.denom, "1");
    assert_hex("normalize_2_coord0", &x2.coord[0], "6");
    assert_hex("normalize_2_coord1", &x2.coord[1], "-3");
    assert_hex("normalize_2_coord2", &x2.coord[2], "0");
    assert_hex("normalize_2_coord3", &x2.coord[3], "32");
}

// -----------------------------------------------------------------------
// Section 2: HNF
// -----------------------------------------------------------------------

#[test]
fn crossval_hnf() {
    let mut lat = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = BigInt::zero();
        }
    }
    lat.basis[0][0] = BigInt::from(1);
    lat.basis[0][3] = BigInt::from(-1);
    lat.basis[1][1] = BigInt::from(-2);
    lat.basis[2][2] = BigInt::from(1);
    lat.basis[2][1] = BigInt::from(1);
    lat.basis[3][3] = BigInt::from(-3);
    lat.denom = BigInt::from(6);
    let result = quat_lattice_hnf(&lat);
    assert_hex("hnf_1_denom", &result.denom, "6");
    assert_hex("hnf_1_basis_00", &result.basis[0][0], "1");
    assert_hex("hnf_1_basis_01", &result.basis[0][1], "0");
    assert_hex("hnf_1_basis_02", &result.basis[0][2], "0");
    assert_hex("hnf_1_basis_03", &result.basis[0][3], "0");
    assert_hex("hnf_1_basis_10", &result.basis[1][0], "0");
    assert_hex("hnf_1_basis_11", &result.basis[1][1], "2");
    assert_hex("hnf_1_basis_12", &result.basis[1][2], "0");
    assert_hex("hnf_1_basis_13", &result.basis[1][3], "0");
    assert_hex("hnf_1_basis_20", &result.basis[2][0], "0");
    assert_hex("hnf_1_basis_21", &result.basis[2][1], "0");
    assert_hex("hnf_1_basis_22", &result.basis[2][2], "1");
    assert_hex("hnf_1_basis_23", &result.basis[2][3], "0");
    assert_hex("hnf_1_basis_30", &result.basis[3][0], "0");
    assert_hex("hnf_1_basis_31", &result.basis[3][1], "0");
    assert_hex("hnf_1_basis_32", &result.basis[3][2], "0");
    assert_hex("hnf_1_basis_33", &result.basis[3][3], "3");
}

// -----------------------------------------------------------------------
// Section 3: lattice
// -----------------------------------------------------------------------

fn make_lat_pair() -> (QuatLattice, QuatLattice) {
    let mut lat1 = QuatLattice::default();
    let mut lat2 = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat1.basis[i][j] = BigInt::zero();
            lat2.basis[i][j] = BigInt::zero();
        }
    }
    lat1.basis[0][0] = BigInt::from(4);
    lat1.basis[0][2] = BigInt::from(3);
    lat2.basis[0][0] = BigInt::from(1);
    lat2.basis[0][3] = BigInt::from(-1);
    lat1.basis[1][1] = BigInt::from(5);
    lat2.basis[1][1] = BigInt::from(-2);
    lat1.basis[2][2] = BigInt::from(3);
    lat2.basis[2][2] = BigInt::from(1);
    lat2.basis[2][1] = BigInt::from(1);
    lat1.basis[3][3] = BigInt::from(7);
    lat2.basis[3][3] = BigInt::from(-3);
    lat1.denom = BigInt::from(4);
    lat2.denom = BigInt::from(6);
    (lat1, lat2)
}

#[test]
fn crossval_lattice_add() {
    let (lat1, lat2) = make_lat_pair();
    let sum = quat_lattice_add(&lat1, &lat2);
    assert_hex("lattice_add_denom", &sum.denom, "c");
    assert_hex("lattice_add_basis_00", &sum.basis[0][0], "2");
    assert_hex("lattice_add_basis_02", &sum.basis[0][2], "1");
    assert_hex("lattice_add_basis_11", &sum.basis[1][1], "1");
    assert_hex("lattice_add_basis_22", &sum.basis[2][2], "1");
    assert_hex("lattice_add_basis_33", &sum.basis[3][3], "3");
    // off-diagonal zeros
    assert_hex("lattice_add_basis_01", &sum.basis[0][1], "0");
    assert_hex("lattice_add_basis_10", &sum.basis[1][0], "0");
}

#[test]
fn crossval_lattice_intersect() {
    let (lat1, lat2) = make_lat_pair();
    let lat1h = quat_lattice_hnf(&lat1);
    let lat2h = quat_lattice_hnf(&lat2);
    let inter = quat_lattice_intersect(&lat1h, &lat2h);
    assert_hex("lattice_inter_denom", &inter.denom, "2");
    assert_hex("lattice_inter_basis_00", &inter.basis[0][0], "2");
    assert_hex("lattice_inter_basis_02", &inter.basis[0][2], "1");
    assert_hex("lattice_inter_basis_11", &inter.basis[1][1], "a");
    assert_hex("lattice_inter_basis_22", &inter.basis[2][2], "3");
    assert_hex("lattice_inter_basis_33", &inter.basis[3][3], "7");
}

#[test]
fn crossval_lattice_index() {
    let mut sublat = QuatLattice::default();
    let mut overlat = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            sublat.basis[i][j] = BigInt::zero();
            overlat.basis[i][j] = BigInt::zero();
        }
    }
    overlat.basis[0][0] = BigInt::one();
    overlat.basis[1][1] = BigInt::one();
    overlat.basis[2][2] = BigInt::one();
    overlat.basis[3][3] = BigInt::one();
    overlat.denom = BigInt::from(2);

    sublat.basis[0][0] = BigInt::from(2);
    sublat.basis[0][2] = BigInt::from(1);
    sublat.basis[1][1] = BigInt::from(4);
    sublat.basis[1][2] = BigInt::from(2);
    sublat.basis[1][3] = BigInt::from(3);
    sublat.basis[2][2] = BigInt::from(1);
    sublat.basis[3][3] = BigInt::from(1);
    sublat.denom = BigInt::from(2);

    let idx = quat_lattice_index(&sublat, &overlat);
    assert_hex("lattice_index", &idx, "8");
}

#[test]
fn crossval_lattice_contains() {
    let mut lat = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = BigInt::zero();
        }
    }
    lat.basis[0][0] = BigInt::from(1);
    lat.basis[0][3] = BigInt::from(-1);
    lat.basis[1][1] = BigInt::from(-2);
    lat.basis[2][2] = BigInt::from(1);
    lat.basis[2][1] = BigInt::from(1);
    lat.basis[3][3] = BigInt::from(-3);
    lat.denom = BigInt::from(6);
    let lat = quat_lattice_hnf(&lat);

    let x = quat_alg_elem_copy_ibz(
        &BigInt::from(3),
        &BigInt::from(1),
        &BigInt::from(-2),
        &BigInt::from(26),
        &BigInt::from(9),
    );
    let coord = quat_lattice_contains(&lat, &x);
    assert!(coord.is_some(), "lattice_contains should succeed");
    let coord = coord.unwrap();
    assert_hex("lattice_contains_coord_0", &coord[0], "2");
    assert_hex("lattice_contains_coord_1", &coord[1], "-2");
    assert_hex("lattice_contains_coord_2", &coord[2], "34");
    assert_hex("lattice_contains_coord_3", &coord[3], "6");
}

#[test]
fn crossval_lattice_conjugate() {
    let mut lat = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = BigInt::zero();
        }
    }
    lat.basis[0][0] = BigInt::from(4);
    lat.basis[0][3] = BigInt::from(1);
    lat.basis[1][1] = BigInt::from(-2);
    lat.basis[2][2] = BigInt::from(-1);
    lat.basis[2][1] = BigInt::from(-1);
    lat.basis[3][3] = BigInt::from(-3);
    lat.denom = BigInt::from(6);
    let lat = quat_lattice_hnf(&lat);
    let conj = quat_lattice_conjugate_without_hnf(&lat);
    let conj = quat_lattice_hnf(&conj);
    assert_hex("lattice_conj_denom", &conj.denom, "6");
    assert_hex("lattice_conj_basis_00", &conj.basis[0][0], "4");
    assert_hex("lattice_conj_basis_03", &conj.basis[0][3], "1");
    assert_hex("lattice_conj_basis_11", &conj.basis[1][1], "2");
    assert_hex("lattice_conj_basis_22", &conj.basis[2][2], "1");
    assert_hex("lattice_conj_basis_33", &conj.basis[3][3], "3");
}

#[test]
fn crossval_lattice_dual() {
    let mut lat = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = BigInt::zero();
        }
    }
    lat.basis[0][0] = BigInt::from(1);
    lat.basis[0][3] = BigInt::from(-1);
    lat.basis[1][1] = BigInt::from(-2);
    lat.basis[2][2] = BigInt::from(1);
    lat.basis[2][1] = BigInt::from(1);
    lat.basis[3][3] = BigInt::from(-3);
    lat.denom = BigInt::from(6);
    let lat = quat_lattice_hnf(&lat);
    let dual = quat_lattice_dual_without_hnf(&lat);
    let dual = quat_lattice_hnf(&dual);
    assert_hex("lattice_dual_denom", &dual.denom, "1");
    assert_hex("lattice_dual_basis_00", &dual.basis[0][0], "6");
    assert_hex("lattice_dual_basis_11", &dual.basis[1][1], "3");
    assert_hex("lattice_dual_basis_22", &dual.basis[2][2], "6");
    assert_hex("lattice_dual_basis_33", &dual.basis[3][3], "2");
}

#[test]
fn crossval_lattice_mul() {
    let alg = quat_alg_init_set_ui(19);
    let mut lat1 = QuatLattice::default();
    let mut lat2 = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            lat1.basis[i][j] = BigInt::zero();
            lat2.basis[i][j] = BigInt::zero();
        }
    }
    lat1.basis[0][0] = BigInt::from(44);
    lat1.basis[0][2] = BigInt::from(3);
    lat1.basis[0][3] = BigInt::from(32);
    lat2.basis[0][0] = BigInt::from(1);
    lat1.basis[1][1] = BigInt::from(5);
    lat2.basis[1][1] = BigInt::from(2);
    lat1.basis[2][2] = BigInt::from(3);
    lat2.basis[2][2] = BigInt::from(1);
    lat1.basis[3][3] = BigInt::from(1);
    lat2.basis[3][3] = BigInt::from(3);
    lat1.denom = BigInt::from(4);
    lat2.denom = BigInt::from(6);
    let prod = quat_lattice_mul(&lat1, &lat2, &alg);
    assert_hex("lattice_mul_denom", &prod.denom, "18");
    assert_hex("lattice_mul_basis_00", &prod.basis[0][0], "1");
    assert_hex("lattice_mul_basis_11", &prod.basis[1][1], "1");
    assert_hex("lattice_mul_basis_22", &prod.basis[2][2], "1");
    assert_hex("lattice_mul_basis_33", &prod.basis[3][3], "1");
    // Off-diagonal zeros
    assert_hex("lattice_mul_basis_01", &prod.basis[0][1], "0");
    assert_hex("lattice_mul_basis_10", &prod.basis[1][0], "0");
}

#[test]
fn crossval_lattice_gram() {
    let alg = quat_alg_init_set_ui(103);
    let o0 = make_o0();
    let gram = quat_lattice_gram(&o0, &alg);
    assert_hex("lattice_gram_00", &gram[0][0], "8");
    assert_hex("lattice_gram_01", &gram[0][1], "0");
    assert_hex("lattice_gram_02", &gram[0][2], "0");
    assert_hex("lattice_gram_03", &gram[0][3], "4");
    assert_hex("lattice_gram_10", &gram[1][0], "0");
    assert_hex("lattice_gram_11", &gram[1][1], "8");
    assert_hex("lattice_gram_12", &gram[1][2], "4");
    assert_hex("lattice_gram_13", &gram[1][3], "0");
    assert_hex("lattice_gram_20", &gram[2][0], "0");
    assert_hex("lattice_gram_21", &gram[2][1], "4");
    assert_hex("lattice_gram_22", &gram[2][2], "d0");
    assert_hex("lattice_gram_23", &gram[2][3], "0");
    assert_hex("lattice_gram_30", &gram[3][0], "4");
    assert_hex("lattice_gram_31", &gram[3][1], "0");
    assert_hex("lattice_gram_32", &gram[3][2], "0");
    assert_hex("lattice_gram_33", &gram[3][3], "d0");
}

// -----------------------------------------------------------------------
// Section 4: ideal
// -----------------------------------------------------------------------

#[test]
fn crossval_lideal_create_principal() {
    let alg = quat_alg_init_set_ui(367);
    let order = make_o0();
    let gamma = quat_alg_elem_copy_ibz(
        &BigInt::from(1),
        &BigInt::from(219),
        &BigInt::from(200),
        &BigInt::from(78),
        &BigInt::from(-1),
    );
    let ideal = quat_lideal_create_principal(&gamma, &order, &alg);
    assert_hex("principal_norm", &ideal.norm, "236b04");
    assert_hex("principal_lat_denom", &ideal.lattice.denom, "1");
    assert_hex(
        "principal_lat_basis_00",
        &ideal.lattice.basis[0][0],
        "11b582",
    );
    assert_hex("principal_lat_basis_01", &ideal.lattice.basis[0][1], "0");
    assert_hex(
        "principal_lat_basis_02",
        &ideal.lattice.basis[0][2],
        "4bb6e",
    );
    assert_hex(
        "principal_lat_basis_03",
        &ideal.lattice.basis[0][3],
        "eec81",
    );
    assert_hex("principal_lat_basis_10", &ideal.lattice.basis[1][0], "0");
    assert_hex(
        "principal_lat_basis_11",
        &ideal.lattice.basis[1][1],
        "11b582",
    );
    assert_hex(
        "principal_lat_basis_12",
        &ideal.lattice.basis[1][2],
        "2c901",
    );
    assert_hex(
        "principal_lat_basis_13",
        &ideal.lattice.basis[1][3],
        "4bb6e",
    );
    assert_hex("principal_lat_basis_20", &ideal.lattice.basis[2][0], "0");
    assert_hex("principal_lat_basis_21", &ideal.lattice.basis[2][1], "0");
    assert_hex("principal_lat_basis_22", &ideal.lattice.basis[2][2], "1");
    assert_hex("principal_lat_basis_23", &ideal.lattice.basis[2][3], "0");
    assert_hex("principal_lat_basis_30", &ideal.lattice.basis[3][0], "0");
    assert_hex("principal_lat_basis_31", &ideal.lattice.basis[3][1], "0");
    assert_hex("principal_lat_basis_32", &ideal.lattice.basis[3][2], "0");
    assert_hex("principal_lat_basis_33", &ideal.lattice.basis[3][3], "1");
}

#[test]
fn crossval_lideal_create() {
    let alg = quat_alg_init_set_ui(367);
    let order = make_o0();
    let gamma = quat_alg_elem_copy_ibz(
        &BigInt::from(1),
        &BigInt::from(219),
        &BigInt::from(200),
        &BigInt::from(78),
        &BigInt::from(-1),
    );
    let n = BigInt::from(31);
    let ideal = quat_lideal_create(&gamma, &n, &order, &alg).unwrap();
    assert_hex("create_norm", &ideal.norm, "1f");
    assert_hex("create_lat_denom", &ideal.lattice.denom, "2");
    assert_hex("create_lat_basis_00", &ideal.lattice.basis[0][0], "3e");
    assert_hex("create_lat_basis_01", &ideal.lattice.basis[0][1], "0");
    assert_hex("create_lat_basis_02", &ideal.lattice.basis[0][2], "2");
    assert_hex("create_lat_basis_03", &ideal.lattice.basis[0][3], "3d");
    assert_hex("create_lat_basis_10", &ideal.lattice.basis[1][0], "0");
    assert_hex("create_lat_basis_11", &ideal.lattice.basis[1][1], "3e");
    assert_hex("create_lat_basis_12", &ideal.lattice.basis[1][2], "1");
    assert_hex("create_lat_basis_13", &ideal.lattice.basis[1][3], "2");
    assert_hex("create_lat_basis_22", &ideal.lattice.basis[2][2], "1");
    assert_hex("create_lat_basis_33", &ideal.lattice.basis[3][3], "1");
}

#[test]
fn crossval_lideal_add_inter_equals() {
    let alg = quat_alg_init_set_ui(103);
    let order = make_o0();

    let gen1 = quat_alg_elem_copy_ibz(
        &BigInt::from(1),
        &BigInt::from(3),
        &BigInt::from(5),
        &BigInt::from(7),
        &BigInt::from(11),
    );
    let n1 = BigInt::from(17);
    let lideal1 = quat_lideal_create(&gen1, &n1, &order, &alg).unwrap();

    // self-intersection should be stable
    let inter = quat_lideal_inter(&lideal1, &lideal1);
    assert!(
        quat_lideal_equals(&inter, &lideal1),
        "self-intersection should be stable"
    );
    assert_hex("selfinter_norm", &inter.norm, "11");
}

#[test]
fn crossval_order_discriminant_maximal() {
    let alg = quat_alg_init_set_ui(43);
    let o0 = make_o0();
    let (ok, disc) = quat_order_discriminant(&o0, &alg);
    assert!(ok);
    assert_hex("disc_O0", &disc, "2b");
    assert!(quat_order_is_maximal(&o0, &alg), "O0 should be maximal");

    // Z^4 is not maximal
    let mut z4 = QuatLattice::default();
    for i in 0..4 {
        for j in 0..4 {
            z4.basis[i][j] = if i == j {
                BigInt::one()
            } else {
                BigInt::zero()
            };
        }
    }
    z4.denom = BigInt::one();
    assert!(
        !quat_order_is_maximal(&z4, &alg),
        "Z^4 should not be maximal"
    );
}

fn make_order_p19() -> QuatLattice {
    let mut order = QuatLattice::default();
    order.basis[0][0] = BigInt::from(4);
    order.basis[0][1] = BigInt::from(0);
    order.basis[0][2] = BigInt::from(2);
    order.basis[0][3] = BigInt::from(2);
    order.basis[1][0] = BigInt::from(0);
    order.basis[1][1] = BigInt::from(8);
    order.basis[1][2] = BigInt::from(4);
    order.basis[1][3] = BigInt::from(3);
    order.basis[2][0] = BigInt::from(0);
    order.basis[2][1] = BigInt::from(0);
    order.basis[2][2] = BigInt::from(2);
    order.basis[2][3] = BigInt::from(0);
    order.basis[3][0] = BigInt::from(0);
    order.basis[3][1] = BigInt::from(0);
    order.basis[3][2] = BigInt::from(0);
    order.basis[3][3] = BigInt::from(1);
    order.denom = BigInt::from(4);
    order
}

#[test]
fn crossval_lideal_right_order() {
    let alg = quat_alg_init_set_ui(19);
    let order = make_order_p19();
    let gen = quat_alg_elem_set(1, 3, 3, 0, 1);
    let norm = BigInt::from(15);
    let lideal = quat_lideal_create(&gen, &norm, &order, &alg).unwrap();
    let rorder = quat_lideal_right_order(&lideal, &alg);

    assert_hex("rorder_denom", &rorder.denom, "4");
    assert_hex("rorder_basis_00", &rorder.basis[0][0], "4");
    assert_hex("rorder_basis_02", &rorder.basis[0][2], "2");
    assert_hex("rorder_basis_03", &rorder.basis[0][3], "2");
    assert_hex("rorder_basis_11", &rorder.basis[1][1], "8");
    assert_hex("rorder_basis_12", &rorder.basis[1][2], "4");
    assert_hex("rorder_basis_13", &rorder.basis[1][3], "3");
    assert_hex("rorder_basis_22", &rorder.basis[2][2], "2");
    assert_hex("rorder_basis_33", &rorder.basis[3][3], "1");

    assert!(
        ibz_mat_4x4_is_hnf(&rorder.basis),
        "right order should be in HNF"
    );

    let one = quat_alg_elem_set(1, 1, 0, 0, 0);
    assert!(
        quat_lattice_contains(&rorder, &one).is_some(),
        "right order should contain 1"
    );
}

#[test]
fn crossval_lideal_inverse_lattice() {
    let alg = quat_alg_init_set_ui(19);
    let order = make_order_p19();
    let gen = quat_alg_elem_set(1, 2, 3, 0, 1);
    let norm = BigInt::from(15);
    let lideal = quat_lideal_create(&gen, &norm, &order, &alg).unwrap();
    let inv = quat_lideal_inverse_lattice_without_hnf(&lideal);
    let prod = quat_lattice_mul(&lideal.lattice, &inv, &alg);
    assert!(
        quat_lattice_equal(&prod, &order),
        "I * Ibar/N should equal the parent order"
    );
}
