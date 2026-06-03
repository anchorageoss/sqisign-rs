//! Cross-validation tests for sqisign-theta Groups 6 and 7 (splitting, extraction, chain).
//!
//! Compares Rust output byte-for-byte against the reference implementation.
//! Reference output captured from tools/c-validate/theta_chain_cval.

use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level1;
use sqisign_verify::theta::splitting::{
    splitting_compute, theta_point_to_montgomery_point, theta_product_structure_to_elliptic_product,
};
use sqisign_verify::theta::theta_structure::is_product_theta_point;
use sqisign_verify::theta::{ThetaPoint, ThetaStructure};

type L1 = Level1;

fn fp2_hex(val: &Fp2<L1>) -> String {
    val.encode().iter().map(|b| format!("{:02x}", b)).collect()
}

fn make_product_null_point() -> ThetaPoint<L1> {
    ThetaPoint {
        x: Fp2::from_small(1),
        y: Fp2::from_small(2),
        z: Fp2::from_small(3),
        t: Fp2::from_small(6),
    }
}

fn make_product_theta_structure() -> ThetaStructure<L1> {
    ThetaStructure {
        null_point: make_product_null_point(),
        precomputation: false,
        ..ThetaStructure::default()
    }
}

// --- Section 1: splitting_compute on product theta structure ---

#[test]
fn test_splitting_compute() {
    let prod = make_product_theta_structure();

    // Verify is_product
    assert!(bool::from(is_product_theta_point(&prod.null_point)));

    let split = splitting_compute(&prod, -1, false, None);
    assert!(split.is_some(), "splitting_compute should succeed");

    let split = split.unwrap();

    // The identity transform (index 9) should be selected, leaving the null point unchanged
    assert_eq!(
        fp2_hex(&split.b.null_point.x),
        "01000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&split.b.null_point.y),
        "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&split.b.null_point.z),
        "03000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&split.b.null_point.t),
        "06000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );

    assert!(bool::from(is_product_theta_point(&split.b.null_point)));
}

// --- Section 2: theta_product_structure_to_elliptic_product ---

#[test]
fn test_theta_product_structure_to_elliptic_product() {
    let prod = make_product_theta_structure();

    let result = theta_product_structure_to_elliptic_product(&prod);
    assert!(result.is_some(), "product_to_elliptic should succeed");

    let (e1, e2) = result.unwrap();

    assert_eq!(
        fp2_hex(&e1.a),
        "5bffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&e1.c),
        "afffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&e2.a),
        "ddffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&e2.c),
        "f0ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
}

// --- Section 3: theta_point_to_montgomery_point ---

#[test]
fn test_theta_point_to_montgomery_point() {
    let prod = make_product_theta_structure();

    // Product point: (5, 7, 10, 14), satisfies 5*14 = 7*10 = 70
    let p = ThetaPoint {
        x: Fp2::from_small(5),
        y: Fp2::from_small(7),
        z: Fp2::from_small(10),
        t: Fp2::from_small(14),
    };

    let result = theta_point_to_montgomery_point(&p, &prod);
    assert!(result.is_some(), "point_to_montgomery should succeed");

    let couple = result.unwrap();

    assert_eq!(
        fp2_hex(&couple.p1.x),
        "19000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p1.z),
        "faffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p2.x),
        "11000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p2.z),
        "fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
}

#[test]
fn test_theta_point_to_montgomery_point_fallback() {
    let prod = make_product_theta_structure();

    // Fallback test: (0, 0, 13, 17), satisfies 0*17 = 0*13 = 0
    let q = ThetaPoint {
        x: Fp2::zero(),
        y: Fp2::zero(),
        z: Fp2::from_small(13),
        t: Fp2::from_small(17),
    };

    let result = theta_point_to_montgomery_point(&q, &prod);
    assert!(
        result.is_some(),
        "point_to_montgomery (fallback) should succeed"
    );

    let couple = result.unwrap();

    assert_eq!(
        fp2_hex(&couple.p1.x),
        "0d000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p1.z),
        "0d000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p2.x),
        "2b000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        fp2_hex(&couple.p2.z),
        "f6ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
}

// --- Section 4: splitting_compute with known zero_index ---

#[test]
fn test_splitting_compute_with_zero_index() {
    let prod = make_product_theta_structure();

    // zero_index=9 (identity transform) should succeed
    let result = splitting_compute(&prod, 9, false, None);
    assert!(
        result.is_some(),
        "splitting with zero_index=9 should succeed"
    );

    // zero_index=0 should fail (wrong index)
    let result = splitting_compute(&prod, 0, false, None);
    assert!(result.is_none(), "splitting with zero_index=0 should fail");
}

// --- Section 5: is_product_theta_point edge cases ---

#[test]
fn test_is_product_theta_point() {
    // (1, 2, 3, 6): 1*6 = 6, 2*3 = 6 → product
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(1),
        y: Fp2::from_small(2),
        z: Fp2::from_small(3),
        t: Fp2::from_small(6),
    };
    assert!(bool::from(is_product_theta_point(&p)));

    // (1, 2, 3, 7): 1*7 = 7, 2*3 = 6 → NOT product
    let q = ThetaPoint::<L1> {
        x: Fp2::from_small(1),
        y: Fp2::from_small(2),
        z: Fp2::from_small(3),
        t: Fp2::from_small(7),
    };
    assert!(!bool::from(is_product_theta_point(&q)));
}
