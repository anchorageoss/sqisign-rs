//! Regression test for Bug 4: Montgomery identity point must be (1:0).
//!
//! The point at infinity on a Montgomery curve in projective (X:Z) coordinates
//! is (1:0), not (0:0). Using (0:0) breaks the Montgomery ladder because
//! xDBL((0:0)) = (0:0) instead of the correct identity behavior.

use sqisign_verify::ec::EcPoint;
use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level1;

type L1 = Level1;

#[test]
fn regression_identity_is_one_zero() {
    let id = EcPoint::<L1>::identity();

    // X must be 1, not 0
    assert!(
        !bool::from(id.x.ct_is_zero()),
        "identity X coordinate must be nonzero"
    );
    assert!(
        bool::from(id.x.ct_equal(&Fp2::<L1>::one())),
        "identity X coordinate must be 1"
    );

    // Z must be 0
    assert!(
        bool::from(id.z.ct_is_zero()),
        "identity Z coordinate must be 0"
    );
}

#[test]
fn regression_identity_detected_as_infinity() {
    let id = EcPoint::<L1>::identity();
    assert!(
        bool::from(id.is_zero()),
        "identity must be detected as the point at infinity"
    );
}

#[test]
fn regression_default_is_identity() {
    let def = EcPoint::<L1>::default();
    let id = EcPoint::<L1>::identity();
    assert!(
        bool::from(def.ct_equal(&id)),
        "Default::default() must equal identity()"
    );
}

#[test]
fn regression_identity_x_is_nonzero() {
    // (0:0) is not a valid projective point. The identity constructor
    // must return (1:0), not (0:0), because xDBL((0:0)) = (0:0) forever
    // while xDBL((1:0)) correctly produces the identity.
    let id = EcPoint::<L1>::identity();
    assert!(
        !bool::from(id.x.ct_is_zero()),
        "identity X must not be zero, (0:0) breaks the Montgomery ladder"
    );
}

#[test]
fn regression_identity_projective_equality() {
    // Any (c:0) with c != 0 represents the same point at infinity as (1:0)
    let id = EcPoint::<L1>::identity();
    let scaled = EcPoint::<L1>::new(Fp2::from_small(42), Fp2::zero());
    assert!(
        bool::from(id.ct_equal(&scaled)),
        "(1:0) and (42:0) must be projectively equal"
    );
}
