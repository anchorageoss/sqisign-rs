//! Algebraic-property tests for Level 3 Fp arithmetic.

mod common;

use common::{DetRng, ITER};
use sqisign_verify::fp::Fp;
use sqisign_verify::params::Level3;

fn eq(a: &Fp<Level3>, b: &Fp<Level3>) -> bool {
    bool::from(a.ct_equal(b))
}

fn is_zero(a: &Fp<Level3>) -> bool {
    bool::from(a.ct_is_zero())
}

#[test]
fn l3_gfp_encode_kat_zero_one_two() {
    let expected_zero = [0u8; 48];
    let zero = Fp::<Level3>::zero();
    let bytes = zero.encode();
    assert_eq!(&bytes[..], &expected_zero[..]);

    let mut expected_one = [0u8; 48];
    expected_one[0] = 1;
    let one = Fp::<Level3>::one();
    let bytes = one.encode();
    assert_eq!(&bytes[..], &expected_one[..]);

    let mut expected_two = [0u8; 48];
    expected_two[0] = 2;
    let two = Fp::<Level3>::from_small(2);
    let bytes = two.encode();
    assert_eq!(&bytes[..], &expected_two[..]);
}

#[test]
fn l3_gfp_encode_kat_p_minus_one() {
    // p = 65 * 2^376 - 1 => p - 1 = 65 * 2^376 - 2
    // In LE bytes: [0xFE, 0xFF x 46, 0x40]
    let mut expected = [0xFFu8; 48];
    expected[0] = 0xFE;
    expected[47] = 0x40;
    let p_minus_one = Fp::<Level3>::zero().sub(&Fp::<Level3>::one());
    let bytes = p_minus_one.encode();
    assert_eq!(&bytes[..], &expected[..]);
}

#[test]
fn l3_gfp_equality() {
    let mut rng = DetRng::new(b"l3_gfp_equality");
    let one = Fp::<Level3>::one();
    let zero = Fp::<Level3>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = a.add(&one);
        assert!(eq(&a, &a));
        assert!(!eq(&a, &b));
        assert!(is_zero(&zero));
        assert!(!is_zero(&one));
    }
}

#[test]
fn l3_gfp_addition() {
    let mut rng = DetRng::new(b"l3_gfp_addition");
    let zero = Fp::<Level3>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = rng.random_fp_level3();
        let c = rng.random_fp_level3();

        let lhs = a.add(&b).add(&c);
        let rhs = a.add(&b.add(&c));
        assert!(eq(&lhs, &rhs), "associativity");

        assert!(eq(&a.add(&b), &b.add(&a)), "commutativity");
        assert!(eq(&a.add(&zero), &a), "identity");
        assert!(is_zero(&a.add(&a.neg())), "inverse");
    }
}

#[test]
fn l3_gfp_subtraction() {
    let mut rng = DetRng::new(b"l3_gfp_subtraction");
    let zero = Fp::<Level3>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = rng.random_fp_level3();
        let c = rng.random_fp_level3();

        let lhs = a.sub(&b).sub(&c);
        let rhs = a.sub(&b.add(&c));
        assert!(eq(&lhs, &rhs));

        assert!(eq(&a.sub(&b), &b.sub(&a).neg()));
        assert!(eq(&a.sub(&zero), &a));
        assert!(is_zero(&a.sub(&a)));
    }
}

#[test]
fn l3_gfp_multiplication() {
    let mut rng = DetRng::new(b"l3_gfp_multiplication");
    let one = Fp::<Level3>::one();
    let zero = Fp::<Level3>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = rng.random_fp_level3();
        let c = rng.random_fp_level3();

        let lhs = a.mul(&b).mul(&c);
        let rhs = a.mul(&b.mul(&c));
        assert!(eq(&lhs, &rhs), "associativity");

        let lhs = a.mul(&b.add(&c));
        let rhs = a.mul(&b).add(&a.mul(&c));
        assert!(eq(&lhs, &rhs), "distributivity");

        assert!(eq(&a.mul(&b), &b.mul(&a)), "commutativity");
        assert!(eq(&a.mul(&one), &a), "identity");
        assert!(is_zero(&a.mul(&zero)), "zero");
    }
}

#[test]
fn l3_gfp_squaring() {
    let mut rng = DetRng::new(b"l3_gfp_squaring");

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        assert!(eq(&a.sqr(), &a.mul(&a)));
    }
    assert!(is_zero(&Fp::<Level3>::zero().sqr()));
}

#[test]
fn l3_gfp_inversion() {
    let mut rng = DetRng::new(b"l3_gfp_inversion");
    let one = Fp::<Level3>::one();
    let zero = Fp::<Level3>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let inv = a.inv();
        assert!(eq(&a.mul(&inv), &one));
    }

    assert!(eq(&zero.inv(), &zero));
}

#[test]
fn l3_gfp_sqrt_and_is_square() {
    let mut rng = DetRng::new(b"l3_gfp_sqrt");

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let c = a.sqr();
        assert!(bool::from(c.is_square()), "a^2 is always a QR");

        let r = c.sqrt();
        let neg_r = r.neg();
        assert!(eq(&a, &r) || eq(&a, &neg_r));
    }
}

#[test]
fn l3_gfp_encode_decode_roundtrip() {
    let mut rng = DetRng::new(b"l3_gfp_encode_decode");
    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let bytes = a.encode();
        let b = Fp::<Level3>::decode(bytes.as_ref()).expect("encoded form must decode");
        assert!(eq(&a, &b));
    }
}

#[test]
fn l3_gfp_decode_rejects_out_of_range() {
    let bytes = [0xffu8; 48];
    assert!(Fp::<Level3>::decode(&bytes).is_none());
}

#[test]
fn l3_gfp_half() {
    let mut rng = DetRng::new(b"l3_gfp_half");

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = a.add(&a);
        let c = b.half();
        assert!(eq(&a, &c));
    }
}

#[test]
fn l3_gfp_div3() {
    let mut rng = DetRng::new(b"l3_gfp_div3");

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = a.add(&a).add(&a);
        let c = b.div3();
        assert!(eq(&a, &c));
    }
    assert!(is_zero(&Fp::<Level3>::zero().div3()));
}

#[test]
fn l3_gfp_mul_small() {
    let mut rng = DetRng::new(b"l3_gfp_mul_small");

    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let val = rng.random_u32();
        let b = a.mul_small(val);
        let c = Fp::<Level3>::from_small(val as u64);
        let d = a.mul(&c);
        assert!(eq(&b, &d));
    }
}

#[test]
fn l3_gfp_select_and_cswap() {
    let mut rng = DetRng::new(b"l3_gfp_select");
    for _ in 0..ITER {
        let a = rng.random_fp_level3();
        let b = rng.random_fp_level3();

        assert!(eq(
            &Fp::<Level3>::select(&a, &b, subtle::Choice::from(0)),
            &a
        ));
        assert!(eq(
            &Fp::<Level3>::select(&a, &b, subtle::Choice::from(1)),
            &b
        ));

        let mut x = a.clone();
        let mut y = b.clone();
        x.cswap(&mut y, subtle::Choice::from(0));
        assert!(eq(&x, &a));
        assert!(eq(&y, &b));
        x.cswap(&mut y, subtle::Choice::from(1));
        assert!(eq(&x, &b));
        assert!(eq(&y, &a));
    }
}
