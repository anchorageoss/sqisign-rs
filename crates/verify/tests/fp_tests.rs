//! Algebraic-property tests for Fp arithmetic.
//!
//! Each test exercises one identity over the field (commutativity,
//! associativity, distributivity, inversion round-trip, sqrt squaring
//! back, etc.) over many random inputs drawn from a deterministic
//! SHAKE256 stream, so results are reproducible. A handful of
//! byte-level known-answer tests near the bottom pin down the wire
//! format independently of the algebraic properties.

mod common;

use common::{DetRng, ITER};
use sqisign_verify::fp::Fp;
use sqisign_verify::params::Level1;

/// Convert the constant-time `Choice` returned by [`Fp::ct_equal`]
/// into a plain bool for use in test assertions.
fn eq(a: &Fp<Level1>, b: &Fp<Level1>) -> bool {
    bool::from(a.ct_equal(b))
}

/// Convert the constant-time `Choice` returned by [`Fp::ct_is_zero`]
/// into a plain bool for use in test assertions.
fn is_zero(a: &Fp<Level1>) -> bool {
    bool::from(a.ct_is_zero())
}

// ---------------------------------------------------------------------
// Equality
// ---------------------------------------------------------------------

#[test]
fn gfp_equality() {
    let mut rng = DetRng::new(b"gfp_equality");
    let one = Fp::<Level1>::one();
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let b = a.add(&one);
        let c = Fp::<Level1>::zero();

        assert!(eq(&a, &a));
        assert!(!eq(&a, &b));
        assert!(eq(&c, &zero));
        assert!(is_zero(&zero));
        assert!(!is_zero(&one));
    }
}

// ---------------------------------------------------------------------
// Multiplication by a small integer
// ---------------------------------------------------------------------

#[test]
fn gfp_mul_small() {
    let mut rng = DetRng::new(b"gfp_mul_small");

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let val = rng.random_u32();

        let b = a.mul_small(val);
        let c = Fp::<Level1>::from_small(val as u64);
        let d = a.mul(&c);

        assert!(eq(&b, &d));
    }
}

// ---------------------------------------------------------------------
// Half (division by 2)
// ---------------------------------------------------------------------

#[test]
fn gfp_half() {
    let mut rng = DetRng::new(b"gfp_half");

    for _ in 0..ITER {
        let a = rng.random_fp_level1();

        // half(a+a) == a
        let b = a.add(&a);
        let c = b.half();
        assert!(eq(&a, &c));

        // a == half(a) + half(a)
        let a = rng.random_fp_level1();
        let b = a.half();
        let c = b.add(&b);
        assert!(eq(&a, &c));
    }
}

// ---------------------------------------------------------------------
// Division by 3
// ---------------------------------------------------------------------

#[test]
fn gfp_div3() {
    let mut rng = DetRng::new(b"gfp_div3");

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        // div3(a+a+a) == a
        let b = a.add(&a).add(&a);
        let c = b.div3();
        assert!(eq(&a, &c));
    }

    // 0/3 == 0
    let zero = Fp::<Level1>::zero();
    let d = zero.div3();
    assert!(is_zero(&d));
}

// ---------------------------------------------------------------------
// Construction from a small integer
// ---------------------------------------------------------------------

#[test]
fn gfp_set_small() {
    let one = Fp::<Level1>::one();
    let b = one.add(&one);
    let two = Fp::<Level1>::from_small(2);
    assert!(eq(&b, &two));

    let four_via_add = one.add(&one).add(&one).add(&one);
    let four = Fp::<Level1>::from_small(4);
    assert!(eq(&four_via_add, &four));
}

// ---------------------------------------------------------------------
// Addition (associative, commutative, zero, negation)
// ---------------------------------------------------------------------

#[test]
fn gfp_addition() {
    let mut rng = DetRng::new(b"gfp_addition");
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let b = rng.random_fp_level1();
        let c = rng.random_fp_level1();

        // (a+b)+c == a+(b+c)
        let lhs = a.add(&b).add(&c);
        let rhs = a.add(&b.add(&c));
        assert!(eq(&lhs, &rhs));

        // a+b == b+a
        assert!(eq(&a.add(&b), &b.add(&a)));

        // a + 0 == a
        assert!(eq(&a.add(&zero), &a));

        // a + (-a) == 0
        let na = a.neg();
        assert!(is_zero(&a.add(&na)));
    }
}

// ---------------------------------------------------------------------
// Subtraction
// ---------------------------------------------------------------------

#[test]
fn gfp_subtraction() {
    let mut rng = DetRng::new(b"gfp_subtraction");
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let b = rng.random_fp_level1();
        let c = rng.random_fp_level1();

        // (a-b)-c == a-(b+c)
        let lhs = a.sub(&b).sub(&c);
        let rhs = a.sub(&b.add(&c));
        assert!(eq(&lhs, &rhs));

        // a-b == -(b-a)
        let d = a.sub(&b);
        let e = b.sub(&a).neg();
        assert!(eq(&d, &e));

        // a - 0 == a
        assert!(eq(&a.sub(&zero), &a));

        // a - a == 0
        assert!(is_zero(&a.sub(&a)));
    }
}

// ---------------------------------------------------------------------
// Multiplication
// ---------------------------------------------------------------------

#[test]
fn gfp_multiplication() {
    let mut rng = DetRng::new(b"gfp_multiplication");
    let one = Fp::<Level1>::one();
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let b = rng.random_fp_level1();
        let c = rng.random_fp_level1();

        // (a*b)*c == a*(b*c)
        let lhs = a.mul(&b).mul(&c);
        let rhs = a.mul(&b.mul(&c));
        assert!(eq(&lhs, &rhs));

        // a*(b+c) == a*b + a*c
        let lhs = a.mul(&b.add(&c));
        let rhs = a.mul(&b).add(&a.mul(&c));
        assert!(eq(&lhs, &rhs));

        // a*b == b*a
        assert!(eq(&a.mul(&b), &b.mul(&a)));

        // a*1 == a
        assert!(eq(&a.mul(&one), &a));

        // a*0 == 0
        assert!(is_zero(&a.mul(&zero)));
    }
}

// ---------------------------------------------------------------------
// Squaring
// ---------------------------------------------------------------------

#[test]
fn gfp_squaring() {
    let mut rng = DetRng::new(b"gfp_squaring");
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        // sqr(a) == a*a
        assert!(eq(&a.sqr(), &a.mul(&a)));
    }

    // sqr(0) == 0
    assert!(is_zero(&zero.sqr()));
}

// ---------------------------------------------------------------------
// Inversion
// ---------------------------------------------------------------------

#[test]
fn gfp_inversion() {
    let mut rng = DetRng::new(b"gfp_inversion");
    let one = Fp::<Level1>::one();
    let zero = Fp::<Level1>::zero();

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let inv = a.inv();
        let c = a.mul(&inv);
        assert!(eq(&c, &one));
    }

    // 0^{-1} returns 0 by convention (no panic on zero input).
    let inv_zero = zero.inv();
    assert!(eq(&inv_zero, &zero));
}

// ---------------------------------------------------------------------
// Square root and quadratic residue detection
// ---------------------------------------------------------------------

#[test]
fn gfp_sqrt_and_is_square() {
    let mut rng = DetRng::new(b"gfp_sqrt");

    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let c = a.sqr();
        assert!(bool::from(c.is_square()), "a^2 is always a QR");

        let r = c.sqrt();
        let neg_r = r.neg();
        // Either r == a or r == -a (sqrt is well-defined up to sign).
        assert!(eq(&a, &r) || eq(&a, &neg_r));
    }
}

// ---------------------------------------------------------------------
// Encode / decode round-trip.
// ---------------------------------------------------------------------

#[test]
fn gfp_encode_decode_roundtrip() {
    let mut rng = DetRng::new(b"gfp_encode_decode");
    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let bytes = a.encode();
        let b = Fp::<Level1>::decode(bytes.as_ref()).expect("encoded form must decode");
        assert!(eq(&a, &b));
    }
}

// ---------------------------------------------------------------------
// decode rejects out-of-range input. The largest canonical value below
// p is p-1; an all-0xFF byte string represents 2^256 - 1 > p.
// ---------------------------------------------------------------------

#[test]
fn gfp_decode_rejects_out_of_range() {
    let bytes = [0xffu8; 32];
    assert!(Fp::<Level1>::decode(&bytes).is_none());
}

// ---------------------------------------------------------------------
// Byte-level known-answer tests pinning down the canonical wire format.
// Encoding 0, 1, 2 must produce specific bytes; decoding those bytes
// must round-trip. These catch any drift in Montgomery handling or
// endianness that the algebraic-identity tests would not notice.
// ---------------------------------------------------------------------

/// Compare the encoded form of `a` against the expected 32 bytes.
fn assert_encoded(a: &Fp<Level1>, expected: &[u8; 32]) {
    let bytes = a.encode();
    let got: &[u8] = &bytes[..];
    assert_eq!(got, expected as &[u8]);
}

#[test]
fn gfp_encode_kat_zero_one_two() {
    assert_encoded(&Fp::<Level1>::zero(), &[0u8; 32]);

    let mut expected_one = [0u8; 32];
    expected_one[0] = 1;
    assert_encoded(&Fp::<Level1>::one(), &expected_one);

    let mut expected_two = [0u8; 32];
    expected_two[0] = 2;
    assert_encoded(&Fp::<Level1>::from_small(2), &expected_two);

    let mut expected_42 = [0u8; 32];
    expected_42[0] = 42;
    assert_encoded(&Fp::<Level1>::from_small(42), &expected_42);
}

/// `p - 1` should encode to canonical bytes `[0xFE, 0xFF, ..., 0xFF, 0x04]`
/// (LE encoding of `5*2^248 - 2`).
#[test]
fn gfp_encode_kat_p_minus_one() {
    let mut expected = [0xFFu8; 32];
    expected[0] = 0xFE;
    expected[31] = 0x04;
    let p_minus_one = Fp::<Level1>::zero().sub(&Fp::<Level1>::one());
    assert_encoded(&p_minus_one, &expected);
}

/// `decode_reduce` on the all-ones 32-byte string equals the canonical
/// representation of `(2^256 - 1) mod p`. We have:
///   2^256 = 256 * 2^248 = 51 * (5*2^248) + 2^248 = 51 * (p + 1) + 2^248
/// so 2^256 - 1 == 2^248 + 50 (mod p).
/// In LE-canonical bytes that is byte[0]=50=0x32, byte[31]=0x01, rest 0.
#[test]
fn gfp_decode_reduce_kat_all_ones() {
    let bytes = [0xFFu8; 32];
    let a = Fp::<Level1>::decode_reduce(&bytes);
    let mut expected = [0u8; 32];
    expected[0] = 0x32;
    expected[31] = 0x01;
    assert_encoded(&a, &expected);
}

// ---------------------------------------------------------------------
// Conditional select / swap.
// ---------------------------------------------------------------------

#[test]
fn gfp_select_and_cswap() {
    let mut rng = DetRng::new(b"gfp_select");
    for _ in 0..ITER {
        let a = rng.random_fp_level1();
        let b = rng.random_fp_level1();

        assert!(eq(
            &Fp::<Level1>::select(&a, &b, subtle::Choice::from(0)),
            &a
        ));
        assert!(eq(
            &Fp::<Level1>::select(&a, &b, subtle::Choice::from(1)),
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
