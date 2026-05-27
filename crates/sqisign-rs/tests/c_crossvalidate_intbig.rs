//! Cross-validation test: compare Rust intbig operations byte-for-byte
//! against the reference output.
//!
//! The expected values come from running `tools/c-validate/intbig_cval`
//! which links against the reference GMP-based intbig implementation.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use sqisign_rs::quaternion::intbig::*;

fn from_hex(s: &str) -> Ibz {
    ibz_set_from_str(s, 16).unwrap()
}

fn to_hex(x: &Ibz) -> String {
    ibz_to_str(x, 16)
}

fn assert_hex_eq(label: &str, got: &Ibz, expected_hex: &str) {
    let got_hex = to_hex(got);
    assert_eq!(
        got_hex, expected_hex,
        "{}: got {} expected {}",
        label, got_hex, expected_hex
    );
}

// -----------------------------------------------------------------------
// Section 1: Basic arithmetic
// -----------------------------------------------------------------------

#[test]
fn crossval_basic_arithmetic() {
    let a = from_hex("deadbeef12345678cafebabe");
    let b = from_hex("1111111122222222333333334444444455555555");

    // add
    let c = &a + &b;
    assert_hex_eq("add", &c, "111111112222222311e0f22256789abd20541013");

    // sub
    let c = &a - &b;
    assert_hex_eq("sub", &c, "-111111112222222154857444320fedcb8a569a97");

    // mul
    let c = &a * &b;
    assert_hex_eq(
        "mul",
        &c,
        "ed8620ffefe016d61dbc99c2ea89668fb7563356f2fe3a2e554de8f11ab1716",
    );

    // neg
    let c = -&a;
    assert_hex_eq("neg_a", &c, "-deadbeef12345678cafebabe");

    // abs of negative
    let c = c.abs();
    assert_hex_eq("abs_neg_a", &c, "deadbeef12345678cafebabe");
}

// -----------------------------------------------------------------------
// Section 2: Division
// -----------------------------------------------------------------------

#[test]
fn crossval_division() {
    // Truncated division
    let a = from_hex("aaaaaaaabbbbbbbbccccccccdddddddd");
    let b = from_hex("1111111122222222");
    let (q, r) = ibz_div(&a, &b);
    assert_hex_eq("tdiv_q", &q, "9fffffff70000001d");
    assert_hex_eq("tdiv_r", &r, "1111110e00000003");

    // Truncated division with negative dividend
    let neg_a = -&a;
    let (q, r) = ibz_div(&neg_a, &b);
    assert_hex_eq("tdiv_negdiv_q", &q, "-9fffffff70000001d");
    assert_hex_eq("tdiv_negdiv_r", &r, "-1111110e00000003");

    // Floor division with negative dividend
    let (q, r) = ibz_div_floor(&neg_a, &b);
    assert_hex_eq("fdiv_negdiv_q", &q, "-9fffffff70000001e");
    assert_hex_eq("fdiv_negdiv_r", &r, "32222221f");

    // ibz_mod: always non-negative
    let a = from_hex("-deadbeefcafebabe");
    let b = from_hex("1234567890abcdef");
    let r = ibz_mod(&a, &b);
    assert_hex_eq("mod_neg", &r, "dfaa52f8dbaba65");

    // div_2exp
    let a = from_hex("ffffffffffffffffffffffffffffffff");
    let q = ibz_div_2exp(&a, 64);
    assert_hex_eq("div2exp_64", &q, "ffffffffffffffff");

    // div_2exp negative
    let neg_a = -&a;
    let q = ibz_div_2exp(&neg_a, 17);
    assert_hex_eq("div2exp_neg_17", &q, "-7fffffffffffffffffffffffffff");
}

// -----------------------------------------------------------------------
// Section 3: Number theory
// -----------------------------------------------------------------------

#[test]
fn crossval_number_theory() {
    // GCD
    let a = from_hex("3b9aca00");
    let b = from_hex("e8d4a51000");
    let g = ibz_gcd(&a, &b);
    assert_hex_eq("gcd", &g, "3b9aca00");

    // pow
    let a = from_hex("ff");
    let c = ibz_pow(&a, 7);
    assert_hex_eq("pow_ff_7", &c, "f914dd22eb06ff");

    // pow_mod
    let a = from_hex("deadbeef");
    let b = from_hex("1234567890");
    let m = from_hex("ffffffffffffffc5");
    let q = ibz_pow_mod(&a, &b, &m);
    assert_hex_eq("powmod", &q, "a161b2252a1780df");

    // invmod
    let a = from_hex("deadbeef12345678");
    let m = from_hex("ffffffffffffffc5");
    let inv = ibz_invmod(&a, &m).expect("invmod should succeed");
    assert_hex_eq("invmod", &inv, "9258ddc1f43bcea2");

    // Legendre
    let p = from_hex("ffffffffffffffc5");
    assert_eq!(ibz_legendre(&BigInt::from(3), &p), -1, "legendre_3");
    assert_eq!(ibz_legendre(&BigInt::from(2), &p), -1, "legendre_2");

    // two_adic
    let a = from_hex("abcdef0000000000");
    assert_eq!(ibz_two_adic(&a), 40, "two_adic");

    // sqrt_floor
    let a = from_hex("10000000000000000"); // 2^64
    let s = ibz_sqrt_floor(&a);
    assert_hex_eq("sqrt_floor", &s, "100000000");

    // sqrt (perfect square)
    let base = from_hex("1000000000000"); // 2^48
    let a = &base * &base; // 2^96
    let s = ibz_sqrt(&a).expect("should be perfect square");
    assert_hex_eq("sqrt_perfect", &s, "1000000000000");
}

#[test]
fn crossval_sqrt_mod_p() {
    // p ≡ 1 mod 8 (Tonelli-Shanks): p = 2^64 - 59
    let p = from_hex("ffffffffffffffc5");
    let a = BigInt::from(9);
    let s = ibz_sqrt_mod_p(&a, &p).expect("sqrt should exist");
    let sq = ibz_pow_mod(&s, &BigInt::from(2), &p);
    assert_hex_eq("sqrtmodp_9_sq", &sq, "9");

    // p ≡ 3 mod 4: p = 31
    let p = from_hex("1f");
    let a = BigInt::from(4);
    let s = ibz_sqrt_mod_p(&a, &p).expect("sqrt should exist");
    let sq = ibz_pow_mod(&s, &BigInt::from(2), &p);
    assert_hex_eq("sqrtmodp_p3m4_sq", &sq, "4");

    // p ≡ 5 mod 8: p = 13
    let p = from_hex("d");
    let a = BigInt::from(4);
    let s = ibz_sqrt_mod_p(&a, &p).expect("sqrt should exist");
    let sq = ibz_pow_mod(&s, &BigInt::from(2), &p);
    assert_hex_eq("sqrtmodp_p5m8_sq", &sq, "4");
}

// -----------------------------------------------------------------------
// Section 4: Digit conversion
// -----------------------------------------------------------------------

#[test]
fn crossval_digit_conversion() {
    // copy_digits: [0x0002, 0x0001] → 0x10000000000000002
    let val = ibz_copy_digits(&[0x0000000000000002u64, 0x0000000000000001u64]);
    assert_hex_eq("copy_digits_2_1", &val, "10000000000000002");

    // to_digits round-trip
    let mut out = [0u64; 2];
    ibz_to_digits(&val, &mut out);
    assert_eq!(out[0], 0x0000000000000002, "to_digits_0");
    assert_eq!(out[1], 0x0000000000000001, "to_digits_1");
}

// -----------------------------------------------------------------------
// Section 5: Comparison / predicates
// -----------------------------------------------------------------------

#[test]
fn crossval_comparison() {
    assert!(BigInt::from(0).is_zero(), "is_zero_0");
    assert!(BigInt::from(1).is_one(), "is_one_1");
    assert!(!BigInt::from(1).is_even(), "is_even_1");
    assert!(BigInt::from(1).is_odd(), "is_odd_1");
    assert!(BigInt::from(42).is_even(), "is_even_42");
    assert!(!BigInt::from(42).is_odd(), "is_odd_42");
    assert_eq!(ibz_bitsize(&BigInt::from(42)), 6, "bitsize_42");

    let a = from_hex("deadbeef12345678");
    assert_eq!(ibz_get(&a), 0x12345678i32, "get_lo32");

    assert_eq!(
        BigInt::from(7).cmp(&BigInt::from(7)),
        std::cmp::Ordering::Equal,
        "cmp_eq"
    );
    assert!(BigInt::from(3) < BigInt::from(7), "cmp_lt");
    assert!(BigInt::from(10) > BigInt::from(7), "cmp_gt");

    let a = ibz_set_from_str("2113309833171849999003363", 10).unwrap();
    assert_eq!(ibz_mod_ui(&a, 3), 0, "mod_ui_3");
    assert_eq!(ibz_mod_ui(&a, 2), 1, "mod_ui_2");

    assert!(
        ibz_divides(&BigInt::from(12), &BigInt::from(3)),
        "divides_12_3"
    );
    assert!(
        !ibz_divides(&BigInt::from(12), &BigInt::from(5)),
        "divides_12_5"
    );

    assert!(ibz_probab_prime(&BigInt::from(17), 25) > 0, "prime_17");
    assert_eq!(ibz_probab_prime(&BigInt::from(15), 25), 0, "prime_15");
}
