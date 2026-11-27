//! Item 6 of paper-shortiso (the part that is arithmetic): the table of
//! `(D, m, h(-D))` for the three SQIsign primes and the six MIKE primes,
//! recomputed from the primes by Corollary B, and the class polynomials the
//! other checks rely on, validated against values known in closed form.

mod shortiso_common;

use num_bigint::{BigInt, BigUint};
use shortiso_common::classpoly::{class_number, class_polynomial, reduced_forms};
use shortiso_common::field::is_probable_prime;
use shortiso_common::{class_polynomial_checked, cor_b, j_class_number_one, pow2};

fn ints(v: &[i128]) -> Vec<BigInt> {
    v.iter().map(|&x| BigInt::from(x)).collect()
}

/// The class polynomial code against the literature: the thirteen orders
/// of class number one, and `H_-15`, `H_-20`, `H_-23`.
#[test]
fn class_polynomials_known_values() {
    for d in [3u64, 4, 7, 8, 11, 12, 16, 19, 27, 28, 43, 67, 163] {
        let h = class_polynomial(d);
        let j = j_class_number_one(d).unwrap();
        assert_eq!(h, vec![-j, BigInt::from(1)], "H_-{d}");
    }
    assert_eq!(
        class_polynomial(15),
        ints(&[-121287375, 191025, 1]),
        "H_-15"
    );
    assert_eq!(
        class_polynomial(20),
        ints(&[-681472000, -1264000, 1]),
        "H_-20"
    );
    assert_eq!(
        class_polynomial(23),
        ints(&[12771880859375, -5151296875, 3491750, 1]),
        "H_-23"
    );
    // class numbers of a few discriminants, from tables
    for (d, h) in [
        (23u64, 3usize),
        (47, 5),
        (71, 7),
        (84, 4),
        (163, 1),
        (164, 8),
        (68, 4),
    ] {
        assert_eq!(class_number(d), h, "h(-{d})");
    }
    println!("INFO|class polynomial code agrees with the 13 class-number-one j-invariants and with H_-15, H_-20, H_-23");
}

/// `(c, f, D, e)` with `p = c 2^f - 1`, `m = 2^e`: the paper's table.
const PAPER_TABLE: [(&str, u64, u32, u64, u32); 9] = [
    ("SQIsign I", 3, 324, 3, 162),
    ("SQIsign III", 27, 500, 27, 250),
    ("SQIsign V", 17, 664, 68, 331),
    ("MIKE conservative I", 117, 374, 468, 186),
    ("MIKE conservative III", 77, 566, 308, 282),
    ("MIKE conservative V", 41, 758, 164, 378),
    ("MIKE fast I", 633, 308, 2532, 153),
    ("MIKE fast III", 593, 474, 2372, 236),
    ("MIKE fast V", 317, 628, 1268, 313),
];

/// The `(D, m)` of the paper follow from the primes by Corollary B, the
/// identity `p = D m^2 - 1` holds, and `h(-D)` is the number of reduced
/// primitive forms; the class polynomials have that degree and are the
/// checked-in records.
#[test]
fn d_m_h_table() {
    for (name, c, f, d, e) in PAPER_TABLE {
        let p = BigUint::from(c) * pow2(f) - 1u32;
        assert!(is_probable_prime(&p), "{name}: prime");
        assert_eq!(cor_b(c, f), (d, e), "{name}: (D, m) by Corollary B");
        assert_eq!(
            BigUint::from(d) * pow2(2 * e) - 1u32,
            p,
            "{name}: p = D m^2 - 1"
        );
        let forms = reduced_forms(d);
        let h = forms.len();
        // the conductor of the order: the largest g with -D / g^2 a discriminant
        let mut cond = 1u64;
        for g in 2..=d {
            if d % (g * g) == 0 && matches!((d / (g * g)) % 4, 0 | 3) {
                cond = g;
            }
        }
        let degree = if h == 1 {
            assert!(j_class_number_one(d).is_some());
            1
        } else {
            class_polynomial_checked(d).len() - 1
        };
        assert_eq!(degree, h, "{name}: deg H_-{d} = h(-{d})");
        println!(
            "ROW|table|{name}|{c}*2^{f}-1|{} bits|D={d}|m=2^{e}|h(-{d})={h}|conductor {cond}|PASS",
            p.bits()
        );
    }
}
