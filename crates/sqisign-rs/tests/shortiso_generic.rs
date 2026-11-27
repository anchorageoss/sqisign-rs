//! Item 3 of paper-shortiso: Corollaries A and B of "Shortest isogenies
//! between supersingular elliptic curves with complex multiplication" on
//! primes the crate has no field for, over the generic `F_{p^2}` of
//! `shortiso_common`. For every prime `p = D m^2 - 1` listed below the
//! kernel `G` of Corollary A is built from `iota` and `pi` evaluated on
//! points, shown cyclic of order exactly `m`, and the isogeny `E -> E/G`
//! is walked by Velu's formulas; the endpoint must be the curve with
//! complex multiplication by the order of discriminant `-D`.
//!
//! Each case prints one `ROW` line; `paper-shortiso/run.sh` collects them.

mod shortiso_common;

use num_bigint::BigUint;
use num_traits::One;
use shortiso_common::curve::{factor_smooth, torsion_basis, Curve};
use shortiso_common::field::{is_probable_prime, Field};
use shortiso_common::{
    chain_facts, cor_b, corollary_a_kernel, endpoint_check, pow2, steps_ascending, walk,
};

fn is_square(d: u64) -> bool {
    let r = (d as f64).sqrt().round() as u64;
    r * r == d
}

/// One prime `p = D m^2 - 1`: every assertion of the chain check.
fn run_case(group: &str, label: &str, d: u64, m: &BigUint) {
    let start = std::time::Instant::now();
    let p = BigUint::from(d) * m * m - 1u32;
    assert!(is_probable_prime(&p), "{label}: p is prime");
    assert!(d % 4 == 0 || d % 4 == 3, "{label}: D = 0 or 3 mod 4");
    assert!(!is_square(d), "{label}: D is not a square");
    let f = Field::new(&p);
    let e0 = Curve::e0(&f);
    let two_m = m * 2u32;
    assert!(((&p + 1u32) % &two_m) == BigUint::from(0u32), "2m | p + 1");
    let fac_m = factor_smooth(m, 1000).expect("m is smooth");
    let fac_2m = factor_smooth(&two_m, 1000).unwrap();

    // the kernel from two different bases of E[2m]
    let basis_a = torsion_basis(&e0, &f, &two_m, &fac_2m, 0);
    let ker_a = corollary_a_kernel(&f, d, m, &fac_m, &basis_a);
    let basis_b = torsion_basis(&e0, &f, &two_m, &fac_2m, 1);
    let ker_b = corollary_a_kernel(&f, d, m, &fac_m, &basis_b);

    // first walk: small primes first
    let steps = steps_ascending(&fac_m);
    let (js, _) = walk(&f, &ker_a.gen, &steps);
    // second walk: the other basis, the steps in the opposite order
    let mut steps_rev = steps.clone();
    steps_rev.reverse();
    let (js_b, _) = walk(&f, &ker_b.gen, &steps_rev);
    let j_end = js.last().unwrap();
    assert_eq!(
        j_end,
        js_b.last().unwrap(),
        "{label}: both walks end on the same j"
    );
    if fac_m.len() == 1 {
        assert_eq!(js, js_b, "{label}: the chain does not depend on the basis");
    }

    let facts = chain_facts(&f, &js);
    let facts_b = chain_facts(&f, &js_b);
    let (endpoint_ok, endpoint) = endpoint_check(&f, d, j_end);
    assert!(endpoint_ok, "{label}: endpoint is not {endpoint}");
    // D = 3, m a power of two: j = 54000 one step before j = 0
    let mut extra = String::new();
    if d == 3 && fac_m.len() == 1 && fac_m[0].0 == 2 {
        assert_eq!(
            js[js.len() - 2],
            f.small(54000),
            "{label}: j = 54000 one step before the end"
        );
        extra = "; j = 54000 at step n-1".into();
    }
    assert!(
        facts.all_in_fp && facts_b.all_in_fp,
        "{label}: a j-invariant of the chain is not in F_p"
    );
    assert!(
        facts.no_repeat && facts_b.no_repeat,
        "{label}: a j-invariant repeats in the chain"
    );
    let m_str = if fac_m.len() == 1 && fac_m[0].0 == 2 {
        format!("2^{}", fac_m[0].1)
    } else {
        fac_m
            .iter()
            .map(|&(l, v)| {
                if v == 1 {
                    format!("{l}")
                } else {
                    format!("{l}^{v}")
                }
            })
            .collect::<Vec<_>>()
            .join("*")
    };
    println!(
        "ROW|{group}|{label}|{} bits|D={d}|m={m_str}|steps={}|{endpoint}{extra}|all j in F_p, none repeats|PASS|{:.1}s",
        p.bits(),
        steps.len(),
        start.elapsed().as_secs_f64()
    );
}

/// `p = c 2^f - 1` with `(D, m)` from Corollary B.
fn run_cor_b(group: &str, c: u64, f: u32) {
    let (d, e) = cor_b(c, f);
    let m = pow2(e);
    let p = BigUint::from(c) * pow2(f) - 1u32;
    assert_eq!(
        BigUint::from(d) * &m * &m - 1u32,
        p,
        "Corollary B: p = D m^2 - 1"
    );
    run_case(group, &format!("{c}*2^{f}-1"), d, &m);
}

/// The six MIKE primes with the paper's `(D, m)`.
#[test]
fn mike_primes() {
    let table: [(u64, u32, u64, u32); 6] = [
        (117, 374, 468, 186),
        (77, 566, 308, 282),
        (41, 758, 164, 378),
        (633, 308, 2532, 153),
        (593, 474, 2372, 236),
        (317, 628, 1268, 313),
    ];
    std::thread::scope(|s| {
        for (c, f, d, e) in table {
            s.spawn(move || {
                assert_eq!(cor_b(c, f), (d, e), "{c}*2^{f}-1: the paper's (D, m)");
                run_cor_b("mike", c, f);
            });
        }
    });
}

/// Every prime `3 * 2^f - 1`, `f` even, `f < 500`: `D = 3`, `m = 2^(f/2)`,
/// endpoint `j = 0`, through `j = 54000`.
#[test]
fn three_times_power_of_two() {
    let fs: Vec<u32> = (2..500)
        .step_by(2)
        .filter(|&f| is_probable_prime(&(BigUint::from(3u32) * pow2(f) - 1u32)))
        .collect();
    println!("INFO|3*2^f-1 prime for even f < 500: f in {fs:?}");
    assert!(fs.contains(&324), "the level I prime is in the list");
    std::thread::scope(|s| {
        for &f in &fs {
            s.spawn(move || {
                assert_eq!(cor_b(3, f), (3, f / 2));
                run_cor_b("3*2^f-1", 3, f);
            });
        }
    });
}

/// The Mersenne primes `2^f - 1`, `5 <= f <= 607`: `D = 8`,
/// `m = 2^((f-3)/2)`, endpoint `j = 8000`.
#[test]
fn mersenne_primes() {
    let fs: Vec<u32> = (5..=607)
        .filter(|&f| is_probable_prime(&(pow2(f) - 1u32)))
        .collect();
    println!("INFO|2^f-1 prime for 5 <= f <= 607: f in {fs:?}");
    assert_eq!(fs, [5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607]);
    std::thread::scope(|s| {
        for &f in &fs {
            s.spawn(move || {
                assert_eq!(cor_b(1, f), (8, (f - 3) / 2));
                run_cor_b("mersenne", 1, f);
            });
        }
    });
}

/// The three SQIsign primes, generically (the crate's own machinery is in
/// `shortiso_sqisign.rs`, which also compares the two).
#[test]
fn sqisign_primes_generic() {
    std::thread::scope(|s| {
        for (c, f, d, e) in [
            (3u64, 324u32, 3u64, 162u32),
            (27, 500, 27, 250),
            (17, 664, 68, 331),
        ] {
            s.spawn(move || {
                assert_eq!(cor_b(c, f), (d, e));
                run_cor_b("sqisign", c, f);
            });
        }
    });
}

/// The integers `2^a 3^b 5^c 7^d 11^e` below `2^130`, sorted.
fn smooth_numbers() -> Vec<BigUint> {
    let limit = pow2(130);
    let mut out = vec![BigUint::one()];
    for l in [2u32, 3, 5, 7, 11] {
        let mut next = Vec::new();
        for x in &out {
            let mut y = x.clone();
            while y < limit {
                next.push(y.clone());
                y *= l;
            }
        }
        out = next;
    }
    out.sort();
    out
}

/// Primes `D m^2 - 1` with `m` composite: for each `D` of class number one
/// in the list, the least 11-smooth `m` with an odd prime factor in each of
/// four size bands for which `D m^2 - 1` is prime.
#[test]
fn class_number_one_composite_m() {
    let smooth = smooth_numbers();
    let bands = [
        BigUint::from(2u32),
        BigUint::from(1000u32),
        BigUint::from(10u32).pow(12),
        BigUint::from(10u32).pow(30),
    ];
    let mut cases = Vec::new();
    for d in [3u64, 7, 8, 11, 12, 19, 27, 28, 43, 67, 163] {
        for lo in &bands {
            let m = smooth
                .iter()
                .find(|m| {
                    *m >= lo
                        && (d % 2 == 0 || (*m % 2u32) == BigUint::from(0u32))
                        && m.trailing_zeros().unwrap_or(0) < m.bits() - 1
                        && is_probable_prime(&(BigUint::from(d) * *m * *m - 1u32))
                })
                .expect("a prime in the band")
                .clone();
            cases.push((d, m));
        }
    }
    std::thread::scope(|s| {
        for (d, m) in &cases {
            s.spawn(move || run_case("D*m^2-1", &format!("{d}*{m}^2-1"), *d, m));
        }
    });
}
