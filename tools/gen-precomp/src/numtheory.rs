//! Integer routines used by the level-constant formulas: primality, next
//! prime, integer square root, and the two Appendix B searches.

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

const SMALL_PRIMES: [u32; 25] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

/// Miller-Rabin with the first 25 primes as bases, after trial division.
/// Deterministic; more than enough for generating constants that are then
/// checked against the reference.
pub fn is_probable_prime(n: &BigUint) -> bool {
    if n < &BigUint::from(2u32) {
        return false;
    }
    for &q in &SMALL_PRIMES {
        let q = BigUint::from(q);
        if n == &q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let one = BigUint::one();
    let n_minus_1 = n - &one;
    let s = n_minus_1.trailing_zeros().expect("n > 1");
    let d = &n_minus_1 >> s;
    'bases: for &a in &SMALL_PRIMES {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x == one || x == n_minus_1 {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % n;
            if x == n_minus_1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// The smallest prime strictly greater than `n`.
pub fn next_prime(n: &BigUint) -> BigUint {
    let mut c = n + BigUint::one();
    while !is_probable_prime(&c) {
        c += BigUint::one();
    }
    c
}

pub fn isqrt(n: &BigUint) -> BigUint {
    n.sqrt()
}

/// Appendix B, Algorithm B.2 (EIBox): the least `b` such that
/// `(1 - 1/(2 log2 p))^((b+1)^4) < 2^-lambda`.
pub fn ei_box(log2p: u64, lambda: u64) -> u64 {
    let b = (1.0 - 1.0 / (2.0 * log2p as f64)).log2();
    let b = -(lambda as f64) / b;
    let b = b.powf(0.25) - 1.0;
    b.ceil() as u64
}

/// Appendix B, Algorithm B.3 (RICofactor), as `precompute_quaternion_data.sage`
/// implements it: the smallest `m >= log2 p - e_rsp + 2` such that
/// `round(floor(2^(e_rsp + m) / p) * log2(1 - 1/(2 (e_rsp + m)))) <= -lambda`,
/// then the next prime after `2^m`.
pub fn ri_cofactor(p: &BigUint, log2p: u64, e_rsp: u64, lambda: u64) -> BigUint {
    let mut m = log2p - e_rsp + 2;
    loop {
        let n = ((BigUint::one() << (e_rsp + m)) / p)
            .to_f64()
            .expect("small quotient");
        let log_m = (e_rsp + m) as f64;
        let val = (n * (1.0 - 1.0 / (2.0 * log_m)).log2()).round();
        if val <= -(lambda as f64) {
            break;
        }
        m += 1;
    }
    next_prime(&(BigUint::one() << m))
}
