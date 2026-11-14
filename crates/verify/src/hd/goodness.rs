//! The goodness test on the response degree `q`, repeated by the verifier
//! before Cornacchia's algorithm and applied by the signer to every
//! candidate: `q < 2^e` and `2^e − q` a probable prime congruent to 1
//! modulo 4 (so `q` is odd and `q ≡ 3 mod 4`). The test is deterministic:
//! forty Miller-Rabin rounds with the first forty primes as bases. Variable
//! time; `q` is public.

use super::uint::{Monty, U6};

/// The first forty primes, the Miller-Rabin bases.
pub const MR_BASES: [u64; 40] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
];

/// Whether `n` is a probable prime: trial division by [`MR_BASES`], then
/// forty Miller-Rabin rounds with those bases.
pub fn is_probable_prime(n: &U6) -> bool {
    if *n == U6::from_u64(2) {
        return true;
    }
    if !n.is_odd() || *n == U6::ONE {
        return false;
    }
    // trial division by the bases: exact for n below 173^2, and it removes
    // most composites before any exponentiation
    for &b in MR_BASES.iter() {
        if n.rem_u64(b) == 0 {
            return *n == U6::from_u64(b);
        }
    }
    let nm1 = n.wrapping_sub(&U6::ONE);
    let s = nm1.trailing_zeros();
    let d = nm1.shr(s);
    let mt = Monty::new(n);
    let one = mt.to(&U6::ONE);
    let minus_one = n.wrapping_sub(&one);
    'bases: for &b in MR_BASES.iter() {
        let mut x = mt.pow(&mt.to(&U6::from_u64(b)), &d);
        if x == one || x == minus_one {
            continue;
        }
        for _ in 1..s {
            x = mt.mul(&x, &x);
            if x == minus_one {
                continue 'bases;
            }
            if x == one {
                return false;
            }
        }
        return false;
    }
    true
}

/// Whether `q` is a good degree for the exponent `e`: `q < 2^e` and
/// `2^e − q` a probable prime congruent to 1 modulo 4.
pub fn is_good_degree(e: u32, q: &U6) -> bool {
    let two_e = U6::pow2(e);
    if !q.lt(&two_e) {
        return false;
    }
    let n = two_e.wrapping_sub(q);
    if n.0[0] & 3 != 1 {
        return false;
    }
    is_probable_prime(&n)
}

#[cfg(test)]
mod tests {
    use super::{is_good_degree, is_probable_prime, U6};

    #[test]
    fn small_numbers_against_trial_division() {
        let trial = |n: u64| n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
        for n in 0..20_000u64 {
            assert_eq!(is_probable_prime(&U6::from_u64(n)), trial(n), "{n}");
        }
        // squares and products of the bases' neighbours
        for n in [173u64 * 173, 173 * 179, 179 * 179, 179 * 181, 65537 * 65539] {
            assert!(!is_probable_prime(&U6::from_u64(n)), "{n}");
        }
    }

    #[test]
    fn known_primes_and_composites_at_the_three_sizes() {
        // 2^e - c prime for these c (checked with an independent tool)
        for (e, c) in [(174u32, 3u64), (264, 275), (346, 363)] {
            let p = U6::pow2(e).wrapping_sub(&U6::from_u64(c));
            assert!(is_probable_prime(&p), "2^{e} - {c}");
            // a strong pseudoprime candidate: the square of a prime is composite
            let half = U6::pow2(e / 2).wrapping_sub(&U6::from_u64(1));
            assert!(!is_probable_prime(&half.wrapping_mul(&half)));
        }
        // Carmichael numbers
        for n in [561u64, 41041, 825265, 321197185, 5394826801, 232250619601] {
            assert!(!is_probable_prime(&U6::from_u64(n)), "{n}");
        }
    }

    #[test]
    fn goodness_rejects_what_the_paper_lists() {
        let e = 174;
        // q >= 2^e
        assert!(!is_good_degree(e, &U6::pow2(e)));
        // q even: 2^e - q even
        assert!(!is_good_degree(e, &U6::from_u64(2)));
        // q = 1 mod 4: 2^e - q = 3 mod 4
        assert!(!is_good_degree(e, &U6::from_u64(1)));
        // q = 3: 2^174 - 3 is prime and 1 mod 4
        assert!(is_good_degree(e, &U6::from_u64(3)));
        // q = 7: 2^174 - 7 = 1 mod 4 and composite
        assert!(!is_good_degree(e, &U6::from_u64(7)));
    }
}
