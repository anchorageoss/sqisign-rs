use super::intbig::{
    ibz_div, ibz_gcd, ibz_pow, ibz_probab_prime, ibz_rand_interval, ibz_sqrt, ibz_sqrt_mod_p, Ibz,
};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use rand::Rng;

/// Generate a random prime of the specified bit size by rejection sampling.
///
/// If `is3mod4` is true, the prime satisfies `p ≡ 3 (mod 4)`.
/// Uses `iters` Miller-Rabin iterations for primality testing.
///
/// Candidates are constructed as `2r + 1` (odd) when `is3mod4` is false,
/// or as `4r + 3` (≡ 3 mod 4) when `is3mod4` is true, where `r` is drawn
/// uniformly at random from an interval that ensures the result has the
/// requested bit size.
pub fn ibz_generate_random_prime(
    rng: &mut impl Rng,
    is3mod4: bool,
    bitsize: u32,
    iters: u32,
) -> Ibz {
    assert!(bitsize != 0);

    let shift = if is3mod4 { 1u32 } else { 0u32 };
    let two_pow = ibz_pow(&BigInt::from(2), (bitsize - 1) - shift);
    let two_powp = ibz_pow(&BigInt::from(2), bitsize - shift);

    loop {
        let mut p = ibz_rand_interval(rng, &two_pow, &two_powp);
        // p = 2 * p
        p = &p + &p;
        if is3mod4 {
            // p = 4 * original_p (double again)
            p = &p + &p;
            // p = 4 * original_p + 2
            p = &p + &BigInt::from(2);
        }
        // p = 2*r+1 (or 4*r+3 if is3mod4)
        p = &p + &BigInt::from(1);

        if ibz_probab_prime(&p, iters) > 0 {
            return p;
        }
    }
}

/// Cornacchia's algorithm: solve `x² + n·y² = p` for a prime `p`.
///
/// Given a prime `p` and a positive integer `n`, attempts to find integers
/// `(x, y)` such that `x² + n·y² = p`. Uses the Euclidean algorithm on
/// a square root of `-n mod p` to reduce the problem to checking whether
/// the remainder yields a valid decomposition.
///
/// Returns `Some((x, y))` if a solution exists, `None` otherwise.
pub fn ibz_cornacchia_prime(n: &Ibz, p: &Ibz) -> Option<(Ibz, Ibz)> {
    let two = BigInt::from(2);

    // Special case: p = 2
    if p == &two {
        if n.is_one() {
            return Some((BigInt::from(1), BigInt::from(1)));
        }
        return None;
    }

    // Special case: p = n
    if p == n {
        return Some((BigInt::zero(), BigInt::from(1)));
    }

    // Test coprimality
    let g = ibz_gcd(p, n);
    if !g.is_one() {
        return None;
    }

    // Compute sqrt(-n mod p)
    let neg_n = -n;
    let r2 = ibz_sqrt_mod_p(&neg_n, p)?;

    // Euclidean algorithm loop.
    // Initialize: r2 = sqrt(-n) mod p, r1 = p
    // We want to find the first remainder r0 such that r0^2 < p.
    let mut r2 = r2;
    let mut r1 = p.clone();
    let mut prod = p.clone(); // Initialize prod = p to enter loop

    let mut r0;
    while prod >= *p {
        let (_q, remainder) = ibz_div(&r2, &r1);
        r0 = remainder;
        prod = &r0 * &r0;
        r2 = r1;
        r1 = r0;
    }

    // Loop exits when r1^2 < p.
    // Then a = p - r1^2; if n | a and a/n is a perfect square, return (r1, sqrt(a/n)).

    // a = p - r1^2
    let a = p - &prod;
    let (q, rem) = ibz_div(&a, n);
    if !rem.is_zero() {
        return None;
    }
    let y = ibz_sqrt(&q)?;

    let x = r1;

    // Verify: x^2 + n * y^2 == p
    let check = &(&x * &x) + &(n * &y * &y);
    if &check != p {
        return None;
    }

    Some((x, y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::intbig::{ibz_bitsize, ibz_mod_ui, ibz_probab_prime};
    use num_bigint::BigInt;
    use rand::SeedableRng;

    #[test]
    fn test_generate_random_prime() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        // 20-bit prime, p ≡ 3 (mod 4)
        let p = ibz_generate_random_prime(&mut rng, true, 20, 30);
        assert!(ibz_probab_prime(&p, 20) > 0, "should be prime");
        assert!(ibz_bitsize(&p) >= 20, "should have at least 20 bits");
        assert_eq!(ibz_mod_ui(&p, 4), 3, "should be 3 mod 4");

        // 30-bit prime, no constraint
        let p = ibz_generate_random_prime(&mut rng, false, 30, 30);
        assert!(ibz_probab_prime(&p, 20) > 0);
        assert!(ibz_bitsize(&p) >= 30);

        // 30-bit prime, p ≡ 3 (mod 4)
        let p = ibz_generate_random_prime(&mut rng, true, 30, 30);
        assert!(ibz_probab_prime(&p, 20) > 0);
        assert!(ibz_bitsize(&p) >= 30);
        assert_eq!(ibz_mod_ui(&p, 4), 3);
    }

    #[test]
    fn test_cornacchia_prime() {
        // n=1, p=5: 1^2 + 1*2^2 = 5
        let n = BigInt::from(1);
        let p = BigInt::from(5);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("should find solution");
        assert_eq!(&x * &x + &n * &y * &y, p);

        // n=1, p=2: 1^2 + 1*1^2 = 2
        let p = BigInt::from(2);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("should find solution");
        assert_eq!(&x * &x + &n * &y * &y, p);

        // n=1, p=41
        let p = BigInt::from(41);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("should find solution");
        assert_eq!(&x * &x + &n * &y * &y, p);

        // n=2, p=3: 1^2 + 2*1^2 = 3
        let n = BigInt::from(2);
        let p = BigInt::from(3);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("should find solution");
        assert_eq!(&x * &x + &n * &y * &y, p);

        // n=3, p=7: 2^2 + 3*1^2 = 7
        let n = BigInt::from(3);
        let p = BigInt::from(7);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("should find solution");
        assert_eq!(&x * &x + &n * &y * &y, p);

        // No solutions:
        // n=1, p=7: 7 ≡ 3 (mod 4), no representation as x²+y²
        let n = BigInt::from(1);
        let p = BigInt::from(7);
        assert!(ibz_cornacchia_prime(&n, &p).is_none());

        // n=1, p=3
        let p = BigInt::from(3);
        assert!(ibz_cornacchia_prime(&n, &p).is_none());

        // n=3, p=5
        let n = BigInt::from(3);
        let p = BigInt::from(5);
        assert!(ibz_cornacchia_prime(&n, &p).is_none());

        // n=3, p=3: special case p=n, should return (0, 1)
        let n = BigInt::from(3);
        let p = BigInt::from(3);
        let (x, y) = ibz_cornacchia_prime(&n, &p).expect("p=n should be solvable");
        assert_eq!(&x * &x + &n * &y * &y, p);
    }
}
