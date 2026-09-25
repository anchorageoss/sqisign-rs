//! Integer helpers of the quaternion module: Cornacchia for primes and
//! random prime generation for tests.

use super::lll::lehmer::dim2_sumofsquares;
use crate::mp::{DefaultDomain, Ibz, Rng};

/// Solve `x^2 + y^2 = p` for a prime `p` (or `p = 1, 2`), with
/// `0 <= x <= y`. `None` if there is no solution (`p = 3 mod 4`) or the
/// square root of `-1` could not be verified.
pub fn cornacchia_prime<const N: usize>(p: &Ibz<N>) -> Option<(Ibz<N>, Ibz<N>)> {
    if *p == Ibz::two() {
        return Some((Ibz::set(1, 2), Ibz::set(1, 2)));
    }
    if *p == Ibz::one() {
        return Some((Ibz::zero(), Ibz::set(1, 2)));
    }
    if (p.get() & 3) != 1 {
        return None;
    }
    let root = Ibz::sqrt_m1_mod_verified(p)?;
    let (mut x, mut y) = dim2_sumofsquares(p, &root);
    let sum = x.mul(&x).add(&y.mul(&y));
    if x > y {
        core::mem::swap(&mut x, &mut y);
    }
    if sum == *p {
        Some((x, y))
    } else {
        None
    }
}

/// A random probable prime of at most `bitsize` bits, `3 mod 4` when
/// requested; `rounds` Miller-Rabin rounds. Mostly for tests. `None` on
/// randomness failure.
pub fn generate_random_prime<const N: usize>(
    is3mod4: bool,
    bitsize: u32,
    rounds: u32,
    rng: &mut impl Rng,
) -> Option<Ibz<N>> {
    debug_assert!(bitsize != 0);
    let two_pow = Ibz::<N>::one().mul_2exp(bitsize - 1 - is3mod4 as u32);
    let two_powp = Ibz::<N>::one().mul_2exp(bitsize - is3mod4 as u32);
    loop {
        let mut p = Ibz::rand_interval(&two_pow, &two_powp, &mut DefaultDomain(&mut *rng))?;
        p = p.mul_2exp(1);
        if is3mod4 {
            p = p.mul_2exp(1);
            p = Ibz::two().add(&p);
        }
        p = Ibz::one().add(&p);
        if p.probab_prime(rounds, rng) {
            return Some(p);
        }
    }
}
