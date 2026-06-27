//!
//! Provides `ibz_probab_prime` and `ibz_sqrt_mod_p` replacements that
//! convert to fixed-precision `crypto-bigint::Uint` at the function
//! boundary, perform all modular arithmetic on the stack with zero
//! heap allocation, then convert back.
//!
//! The primality and sqrt functions are generic over `const N: usize`
//! (the number of 64-bit limbs). Entry points dispatch to the smallest
//! width that fits the input:
//! - `Uint<8>` (512 bits) for values up to 512 bits
//! - `Uint<16>` (1024 bits) for larger values
//! - `num-bigint` fallback for anything exceeding 1024 bits (safety net)

use super::intbig::{ibz_mod, ibz_mod_ui, Ibz};
use alloc::vec;
use crypto_bigint::modular::{MontyForm, MontyParams};
use crypto_bigint::{Limb, Odd, Uint};
use num_bigint::{BigInt, Sign};
use num_traits::Zero;
// Conversion helpers (generic over limb count)

/// Returns true if `n` fits in `N * 64` bits.
fn fits_in<const N: usize>(n: &BigInt) -> bool {
    let (_, bytes) = n.to_bytes_le();
    bytes.len() <= N * 8
}

/// Convert a non-negative `BigInt` to `Uint<N>`. Panics if negative or too large.
fn bigint_to_uint<const N: usize>(n: &BigInt) -> Uint<N> {
    debug_assert!(n.sign() != Sign::Minus, "bigint_to_uint: negative input");
    let (_, bytes) = n.to_bytes_le();
    let mut limb_words = [0u64; N];
    for (i, chunk) in bytes.chunks(8).enumerate() {
        // Callers gate on fits_in::<N>(), so this should never be reached.
        debug_assert!(
            i < N,
            "bigint_to_uint: value too large for Uint<{N}> ({} bytes)",
            bytes.len()
        );
        if i >= N {
            break;
        }
        let mut word_bytes = [0u8; 8];
        word_bytes[..chunk.len()].copy_from_slice(chunk);
        limb_words[i] = u64::from_le_bytes(word_bytes);
    }
    let limbs: [Limb; N] = core::array::from_fn(|i| Limb(limb_words[i]));
    Uint::new(limbs)
}

/// Convert `Uint<N>` back to a non-negative `BigInt`.
fn uint_to_bigint<const N: usize>(n: &Uint<N>) -> BigInt {
    let mut bytes = vec![0u8; N * 8];
    for (i, limb) in n.as_limbs().iter().enumerate() {
        bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb.0.to_le_bytes());
    }
    BigInt::from_bytes_le(Sign::Plus, &bytes)
}

/// Modular exponentiation: `base^exp mod modulus` using crypto-bigint Montgomery form.
fn modpow_uint<const N: usize>(base: &Uint<N>, exp: &Uint<N>, params: &MontyParams<N>) -> Uint<N> {
    let base_monty = MontyForm::new(base, *params);
    let result = base_monty.pow(exp);
    result.retrieve()
}

/// Product of small odd primes: 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29
const SMALL_PRIME_PRODUCT: u64 = 3234846615;

/// Bit (p+1)/2 is set for each odd prime <= 61.
const SMALL_PRIME_MASK: u32 = 0xc96996dc;

/// Probabilistic primality test using BPSW + extra Miller-Rabin rounds.
///
/// 1. Even/small/trial-division fast path
/// 2. BPSW: Miller-Rabin base 2 + Strong Lucas test
/// 3. `reps - 24` extra MR rounds with Euler polynomial bases `j^2+j+41`
///
/// Deterministic, does not use any RNG.
pub fn ibz_probab_prime(n: &Ibz, reps: u32) -> i32 {
    let n_abs = if n.sign() == Sign::Minus {
        -n.clone()
    } else {
        n.clone()
    };
    let (_, n_bytes) = n_abs.to_bytes_le();

    // Even check
    if n_bytes[0] & 1 == 0 {
        return if n_bytes.len() == 1 && n_bytes[0] == 2 {
            2
        } else {
            0
        };
    }

    // Dispatch to the smallest fixed-precision width that fits.
    if fits_in::<8>(&n_abs) {
        ibz_probab_prime_uint::<8>(&n_abs, &n_bytes, reps)
    } else if fits_in::<16>(&n_abs) {
        ibz_probab_prime_uint::<16>(&n_abs, &n_bytes, reps)
    } else {
        ibz_probab_prime_bigint(&n_abs, reps)
    }
}

/// Fixed-precision primality test using `Uint<N>` Montgomery arithmetic.
fn ibz_probab_prime_uint<const N: usize>(n_abs: &BigInt, n_bytes: &[u8], reps: u32) -> i32 {
    let n_u: Uint<N> = bigint_to_uint(n_abs);
    let bits = (N * 64) as u32;

    // Small primes via bitmask (n < 64)
    if n_u < Uint::<N>::from_u8(64) {
        let low = n_bytes[0] as u32;
        return ((SMALL_PRIME_MASK >> (low >> 1)) & 2) as i32;
    }

    // Trial division: gcd(n, 3*5*7*...*29)
    // SMALL_PRIME_PRODUCT is a non-zero constant.
    let rem = n_u.rem_limb(
        crypto_bigint::NonZero::new(Limb(SMALL_PRIME_PRODUCT))
            .expect("invariant: SMALL_PRIME_PRODUCT is a non-zero constant"),
    );
    if gcd_pair(rem.0, SMALL_PRIME_PRODUCT) != 1 {
        return 0;
    }

    // All prime factors >= 31, so if n < 31*31 it must be prime
    if n_u < Uint::<N>::from_u16(31 * 31) {
        return 2;
    }

    // n is guaranteed odd: even values were returned at the top of ibz_probab_prime.
    let n_odd = Odd::new(n_u).expect("invariant: n is odd (even values returned earlier)");
    let params = MontyParams::<N>::new_vartime(n_odd);
    let nm1 = n_u.wrapping_sub(&Uint::<N>::ONE);
    let k = nm1.trailing_zeros();
    let q = nm1.shr(k);

    // BPSW: Miller-Rabin base 2
    let mut is_prime = miller_rabin_uint(&Uint::<N>::from_u8(2), &q, k, &nm1, &params);

    // BPSW: Strong Lucas test
    if is_prime {
        is_prime = strong_lucas_uint(&n_u, &params, bits);
    }

    // Extra MR rounds with Euler polynomial bases j^2+j+41
    let extra_reps = (reps as i32) - 24;
    let mut j: u32 = 0;
    while is_prime && (j as i32) < extra_reps {
        let base_val = (j as u64) * (j as u64) + (j as u64) + 41;
        let base = Uint::<N>::from_u64(base_val);
        if base >= nm1 {
            break;
        }
        is_prime = miller_rabin_uint(&base, &q, k, &nm1, &params);
        j += 1;
    }

    if is_prime {
        1
    } else {
        0
    }
}

/// Fallback primality test for values exceeding all fixed-precision tiers.
/// Uses `num-bigint::BigUint::modpow`, works for any bit width.
fn ibz_probab_prime_bigint(n_abs: &BigInt, reps: u32) -> i32 {
    use num_bigint::BigUint;
    use num_traits::One;

    let (_, n_bytes) = n_abs.to_bytes_le();
    let n_u = BigUint::from_bytes_le(&n_bytes);

    // Small primes via bitmask (n < 64), shouldn't hit for oversized values, but safe
    if n_u < BigUint::from(64u32) {
        let low = n_bytes[0] as u32;
        return ((SMALL_PRIME_MASK >> (low >> 1)) & 2) as i32;
    }

    // Trial division: gcd(n, 3*5*7*...*29)
    let rem_u = &n_u % BigUint::from(SMALL_PRIME_PRODUCT);
    let rem64: u64 = rem_u.try_into().unwrap_or(0);
    if gcd_pair(rem64, SMALL_PRIME_PRODUCT) != 1 {
        return 0;
    }

    // n - 1 = q * 2^k with q odd
    let nm1 = &n_u - BigUint::one();
    let k = nm1.trailing_zeros().unwrap_or(0);
    let q_u = &nm1 >> k;

    // BPSW: Miller-Rabin base 2
    let mut is_prime = miller_rabin_biguint(&BigUint::from(2u32), &q_u, k as u32, &nm1, &n_u);

    // BPSW: Strong Lucas test
    if is_prime {
        is_prime = strong_lucas_biguint(&n_u);
    }

    // Extra MR rounds with Euler polynomial bases j^2+j+41
    let extra_reps = (reps as i32) - 24;
    let mut j: u32 = 0;
    while is_prime && (j as i32) < extra_reps {
        let base_val = (j as u64) * (j as u64) + (j as u64) + 41;
        let base = BigUint::from(base_val);
        if base >= nm1 {
            break;
        }
        is_prime = miller_rabin_biguint(&base, &q_u, k as u32, &nm1, &n_u);
        j += 1;
    }

    if is_prime {
        1
    } else {
        0
    }
}

fn gcd_pair(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Single Miller-Rabin witness test using crypto-bigint Montgomery arithmetic.
fn miller_rabin_uint<const N: usize>(
    a: &Uint<N>,
    q: &Uint<N>,
    k: u32,
    nm1: &Uint<N>,
    params: &MontyParams<N>,
) -> bool {
    let one = MontyForm::one(*params);
    let neg_one = MontyForm::new(nm1, *params);
    let mut y = MontyForm::new(a, *params).pow(q);

    if y == one || y == neg_one {
        return true;
    }
    for _ in 1..k {
        y = y * y;
        if y == neg_one {
            return true;
        }
    }
    false
}

/// Strong Lucas pseudoprime test with Selfridge's Method A parameter selection.
fn strong_lucas_uint<const N: usize>(n: &Uint<N>, params: &MontyParams<N>, bits: u32) -> bool {
    if is_perfect_square_uint(n, bits) {
        return false;
    }

    // Selfridge's Method A: find D in {5, -7, 9, -11, ...} with (D/n) = -1.
    let sqrt_n = isqrt_uint(n, bits);
    let max_d: u64 = if sqrt_n < Uint::<N>::from_u64(u64::MAX) {
        sqrt_n.as_limbs()[0].0.saturating_sub(1)
    } else {
        u64::MAX
    };

    let mut d_abs: u64 = 3;
    loop {
        if d_abs >= max_d {
            return true;
        }
        d_abs += 2;
        // d_abs starts at 5 and increases by 2 each iteration, always non-zero.
        let tl = n
            .rem_limb(
                crypto_bigint::NonZero::new(Limb(d_abs))
                    .expect("invariant: d_abs starts at 5 and increases, always non-zero"),
            )
            .0;
        if tl == 0 {
            return false;
        }
        if jacobi_coprime(tl, d_abs) != 1 {
            break;
        }
    }

    // Selfridge parameter: Q = (1 - D) / 4, with D = ±d_abs.
    let q_val: i64 = if d_abs & 2 != 0 {
        (d_abs >> 2) as i64 + 1
    } else {
        -((d_abs >> 2) as i64)
    };

    // n + 1 = d_odd * 2^b0
    let np1 = n.wrapping_add(&Uint::<N>::ONE);
    let b0 = np1.trailing_zeros();

    let (u_is_zero, mut v_monty, mut qk_monty) = lucas_mod_uint(n, params, q_val, b0, bits);

    if u_is_zero {
        return true;
    }

    // Check doubling steps: while V != 0 and remaining b0 > 1
    let zero = MontyForm::new(&Uint::<N>::ZERO, *params);
    let two = MontyForm::new(&Uint::<N>::from_u8(2), *params);
    let mut remaining = b0;

    while v_monty != zero && remaining > 1 {
        remaining -= 1;
        let v_sq = v_monty * v_monty;
        v_monty = v_sq - two * qk_monty;
        qk_monty = qk_monty * qk_monty;
    }

    remaining > 0
}

/// Jacobi symbol `(a/b)` for coprime `a`, `b` with `b` odd.
fn jacobi_coprime(mut a: u64, mut b: u64) -> i32 {
    debug_assert!(b & 1 == 1);
    debug_assert!(a != 0);

    let mut bit: u64 = 0;
    b >>= 1;

    let mut c = a.trailing_zeros();
    a >>= 1;

    loop {
        a >>= c;
        bit ^= (c as u64) & (b ^ (b >> 1));
        if a < b {
            if a == 0 {
                return if bit & 1 != 0 { -1 } else { 1 };
            }
            bit ^= a & b;
            let new_a = b - a;
            b = a;
            a = new_a;
        } else {
            a -= b;
            debug_assert!(a != 0);
        }

        c = a.trailing_zeros() + 1;
    }
}

/// Integer square root of a `Uint<N>` value (floor).
fn isqrt_uint<const N: usize>(n: &Uint<N>, bits: u32) -> Uint<N> {
    if *n == Uint::<N>::ZERO {
        return Uint::<N>::ZERO;
    }
    let nbits = bits - n.leading_zeros();
    let mut x = Uint::<N>::ONE.shl(nbits.div_ceil(2));
    loop {
        // x starts as 2^(nbits/2) > 0 and remains positive throughout Newton iteration.
        let (q, _) = n.div_rem(
            &crypto_bigint::NonZero::new(x)
                .expect("invariant: Newton iteration x is always positive"),
        );
        let sum = x.wrapping_add(&q);
        let x_new = sum.shr(1);
        if x_new >= x {
            return x;
        }
        x = x_new;
    }
}

fn is_perfect_square_uint<const N: usize>(n: &Uint<N>, bits: u32) -> bool {
    let s = isqrt_uint(n, bits);
    let ss = s.wrapping_mul(&s);
    ss == *n
}

/// Computes the Lucas sequence U_k, V_k, Q^k (mod n) with P=1, Q=q_val.
fn lucas_mod_uint<const N: usize>(
    n: &Uint<N>,
    params: &MontyParams<N>,
    q_val: i64,
    b0: u32,
    bits: u32,
) -> (bool, MontyForm<N>, MontyForm<N>) {
    let one = MontyForm::one(*params);
    let zero_mf = MontyForm::new(&Uint::<N>::ZERO, *params);

    let q_abs = q_val.unsigned_abs();
    let q_u = Uint::<N>::from_u64(q_abs);
    let q_monty = if q_val >= 0 {
        MontyForm::new(&q_u, *params)
    } else {
        let n_minus_q = n.wrapping_sub(&q_u);
        MontyForm::new(&n_minus_q, *params)
    };

    let mut u = one;
    let mut v = one;
    let mut qk = q_monty;

    let n_bits = bits - n.leading_zeros();
    if n_bits < 2 {
        let u_zero = u == zero_mf;
        return (u_zero, v, qk);
    }

    let start = n_bits - 1;
    let mut bs = start as i32 - 1;
    let two = MontyForm::new(&Uint::<N>::from_u8(2), *params);

    while bs >= b0 as i32 {
        u *= v;
        v = v * v - two * qk;
        qk = qk * qk;

        let do_step = bs == b0 as i32 || n.bit(bs as u32) == crypto_bigint::ConstChoice::TRUE;
        if do_step {
            qk *= q_monty;

            let old_u = u;
            let old_v = v;

            let sum = old_u + old_v;
            u = half_monty(sum, n, params);

            let neg_2q = if q_val >= 0 {
                let val = Uint::<N>::from_u64(2 * q_abs);
                let n_minus_val = n.wrapping_sub(&val);
                MontyForm::new(&n_minus_val, *params)
            } else {
                MontyForm::new(&Uint::<N>::from_u64(2 * q_abs), *params)
            };
            v = u + neg_2q * old_u;
        }

        bs -= 1;
    }

    let u_zero = u == zero_mf;
    (u_zero, v, qk)
}

/// Divide a Montgomery form value by 2 (mod n).
fn half_monty<const N: usize>(
    val: MontyForm<N>,
    n: &Uint<N>,
    params: &MontyParams<N>,
) -> MontyForm<N> {
    let v = val.retrieve();
    let result = if v.bit(0) == crypto_bigint::ConstChoice::TRUE {
        v.wrapping_add(n).shr(1)
    } else {
        v.shr(1)
    };
    MontyForm::new(&result, *params)
}

/// Miller-Rabin witness test using num-bigint.
fn miller_rabin_biguint(
    a: &num_bigint::BigUint,
    q: &num_bigint::BigUint,
    k: u32,
    nm1: &num_bigint::BigUint,
    n: &num_bigint::BigUint,
) -> bool {
    use num_bigint::BigUint;
    use num_traits::One;

    let one = BigUint::one();
    let mut y = a.modpow(q, n);

    if y == one || y == *nm1 {
        return true;
    }
    for _ in 1..k {
        y = y.modpow(&BigUint::from(2u32), n);
        if y == *nm1 {
            return true;
        }
    }
    false
}

/// Strong Lucas pseudoprime test using num-bigint. Works for any size.
fn strong_lucas_biguint(n: &num_bigint::BigUint) -> bool {
    use num_bigint::BigUint;
    use num_traits::One;

    let s = n.sqrt();
    if &s * &s == *n {
        return false;
    }

    let n_bi = BigInt::from(n.clone());
    let mut d_abs: u64 = 3;
    loop {
        d_abs += 2;
        let rem = &n_bi % BigInt::from(d_abs);
        let rem_abs = if rem < BigInt::ZERO {
            rem + BigInt::from(d_abs)
        } else {
            rem
        };
        if rem_abs.is_zero() {
            return false;
        }
        let rem_u64: u64 = rem_abs.try_into().unwrap_or(0);
        if jacobi_coprime(rem_u64, d_abs) != 1 {
            break;
        }
        if d_abs > 1_000_000 {
            return true;
        }
    }

    let q_val: i64 = if d_abs & 2 != 0 {
        (d_abs >> 2) as i64 + 1
    } else {
        -((d_abs >> 2) as i64)
    };

    let np1 = n + BigUint::one();
    let b0 = np1.trailing_zeros().unwrap_or(0);

    let (u_is_zero, mut v, mut qk) = lucas_mod_biguint(n, q_val, b0 as u32);

    if u_is_zero {
        return true;
    }

    let zero = BigUint::ZERO;
    let two = BigUint::from(2u32);
    let mut remaining = b0 as u32;

    while v != zero && remaining > 1 {
        remaining -= 1;
        let v_sq = (&v * &v) % n;
        let twoqk = (&two * &qk) % n;
        v = if v_sq >= twoqk {
            (v_sq - &twoqk) % n
        } else {
            (n - &twoqk + v_sq) % n
        };
        qk = (&qk * &qk) % n;
    }

    remaining > 0
}

/// Lucas sequence computation using num-bigint.
fn lucas_mod_biguint(
    n: &num_bigint::BigUint,
    q_val: i64,
    b0: u32,
) -> (bool, num_bigint::BigUint, num_bigint::BigUint) {
    use num_bigint::BigUint;
    use num_traits::One;

    let q_mod_n: BigUint = if q_val >= 0 {
        BigUint::from(q_val as u64) % n
    } else {
        let neg_q = BigUint::from(q_val.unsigned_abs());
        if neg_q <= *n {
            (n - neg_q) % n
        } else {
            let r = neg_q % n;
            if r.is_zero() {
                BigUint::ZERO
            } else {
                n - r
            }
        }
    };

    let one = BigUint::one();
    let two = BigUint::from(2u32);

    let mut u = one.clone();
    let mut v = one.clone();
    let mut qk = q_mod_n.clone();

    let n_bits = n.bits();
    if n_bits < 2 {
        return (u.is_zero(), v, qk);
    }

    let start = n_bits - 1;
    let mut bs = start as i64 - 1;

    while bs >= b0 as i64 {
        u = (&u * &v) % n;

        let v_sq = (&v * &v) % n;
        let twoqk = (&two * &qk) % n;
        v = if v_sq >= twoqk {
            (&v_sq - &twoqk) % n
        } else {
            (n - &twoqk + &v_sq) % n
        };

        qk = (&qk * &qk) % n;

        let do_step = bs == b0 as i64 || n.bit(bs as u64);
        if do_step {
            qk = (&qk * &q_mod_n) % n;

            let old_u = u.clone();
            let old_v = v.clone();

            let sum = (&old_u + &old_v) % n;
            u = biguint_half_mod(&sum, n);

            let neg_2q = if q_val >= 0 {
                let val = BigUint::from(2u64 * q_val.unsigned_abs());
                if val <= *n {
                    (n - &val) % n
                } else {
                    let r = &val % n;
                    if r.is_zero() {
                        BigUint::ZERO
                    } else {
                        n - &r
                    }
                }
            } else {
                BigUint::from(2u64 * q_val.unsigned_abs()) % n
            };
            v = (&u + &((&neg_2q * &old_u) % n)) % n;
        }

        bs -= 1;
    }

    (u.is_zero(), v, qk)
}

/// (val / 2) mod n for BigUint.
fn biguint_half_mod(val: &num_bigint::BigUint, n: &num_bigint::BigUint) -> num_bigint::BigUint {
    if val.bit(0) {
        (val + n) >> 1
    } else {
        val >> 1
    }
}

/// Legendre symbol `(a/p)` using crypto-bigint.
fn legendre_uint<const N: usize>(
    a: &Uint<N>,
    p_minus_1_half: &Uint<N>,
    params: &MontyParams<N>,
) -> i32 {
    let r = modpow_uint(a, p_minus_1_half, params);
    if r == Uint::<N>::ZERO {
        0
    } else if r == Uint::<N>::ONE {
        1
    } else {
        -1
    }
}

/// Square root modulo a prime `p` (Tonelli-Shanks).
///
/// Returns `Some(sqrt)` where `sqrt^2 ≡ a (mod p)`, or `None` if `a` is
/// not a quadratic residue mod `p`.
pub fn ibz_sqrt_mod_p(a: &Ibz, p: &Ibz) -> Option<Ibz> {
    let amod = ibz_mod(a, p);
    if amod.is_zero() {
        return Some(Ibz::zero());
    }

    if fits_in::<8>(p) {
        ibz_sqrt_mod_p_uint::<8>(&amod, p)
    } else if fits_in::<16>(p) {
        ibz_sqrt_mod_p_uint::<16>(&amod, p)
    } else {
        ibz_sqrt_mod_p_biguint(&amod, p)
    }
}

/// Fixed-precision path for sqrt_mod_p using `Uint<N>`.
fn ibz_sqrt_mod_p_uint<const N: usize>(amod: &Ibz, p: &Ibz) -> Option<Ibz> {
    let p_u: Uint<N> = bigint_to_uint(p);
    // p is a prime > 2, so it is odd.
    let p_odd: Odd<Uint<N>> = Option::from(Odd::new(p_u))?;
    let params = MontyParams::<N>::new_vartime(p_odd);
    let a_u: Uint<N> = bigint_to_uint(amod);

    let pm1 = p_u.wrapping_sub(&Uint::<N>::ONE);
    let pm1_half = pm1.shr(1);

    if legendre_uint(&a_u, &pm1_half, &params) != 1 {
        return None;
    }

    let p_mod_4 = ibz_mod_ui(p, 4);
    let p_mod_8 = ibz_mod_ui(p, 8);

    if p_mod_4 == 3 {
        let exp = p_u.wrapping_add(&Uint::<N>::ONE).shr(2);
        let result = modpow_uint(&a_u, &exp, &params);
        Some(uint_to_bigint(&result))
    } else if p_mod_8 == 5 {
        let exp_quarter = pm1.shr(2);
        let test = modpow_uint(&a_u, &exp_quarter, &params);
        if test == Uint::<N>::ONE {
            let exp = p_u.wrapping_add(&Uint::<N>::from_u8(3)).shr(3);
            let result = modpow_uint(&a_u, &exp, &params);
            Some(uint_to_bigint(&result))
        } else {
            let exp = p_u.wrapping_sub(&Uint::<N>::from_u8(5)).shr(3);
            let a_monty = MontyForm::new(&a_u, params);
            let four = Uint::<N>::from_u8(4);
            let a4_u = {
                let a4_monty = MontyForm::new(&four, params) * a_monty;
                a4_monty.retrieve()
            };
            let tmp = modpow_uint(&a4_u, &exp, &params);
            let tmp_monty = MontyForm::new(&tmp, params);
            let two = Uint::<N>::from_u8(2);
            let a2_monty = MontyForm::new(&two, params) * a_monty;
            let result = (a2_monty * tmp_monty).retrieve();
            Some(uint_to_bigint(&result))
        }
    } else {
        // p ≡ 1 (mod 8): Tonelli-Shanks
        let mut q = pm1;
        let mut e = 0u32;
        while q.bit(0) == crypto_bigint::ConstChoice::FALSE {
            q = q.shr(1);
            e += 1;
        }

        let mut qnr = Uint::<N>::from_u8(2);
        while legendre_uint(&qnr, &pm1_half, &params) != -1 {
            qnr = qnr.wrapping_add(&Uint::<N>::ONE);
        }

        let mut z_monty = MontyForm::new(&qnr, params).pow(&q);
        let mut y_monty = MontyForm::new(&a_u, params).pow(&q);
        let q_plus_1_half = q.wrapping_add(&Uint::<N>::ONE).shr(1);
        let mut x_monty = MontyForm::new(&a_u, params).pow(&q_plus_1_half);

        let neg_one_monty = MontyForm::new(&pm1, params);
        let mut exp_u = Uint::<N>::ONE.shl(e - 2);

        for _ in 0..e {
            let b = y_monty.pow(&exp_u);

            if b == neg_one_monty {
                x_monty *= z_monty;
                y_monty = y_monty * z_monty * z_monty;
            }

            z_monty = z_monty * z_monty;
            exp_u = exp_u.shr(1);
        }

        Some(uint_to_bigint(&x_monty.retrieve()))
    }
}

/// Fallback path for sqrt_mod_p when `p` exceeds all fixed-precision tiers.
fn ibz_sqrt_mod_p_biguint(amod: &Ibz, p: &Ibz) -> Option<Ibz> {
    use num_bigint::BigUint;
    use num_traits::One;

    let leg = ibz_legendre(amod, p);
    if leg != 1 {
        return None;
    }

    let p_mod_4 = ibz_mod_ui(p, 4);
    let p_mod_8 = ibz_mod_ui(p, 8);

    let (_, a_bytes) = amod.to_bytes_le();
    let (_, p_bytes) = p.to_bytes_le();
    let a_u = BigUint::from_bytes_le(&a_bytes);
    let p_u = BigUint::from_bytes_le(&p_bytes);

    if p_mod_4 == 3 {
        let exp: BigInt = (p + 1) / 4;
        let (_, exp_bytes) = exp.to_bytes_le();
        let exp_u = BigUint::from_bytes_le(&exp_bytes);
        let result = a_u.modpow(&exp_u, &p_u);
        Some(BigInt::from(result))
    } else if p_mod_8 == 5 {
        let pm1_u = &p_u - BigUint::one();
        let exp_quarter = &pm1_u / BigUint::from(4u32);
        let test = a_u.modpow(&exp_quarter, &p_u);
        if test.is_one() {
            let exp: BigInt = (p + 3) / 8;
            let (_, exp_bytes) = exp.to_bytes_le();
            let exp_u = BigUint::from_bytes_le(&exp_bytes);
            let result = a_u.modpow(&exp_u, &p_u);
            Some(BigInt::from(result))
        } else {
            let exp: BigInt = (p - 5) / 8;
            let (_, exp_bytes) = exp.to_bytes_le();
            let exp_u = BigUint::from_bytes_le(&exp_bytes);
            let a4_u = (&a_u * BigUint::from(4u32)) % &p_u;
            let tmp = a4_u.modpow(&exp_u, &p_u);
            let a2 = (&a_u * BigUint::from(2u32)) % &p_u;
            let result = (&a2 * &tmp) % &p_u;
            Some(BigInt::from(result))
        }
    } else {
        // p ≡ 1 (mod 8): Tonelli-Shanks
        let pm1_u = &p_u - BigUint::one();
        let mut q_u = pm1_u.clone();
        let mut e = 0u32;
        while !q_u.bit(0) {
            q_u >>= 1;
            e += 1;
        }

        let mut qnr = BigInt::from(2);
        while ibz_legendre(&qnr, p) != -1 {
            qnr += 1;
        }
        let (_, qnr_bytes) = qnr.to_bytes_le();
        let qnr_u = BigUint::from_bytes_le(&qnr_bytes);

        let mut z = qnr_u.modpow(&q_u, &p_u);
        let mut y = a_u.modpow(&q_u, &p_u);
        let q_plus_1_half = (&q_u + BigUint::one()) / BigUint::from(2u32);
        let mut x = a_u.modpow(&q_plus_1_half, &p_u);

        let mut exp_u = BigUint::one() << (e - 2) as usize;

        for _ in 0..e {
            let b = y.modpow(&exp_u, &p_u);

            if b == pm1_u {
                x = (&x * &z) % &p_u;
                y = (&y * &z * &z) % &p_u;
            }

            z = (&z * &z) % &p_u;
            exp_u >>= 1;
        }

        Some(BigInt::from(x))
    }
}

/// Legendre symbol `(a/p)`. Returns 1, -1, or 0.
///
/// Precondition: `p` must be an odd prime.
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    let amod = ibz_mod(a, p);

    if fits_in::<8>(p) {
        let p_u: Uint<8> = bigint_to_uint(p);
        let a_u: Uint<8> = bigint_to_uint(&amod);
        // p is an odd prime (precondition).
        let p_odd = Odd::new(p_u).expect("invariant: p is an odd prime");
        let params = MontyParams::<8>::new_vartime(p_odd);
        let pm1_half = p_u.wrapping_sub(&Uint::<8>::ONE).shr(1);
        legendre_uint(&a_u, &pm1_half, &params)
    } else if fits_in::<16>(p) {
        let p_u: Uint<16> = bigint_to_uint(p);
        let a_u: Uint<16> = bigint_to_uint(&amod);
        // p is an odd prime (precondition).
        let p_odd = Odd::new(p_u).expect("invariant: p is an odd prime");
        let params = MontyParams::<16>::new_vartime(p_odd);
        let pm1_half = p_u.wrapping_sub(&Uint::<16>::ONE).shr(1);
        legendre_uint(&a_u, &pm1_half, &params)
    } else {
        use num_bigint::BigUint;
        use num_traits::One;
        let exp: BigInt = (p - 1) / 2;
        let (_, a_bytes) = amod.to_bytes_le();
        let (_, exp_bytes) = exp.to_bytes_le();
        let (_, p_bytes) = p.to_bytes_le();
        let r = BigUint::from_bytes_le(&a_bytes).modpow(
            &BigUint::from_bytes_le(&exp_bytes),
            &BigUint::from_bytes_le(&p_bytes),
        );
        if r.is_zero() {
            0
        } else if r.is_one() {
            1
        } else {
            -1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::intbig as old;

    #[test]
    fn test_probab_prime_matches_old() {
        let primes: &[i64] = &[2, 3, 5, 7, 13, 97, 1009, 1000000007];
        for &p in primes {
            let n = BigInt::from(p);
            let new_res = ibz_probab_prime(&n, 32);
            assert!(new_res > 0, "expected prime for {p}, got {new_res}");
        }
        let composites: &[i64] = &[4, 6, 8, 9, 15, 100, 1000000006];
        for &c in composites {
            let n = BigInt::from(c);
            let new_res = ibz_probab_prime(&n, 32);
            assert_eq!(new_res, 0, "expected composite for {c}, got {new_res}");
        }
    }

    #[test]
    fn test_probab_prime_large() {
        // A known 256-bit prime (secp256k1 order)
        let n = BigInt::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap();
        assert!(ibz_probab_prime(&n, 32) > 0);
        // n-1 is composite
        let nm1 = &n - 1;
        assert_eq!(ibz_probab_prime(&nm1, 32), 0);
    }

    #[test]
    fn test_jacobi_coprime() {
        assert_eq!(jacobi_coprime(2, 3), -1);
        assert_eq!(jacobi_coprime(1, 3), 1);
        assert_eq!(jacobi_coprime(2, 7), 1);
        assert_eq!(jacobi_coprime(3, 7), -1);
    }

    #[test]
    fn test_strong_lucas() {
        for &p in &[5u64, 7, 11, 13, 17, 101, 1009, 10007] {
            let n = Uint::<8>::from_u64(p);
            let n_odd = Odd::new(n).unwrap();
            let params = MontyParams::<8>::new_vartime(n_odd);
            assert!(
                strong_lucas_uint(&n, &params, 512),
                "Strong Lucas failed for prime {p}"
            );
        }
    }

    #[test]
    fn test_sqrt_mod_p_matches_old() {
        // p ≡ 3 (mod 4)
        let p = BigInt::from(23);
        for a in 0..23 {
            let a_bi = BigInt::from(a);
            let new_res = ibz_sqrt_mod_p(&a_bi, &p);
            let old_res = old::ibz_sqrt_mod_p(&a_bi, &p);
            assert_eq!(new_res, old_res, "sqrt_mod_p mismatch for a={a}, p=23");
        }

        // p ≡ 1 (mod 8)
        let p = BigInt::from(41);
        for a in 0..41 {
            let a_bi = BigInt::from(a);
            let new_res = ibz_sqrt_mod_p(&a_bi, &p);
            let old_res = old::ibz_sqrt_mod_p(&a_bi, &p);
            assert_eq!(new_res, old_res, "sqrt_mod_p mismatch for a={a}, p=41");
        }
    }

    #[test]
    fn test_legendre_matches_old() {
        let p = BigInt::from(1000000007);
        for a in 1..100 {
            let a_bi = BigInt::from(a);
            let new_res = ibz_legendre(&a_bi, &p);
            let old_res = old::ibz_legendre(&a_bi, &p);
            assert_eq!(new_res, old_res, "legendre mismatch for a={a}");
        }
    }

    #[test]
    fn bench_probab_prime() {
        // 271-bit odd composite
        let n = BigInt::parse_bytes(
            b"5F2A3B4C5D6E7F8091A2B3C4D5E6F7081920A1B2C3D4E5F60718293A4B5C6D7F",
            16,
        )
        .unwrap()
            | BigInt::from(1);
        let iters = 10000;
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            let _ = ibz_probab_prime(&n, 32);
        }
        let elapsed = t0.elapsed();
        eprintln!(
            "crypto-bigint ibz_probab_prime: {:.2}us/call ({iters} iters in {:.0}ms)",
            elapsed.as_micros() as f64 / iters as f64,
            elapsed.as_millis()
        );
    }

    #[test]
    fn regression_bpsw_pseudoprimes() {
        let spsp2: &[u64] = &[
            2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 80581, 85489,
            88357, 90751,
        ];
        for &n in spsp2 {
            assert_eq!(
                ibz_probab_prime(&BigInt::from(n), 32),
                0,
                "strong psp base 2: {n} should be composite"
            );
        }

        let carmichael: &[u64] = &[561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341];
        for &n in carmichael {
            assert_eq!(
                ibz_probab_prime(&BigInt::from(n), 32),
                0,
                "Carmichael number {n} should be composite"
            );
        }

        let primes: &[u64] = &[
            2053, 3271, 4049, 4673, 8329, 15859, 29347, 42797, 52631, 65287,
        ];
        for &p in primes {
            assert!(
                ibz_probab_prime(&BigInt::from(p), 32) > 0,
                "actual prime {p} should be identified as prime"
            );
        }

        let spsp23: &[u64] = &[1373653, 1530787];
        for &n in spsp23 {
            assert_eq!(
                ibz_probab_prime(&BigInt::from(n), 32),
                0,
                "strong psp bases 2,3: {n} should be composite"
            );
        }

        let p = BigInt::from(1000000007i64);
        assert!(ibz_probab_prime(&(-&p), 32) > 0);

        assert_eq!(ibz_probab_prime(&BigInt::from(0), 32), 0);
        assert_eq!(ibz_probab_prime(&BigInt::from(1), 32), 0);
        assert!(ibz_probab_prime(&BigInt::from(2), 32) > 0);
        assert!(ibz_probab_prime(&BigInt::from(3), 32) > 0);
    }

    #[test]
    fn test_probab_prime_oversized_uses_uint16() {
        // 521-bit prime (exceeds 512-bit Uint<8>, uses Uint<16> tier)
        let mersenne521 = BigInt::parse_bytes(
            b"6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151",
            10,
        ).unwrap();
        assert!(
            !fits_in::<8>(&mersenne521),
            "test value must exceed Uint<8>"
        );
        assert!(
            fits_in::<16>(&mersenne521),
            "test value must fit in Uint<16>"
        );
        assert!(ibz_probab_prime(&mersenne521, 32) > 0, "M521 is prime");

        // M521 - 2 is even, so composite
        assert_eq!(ibz_probab_prime(&(&mersenne521 - 2), 32), 0);

        // M521 + 2 should be composite
        assert_eq!(ibz_probab_prime(&(&mersenne521 + 2), 32), 0);
    }
}
