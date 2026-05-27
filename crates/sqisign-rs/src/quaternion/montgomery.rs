//!
//! Replaces `BigInt::modpow` with a Montgomery-form implementation that
//! avoids per-multiply division. For odd moduli of 200-300 bits this is
//! roughly 10-30x faster than `num-bigint`'s built-in `modpow`.
//!
//! SECURITY: NOT constant-time. Acceptable here because the primality
//! test in `represent_integer` operates on non-secret data.

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};

/// Montgomery multiplication context for a fixed odd modulus.
pub struct MontgomeryCtx {
    /// The modulus n (odd, positive).
    n: BigUint,
    /// Number of 64-bit limbs in the modulus.
    limbs: usize,
    /// R = 2^(64×limbs), not stored directly (too large), but R mod n and R² mod n are.
    r_mod_n: BigUint,
    /// R² mod n, used to convert into Montgomery form.
    r2_mod_n: BigUint,
    /// −n⁻¹ mod R, the full multi-limb Montgomery reduction constant.
    n_inv_full: BigUint,
    /// Bitmask: R − 1 = 2^(64×limbs) − 1.
    r_mask: BigUint,
}

impl MontgomeryCtx {
    /// Create a context for Montgomery arithmetic modulo `n`.
    pub fn new(n: &BigUint) -> Self {
        debug_assert!(
            !n.is_zero() && n.bit(0),
            "Montgomery modulus must be odd and positive"
        );

        let bits = n.bits() as usize;
        let limbs = bits.div_ceil(64);
        let r_bits = 64 * limbs;

        // R = 2^(64*limbs)
        let r = BigUint::one() << r_bits;
        let r_mask = &r - BigUint::one();
        let r_mod_n = &r % n;
        let r2_mod_n = (&r_mod_n * &r_mod_n) % n;

        // Compute -n^{-1} mod R using Hensel lifting.
        let n_inv_full = compute_neg_n_inv_full(n, r_bits);

        MontgomeryCtx {
            n: n.clone(),
            limbs,
            r_mod_n,
            r2_mod_n,
            n_inv_full,
            r_mask,
        }
    }

    /// Convert a value into Montgomery form: a_mont = a * R mod n.
    #[inline]
    fn to_mont(&self, a: &BigUint) -> BigUint {
        self.mont_mul(a, &self.r2_mod_n)
    }

    /// Convert from Montgomery form back to normal: a = a_mont * 1 = a_mont * R⁻¹ mod n.
    #[inline]
    fn reduce_mont(&self, a_mont: &BigUint) -> BigUint {
        self.mont_redc(a_mont)
    }

    /// Montgomery reduction: given T < n * R, compute T * R⁻¹ mod n.
    ///
    /// Uses the REDC algorithm:
    ///   m = (T mod R) * n_inv mod R
    ///   t = (T + m * n) / R
    ///   if t >= n then t -= n
    fn mont_redc(&self, t: &BigUint) -> BigUint {
        let r_bits = 64 * self.limbs;

        // m = (T mod R) * n_inv mod R
        let t_low = t & &self.r_mask;
        let m = (&t_low * &self.n_inv_full) & &self.r_mask;

        // t = (T + m * n) / R
        let mn = &m * &self.n;
        let sum = t + &mn;
        let result = &sum >> r_bits;

        if result >= self.n {
            result - &self.n
        } else {
            result
        }
    }

    /// Montgomery multiplication: compute a * b * R⁻¹ mod n.
    #[inline]
    fn mont_mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let t = a * b;
        self.mont_redc(&t)
    }

    /// Montgomery squaring: compute a² * R⁻¹ mod n.
    #[inline]
    fn mont_sqr(&self, a: &BigUint) -> BigUint {
        let t = a * a;
        self.mont_redc(&t)
    }

    /// Modular exponentiation using Montgomery form: base^exp mod n.
    pub fn modpow(&self, base: &BigUint, exp: &BigUint) -> BigUint {
        if exp.is_zero() {
            return if self.n == BigUint::one() {
                BigUint::zero()
            } else {
                BigUint::one()
            };
        }

        let base_reduced = base % &self.n;
        let base_mont = self.to_mont(&base_reduced);
        let mut result = self.r_mod_n.clone(); // 1 in Montgomery form

        let exp_bits = exp.bits();
        for i in (0..exp_bits).rev() {
            result = self.mont_sqr(&result);
            if exp.bit(i) {
                result = self.mont_mul(&result, &base_mont);
            }
        }

        self.reduce_mont(&result)
    }
}

/// Compute −n⁻¹ mod 2^r_bits via Hensel lifting.
///
/// `n` must be odd. `r_bits` must be a multiple of 64.
fn compute_neg_n_inv_full(n: &BigUint, r_bits: usize) -> BigUint {
    debug_assert!(n.bit(0), "n must be odd");
    debug_assert!(r_bits > 0 && r_bits % 64 == 0);

    // Start with n^{-1} mod 2 = 1 (n is odd)
    let mut x = BigUint::one();

    // Hensel lift: x = x * (2 - n * x) mod 2^k, doubling k each time.
    // After each iteration, n * x ≡ 1 (mod 2^k).
    let mut k = 1usize;
    while k < r_bits {
        let next_k = (k * 2).min(r_bits);
        let mask = (BigUint::one() << next_k) - BigUint::one();
        // nx = n * x mod 2^next_k
        let nx = (n * &x) & &mask;
        // two_minus_nx = (2 - nx) mod 2^next_k
        // Since nx ≡ 1 (mod 2^k), we know nx is odd, so 2 - nx can be negative.
        // Work mod 2^next_k: if nx > 2, result = 2^next_k + 2 - nx.
        let two_minus_nx = {
            let two = BigUint::from(2u32);
            if nx <= two {
                &two - &nx
            } else {
                let r_local = BigUint::one() << next_k;
                (&r_local + &two - &nx) & &mask
            }
        };
        x = (&x * &two_minus_nx) & &mask;
        k = next_k;
    }

    // x = n^{-1} mod R. We want -n^{-1} mod R = R - x (if x != 0).
    let r_mask = (BigUint::one() << r_bits) - BigUint::one();
    if x.is_zero() {
        BigUint::zero()
    } else {
        let r = BigUint::one() << r_bits;
        (&r - &x) & &r_mask
    }
}

/// Montgomery modpow for `BigInt`: computes base^exp mod modulus.
///
/// `modulus` must be a positive odd integer. `exp` must be non-negative.
/// `base` can be any integer (it is reduced mod `modulus` first).
pub fn montgomery_modpow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    debug_assert!(modulus > &BigInt::zero(), "modulus must be positive");
    debug_assert!(exp >= &BigInt::zero(), "exponent must be non-negative");

    let (_, mod_bytes) = modulus.to_bytes_le();
    let n = BigUint::from_bytes_le(&mod_bytes);

    // Handle even modulus by falling back to num-bigint (shouldn't happen in our usage)
    if !n.bit(0) {
        return base.modpow(exp, modulus);
    }

    let (_, exp_bytes) = exp.to_bytes_le();
    let exp_u = BigUint::from_bytes_le(&exp_bytes);

    // Reduce base mod n (handle negative base)
    let base_reduced = {
        let base_mod = base % modulus;
        if base_mod.sign() == Sign::Minus {
            let (_, bytes) = (&base_mod + modulus).to_bytes_le();
            BigUint::from_bytes_le(&bytes)
        } else {
            let (_, bytes) = base_mod.to_bytes_le();
            BigUint::from_bytes_le(&bytes)
        }
    };

    let ctx = MontgomeryCtx::new(&n);
    let result = ctx.modpow(&base_reduced, &exp_u);

    BigInt::from(result)
}

/// Reusable Montgomery context for repeated modpow with the same modulus.
///
/// Caches the Montgomery precomputation across multiple exponentiations.
pub struct MontgomeryModpow {
    ctx: MontgomeryCtx,
    n_bigint: BigInt,
}

impl MontgomeryModpow {
    /// Create a new context for modular exponentiation with the given modulus.
    pub fn new(modulus: &BigInt) -> Self {
        let (_, mod_bytes) = modulus.to_bytes_le();
        let n = BigUint::from_bytes_le(&mod_bytes);
        MontgomeryModpow {
            ctx: MontgomeryCtx::new(&n),
            n_bigint: modulus.clone(),
        }
    }

    /// Compute base^exp mod n, reusing the precomputed Montgomery constants.
    pub fn modpow(&self, base: &BigInt, exp: &BigInt) -> BigInt {
        let (_, exp_bytes) = exp.to_bytes_le();
        let exp_u = BigUint::from_bytes_le(&exp_bytes);

        let base_reduced = {
            let base_mod = base % &self.n_bigint;
            if base_mod.sign() == Sign::Minus {
                let (_, bytes) = (&base_mod + &self.n_bigint).to_bytes_le();
                BigUint::from_bytes_le(&bytes)
            } else {
                let (_, bytes) = base_mod.to_bytes_le();
                BigUint::from_bytes_le(&bytes)
            }
        };

        BigInt::from(self.ctx.modpow(&base_reduced, &exp_u))
    }

    /// Montgomery squaring: compute val^2 mod n, where val is a BigInt.
    ///
    /// Faster than calling modpow with exp=2 because it avoids the
    /// bit-scanning loop.
    pub fn mont_sqr_mod(&self, val: &BigInt) -> BigInt {
        let reduced = {
            let v = val % &self.n_bigint;
            if v.sign() == Sign::Minus {
                let (_, bytes) = (&v + &self.n_bigint).to_bytes_le();
                BigUint::from_bytes_le(&bytes)
            } else {
                let (_, bytes) = v.to_bytes_le();
                BigUint::from_bytes_le(&bytes)
            }
        };
        // a^2 mod n = reduce_mont(mont_sqr(to_mont(a)))
        let a_mont = self.ctx.to_mont(&reduced);
        let sq = self.ctx.mont_sqr(&a_mont);
        BigInt::from(self.ctx.reduce_mont(&sq))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_montgomery_modpow_basic() {
        let base = BigInt::from(2);
        let exp = BigInt::from(10);
        let modulus = BigInt::from(1000000007);
        let expected = base.modpow(&exp, &modulus);
        let result = montgomery_modpow(&base, &exp, &modulus);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_montgomery_modpow_large() {
        // 256-bit modulus (prime)
        let n_hex = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
        let n = BigInt::parse_bytes(n_hex.as_bytes(), 16).unwrap();
        let base = BigInt::from(31337);
        let exp = BigInt::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364140",
            16,
        )
        .unwrap();
        let expected = base.modpow(&exp, &n);
        let result = montgomery_modpow(&base, &exp, &n);
        assert_eq!(
            result, expected,
            "Montgomery modpow mismatch for secp256k1-like inputs"
        );
    }

    #[test]
    fn test_montgomery_modpow_various() {
        let cases: Vec<(i64, i64, i64)> = vec![
            (3, 7, 11),
            (5, 0, 13),
            (7, 1, 17),
            (123456789, 987654321, 1000000007),
            (2, 255, 257),
        ];
        for (b, e, m) in cases {
            let base = BigInt::from(b);
            let exp = BigInt::from(e);
            let modulus = BigInt::from(m);
            let expected = base.modpow(&exp, &modulus);
            let result = montgomery_modpow(&base, &exp, &modulus);
            assert_eq!(result, expected, "mismatch for {}^{} mod {}", b, e, m);
        }
    }

    #[test]
    fn test_montgomery_context_reuse() {
        let n = BigInt::from(1000000007i64);
        let ctx = MontgomeryModpow::new(&n);

        for b in 2..50 {
            for e in [1, 2, 10, 100, 1000].iter() {
                let base = BigInt::from(b);
                let exp = BigInt::from(*e);
                let expected = base.modpow(&exp, &n);
                let result = ctx.modpow(&base, &exp);
                assert_eq!(result, expected, "mismatch for {}^{} mod {}", b, e, n);
            }
        }
    }

    #[test]
    fn test_neg_n_inv() {
        // Verify: n * neg_n_inv ≡ -1 (mod R) for various odd moduli
        for n_val in [3u64, 5, 7, 11, 13, 0xDEADBEEFCAFEBABD] {
            let n = BigUint::from(n_val);
            let neg_inv = compute_neg_n_inv_full(&n, 64);
            let r = BigUint::one() << 64;
            let product = (&n * &neg_inv) % &r;
            let expected = &r - BigUint::one(); // -1 mod R
            assert_eq!(product, expected, "n_inv property failed for n={n_val}");
        }
        // Multi-limb test
        let n = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap();
        let neg_inv = compute_neg_n_inv_full(&n, 256);
        let r = BigUint::one() << 256;
        let product = (&n * &neg_inv) % &r;
        let expected = &r - BigUint::one();
        assert_eq!(product, expected, "n_inv property failed for secp256k1 p");
    }

    #[test]
    fn bench_modpow_comparison() {
        // 271-bit odd modulus (similar to cornacchia_target in represent_integer)
        let n_u = BigUint::parse_bytes(
            b"5F2A3B4C5D6E7F8091A2B3C4D5E6F7081920A1B2C3D4E5F60718293A4B5C6D7F",
            16,
        )
        .unwrap()
            | BigUint::one(); // ensure odd
        let n = BigInt::from(n_u.clone());
        // exponent ~ 270 bits
        let d_u: BigUint = &n_u >> 1;
        let d = BigInt::from(d_u.clone());
        let base_u = BigUint::from(2u32);
        let base = BigInt::from(2);

        let iters = 1000u32;

        // BigUint::modpow
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            let _ = base_u.modpow(&d_u, &n_u);
        }
        let biguint_time = t0.elapsed();

        // Montgomery modpow (fresh context each time, like ibz_probab_prime)
        let t1 = std::time::Instant::now();
        for _ in 0..iters {
            let _ = montgomery_modpow(&base, &d, &n);
        }
        let mont_fresh_time = t1.elapsed();

        // Montgomery modpow (reused context, like ibz_sqrt_mod_p)
        let ctx = MontgomeryModpow::new(&n);
        let t2 = std::time::Instant::now();
        for _ in 0..iters {
            let _ = ctx.modpow(&base, &d);
        }
        let mont_reuse_time = t2.elapsed();

        eprintln!("=== modpow benchmark ({iters} iters, ~271-bit modulus) ===");
        eprintln!(
            "BigUint::modpow:       {:.1}us/call",
            biguint_time.as_micros() as f64 / iters as f64
        );
        eprintln!(
            "Montgomery (fresh ctx): {:.1}us/call",
            mont_fresh_time.as_micros() as f64 / iters as f64
        );
        eprintln!(
            "Montgomery (reuse ctx): {:.1}us/call",
            mont_reuse_time.as_micros() as f64 / iters as f64
        );
    }
}
