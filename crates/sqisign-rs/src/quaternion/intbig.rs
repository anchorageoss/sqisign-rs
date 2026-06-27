//!
//! Provides the `Ibz` type (an alias for `num_bigint::BigInt`) and free
//! functions for arbitrary-precision integer arithmetic. Division follows
//! truncated (C99) semantics (`ibz_div`), while `ibz_mod` always returns
//! a non-negative result.
//!
//! SECURITY: `num-bigint` is NOT constant-time. Signing path only; not used in verification.

use alloc::{
    format,
    string::{String, ToString},
    vec,
    vec::Vec,
};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rand::Rng;

/// Signed arbitrary-precision integer.
pub type Ibz = BigInt;

/// Best-effort zeroization of a BigInt.
///
/// Overwrites the value with zero. The old heap allocation backing the
/// BigInt's internal `Vec<u64>` is freed but NOT scrubbed, `num-bigint`
/// does not expose its backing storage. Use the ZeroizingAllocator (Tier 3)
/// for comprehensive heap scrubbing.
pub fn ibz_zeroize(v: &mut Ibz) {
    // Assign to a new zero, dropping the old allocation.
    *v = Ibz::zero();
}

#[allow(non_snake_case)]
mod constants {
    use super::*;
    lazy_static_ibz! {
        pub static IBZ_ZERO: Ibz = Ibz::zero();
        pub static IBZ_ONE: Ibz = Ibz::one();
        pub static IBZ_TWO: Ibz = BigInt::from(2);
        pub static IBZ_THREE: Ibz = BigInt::from(3);
    }
}
pub use constants::*;

/// Convenience macro for small BigInt constants. `BigInt` is not
/// const-constructible (it owns a `Vec`), so each constant is an accessor
/// function that builds the value on call. The values are tiny (0..3), so this
/// is cheap; it needs no synchronization and no `std`.
#[doc(hidden)]
macro_rules! lazy_static_ibz {
    ($($(#[$meta:meta])* pub static $name:ident : $ty:ty = $init:expr;)*) => {
        $(
            $(#[$meta])*
            pub fn $name() -> $ty {
                $init
            }
        )*
    };
}

// Expand the macro before its first use, Rust requires the macro to be
// defined before the call site in the source file. We place the actual
// definitions here.
use lazy_static_ibz;

// Re-generate the constants with the macro now that it's defined.
// (The block above already expanded them.)

/// Truncated division: `a = q * b + r` where `|r| < |b|` and `sign(r) == sign(a)`.
pub fn ibz_div(a: &Ibz, b: &Ibz) -> (Ibz, Ibz) {
    let (q, r) = a.div_rem(b);
    (q, r)
}

/// Truncated division by `2^exp`: right-shift of `|a|` with sign preserved.
pub fn ibz_div_2exp(a: &Ibz, exp: u32) -> Ibz {
    if a.is_negative() {
        -((-a) >> exp as usize)
    } else {
        a >> exp as usize
    }
}

/// Floor division: `a = q * b + r` where `0 <= r < |b|` when `b > 0`.
pub fn ibz_div_floor(n: &Ibz, d: &Ibz) -> (Ibz, Ibz) {
    n.div_mod_floor(d)
}

/// Non-negative modular reduction: result is always in `[0, |b|)`.
pub fn ibz_mod(a: &Ibz, b: &Ibz) -> Ibz {
    a.mod_floor(&b.abs())
}

/// `n mod d` for small unsigned `d`. Always non-negative.
pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    let d_big = BigInt::from(d);
    let r = n.mod_floor(&d_big);
    r.to_u64().unwrap_or(0)
}

/// Test if `a` is divisible by `b` (`a mod b == 0`).
pub fn ibz_divides(a: &Ibz, b: &Ibz) -> bool {
    if b.is_zero() {
        return a.is_zero();
    }
    (a % b).is_zero()
}

/// `x^e` for unsigned exponent.
pub fn ibz_pow(x: &Ibz, e: u32) -> Ibz {
    num_traits::pow::Pow::pow(x, e as usize)
}

/// `x^e mod m` (non-negative exponent).
pub fn ibz_pow_mod(x: &Ibz, e: &Ibz, m: &Ibz) -> Ibz {
    assert!(
        e.sign() != Sign::Minus,
        "ibz_pow_mod: negative exponent not supported"
    );
    let x_mod = ibz_mod(x, m);
    let (_, x_bytes) = x_mod.to_bytes_le();
    let (_, e_bytes) = e.to_bytes_le();
    let (_, m_bytes) = m.to_bytes_le();
    let x_u = num_bigint::BigUint::from_bytes_le(&x_bytes);
    let e_u = num_bigint::BigUint::from_bytes_le(&e_bytes);
    let m_u = num_bigint::BigUint::from_bytes_le(&m_bytes);
    let r = x_u.modpow(&e_u, &m_u);
    BigInt::from(r)
}

/// Two-adic valuation: index of the lowest set bit.
///
/// Returns 0 for odd numbers, and the number of trailing zero bits
/// otherwise. For zero, returns 0.
pub fn ibz_two_adic(x: &Ibz) -> u32 {
    if x.is_zero() {
        return 0;
    }
    let (_, bytes) = x.to_bytes_le();
    let mut count: u32 = 0;
    for &byte in &bytes {
        if byte == 0 {
            count += 8;
        } else {
            count += byte.trailing_zeros();
            break;
        }
    }
    count
}

/// Bit size (number of bits needed to represent `|x|`).
///
/// Returns 1 for zero.
pub fn ibz_bitsize(a: &Ibz) -> u32 {
    if a.is_zero() {
        return 1;
    }
    a.bits() as u32
}

/// Extract the low 32 bits with the original sign (sign + low 31 bits).
pub fn ibz_get(x: &Ibz) -> i32 {
    let (sign, digits) = x.to_u64_digits();
    let low = if digits.is_empty() { 0u64 } else { digits[0] };
    let magnitude = (low & 0x7FFF_FFFF) as i32;
    match sign {
        Sign::Minus => {
            if low & 0x8000_0000 != 0 && magnitude == 0 {
                i32::MIN
            } else {
                -magnitude
            }
        }
        _ => {
            if low & 0x8000_0000 != 0 {
                low as u32 as i32
            } else {
                magnitude
            }
        }
    }
}

/// Set from a string in the given base (10 or 16). Returns `Some(value)`
/// on success, `None` on parse failure.
pub fn ibz_set_from_str(s: &str, base: u32) -> Option<Ibz> {
    // num-bigint doesn't have a direct "parse with base" on BigInt,
    // but we can use the standard radix parsing
    if base == 10 {
        s.parse::<BigInt>().ok()
    } else if base == 16 {
        let (neg, hex) = if let Some(rest) = s.strip_prefix('-') {
            (true, rest)
        } else {
            (false, s)
        };
        // Pad to even length so byte pairs align correctly
        let hex_padded = if hex.len() % 2 == 1 {
            format!("0{}", hex)
        } else {
            hex.to_string()
        };
        let bytes: Vec<u8> = (0..hex_padded.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&hex_padded[i..i + 2], 16))
            .collect::<Result<Vec<_>, _>>()
            .ok()?;
        let val = BigInt::from_bytes_be(Sign::Plus, &bytes);
        if neg {
            Some(-val)
        } else {
            Some(val)
        }
    } else {
        None
    }
}

/// Convert to string in the given base (10 or 16).
pub fn ibz_to_str(x: &Ibz, base: u32) -> String {
    match base {
        10 => format!("{}", x),
        16 => format!("{:x}", x),
        _ => String::new(),
    }
}

/// Import from little-endian `u64` digit array (unsigned / non-negative).
pub fn ibz_copy_digits(digits: &[u64]) -> Ibz {
    if digits.is_empty() {
        return Ibz::zero();
    }
    // Build big-endian byte array from little-endian u64 limbs
    let mut bytes = Vec::with_capacity(digits.len() * 8);
    for &d in digits.iter().rev() {
        bytes.extend_from_slice(&d.to_be_bytes());
    }
    BigInt::from_bytes_be(Sign::Plus, &bytes)
}

/// Export to little-endian `u64` digit array. Assumes `x >= 0`.
///
/// Writes into `target`; caller must ensure it is large enough.
pub fn ibz_to_digits(x: &Ibz, target: &mut [u64]) {
    for t in target.iter_mut() {
        *t = 0;
    }
    if x.is_zero() {
        return;
    }
    let (_, bytes) = x.to_bytes_le();
    for (i, chunk) in bytes.chunks(8).enumerate() {
        if i >= target.len() {
            break;
        }
        let mut buf = [0u8; 8];
        buf[..chunk.len()].copy_from_slice(chunk);
        target[i] = u64::from_le_bytes(buf);
    }
}

/// BPSW primality test (Miller-Rabin base 2 + strong Lucas).
///
/// Returns 0 if definitely composite, 1 if probably prime, 2 if
/// certainly prime (for small values).
pub fn ibz_probab_prime(n: &Ibz, reps: u32) -> i32 {
    super::fast_modpow::ibz_probab_prime(n, reps)
}

/// GCD of two integers.
pub fn ibz_gcd(a: &Ibz, b: &Ibz) -> Ibz {
    a.gcd(b)
}

/// Extended GCD: returns `(gcd, u, v)` where `gcd = a*u + b*v`.
pub fn ibz_xgcd(a: &Ibz, b: &Ibz) -> (Ibz, Ibz, Ibz) {
    let g = a.extended_gcd(b);
    (g.gcd, g.x, g.y)
}

/// Modular inverse: returns `Some(inv)` where `inv` is in `[0, m)` and
/// `a * inv ≡ 1 (mod m)`, or `None` if no inverse exists.
pub fn ibz_invmod(a: &Ibz, m: &Ibz) -> Option<Ibz> {
    let g = a.extended_gcd(m);
    if !g.gcd.is_one() {
        return None;
    }
    Some(ibz_mod(&g.x, m))
}

/// Legendre symbol `(a/p)`. Returns 1, -1, or 0.
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    super::fast_modpow::ibz_legendre(a, p)
}

/// Integer square root of a perfect square.
///
/// Returns `Some(sqrt)` if `a` is a perfect square, `None` otherwise.
pub fn ibz_sqrt(a: &Ibz) -> Option<Ibz> {
    if a.is_negative() {
        return None;
    }
    let s = a.sqrt();
    if &(&s * &s) == a {
        Some(s)
    } else {
        None
    }
}

/// Floor of the integer square root: `floor(sqrt(a))`.
pub fn ibz_sqrt_floor(a: &Ibz) -> Ibz {
    a.sqrt()
}

/// Square root modulo a prime `p` (Tonelli-Shanks).
///
/// Returns `Some(sqrt)` where `sqrt^2 ≡ a (mod p)`, or `None` if `a` is
/// not a quadratic residue mod `p`.
///
/// Delegates to the `crypto-bigint`-backed implementation in `fast_modpow`.
pub fn ibz_sqrt_mod_p(a: &Ibz, p: &Ibz) -> Option<Ibz> {
    super::fast_modpow::ibz_sqrt_mod_p(a, p)
}

/// Random integer in `[a, b]` (inclusive) via rejection sampling.
///
/// Each attempt reads `ceil(bit_length(b-a) / 8)` bytes from `rng`,
/// loads them little-endian into limbs, masks the top bits, and accepts
/// if the result is `<= b - a`.
pub fn ibz_rand_interval(rng: &mut impl Rng, a: &Ibz, b: &Ibz) -> Ibz {
    let bmina = b - a;

    if bmina.is_zero() {
        return a.clone();
    }

    let len_bits = ibz_bitsize(&bmina) as usize;
    let len_bytes = len_bits.div_ceil(8);
    let sizeof_limb: usize = 8;
    let sizeof_limb_bits = sizeof_limb * 8;
    let len_limbs = len_bytes.div_ceil(sizeof_limb);

    let mask: u64 = if len_bits % sizeof_limb_bits == 0 {
        u64::MAX
    } else {
        (1u64 << (len_bits % sizeof_limb_bits)) - 1
    };

    loop {
        let mut bytes = vec![0u8; len_bytes];
        rng.fill_bytes(&mut bytes);

        let mut limbs = vec![0u64; len_limbs];
        for (i, chunk) in bytes.chunks(8).enumerate() {
            let mut word = [0u8; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            limbs[i] = u64::from_le_bytes(word);
        }

        limbs[len_limbs - 1] &= mask;

        let mut result_bytes = vec![0u8; len_limbs * 8];
        for (i, &limb) in limbs.iter().enumerate() {
            result_bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb.to_le_bytes());
        }
        let r = BigInt::from_bytes_le(Sign::Plus, &result_bytes);

        if r <= bmina {
            return a + r;
        }
    }
}

/// Random integer in `[a, b]` for small non-negative `i32` values.
///
/// Reads 4 bytes per rejection sampling attempt, interprets as
/// little-endian u32, masks to the required bit width, and accepts
/// if `<= b - a`.
pub fn ibz_rand_interval_i(rng: &mut impl Rng, a: i32, b: i32) -> Ibz {
    assert!(a >= 0 && b >= 0 && b > a);
    let diff = (b - a) as u32;
    let bits = 32 - diff.leading_zeros();
    let mask = (1u32 << bits) - 1;

    loop {
        let mut buf = [0u8; 4];
        rng.fill_bytes(&mut buf);
        let rand32 = u32::from_le_bytes(buf) & mask;
        if rand32 <= diff {
            return BigInt::from(rand32 as i32 + a);
        }
    }
}

/// Random integer in `[-m, m]`.
pub fn ibz_rand_interval_minm_m(rng: &mut impl Rng, m: i32) -> Ibz {
    let two_m = BigInt::from(2 * m as i64);
    let r = ibz_rand_interval(rng, &Ibz::zero(), &two_m);
    r - m
}

/// Random integer in `[−2ᵐ − m, 2ᵐ − m]`.
///
/// Samples uniformly from `[−2ᵐ, 2ᵐ]` then subtracts `m`.
pub fn ibz_rand_interval_bits(rng: &mut impl Rng, m: u32) -> Ibz {
    let bound = BigInt::one() << m as usize;
    let low = -&bound;
    let r = ibz_rand_interval(rng, &low, &bound);
    r - BigInt::from(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_str_10(s: &str) -> Ibz {
        ibz_set_from_str(s, 10).unwrap()
    }

    fn from_str_16(s: &str) -> Ibz {
        ibz_set_from_str(s, 16).unwrap()
    }

    #[test]
    fn test_init_set_cmp() {
        let a_init = Ibz::zero();
        assert!(a_init.is_zero());

        let mut a = BigInt::from(1);
        assert!(a.is_one());
        assert!(!a.is_zero());

        let b = Ibz::zero();
        let c = Ibz::zero();
        assert_eq!(b, c);
        assert_ne!(a, c);
        assert!(a.is_odd());
        assert!(!a.is_even());
        assert!(b.is_even());
        assert!(!b.is_odd());

        let mut b = a.clone();
        assert!(a.is_one());
        assert!(b.is_one());
        assert_eq!(a, b);
        assert_ne!(c, b);

        core::mem::swap(&mut b, &mut a);
        // After swap: a=1(old b), b=1(old a), actually same values since both were 1
        // Let's redo with different values
        a = BigInt::from(1);
        b = Ibz::zero();
        let _c = BigInt::from(-1);
        core::mem::swap(&mut a, &mut b);
        assert!(a.is_zero());
        assert!(b.is_one());

        assert_eq!(ibz_bitsize(&BigInt::from(1)), 1);
        assert_eq!(ibz_get(&BigInt::from(1)), 1);
        assert_eq!(ibz_get(&Ibz::zero()), 0);

        let a = BigInt::from(-1);
        assert!(a < Ibz::zero());

        let a = from_str_10("-10000000000000000011111100000001");
        assert!(a < Ibz::zero());
        let b = a.clone();
        assert_eq!(a, b);

        let a = from_str_16("1aaaa00000000000000000123");
        assert!(a > Ibz::zero());
        assert_eq!(ibz_bitsize(&a), 4 * 24 + 1);
        assert_eq!(ibz_get(&a), 16 * 18 + 3);

        let a = from_str_16("deadbeef12345678");
        assert_eq!(ibz_get(&a), 0x12345678);

        // Test INT32_MIN / INT32_MAX
        let a = from_str_10("-2147483648");
        let b = BigInt::from(i32::MIN);
        assert_eq!(ibz_get(&a), -2147483648i32);
        assert_eq!(a, b);

        let a = from_str_10("2147483647");
        let b = BigInt::from(i32::MAX);
        assert_eq!(ibz_get(&a), 2147483647);
        assert_eq!(a, b);
    }

    #[test]
    fn test_add_sub_neg_abs() {
        let a = from_str_16("10000111100002222");
        let b = Ibz::zero();

        // a + 0 == a
        assert_eq!(&a + &b, a);
        assert_eq!(&b + &a, a);

        let a = from_str_16("10000111100002222");
        let b = from_str_16("20000111100002223");
        let c = from_str_16("30000222200004445");
        assert_eq!(&a + &b, c);

        let neg_a = -&a;
        assert_eq!(&neg_a + &(&a + &b), b);
        assert_eq!(-&neg_a, a);

        // sub
        assert_eq!(&b - &(-&a), c);
        assert_eq!(&(&a + &b) - &a, b);

        let d = -&c;
        let neg_a2 = -&a;
        let neg_b = -&b;
        assert_eq!(&d - &neg_a2, neg_b);
        assert_eq!(&neg_a2 - &b, d);

        // abs
        assert_eq!(Ibz::zero().abs(), Ibz::zero());
        assert_eq!(a.abs(), a);
        assert_eq!(neg_a.abs(), a);
    }

    #[test]
    fn test_mul_sqrt() {
        let a = from_str_10("2113309833171849999003363");

        // zero
        assert_eq!(&a * &Ibz::zero(), Ibz::zero());
        assert_eq!(&Ibz::zero() * &a, Ibz::zero());

        // one
        assert_eq!(&a * &BigInt::from(1), a);
        assert_eq!(&BigInt::from(1) * &a, a);

        // -1
        assert_eq!(&a * &BigInt::from(-1), -&a);
        assert_eq!(&BigInt::from(-1) * &a, -&a);

        // larger
        let b = from_str_10("34575345632322576567896");
        let c = from_str_10("73068417910102676801285574959599851857101834248");
        assert_eq!(&a * &b, c);
        assert_eq!(&(-&a) * &(-&b), c);
        assert_eq!(&a * &(-&b), -&c);
        assert_eq!(&(-&a) * &b, -&c);

        // sqrt_floor tests
        let a_abs = a.abs();
        let b_abs = b.abs();
        assert_eq!(ibz_sqrt_floor(&(&a_abs * &a_abs)), a_abs);
        assert_eq!(ibz_sqrt_floor(&(&b_abs * &b_abs)), b_abs);
        let bsq_minus_1 = &b_abs * &b_abs - 1;
        assert_eq!(ibz_sqrt_floor(&bsq_minus_1), &b_abs - 1);
    }

    #[test]
    fn test_div() {
        let a = from_str_10("2113309833171849999003363");

        // Divide by 1
        let (q, r) = ibz_div(&a, &BigInt::from(1));
        assert!(r.is_zero());
        assert_eq!(q, a);
        assert!(ibz_divides(&a, &BigInt::from(1)));

        // Not one, zero remainder
        let b = from_str_10("15678200126527887351125");
        let d = &a * &b;
        assert!(ibz_divides(&d, &b));
        assert!(ibz_divides(&d, &a));
        let (q, r) = ibz_div(&d, &b);
        assert!(r.is_zero());
        assert_eq!(q, a);
        let (q, r) = ibz_div(&d, &a);
        assert!(r.is_zero());
        assert_eq!(q, b);

        // Flipping signs
        let neg_a = -&a;
        let neg_b = -&b;
        let (q, r) = ibz_div(&d, &neg_b);
        assert!(r.is_zero());
        assert_eq!(q, neg_a);
        let (q, r) = ibz_div(&d, &neg_a);
        assert!(r.is_zero());
        assert_eq!(q, neg_b);

        let neg_d = -&d;
        let (q, r) = ibz_div(&neg_d, &neg_b);
        assert!(r.is_zero());
        assert_eq!(q, a);
        let (q, r) = ibz_div(&neg_d, &a);
        assert!(r.is_zero());
        assert_eq!(q, neg_b);

        // Non-zero remainder
        let c = from_str_10("8678205677345432110000");
        let d2 = &d + &c;
        let (q, r) = ibz_div(&d2, &b);
        // Verify: q * b + r == d2
        assert_eq!(&q * &b + &r, d2);
        assert!(r.abs() < b.abs());

        // Verify truncated semantics: remainder has same sign as dividend
        // (or is zero)
        if !r.is_zero() {
            assert!((r.sign() == d2.sign()) || r.is_zero());
        }

        // Flip signs and verify invariant
        let neg_d2 = -&d2;
        let (q, r) = ibz_div(&neg_d2, &b);
        assert_eq!(&q * &b + &r, neg_d2);
        assert!(r.abs() < b.abs());

        let (q, r) = ibz_div(&neg_d2, &neg_b);
        assert_eq!(&q * &neg_b + &r, neg_d2);
        assert!(r.abs() < b.abs());

        let (q, r) = ibz_div(&d2, &neg_b);
        assert_eq!(&q * &neg_b + &r, d2);
        assert!(r.abs() < b.abs());
    }

    #[test]
    fn test_mod() {
        let a = from_str_10("2113309833171849999003363");
        assert_eq!(ibz_mod_ui(&a, 3), 0);
        assert_eq!(ibz_mod_ui(&a, 2), 1);

        let m = from_str_10("2113309833171840000000000");
        let r = ibz_mod(&a, &m);
        assert_eq!(&r + &m, a);

        let neg_a = -&a;
        assert_eq!(ibz_mod_ui(&neg_a, 3), 0);
        assert_eq!(ibz_mod_ui(&neg_a, 2), 1);
        let r = ibz_mod(&neg_a, &m);
        // r - m - m == -a => r = -a + 2m
        assert_eq!(&r - &m - &m, neg_a);
    }

    #[test]
    fn test_pow() {
        let a = from_str_16("aaaaaaaabbbbbbbb2222221111");
        let exp = 10u32;

        let mut manual = BigInt::from(1);
        for _ in 0..exp {
            manual = &manual * &a;
        }
        assert_eq!(ibz_pow(&a, exp), manual);

        // Negative base, even exponent
        let neg_a = -&a;
        assert_eq!(ibz_pow(&neg_a, exp), manual);

        // Odd exponent
        let exp_odd = 9u32;
        let mut manual_odd = BigInt::from(1);
        for _ in 0..exp_odd {
            manual_odd = &manual_odd * &neg_a;
        }
        assert_eq!(ibz_pow(&neg_a, exp_odd), manual_odd);

        // div_2exp
        let shift = 23u32;
        let two_pow = ibz_pow(&BigInt::from(2), shift);
        let (expected, _) = ibz_div(&a, &two_pow);
        assert_eq!(ibz_div_2exp(&a, shift), expected);
    }

    #[test]
    fn test_gcd() {
        let c = from_str_10("25791357069084");
        let a = from_str_10("6173271838293993987767");
        let b = from_str_10("89882267321617266071838286");
        let ac = &a * &c;
        let bc = &b * &c;
        assert_eq!(ibz_gcd(&ac, &bc), c);

        // gcd is always positive
        assert_eq!(ibz_gcd(&(-&ac), &bc), c);
        assert_eq!(ibz_gcd(&(-&ac), &(-&bc)), c);
        assert_eq!(ibz_gcd(&ac, &(-&bc)), c);

        // Different sizes
        assert_eq!(ibz_gcd(&ac, &BigInt::from(2)), BigInt::from(2));
    }

    #[test]
    fn test_invmod() {
        let a = from_str_10("6173271838293993987767");
        let b = from_str_10("89882267321617266071838286");

        let inv = ibz_invmod(&a, &b).unwrap();
        assert!(inv > Ibz::zero());
        assert!(inv < b);
        assert_eq!(ibz_mod(&(&inv * &a), &b), BigInt::from(1));

        // Negative a
        let inv2 = ibz_invmod(&(-&a), &b).unwrap();
        assert!(inv2 > Ibz::zero());
        assert!(inv2 < b);
        assert_eq!(ibz_mod(&(&inv2 * &(-&a)), &b), BigInt::from(1));
    }

    #[test]
    fn test_sqrt_mod_p() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        for prime_n in [67, 103] {
            let two_to_n = BigInt::one() << prime_n;
            let upper = &two_to_n - 1;

            for _ in 0..3 {
                let mut found_p3 = false;
                let mut found_p5 = false;
                let mut found_p1 = false;
                let mut prime_p3 = Ibz::zero();
                let mut prime_p5 = Ibz::zero();
                let mut prime_p1 = Ibz::zero();

                let mut candidate = ibz_rand_interval(&mut rng, &Ibz::zero(), &upper);
                if candidate.is_even() {
                    candidate += 1;
                }

                while !found_p3 || !found_p5 || !found_p1 {
                    candidate += 2;
                    if ibz_probab_prime(&candidate, 25) == 0 {
                        continue;
                    }
                    let m4 = ibz_mod_ui(&candidate, 4);
                    let m8 = ibz_mod_ui(&candidate, 8);
                    if m4 == 3 && !found_p3 {
                        prime_p3 = candidate.clone();
                        found_p3 = true;
                    } else if m8 == 5 && !found_p5 {
                        prime_p5 = candidate.clone();
                        found_p5 = true;
                    } else if m8 == 1 && !found_p1 {
                        prime_p1 = candidate.clone();
                        found_p1 = true;
                    }
                }

                for p in [&prime_p3, &prime_p5, &prime_p1] {
                    let p_minus_1 = p - 1;
                    let a = ibz_rand_interval(&mut rng, &Ibz::zero(), &p_minus_1);
                    let asq = ibz_mod(&(&a * &a), p);
                    let prime_minus_a = p - &a;

                    let sqrt = ibz_sqrt_mod_p(&asq, p).expect("sqrt should exist for a square");
                    assert!(
                        sqrt == a || sqrt == prime_minus_a,
                        "sqrt^2 should equal a^2 mod p"
                    );
                }
            }
        }
    }

    #[test]
    fn test_two_adic() {
        assert_eq!(ibz_two_adic(&BigInt::from(1)), 0);
        assert_eq!(ibz_two_adic(&BigInt::from(2)), 1);
        assert_eq!(ibz_two_adic(&BigInt::from(4)), 2);
        assert_eq!(ibz_two_adic(&BigInt::from(8)), 3);
        assert_eq!(ibz_two_adic(&BigInt::from(12)), 2); // 1100 binary
        assert_eq!(ibz_two_adic(&BigInt::from(0)), 0);
    }

    #[test]
    fn test_rand_interval() {
        let mut rng = rand::thread_rng();
        let low = from_str_16("ffa");
        let high = from_str_16("eeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef");

        for _ in 0..10 {
            let r = ibz_rand_interval(&mut rng, &low, &high);
            assert!(r >= low);
            assert!(r <= high);
        }
    }

    #[test]
    fn test_rand_interval_i() {
        let mut rng = rand::thread_rng();
        for _ in 0..20 {
            let r = ibz_rand_interval_i(&mut rng, 10, 100);
            assert!(r >= BigInt::from(10));
            assert!(r <= BigInt::from(100));
        }
    }

    #[test]
    fn test_rand_interval_minm_m() {
        let mut rng = rand::thread_rng();
        let m = 1000i32;
        for _ in 0..20 {
            let r = ibz_rand_interval_minm_m(&mut rng, m);
            assert!(r >= BigInt::from(-m));
            assert!(r <= BigInt::from(m));
        }
    }

    #[test]
    fn test_copy_digits() {
        let d1: &[u64] = &[0x12345678];
        let val = ibz_copy_digits(d1);
        assert_eq!(ibz_to_str(&val, 16), "12345678");

        let d2: &[u64] = &[2, 1];
        let val = ibz_copy_digits(d2);
        assert_eq!(ibz_to_str(&val, 16), "10000000000000002");
    }

    #[test]
    fn test_to_digits() {
        let val = from_str_16("12345678");
        let mut digits = [0u64; 1];
        ibz_to_digits(&val, &mut digits);
        assert_eq!(digits[0], 0x12345678);

        let val = from_str_16("10000000000000002");
        let mut digits = [0u64; 2];
        ibz_to_digits(&val, &mut digits);
        assert_eq!(digits[0], 2);
        assert_eq!(digits[1], 1);

        // Round-trip
        let val = from_str_10("1617406613339667622221321");
        let nbits = ibz_bitsize(&val) as usize;
        let ndigits = nbits.div_ceil(64);
        let mut digits = vec![0u64; ndigits];
        ibz_to_digits(&val, &mut digits);
        let recovered = ibz_copy_digits(&digits);
        assert_eq!(val, recovered);

        // Zero
        let mut digits = [0u64; 1];
        ibz_to_digits(&Ibz::zero(), &mut digits);
        assert_eq!(digits[0], 0);
    }

    #[test]
    fn test_probab_prime() {
        assert!(ibz_probab_prime(&BigInt::from(2), 25) > 0);
        assert!(ibz_probab_prime(&BigInt::from(3), 25) > 0);
        assert!(ibz_probab_prime(&BigInt::from(5), 25) > 0);
        assert!(ibz_probab_prime(&BigInt::from(7), 25) > 0);
        assert_eq!(ibz_probab_prime(&BigInt::from(4), 25), 0);
        assert_eq!(ibz_probab_prime(&BigInt::from(9), 25), 0);
        assert_eq!(ibz_probab_prime(&BigInt::from(1), 25), 0);

        // Large prime
        let p = from_str_10("170141183460469231731687303715884105727"); // 2^127 - 1 (Mersenne prime)
        assert!(ibz_probab_prime(&p, 25) > 0);

        // Large composite
        let c = &p * &BigInt::from(3);
        assert_eq!(ibz_probab_prime(&c, 25), 0);
    }

    #[test]
    fn test_legendre() {
        let p = BigInt::from(7);
        // 1 is QR mod 7
        assert_eq!(ibz_legendre(&BigInt::from(1), &p), 1);
        // 2 is QR mod 7 (3^2 = 2 mod 7)
        assert_eq!(ibz_legendre(&BigInt::from(2), &p), 1);
        // 3 is not QR mod 7
        assert_eq!(ibz_legendre(&BigInt::from(3), &p), -1);
        // 0 mod p
        assert_eq!(ibz_legendre(&BigInt::from(7), &p), 0);
    }

    #[test]
    fn test_divides() {
        assert!(ibz_divides(&BigInt::from(12), &BigInt::from(3)));
        assert!(ibz_divides(&BigInt::from(12), &BigInt::from(4)));
        assert!(!ibz_divides(&BigInt::from(12), &BigInt::from(5)));
        assert!(ibz_divides(&BigInt::from(0), &BigInt::from(5)));
    }

    /// RNG wrapper that counts bytes consumed and fills with a fixed value.
    struct ByteCountingRng {
        count: usize,
        fill_value: u8,
    }

    impl ByteCountingRng {
        fn new(fill_value: u8) -> Self {
            Self {
                count: 0,
                fill_value,
            }
        }
    }

    impl rand::RngCore for ByteCountingRng {
        fn next_u32(&mut self) -> u32 {
            let mut buf = [0u8; 4];
            self.fill_bytes(&mut buf);
            u32::from_le_bytes(buf)
        }
        fn next_u64(&mut self) -> u64 {
            let mut buf = [0u8; 8];
            self.fill_bytes(&mut buf);
            u64::from_le_bytes(buf)
        }
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            self.count += dest.len();
            for b in dest.iter_mut() {
                *b = self.fill_value;
            }
        }
        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
            self.fill_bytes(dest);
            Ok(())
        }
    }

    /// Bug 1: ibz_rand_interval must consume exactly ceil(bitsize/8) bytes
    /// per rejection-sampling attempt.
    #[test]
    fn regression_rand_interval_byte_consumption() {
        // fill_value=0 → r=0 which is always <= b-a, so accepted on first try
        let mut rng = ByteCountingRng::new(0);

        // a=0, b=2^255-1: bmina has 255 bits → ceil(255/8) = 32 bytes
        let a = Ibz::zero();
        let b = (BigInt::one() << 255) - 1;
        let _ = ibz_rand_interval(&mut rng, &a, &b);
        assert_eq!(rng.count, 32, "255-bit range should consume 32 bytes");

        // Reset and test 256-bit range
        rng.count = 0;
        let b = (BigInt::one() << 256) - 1;
        let _ = ibz_rand_interval(&mut rng, &a, &b);
        assert_eq!(rng.count, 32, "256-bit range should consume 32 bytes");

        // Reset and test 257-bit range (2^256, which has 257 bits)
        rng.count = 0;
        let b = BigInt::one() << 256;
        let _ = ibz_rand_interval(&mut rng, &a, &b);
        assert_eq!(rng.count, 33, "257-bit range should consume 33 bytes");

        // Reset and test small range: 64 bits
        rng.count = 0;
        let b = (BigInt::one() << 64) - 1;
        let _ = ibz_rand_interval(&mut rng, &a, &b);
        assert_eq!(rng.count, 8, "64-bit range should consume 8 bytes");

        // Reset and test boundary: 65 bits
        rng.count = 0;
        let b = BigInt::one() << 64;
        let _ = ibz_rand_interval(&mut rng, &a, &b);
        assert_eq!(rng.count, 9, "65-bit range should consume 9 bytes");

        // Single byte: range with 1 bit
        rng.count = 0;
        let _ = ibz_rand_interval(&mut rng, &Ibz::zero(), &BigInt::from(1));
        assert_eq!(rng.count, 1, "1-bit range should consume 1 byte");
    }

    /// ibz_probab_prime is deterministic (does not consume RNG).
    #[test]
    fn regression_probab_prime_no_rng_consumption() {
        // Small primes
        assert!(ibz_probab_prime(&BigInt::from(2), 32) > 0);
        assert!(ibz_probab_prime(&BigInt::from(3), 32) > 0);
        assert!(ibz_probab_prime(&BigInt::from(97), 32) > 0);

        // Small composites
        assert_eq!(ibz_probab_prime(&BigInt::from(4), 32), 0);
        assert_eq!(ibz_probab_prime(&BigInt::from(100), 32), 0);

        // Large prime (Mersenne prime 2^127 - 1)
        let p = from_str_10("170141183460469231731687303715884105727");
        assert!(ibz_probab_prime(&p, 32) > 0);

        // Large composite
        let c = &p * &BigInt::from(7);
        assert_eq!(ibz_probab_prime(&c, 32), 0);

        // 256-bit prime (secp256k1 order)
        let p256 = BigInt::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap();
        assert!(ibz_probab_prime(&p256, 32) > 0);
    }

    /// ibz_get must return the low 32 bits (sign + low 31 bits) for
    /// values of any magnitude.
    #[test]
    fn regression_ibz_get_large_values() {
        // Boundary: i32::MAX
        assert_eq!(ibz_get(&BigInt::from(i32::MAX)), i32::MAX);
        // Boundary: i32::MIN
        assert_eq!(ibz_get(&BigInt::from(i32::MIN)), i32::MIN);
        // Zero
        assert_eq!(ibz_get(&Ibz::zero()), 0);

        // Large positive with small low bits
        let v = (BigInt::one() << 128) + BigInt::from(42);
        assert_eq!(ibz_get(&v), 42);

        // Large positive with zero low bits
        let v = BigInt::one() << 128;
        assert_eq!(ibz_get(&v), 0);

        // Large positive with bit 31 set in low word
        let v = (BigInt::one() << 128) + BigInt::from(0xDEADBEEFu32);
        assert_eq!(ibz_get(&v), 0xDEADBEEFu32 as i32);

        // Large negative with small low bits
        let big128: BigInt = BigInt::one() << 128;
        let v = -(big128.clone() + BigInt::from(123));
        assert_eq!(ibz_get(&v), -123);

        // Large negative with bit 31 set (0x80000000 in low word → i32::MIN)
        let v = -(big128 + BigInt::from(0x80000000u32));
        assert_eq!(ibz_get(&v), i32::MIN);

        // Value just above u32::MAX: 2^32 + 1
        // low 64-bit limb = 0x1_0000_0001, low 32 bits = 1
        let v = BigInt::from(0x1_0000_0001u64);
        assert_eq!(ibz_get(&v), 1);

        // Exactly 2^31 (positive, bit 31 set)
        let v = BigInt::from(0x80000000u32);
        assert_eq!(ibz_get(&v), i32::MIN); // 0x80000000 as i32

        // Exactly 2^32 (positive, low 32 bits = 0)
        let v = BigInt::from(0x100000000u64);
        assert_eq!(ibz_get(&v), 0);

        // Multi-limb: 2^200 + 0x7FFF_FFFF (max positive from low 31 bits)
        let v = (BigInt::one() << 200) + BigInt::from(0x7FFFFFFFu32);
        assert_eq!(ibz_get(&v), 0x7FFFFFFF);
    }
}
