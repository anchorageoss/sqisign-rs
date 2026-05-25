//!
//! A DPE value represents `mantissa × 2^exponent` where the mantissa is an
//! `f64` normalized to `[0.5, 1.0)` (or zero). This extends the dynamic
//! range of `f64` far beyond its native 11-bit exponent while preserving
//! 53 bits of significand precision.
//!
//! Used exclusively by the LLL lattice reduction algorithm for approximate
//! Gram-Schmidt orthogonalization.

use super::intbig::Ibz;
use num_traits::{ToPrimitive, Zero};

/// Number of significant bits in the mantissa (IEEE 754 double precision).
const DPE_BITSIZE: i32 = 53;
/// Sentinel exponent representing zero or underflow.
const DPE_EXPMIN: i32 = i32::MIN;
/// `2^53`: scale factor for integer-mantissa ↔ float-mantissa conversion.
const DPE_2_POW_BITSIZE: f64 = 9007199254740992.0;

/// Double-Precision Extended floating-point value.
///
/// Represents `mantissa × 2^exponent` with the invariant that
/// `mantissa` is in `[0.5, 1.0)` after normalization (or zero).
#[derive(Clone, Copy, Debug)]
pub struct Dpe {
    pub mantissa: f64,
    pub exponent: i32,
}

impl Default for Dpe {
    fn default() -> Self {
        Self {
            mantissa: 0.0,
            exponent: DPE_EXPMIN,
        }
    }
}

fn dpe_normalize(x: &mut Dpe) {
    if x.mantissa == 0.0 || !x.mantissa.is_finite() {
        if x.mantissa == 0.0 {
            x.exponent = DPE_EXPMIN;
        }
    } else {
        let (m, e) = frexp(x.mantissa);
        x.mantissa = m;
        x.exponent = x.exponent.wrapping_add(e);
    }
}

/// Decompose `x` into `(m, e)` such that `x = m × 2^e` with `0.5 <= |m| < 1`.
fn frexp(x: f64) -> (f64, i32) {
    if x == 0.0 || !x.is_finite() {
        return (x, 0);
    }
    let bits = x.to_bits();
    let exp_bits = ((bits >> 52) & 0x7FF) as i32;
    let e = exp_bits - 1022;
    // Set exponent to -1 (biased: 1022) to get mantissa in [0.5, 1.0)
    let new_bits = (bits & 0x800F_FFFF_FFFF_FFFF) | (1022u64 << 52);
    (f64::from_bits(new_bits), e)
}

/// Reconstruct `m × 2^e`.
fn ldexp(m: f64, e: i32) -> f64 {
    m * (2.0f64).powi(e)
}

fn dpe_scale(d: f64, s: i32) -> f64 {
    // s is in range (-DPE_BITSIZE, 0] and 0.5 <= |d| < 1
    ldexp(d, s)
}

pub fn dpe_set(y: &Dpe) -> Dpe {
    *y
}

pub fn dpe_set_d(y: f64) -> Dpe {
    let mut x = Dpe {
        mantissa: y,
        exponent: 0,
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_set_si(y: i64) -> Dpe {
    let mut x = Dpe {
        mantissa: y as f64,
        exponent: 0,
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_set_ui(y: u64) -> Dpe {
    let mut x = Dpe {
        mantissa: y as f64,
        exponent: 0,
    };
    dpe_normalize(&mut x);
    x
}

/// Set from a big integer, extracting mantissa and exponent.
///
/// Extracts the top ~53 significant bits as an f64 mantissa in [0.5, 1.0)
/// and the exponent separately, so arbitrarily large integers are handled
/// without overflow.
pub fn dpe_set_z(y: &Ibz) -> Dpe {
    if y.is_zero() {
        return Dpe {
            mantissa: 0.0,
            exponent: DPE_EXPMIN,
        };
    }
    let bits = super::intbig::ibz_bitsize(y) as i32;
    if bits <= DPE_BITSIZE {
        let d = y.to_f64().unwrap_or(0.0);
        let (m, e) = frexp(d);
        return Dpe {
            mantissa: m,
            exponent: e,
        };
    }
    let sign = if y < &Ibz::zero() { -1i8 } else { 1 };
    let abs_y = if sign < 0 { -y } else { y.clone() };
    let shift = bits - DPE_BITSIZE;
    let top_bits = &abs_y >> (shift as usize);
    let d = top_bits.to_f64().unwrap_or(0.0);
    let (m, e) = frexp(d);
    Dpe {
        mantissa: if sign < 0 { -m } else { m },
        exponent: e + shift,
    }
}

/// Convert DPE to big integer, rounded to nearest.
pub fn dpe_get_z(y: &Dpe) -> Ibz {
    let ey = y.exponent;
    if ey >= DPE_BITSIZE {
        // y is an integer
        let d = y.mantissa * DPE_2_POW_BITSIZE;
        let d_int = d as i64;
        let base = Ibz::from(d_int);
        let shift = (ey - DPE_BITSIZE) as usize;
        base << shift
    } else if ey < 0 {
        // |y| < 1/2
        Ibz::zero()
    } else {
        let d = ldexp(y.mantissa, ey);
        Ibz::from(d.round() as i64)
    }
}

pub fn dpe_get_d(x: &Dpe) -> f64 {
    ldexp(x.mantissa, x.exponent)
}

pub fn dpe_get_si(x: &Dpe) -> i64 {
    ldexp(x.mantissa, x.exponent) as i64
}

pub fn dpe_neg(y: &Dpe) -> Dpe {
    Dpe {
        mantissa: -y.mantissa,
        exponent: y.exponent,
    }
}

pub fn dpe_abs(y: &Dpe) -> Dpe {
    Dpe {
        mantissa: y.mantissa.abs(),
        exponent: y.exponent,
    }
}

pub fn dpe_zero_p(x: &Dpe) -> bool {
    x.mantissa == 0.0
}

pub fn dpe_sign(x: &Dpe) -> i32 {
    if x.mantissa < 0.0 {
        -1
    } else if x.mantissa > 0.0 {
        1
    } else {
        0
    }
}

pub fn dpe_add(y: &Dpe, z: &Dpe) -> Dpe {
    if (y.exponent as i64) > (z.exponent as i64) + (DPE_BITSIZE as i64) {
        return *y;
    } else if (z.exponent as i64) > (y.exponent as i64) + (DPE_BITSIZE as i64) {
        return *z;
    }

    let d = y.exponent.wrapping_sub(z.exponent);
    let mut x = if d >= 0 {
        Dpe {
            mantissa: y.mantissa + dpe_scale(z.mantissa, -d),
            exponent: y.exponent,
        }
    } else {
        Dpe {
            mantissa: z.mantissa + dpe_scale(y.mantissa, d),
            exponent: z.exponent,
        }
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_sub(y: &Dpe, z: &Dpe) -> Dpe {
    if (y.exponent as i64) > (z.exponent as i64) + (DPE_BITSIZE as i64) {
        return *y;
    } else if (z.exponent as i64) > (y.exponent as i64) + (DPE_BITSIZE as i64) {
        return dpe_neg(z);
    }

    let d = y.exponent.wrapping_sub(z.exponent);
    let mut x = if d >= 0 {
        Dpe {
            mantissa: y.mantissa - dpe_scale(z.mantissa, -d),
            exponent: y.exponent,
        }
    } else {
        Dpe {
            mantissa: dpe_scale(y.mantissa, d) - z.mantissa,
            exponent: z.exponent,
        }
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_mul(y: &Dpe, z: &Dpe) -> Dpe {
    let mut x = Dpe {
        mantissa: y.mantissa * z.mantissa,
        exponent: y.exponent.wrapping_add(z.exponent),
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_div(y: &Dpe, z: &Dpe) -> Dpe {
    let mut x = Dpe {
        mantissa: y.mantissa / z.mantissa,
        exponent: y.exponent.wrapping_sub(z.exponent),
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_sqrt(y: &Dpe) -> Dpe {
    let ey = y.exponent;
    let mut x = if ey % 2 != 0 {
        Dpe {
            mantissa: (0.5 * y.mantissa).sqrt(),
            exponent: (ey + 1) / 2,
        }
    } else {
        Dpe {
            mantissa: y.mantissa.sqrt(),
            exponent: ey / 2,
        }
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_mul_ui(y: &Dpe, z: u64) -> Dpe {
    let mut x = Dpe {
        mantissa: y.mantissa * (z as f64),
        exponent: y.exponent,
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_div_ui(y: &Dpe, z: u64) -> Dpe {
    let mut x = Dpe {
        mantissa: y.mantissa / (z as f64),
        exponent: y.exponent,
    };
    dpe_normalize(&mut x);
    x
}

pub fn dpe_mul_2exp(y: &Dpe, e: u64) -> Dpe {
    Dpe {
        mantissa: y.mantissa,
        exponent: y.exponent.wrapping_add(e as i32),
    }
}

pub fn dpe_div_2exp(y: &Dpe, e: u64) -> Dpe {
    Dpe {
        mantissa: y.mantissa,
        exponent: y.exponent.wrapping_sub(e as i32),
    }
}

pub fn dpe_cmp(x: &Dpe, y: &Dpe) -> i32 {
    let sx = dpe_sign(x);
    let d = sx - dpe_sign(y);
    if d != 0 {
        return d;
    }
    if x.exponent > y.exponent {
        return if sx > 0 { 1 } else { -1 };
    } else if y.exponent > x.exponent {
        return if sx > 0 { -1 } else { 1 };
    }
    if x.mantissa < y.mantissa {
        -1
    } else if x.mantissa > y.mantissa {
        1
    } else {
        0
    }
}

pub fn dpe_cmp_d(x: &Dpe, d: f64) -> i32 {
    let y = dpe_set_d(d);
    dpe_cmp(x, &y)
}

pub fn dpe_cmp_ui(x: &Dpe, d: u64) -> i32 {
    let y = dpe_set_ui(d);
    dpe_cmp(x, &y)
}

pub fn dpe_cmp_si(x: &Dpe, d: i64) -> i32 {
    let y = dpe_set_si(d);
    dpe_cmp(x, &y)
}

pub fn dpe_round(y: &Dpe) -> Dpe {
    if y.exponent < 0 {
        dpe_set_ui(0)
    } else if y.exponent >= DPE_BITSIZE {
        *y
    } else {
        let d = ldexp(y.mantissa, y.exponent);
        dpe_set_d(d.round())
    }
}

pub fn dpe_floor(y: &Dpe) -> Dpe {
    if y.exponent <= 0 {
        if dpe_sign(y) >= 0 {
            dpe_set_ui(0)
        } else {
            dpe_set_si(-1)
        }
    } else if y.exponent >= DPE_BITSIZE {
        *y
    } else {
        let d = ldexp(y.mantissa, y.exponent);
        dpe_set_d(d.floor())
    }
}

pub fn dpe_ceil(y: &Dpe) -> Dpe {
    if y.exponent <= 0 {
        if dpe_sign(y) > 0 {
            dpe_set_ui(1)
        } else {
            dpe_set_si(0)
        }
    } else if y.exponent >= DPE_BITSIZE {
        *y
    } else {
        let d = ldexp(y.mantissa, y.exponent);
        dpe_set_d(d.ceil())
    }
}

pub fn dpe_frac(y: &Dpe) -> Dpe {
    if y.exponent <= 0 {
        *y
    } else if y.exponent >= DPE_BITSIZE {
        dpe_set_ui(0)
    } else {
        let d = ldexp(y.mantissa, y.exponent);
        dpe_set_d(d - d.trunc())
    }
}

pub fn dpe_swap(x: &mut Dpe, y: &mut Dpe) {
    std::mem::swap(x, y);
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn test_dpe_set_d_and_get() {
        let x = dpe_set_d(1.0);
        assert_eq!(dpe_get_d(&x), 1.0);
        assert!(x.mantissa >= 0.5 && x.mantissa < 1.0);

        let x = dpe_set_d(0.0);
        assert_eq!(dpe_get_d(&x), 0.0);
        assert!(dpe_zero_p(&x));

        let x = dpe_set_d(-3.5);
        assert!((dpe_get_d(&x) - (-3.5)).abs() < 1e-15);
    }

    #[test]
    fn test_dpe_arithmetic() {
        let a = dpe_set_d(2.5);
        let b = dpe_set_d(1.5);

        let sum = dpe_add(&a, &b);
        assert!((dpe_get_d(&sum) - 4.0).abs() < 1e-15);

        let diff = dpe_sub(&a, &b);
        assert!((dpe_get_d(&diff) - 1.0).abs() < 1e-15);

        let prod = dpe_mul(&a, &b);
        assert!((dpe_get_d(&prod) - 3.75).abs() < 1e-15);

        let quot = dpe_div(&a, &b);
        assert!((dpe_get_d(&quot) - (2.5 / 1.5)).abs() < 1e-15);
    }

    #[test]
    fn test_dpe_cmp() {
        let a = dpe_set_d(2.5);
        let b = dpe_set_d(1.5);
        assert!(dpe_cmp(&a, &b) > 0);
        assert!(dpe_cmp(&b, &a) < 0);
        assert_eq!(dpe_cmp(&a, &a), 0);

        assert!(dpe_cmp_d(&a, 2.0) > 0);
        assert!(dpe_cmp_d(&a, 3.0) < 0);
    }

    #[test]
    fn test_dpe_round() {
        let x = dpe_set_d(2.7);
        let r = dpe_round(&x);
        assert!((dpe_get_d(&r) - 3.0).abs() < 1e-15);

        let x = dpe_set_d(2.3);
        let r = dpe_round(&x);
        assert!((dpe_get_d(&r) - 2.0).abs() < 1e-15);

        let x = dpe_set_d(-0.3);
        let r = dpe_round(&x);
        assert!((dpe_get_d(&r) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_dpe_set_z_get_z() {
        let z = BigInt::from(42);
        let x = dpe_set_z(&z);
        let back = dpe_get_z(&x);
        assert_eq!(back, BigInt::from(42));

        let z = BigInt::from(-100);
        let x = dpe_set_z(&z);
        let back = dpe_get_z(&x);
        assert_eq!(back, BigInt::from(-100));

        let z = BigInt::from(0);
        let x = dpe_set_z(&z);
        assert!(dpe_zero_p(&x));
    }

    #[test]
    fn test_dpe_large_z() {
        let z = BigInt::from(1) << 100;
        let x = dpe_set_z(&z);
        assert!(x.exponent > 100);
        let back = dpe_get_z(&x);
        // Should be approximately equal (within f64 precision)
        let diff = &back - &z;
        let ratio = diff.to_f64().unwrap_or(1.0) / z.to_f64().unwrap_or(1.0);
        assert!(ratio.abs() < 1e-10);
    }

    #[test]
    fn test_dpe_neg_abs() {
        let x = dpe_set_d(3.0);
        let n = dpe_neg(&x);
        assert!((dpe_get_d(&n) - (-3.0)).abs() < 1e-15);

        let a = dpe_abs(&n);
        assert!((dpe_get_d(&a) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_dpe_sqrt() {
        let x = dpe_set_d(4.0);
        let s = dpe_sqrt(&x);
        assert!((dpe_get_d(&s) - 2.0).abs() < 1e-15);

        let x = dpe_set_d(2.0);
        let s = dpe_sqrt(&x);
        assert!((dpe_get_d(&s) - std::f64::consts::SQRT_2).abs() < 1e-15);
    }
}
