//!
//! A rational number is represented as a pair `(numerator, denominator)` of
//! big integers. Operations used by the LLL verification code for exact
//! size-reduction checks.

use super::intbig::{ibz_div, ibz_gcd, ibz_mod, Ibz};
use num_traits::{One, Zero};
use std::ops::{Index, IndexMut};

/// Rational number: `numerator / denominator`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ibq {
    pub num: Ibz,
    pub den: Ibz,
}

impl Default for Ibq {
    fn default() -> Self {
        Self {
            num: Ibz::zero(),
            den: Ibz::one(),
        }
    }
}

/// 4-element vector of rationals.
#[derive(Clone, Debug, Default)]
pub struct IbqVec4(pub [Ibq; 4]);

impl Index<usize> for IbqVec4 {
    type Output = Ibq;
    fn index(&self, i: usize) -> &Ibq {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbqVec4 {
    fn index_mut(&mut self, i: usize) -> &mut Ibq {
        &mut self.0[i]
    }
}

/// 4x4 matrix of rationals.
#[derive(Clone, Debug, Default)]
pub struct IbqMat4x4(pub [IbqVec4; 4]);

impl Index<usize> for IbqMat4x4 {
    type Output = IbqVec4;
    fn index(&self, i: usize) -> &IbqVec4 {
        &self.0[i]
    }
}

impl IndexMut<usize> for IbqMat4x4 {
    fn index_mut(&mut self, i: usize) -> &mut IbqVec4 {
        &mut self.0[i]
    }
}

/// Reduce a rational by dividing numerator and denominator by their GCD.
pub fn ibq_reduce(x: &mut Ibq) {
    let g = ibz_gcd(&x.num, &x.den);
    if !g.is_zero() && !g.is_one() {
        let (q1, r1) = ibz_div(&x.num, &g);
        debug_assert!(r1.is_zero());
        x.num = q1;
        let (q2, r2) = ibz_div(&x.den, &g);
        debug_assert!(r2.is_zero());
        x.den = q2;
    }
}

/// `sum = a + b` (unreduced).
pub fn ibq_add(a: &Ibq, b: &Ibq) -> Ibq {
    // a.num/a.den + b.num/b.den = (a.num*b.den + b.num*a.den) / (a.den*b.den)
    let num = &a.num * &b.den + &b.num * &a.den;
    let den = &a.den * &b.den;
    Ibq { num, den }
}

/// `neg = -x`.
pub fn ibq_neg(x: &Ibq) -> Ibq {
    Ibq {
        num: -&x.num,
        den: x.den.clone(),
    }
}

/// `diff = a - b`.
pub fn ibq_sub(a: &Ibq, b: &Ibq) -> Ibq {
    let neg_b = ibq_neg(b);
    ibq_add(a, &neg_b)
}

/// `abs = |x|`.
pub fn ibq_abs(x: &Ibq) -> Ibq {
    let neg = ibq_neg(x);
    if ibq_cmp(x, &neg) < 0 {
        neg
    } else {
        x.clone()
    }
}

/// `prod = a * b` (unreduced).
pub fn ibq_mul(a: &Ibq, b: &Ibq) -> Ibq {
    Ibq {
        num: &a.num * &b.num,
        den: &a.den * &b.den,
    }
}

/// `inv = 1/x`. Returns `None` if `x` is zero.
pub fn ibq_inv(x: &Ibq) -> Option<Ibq> {
    if ibq_is_zero(x) {
        None
    } else {
        Some(Ibq {
            num: x.den.clone(),
            den: x.num.clone(),
        })
    }
}

/// Compare `a` and `b`. Returns positive if `a > b`, negative if `a < b`, 0 if equal.
pub fn ibq_cmp(a: &Ibq, b: &Ibq) -> i32 {
    // a.num/a.den vs b.num/b.den
    // Cross-multiply: a.num*b.den vs b.num*a.den
    // But must account for signs of denominators
    let mut x = &a.num * &b.den;
    let mut y = &b.num * &a.den;

    // If a.den > 0, the inequality direction is preserved; if < 0, reversed.
    // Negate both cross-products when den > 0 (flips comparison direction):
    if a.den > Ibz::zero() {
        x = -x;
        y = -y;
    }
    if b.den > Ibz::zero() {
        x = -x;
        y = -y;
    }

    if x < y {
        -1
    } else if x > y {
        1
    } else {
        0
    }
}

/// Test if `x` is zero.
pub fn ibq_is_zero(x: &Ibq) -> bool {
    x.num.is_zero()
}

/// Test if `x` is one.
pub fn ibq_is_one(x: &Ibq) -> bool {
    x.num == x.den
}

/// Set `q = a/b`. Returns false if `b` is zero.
pub fn ibq_set(a: &Ibz, b: &Ibz) -> Option<Ibq> {
    if b.is_zero() {
        None
    } else {
        Some(Ibq {
            num: a.clone(),
            den: b.clone(),
        })
    }
}

/// Copy a rational.
pub fn ibq_copy(value: &Ibq) -> Ibq {
    value.clone()
}

/// Check if `q` is an integer (denominator divides numerator).
pub fn ibq_is_ibz(q: &Ibq) -> bool {
    let r = ibz_mod(&q.num, &q.den);
    r.is_zero()
}

/// Convert rational to integer if it is one. Returns `None` if not an integer.
pub fn ibq_to_ibz(q: &Ibq) -> Option<Ibz> {
    let (quotient, rem) = ibz_div(&q.num, &q.den);
    if rem.is_zero() {
        Some(quotient)
    } else {
        None
    }
}

/// Set an `IbqVec4` from four integer values.
pub fn ibq_vec_4_copy_ibz(c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) -> IbqVec4 {
    let one = Ibz::one();
    IbqVec4([
        Ibq {
            num: c0.clone(),
            den: one.clone(),
        },
        Ibq {
            num: c1.clone(),
            den: one.clone(),
        },
        Ibq {
            num: c2.clone(),
            den: one.clone(),
        },
        Ibq {
            num: c3.clone(),
            den: one,
        },
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn test_ibq_basics() {
        let q = ibq_set(&BigInt::from(123), &BigInt::from(-123)).unwrap();
        assert!(!ibq_is_one(&q));
        assert!(ibq_is_ibz(&q));
        let z = ibq_to_ibz(&q).unwrap();
        assert_eq!(z, BigInt::from(-1));
    }

    #[test]
    fn test_ibq_is_one() {
        let q = ibq_set(&BigInt::from(123), &BigInt::from(123)).unwrap();
        assert!(ibq_is_one(&q));
    }

    #[test]
    fn test_ibq_is_zero() {
        let q = ibq_set(&BigInt::from(0), &BigInt::from(123)).unwrap();
        assert!(ibq_is_zero(&q));
        assert!(ibq_is_ibz(&q));
        let z = ibq_to_ibz(&q).unwrap();
        assert!(z.is_zero());
    }

    #[test]
    fn test_ibq_add_mul() {
        let a = ibq_set(&BigInt::from(1), &BigInt::from(3)).unwrap();
        let b = ibq_set(&BigInt::from(1), &BigInt::from(6)).unwrap();
        let sum = ibq_add(&a, &b);
        // 1/3 + 1/6 = 3/18 = 1/2
        let half = ibq_set(&BigInt::from(1), &BigInt::from(2)).unwrap();
        assert_eq!(ibq_cmp(&sum, &half), 0);
    }

    #[test]
    fn test_ibq_cmp() {
        let a = ibq_set(&BigInt::from(1), &BigInt::from(3)).unwrap();
        let b = ibq_set(&BigInt::from(1), &BigInt::from(2)).unwrap();
        assert!(ibq_cmp(&a, &b) < 0);
        assert!(ibq_cmp(&b, &a) > 0);
        assert_eq!(ibq_cmp(&a, &a), 0);
    }

    #[test]
    fn test_ibq_inv() {
        let a = ibq_set(&BigInt::from(3), &BigInt::from(7)).unwrap();
        let inv = ibq_inv(&a).unwrap();
        let prod = ibq_mul(&a, &inv);
        assert!(ibq_is_one(&prod));
    }

    #[test]
    fn test_ibq_zero_inv() {
        let zero = Ibq::default();
        assert!(ibq_inv(&zero).is_none());
    }

    #[test]
    fn test_ibq_vec_4_copy_ibz() {
        let v = ibq_vec_4_copy_ibz(
            &BigInt::from(2),
            &BigInt::from(3),
            &BigInt::from(4),
            &BigInt::from(5),
        );
        for i in 0..4 {
            let z = ibq_to_ibz(&v[i]).unwrap();
            assert_eq!(z, BigInt::from(i as i32 + 2));
        }
    }
}
