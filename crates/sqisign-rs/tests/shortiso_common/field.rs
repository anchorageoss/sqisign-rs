//! `F_p` and `F_{p^2} = F_p[i] / (i^2 + 1)` over an arbitrary prime
//! `p = 3 mod 4`, on `num-bigint`. Test-only: no performance work, nothing
//! constant time.

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};

/// An element `re + im i` of `F_{p^2}`, both parts reduced.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct F2 {
    pub re: BigUint,
    pub im: BigUint,
}

/// The field `F_{p^2}`, `p = 3 mod 4`: every operation goes through it.
#[derive(Clone, Debug)]
pub struct Field {
    pub p: BigUint,
}

impl Field {
    pub fn new(p: &BigUint) -> Self {
        assert_eq!(p % 4u32, BigUint::from(3u32), "p = 3 mod 4");
        Self { p: p.clone() }
    }

    pub fn zero(&self) -> F2 {
        F2 {
            re: BigUint::zero(),
            im: BigUint::zero(),
        }
    }

    pub fn one(&self) -> F2 {
        F2 {
            re: BigUint::one(),
            im: BigUint::zero(),
        }
    }

    /// The image of an integer in `F_p`.
    pub fn int(&self, v: &BigInt) -> F2 {
        let p = BigInt::from_biguint(Sign::Plus, self.p.clone());
        let r = v.mod_floor(&p);
        F2 {
            re: r.to_biguint().expect("non-negative residue"),
            im: BigUint::zero(),
        }
    }

    pub fn small(&self, v: i64) -> F2 {
        self.int(&BigInt::from(v))
    }

    pub fn elem(&self, re: &BigUint, im: &BigUint) -> F2 {
        F2 {
            re: re % &self.p,
            im: im % &self.p,
        }
    }

    fn addp(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let s = a + b;
        if s >= self.p {
            s - &self.p
        } else {
            s
        }
    }

    fn subp(&self, a: &BigUint, b: &BigUint) -> BigUint {
        if a >= b {
            a - b
        } else {
            &self.p - b + a
        }
    }

    fn negp(&self, a: &BigUint) -> BigUint {
        if a.is_zero() {
            BigUint::zero()
        } else {
            &self.p - a
        }
    }

    pub fn add(&self, a: &F2, b: &F2) -> F2 {
        F2 {
            re: self.addp(&a.re, &b.re),
            im: self.addp(&a.im, &b.im),
        }
    }

    pub fn sub(&self, a: &F2, b: &F2) -> F2 {
        F2 {
            re: self.subp(&a.re, &b.re),
            im: self.subp(&a.im, &b.im),
        }
    }

    pub fn neg(&self, a: &F2) -> F2 {
        F2 {
            re: self.negp(&a.re),
            im: self.negp(&a.im),
        }
    }

    pub fn mul(&self, a: &F2, b: &F2) -> F2 {
        let ac = &a.re * &b.re;
        let bd = &a.im * &b.im;
        let cross = (&a.re + &a.im) * (&b.re + &b.im);
        // re = ac - bd, im = cross - ac - bd
        let re = self.subp(&(&ac % &self.p), &(&bd % &self.p));
        let im = (cross - &ac - &bd) % &self.p;
        F2 { re, im }
    }

    pub fn sqr(&self, a: &F2) -> F2 {
        self.mul(a, a)
    }

    /// Multiplication by a small integer.
    pub fn muls(&self, a: &F2, k: i64) -> F2 {
        self.mul(a, &self.small(k))
    }

    /// Multiplication by `i`.
    pub fn mul_i(&self, a: &F2) -> F2 {
        F2 {
            re: self.negp(&a.im),
            im: a.re.clone(),
        }
    }

    /// The `p`-power Frobenius, complex conjugation.
    pub fn frob(&self, a: &F2) -> F2 {
        F2 {
            re: a.re.clone(),
            im: self.negp(&a.im),
        }
    }

    pub fn is_zero(&self, a: &F2) -> bool {
        a.re.is_zero() && a.im.is_zero()
    }

    /// Whether `a` lies in `F_p`.
    pub fn in_fp(&self, a: &F2) -> bool {
        a.im.is_zero()
    }

    fn invp(&self, a: &BigUint) -> BigUint {
        assert!(!a.is_zero(), "inverse of zero");
        a.modpow(&(&self.p - 2u32), &self.p)
    }

    pub fn inv(&self, a: &F2) -> F2 {
        let n = (&a.re * &a.re + &a.im * &a.im) % &self.p;
        let ni = self.invp(&n);
        F2 {
            re: (&a.re * &ni) % &self.p,
            im: (self.negp(&a.im) * &ni) % &self.p,
        }
    }

    pub fn pow(&self, a: &F2, e: &BigUint) -> F2 {
        let mut r = self.one();
        for i in (0..e.bits()).rev() {
            r = self.sqr(&r);
            if e.bit(i) {
                r = self.mul(&r, a);
            }
        }
        r
    }

    /// A square root in `F_{p^2}`, or `None` for a non-square (the
    /// `p = 3 mod 4` algorithm of Adj and Rodriguez-Henriquez), verified.
    pub fn sqrt(&self, a: &F2) -> Option<F2> {
        if self.is_zero(a) {
            return Some(self.zero());
        }
        let e1 = (&self.p - 3u32) >> 2;
        let a1 = self.pow(a, &e1);
        let x0 = self.mul(&a1, a);
        let alpha = self.mul(&a1, &x0);
        let minus_one = self.small(-1);
        let x = if alpha == minus_one {
            self.mul_i(&x0)
        } else {
            let e2 = (&self.p - 1u32) >> 1;
            let b = self.pow(&self.add(&self.one(), &alpha), &e2);
            self.mul(&b, &x0)
        };
        if self.sqr(&x) == *a {
            Some(x)
        } else {
            None
        }
    }

    /// Horner evaluation of an integer polynomial (coefficients low to high).
    pub fn eval_poly(&self, coeffs: &[BigInt], x: &F2) -> F2 {
        let mut acc = self.zero();
        for c in coeffs.iter().rev() {
            acc = self.add(&self.mul(&acc, x), &self.int(c));
        }
        acc
    }

    pub fn hex(&self, a: &F2) -> String {
        if a.im.is_zero() {
            format!("0x{:x}", a.re)
        } else {
            format!("0x{:x} + 0x{:x}*i", a.re, a.im)
        }
    }
}

/// Strong probable-prime test to the first 24 prime bases, after trial
/// division. Deterministic; a composite passes with negligible probability
/// only.
pub fn is_probable_prime(n: &BigUint) -> bool {
    const SMALL: [u32; 24] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    ];
    if *n < BigUint::from(2u32) {
        return false;
    }
    for q in SMALL {
        let q = BigUint::from(q);
        if *n == q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let one = BigUint::one();
    let nm1 = n - &one;
    let s = nm1.trailing_zeros().expect("n - 1 non-zero");
    let d = &nm1 >> s;
    'base: for a in SMALL {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x == one || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % n;
            if x == nm1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}
