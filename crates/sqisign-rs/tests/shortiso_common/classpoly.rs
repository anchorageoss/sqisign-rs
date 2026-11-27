//! Class numbers and Hilbert (ring) class polynomials of imaginary
//! quadratic orders, computed from scratch in fixed-point arithmetic on
//! `num-bigint`: the reduced primitive forms `(a, b, c)` of discriminant
//! `-D`, the values `j((-b + sqrt(-D)) / 2a)` from the `q`-expansions of
//! `E4` and `Delta`, the product of `X - j`, rounded to integers with the
//! rounding error checked.

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// The reduced primitive positive definite forms of discriminant `-d`.
pub fn reduced_forms(d: u64) -> Vec<(i64, i64, i64)> {
    assert!(d % 4 == 0 || d % 4 == 3, "-D is a discriminant");
    let d = d as i64;
    let mut out = Vec::new();
    let mut a = 1i64;
    while 3 * a * a <= d {
        for b in -a + 1..=a {
            if (b * b + d) % (4 * a) != 0 {
                continue;
            }
            let c = (b * b + d) / (4 * a);
            if c < a || (c == a && b < 0) {
                continue;
            }
            if a.gcd(&b).gcd(&c) != 1 {
                continue;
            }
            out.push((a, b, c));
        }
        a += 1;
    }
    out
}

/// The class number `h(-D)` of the order of discriminant `-D`.
pub fn class_number(d: u64) -> usize {
    reduced_forms(d).len()
}

/// Fixed-point numbers: integers scaled by `2^prec`.
struct Fx {
    prec: u32,
}

type C = (BigInt, BigInt);

impl Fx {
    fn one(&self) -> BigInt {
        BigInt::one() << self.prec
    }

    fn mul(&self, a: &BigInt, b: &BigInt) -> BigInt {
        (a * b) >> self.prec
    }

    fn cmul(&self, a: &C, b: &C) -> C {
        (
            (&a.0 * &b.0 - &a.1 * &b.1) >> self.prec,
            (&a.0 * &b.1 + &a.1 * &b.0) >> self.prec,
        )
    }

    fn cdiv(&self, a: &C, b: &C) -> C {
        let n = &b.0 * &b.0 + &b.1 * &b.1;
        let re = &a.0 * &b.0 + &a.1 * &b.1;
        let im = &a.1 * &b.0 - &a.0 * &b.1;
        ((re << self.prec) / &n, (im << self.prec) / &n)
    }

    /// `atan(1/n)` by its Taylor series.
    fn atan_inv(&self, n: u32) -> BigInt {
        let n2 = BigInt::from(n) * BigInt::from(n);
        let mut term = self.one() / BigInt::from(n);
        let mut sum = term.clone();
        let mut k = 1u32;
        while !term.is_zero() {
            term = -term / &n2;
            sum += &term / BigInt::from(2 * k + 1);
            k += 1;
        }
        sum
    }

    /// Machin: `pi = 16 atan(1/5) - 4 atan(1/239)`.
    fn pi(&self) -> BigInt {
        self.atan_inv(5) * 16 - self.atan_inv(239) * 4
    }

    /// `exp(z)` for complex `z`: scaling by `2^-k`, Taylor, `k` squarings.
    fn cexp(&self, z: &C) -> C {
        let mag = z.0.abs().max(z.1.abs());
        let k = (mag.bits() as i64 - self.prec as i64 + 16).max(0) as u32;
        let zs: C = (&z.0 >> k, &z.1 >> k);
        let mut term: C = (self.one(), BigInt::zero());
        let mut sum = term.clone();
        let mut n = 1u32;
        while !(term.0.is_zero() && term.1.is_zero()) {
            term = self.cmul(&term, &zs);
            term = (&term.0 / BigInt::from(n), &term.1 / BigInt::from(n));
            sum = (&sum.0 + &term.0, &sum.1 + &term.1);
            n += 1;
        }
        for _ in 0..k {
            sum = self.cmul(&sum, &sum);
        }
        sum
    }
}

/// `sigma_3(n)` for `n <= max`.
fn sigma3(max: usize) -> Vec<BigInt> {
    let mut s = vec![BigInt::zero(); max + 1];
    for d in 1..=max {
        let d3 = BigInt::from(d).pow(3);
        let mut m = d;
        while m <= max {
            s[m] += &d3;
            m += d;
        }
    }
    s
}

/// The class polynomial `H_{-D}(X)` of the order of discriminant `-D`:
/// integer coefficients, low to high, monic of degree `h(-D)`.
///
/// Working precision is twice the estimated size of the largest
/// coefficient plus 320 bits; the function panics unless every coefficient
/// is within `2^-64` of an integer with imaginary part below `2^-64`.
pub fn class_polynomial(d: u64) -> Vec<BigInt> {
    let forms = reduced_forms(d);
    let sqrt_d = (d as f64).sqrt();
    // log2 of the product of the |j|: about pi sqrt(D) sum 1/a / ln 2
    let bits: f64 = forms
        .iter()
        .map(|&(a, _, _)| {
            std::f64::consts::PI * sqrt_d / (a as f64) / std::f64::consts::LN_2 + 12.0
        })
        .sum();
    let prec = 2 * bits.ceil() as u32 + 320;
    let fx = Fx { prec };
    let pi = fx.pi();
    let sqrt_d_fx = BigInt::from_biguint(
        Sign::Plus,
        (BigInt::from(d) << (2 * prec)).to_biguint().unwrap().sqrt(),
    );
    let pi_sqrt_d = fx.mul(&pi, &sqrt_d_fx);

    // number of q-terms for the worst form: |q| = exp(-pi sqrt(D) / a)
    let mut poly: Vec<C> = vec![(fx.one(), BigInt::zero())];
    for &(a, b, _) in &forms {
        // 2 pi i tau = -pi sqrt(D) / a - i pi b / a
        let z: C = (
            -(&pi_sqrt_d / BigInt::from(a)),
            -(&pi * BigInt::from(b)) / BigInt::from(a),
        );
        let q = fx.cexp(&z);
        let log2_q = std::f64::consts::PI * sqrt_d / (a as f64) / std::f64::consts::LN_2;
        let nterms = (prec as f64 / log2_q).ceil() as usize + 4;
        let s3 = sigma3(nterms);
        // powers of q
        let mut qp: Vec<C> = Vec::with_capacity(nterms + 1);
        qp.push((fx.one(), BigInt::zero()));
        for n in 1..=nterms {
            let next = fx.cmul(&qp[n - 1], &q);
            qp.push(next);
        }
        // E4 = 1 + 240 sum sigma_3(n) q^n
        let mut e4: C = (fx.one(), BigInt::zero());
        for n in 1..=nterms {
            e4.0 += &qp[n].0 * &s3[n] * 240;
            e4.1 += &qp[n].1 * &s3[n] * 240;
        }
        // eta^24 / q = (sum_{n in Z} (-1)^n q^(n (3n - 1) / 2))^24
        let mut euler: C = (fx.one(), BigInt::zero());
        let mut n = 1usize;
        loop {
            let e1 = n * (3 * n - 1) / 2;
            let e2 = n * (3 * n + 1) / 2;
            if e1 > nterms {
                break;
            }
            let sign = if n % 2 == 0 { 1 } else { -1 };
            euler.0 += &qp[e1].0 * sign;
            euler.1 += &qp[e1].1 * sign;
            if e2 <= nterms {
                euler.0 += &qp[e2].0 * sign;
                euler.1 += &qp[e2].1 * sign;
            }
            n += 1;
        }
        let mut delta = q.clone();
        let e2 = fx.cmul(&euler, &euler);
        let e4p = fx.cmul(&e2, &e2);
        let e8 = fx.cmul(&e4p, &e4p);
        let e24 = fx.cmul(&fx.cmul(&e8, &e8), &e8);
        delta = fx.cmul(&delta, &e24);
        let e4_cubed = fx.cmul(&fx.cmul(&e4, &e4), &e4);
        let j = fx.cdiv(&e4_cubed, &delta);
        // poly *= (X - j)
        let mut next: Vec<C> = vec![(BigInt::zero(), BigInt::zero()); poly.len() + 1];
        for (i, c) in poly.iter().enumerate() {
            next[i + 1].0 += &c.0;
            next[i + 1].1 += &c.1;
            let m = fx.cmul(c, &j);
            next[i].0 -= m.0;
            next[i].1 -= m.1;
        }
        poly = next;
    }
    let half = BigInt::one() << (prec - 1);
    let tol = BigInt::one() << (prec - 64);
    poly.iter()
        .map(|(re, im)| {
            let r = (re + &half) >> prec;
            let err = (re - (&r << prec)).abs();
            assert!(
                err < tol && im.abs() < tol,
                "H_-{d}: a coefficient is not an integer to 64 fractional bits"
            );
            r
        })
        .collect()
}

/// The polynomial as text: one decimal coefficient per line, constant term
/// first, after a header line.
pub fn to_text(d: u64, poly: &[BigInt]) -> String {
    let mut s = format!(
        "# H_-{d}(X), class polynomial of the order of discriminant -{d}; degree {}; \
         coefficients of X^0, X^1, ... one per line\n",
        poly.len() - 1
    );
    for c in poly {
        s.push_str(&c.to_string());
        s.push('\n');
    }
    s
}

pub fn from_text(s: &str) -> Vec<BigInt> {
    s.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.trim().parse().expect("decimal coefficient"))
        .collect()
}
