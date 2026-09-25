//! The x86-64 field backends (`fp/<prime>_asm.rs`) against big-integer
//! arithmetic, operation by operation, on random inputs, for the three
//! round-3 primes, under both run-time paths: the saturated representation
//! must agree with the canonical encoding on every operation, the `Fp2`
//! layer built on the backend must agree with complex arithmetic modulo
//! `p`, and the assembly and the portable kernels must agree limb for limb
//! (they run the same schedule). On other architectures the test is empty.
//!
//! The kernel choice is process-wide, so the tests serialise on a mutex
//! and restore the detection afterwards.

#![cfg(target_arch = "x86_64")]

use num_bigint::BigUint;
use num_traits::{One, Zero};
use sqisign_verify::fp::{
    arithmetic_backend, force_auto_arithmetic, force_portable_arithmetic, Fp, Fp2, FpBackend,
};
use sqisign_verify::params::{P324_3, P500_27, P664_17};
use std::sync::Mutex;

static MODE: Mutex<()> = Mutex::new(());

/// A small deterministic generator (xorshift128+), enough for test inputs.
struct Xs(u64, u64);
impl Xs {
    fn next(&mut self) -> u64 {
        let mut s1 = self.0;
        let s0 = self.1;
        self.0 = s0;
        s1 ^= s1 << 23;
        self.1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
        self.1.wrapping_add(s0)
    }
    fn big(&mut self, modulus: &BigUint) -> BigUint {
        let mut bytes = [0u8; 96];
        for c in bytes.chunks_mut(8) {
            c.copy_from_slice(&self.next().to_le_bytes());
        }
        BigUint::from_bytes_le(&bytes) % modulus
    }
}

struct Prime {
    p: BigUint,
    bytes: usize,
}

fn prime(c: u32, e: u32, bytes: usize) -> Prime {
    Prime {
        p: (BigUint::from(c) << e) - BigUint::one(),
        bytes,
    }
}

fn to_bytes(x: &BigUint, n: usize) -> Vec<u8> {
    let mut out = x.to_bytes_le();
    out.resize(n, 0);
    out
}

fn fp<L: FpBackend>(pr: &Prime, x: &BigUint) -> Fp<L> {
    Fp::<L>::decode(&to_bytes(x, pr.bytes)).expect("canonical")
}

fn big<L: FpBackend>(x: &Fp<L>) -> BigUint {
    BigUint::from_bytes_le(&x.encode())
}

fn legendre_is_square(x: &BigUint, p: &BigUint) -> bool {
    if x.is_zero() {
        return true;
    }
    let e = (p - BigUint::one()) >> 1u32;
    x.modpow(&e, p).is_one()
}

fn fp_operations<L: FpBackend>(pr: &Prime, rounds: usize) {
    let p = &pr.p;
    let mut r = Xs(0x9e3779b97f4a7c15, 0xd1b54a32d192ed03);
    let two_inv = (p + BigUint::one()) >> 1u32;
    let three_inv = BigUint::from(3u32).modpow(&(p - BigUint::from(2u32)), p);
    for i in 0..rounds {
        let a = r.big(p);
        let b = r.big(p);
        let (fa, fb) = (fp::<L>(pr, &a), fp::<L>(pr, &b));
        assert_eq!(big(&fa), a, "decode/encode round trip");
        assert_eq!(big(&fa.add(&fb)), (&a + &b) % p, "add #{i}");
        assert_eq!(big(&fa.sub(&fb)), (&a + p - &b) % p, "sub #{i}");
        assert_eq!(big(&fa.neg()), (p - &a) % p, "neg #{i}");
        assert_eq!(big(&fa.mul(&fb)), (&a * &b) % p, "mul #{i}");
        assert_eq!(big(&fa.sqr()), (&a * &a) % p, "sqr #{i}");
        assert_eq!(big(&fa.half()), (&a * &two_inv) % p, "half #{i}");
        assert_eq!(big(&fa.div3()), (&a * &three_inv) % p, "div3 #{i}");
        let small = (r.next() >> 32) as u32;
        assert_eq!(
            big(&fa.mul_small(small)),
            (&a * BigUint::from(small)) % p,
            "mul_small #{i}"
        );
        if !a.is_zero() {
            let inv = a.modpow(&(p - BigUint::from(2u32)), p);
            assert_eq!(big(&fa.inv()), inv, "inv #{i}");
        }
        assert_eq!(
            bool::from(fa.is_square()),
            legendre_is_square(&a, p),
            "is_square #{i}"
        );
        let sq = (&a * &a) % p;
        let root = big(&fp::<L>(pr, &sq).sqrt());
        assert!(root == a || root == (p - &a) % p, "sqrt #{i}");
        assert!(
            bool::from(fa.add(&fb).sub(&fb).ct_equal(&fa)),
            "ct_equal #{i}"
        );
        assert!(bool::from(fa.sub(&fa).ct_is_zero()), "ct_is_zero #{i}");
        let lazy = fa.add(&fb).add(&fa).sub(&fb);
        assert_eq!(
            big(&lazy),
            (BigUint::from(2u32) * &a) % p,
            "lazy encode #{i}"
        );
        // long lazy chains stay within the representation's bound
        let mut acc = fa.clone();
        for _ in 0..64 {
            acc = acc.add(&fb).add(&fa).sub(&fb);
        }
        assert_eq!(
            big(&acc),
            (BigUint::from(65u32) * &a) % p,
            "lazy chain #{i}"
        );
    }
    assert!(bool::from(Fp::<L>::zero().ct_is_zero()));
    assert_eq!(big(&Fp::<L>::one()), BigUint::one());
    assert_eq!(big(&Fp::<L>::from_small(12345)), BigUint::from(12345u32));
    assert!(bool::from(Fp::<L>::zero().inv().ct_is_zero()));
    assert!(bool::from(Fp::<L>::zero().is_square()));
    assert!(Fp::<L>::decode(&to_bytes(p, pr.bytes)).is_none());
    assert!(Fp::<L>::decode(&to_bytes(&(p + BigUint::from(7u32)), pr.bytes)).is_none());
    // decode_reduce
    let chunk = pr.bytes - 1;
    let mut r = Xs(0x1234567887654321, 0xfedcba9876543210);
    for len in [
        0usize,
        1,
        7,
        chunk,
        pr.bytes,
        pr.bytes + 7,
        2 * chunk,
        2 * chunk + 1,
        3 * chunk + 5,
    ] {
        let mut bytes = vec![0u8; len];
        for b in bytes.iter_mut() {
            *b = r.next() as u8;
        }
        let want = BigUint::from_bytes_le(&bytes) % p;
        let got = big(&Fp::<L>::decode_reduce(&bytes));
        // the portable rule: a trailing full-length chunk not below p is zeroed
        let rem = len % chunk;
        let tail_ok = rem != pr.bytes || BigUint::from_bytes_le(&bytes[len - rem..]) < *p;
        if tail_ok {
            assert_eq!(got, want, "decode_reduce of {len} bytes");
        }
    }
}

fn fp2_operations<L: FpBackend>(pr: &Prime, rounds: usize) {
    let p = &pr.p;
    let mut r = Xs(0x0123456789abcdef, 0x0fedcba987654321);
    let modp = |x: &BigUint| x % p;
    let sub = |x: &BigUint, y: &BigUint| (x + p - y) % p;
    for i in 0..rounds {
        let (a0, a1, b0, b1) = (r.big(p), r.big(p), r.big(p), r.big(p));
        let fa = Fp2 {
            re: fp::<L>(pr, &a0),
            im: fp::<L>(pr, &a1),
        };
        let fb = Fp2 {
            re: fp::<L>(pr, &b0),
            im: fp::<L>(pr, &b1),
        };
        let m = fa.mul(&fb);
        assert_eq!(
            big(&m.re),
            sub(&modp(&(&a0 * &b0)), &modp(&(&a1 * &b1))),
            "fp2 mul re #{i}"
        );
        assert_eq!(
            big(&m.im),
            modp(&(&a0 * &b1 + &a1 * &b0)),
            "fp2 mul im #{i}"
        );
        let s = fa.sqr();
        assert_eq!(
            big(&s.re),
            sub(&modp(&(&a0 * &a0)), &modp(&(&a1 * &a1))),
            "fp2 sqr re #{i}"
        );
        assert_eq!(
            big(&s.im),
            modp(&(BigUint::from(2u32) * &a0 * &a1)),
            "fp2 sqr im #{i}"
        );
        let one = fa.mul(&fa.inv());
        assert_eq!(big(&one.re), BigUint::one(), "fp2 inv re #{i}");
        assert!(big(&one.im).is_zero(), "fp2 inv im #{i}");
        let root = s.sqrt();
        assert!(bool::from(root.sqr().ct_equal(&s)), "fp2 sqrt #{i}");
        assert!(bool::from(s.is_square()), "fp2 is_square #{i}");
        let n = fa.mul(&fa.frob());
        assert_eq!(big(&n.re), modp(&(&a0 * &a0 + &a1 * &a1)), "norm #{i}");
        assert!(big(&n.im).is_zero(), "norm im #{i}");
        // the lazy inputs of the fused squaring: a sum of many terms squared
        let mut acc = fa.clone();
        for _ in 0..32 {
            acc = acc.add(&fa);
        }
        let sq = acc.sqr();
        let (c0, c1) = (
            modp(&(BigUint::from(33u32) * &a0)),
            modp(&(BigUint::from(33u32) * &a1)),
        );
        assert_eq!(
            big(&sq.re),
            sub(&modp(&(&c0 * &c0)), &modp(&(&c1 * &c1))),
            "lazy fp2 sqr #{i}"
        );
    }
}

/// The raw limbs of every operation's result on random inputs, including
/// lazily reduced ones, for the limb-for-limb comparison of the two paths.
fn raw_results<L: FpBackend>(pr: &Prime, rounds: usize) -> Vec<Vec<u64>> {
    let p = &pr.p;
    let mut r = Xs(0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15);
    let mut out = Vec::new();
    let mut push = |x: &Fp<L>| out.push(x.as_limbs().to_vec());
    for _ in 0..rounds {
        let (a, b) = (r.big(p), r.big(p));
        let (fa, fb) = (fp::<L>(pr, &a), fp::<L>(pr, &b));
        let lazy = fa.add(&fb).add(&fa);
        push(&fa.mul(&fb));
        push(&fa.sqr());
        push(&lazy.mul(&fb));
        push(&lazy.sqr());
        push(&fa.inv());
        push(&fa.sqrt());
        push(&fa.mul_small(0xdead_beef));
        let x = Fp2 {
            re: fa.clone(),
            im: fb.clone(),
        };
        let y = Fp2 {
            re: lazy.clone(),
            im: fa.sub(&fb),
        };
        let m = x.mul(&y);
        push(&m.re);
        push(&m.im);
        let s = y.sqr();
        push(&s.re);
        push(&s.im);
        let i = x.inv();
        push(&i.re);
        push(&i.im);
    }
    out
}

fn both_paths<L: FpBackend>(pr: &Prime) {
    let _guard = MODE.lock().unwrap_or_else(|e| e.into_inner());
    force_auto_arithmetic();
    let detected = arithmetic_backend();
    fp_operations::<L>(pr, 300);
    fp2_operations::<L>(pr, 150);
    let fast = raw_results::<L>(pr, 64);
    force_portable_arithmetic();
    assert!(arithmetic_backend().starts_with("portable"));
    fp_operations::<L>(pr, 300);
    fp2_operations::<L>(pr, 150);
    let portable = raw_results::<L>(pr, 64);
    force_auto_arithmetic();
    assert_eq!(
        fast, portable,
        "the two paths differ in a limb ({detected})"
    );
    eprintln!(
        "{} bits: detected {detected}; portable path agrees limb for limb",
        pr.p.bits()
    );
}

#[test]
fn p324_3_matches_big_integers_on_both_paths() {
    both_paths::<P324_3>(&prime(3, 324, 41));
}

#[test]
fn p500_27_matches_big_integers_on_both_paths() {
    both_paths::<P500_27>(&prime(27, 500, 64));
}

#[test]
fn p664_17_matches_big_integers_on_both_paths() {
    both_paths::<P664_17>(&prime(17, 664, 84));
}
