//! Algebraic-property tests for Level 5 Fp2 arithmetic.

mod common;

use common::{DetRng, ITER};
use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level5;

fn eq(a: &Fp2<Level5>, b: &Fp2<Level5>) -> bool {
    bool::from(a.ct_equal(b))
}

fn is_zero(a: &Fp2<Level5>) -> bool {
    bool::from(a.ct_is_zero())
}

#[test]
fn l5_gfp2_addition() {
    let mut rng = DetRng::new(b"l5_gfp2_addition");
    let zero = Fp2::<Level5>::zero();
    let one = Fp2::<Level5>::one();

    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let b = rng.random_fp2_level5();
        let c = rng.random_fp2_level5();

        assert!(eq(&a.add(&b).add(&c), &a.add(&b.add(&c))));
        assert!(eq(&a.add(&b), &b.add(&a)));
        assert!(eq(&a.add(&zero), &a));
        assert!(is_zero(&a.add(&a.neg())));
        assert!(eq(&a.add(&one), &a.add_one()));
    }
}

#[test]
fn l5_gfp2_subtraction() {
    let mut rng = DetRng::new(b"l5_gfp2_subtraction");
    let zero = Fp2::<Level5>::zero();
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let b = rng.random_fp2_level5();
        let c = rng.random_fp2_level5();

        assert!(eq(&a.sub(&b).sub(&c), &a.sub(&b.add(&c))));
        assert!(eq(&a.sub(&b), &b.sub(&a).neg()));
        assert!(eq(&a.sub(&zero), &a));
        assert!(is_zero(&a.sub(&a)));
    }
}

#[test]
fn l5_gfp2_multiplication() {
    let mut rng = DetRng::new(b"l5_gfp2_multiplication");
    let zero = Fp2::<Level5>::zero();
    let one = Fp2::<Level5>::one();

    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let b = rng.random_fp2_level5();
        let c = rng.random_fp2_level5();

        assert!(eq(&a.mul(&b).mul(&c), &a.mul(&b.mul(&c))));
        assert!(eq(&a.mul(&b.add(&c)), &a.mul(&b).add(&a.mul(&c))));
        assert!(eq(&a.mul(&b), &b.mul(&a)));
        assert!(eq(&a.mul(&one), &a));
        assert!(is_zero(&a.mul(&zero)));
    }
}

#[test]
fn l5_gfp2_squaring() {
    let mut rng = DetRng::new(b"l5_gfp2_squaring");
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        assert!(eq(&a.sqr(), &a.mul(&a)));
    }
    assert!(is_zero(&Fp2::<Level5>::zero().sqr()));
}

#[test]
fn l5_gfp2_inversion() {
    let mut rng = DetRng::new(b"l5_gfp2_inversion");
    let one = Fp2::<Level5>::one();
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let inv = a.inv();
        assert!(eq(&a.mul(&inv), &one));
    }
}

#[test]
fn l5_gfp2_sqrt_and_is_square() {
    let mut rng = DetRng::new(b"l5_gfp2_sqrt");
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let c = a.sqr();
        assert!(bool::from(c.is_square()));

        let mut c_clone = c.clone();
        assert!(bool::from(c_clone.sqrt_verify()));

        let r = c.sqrt();
        let neg_r = r.neg();
        assert!(eq(&a, &r) || eq(&a, &neg_r));
    }
}

#[test]
fn l5_gfp2_encode_decode_roundtrip() {
    let mut rng = DetRng::new(b"l5_gfp2_encode_decode");
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let bytes = a.encode();
        let b = Fp2::<Level5>::decode(bytes.as_ref()).expect("encoded form must decode");
        assert!(eq(&a, &b));
    }
}

#[test]
fn l5_gfp2_mul_small() {
    let mut rng = DetRng::new(b"l5_gfp2_mul_small");
    for _ in 0..ITER {
        let a = rng.random_fp2_level5();
        let v = rng.random_u32();
        let b = a.mul_small(v);
        let c = Fp2::<Level5>::from_small(v as u64);
        let d = a.mul(&c);
        assert!(eq(&b, &d));
    }
}

mod alloc_like {
    pub use std::vec::Vec;
}

#[test]
fn l5_gfp2_batched_inv() {
    let mut rng = DetRng::new(b"l5_gfp2_batched_inv");
    let one = Fp2::<Level5>::one();

    for &len in &[1usize, 2, 4, 8, 11] {
        let mut xs: alloc_like::Vec<Fp2<Level5>> =
            (0..len).map(|_| rng.random_fp2_level5()).collect();
        let mut t1: alloc_like::Vec<Fp2<Level5>> =
            (0..len).map(|_| Fp2::<Level5>::zero()).collect();
        let mut t2: alloc_like::Vec<Fp2<Level5>> =
            (0..len).map(|_| Fp2::<Level5>::zero()).collect();

        let expected: alloc_like::Vec<Fp2<Level5>> = xs.iter().map(|x| x.inv()).collect();
        let originals = xs.clone();

        Fp2::<Level5>::batched_inv(&mut xs, &mut t1, &mut t2);

        for (i, (got, exp)) in xs.iter().zip(expected.iter()).enumerate() {
            assert!(eq(got, exp), "batched_inv mismatch at i={i} len={len}");
        }

        for i in 0..len {
            let prod = xs[i].mul(&originals[i]);
            assert!(eq(&prod, &one), "x[i] * x[i]^-1 != 1 at i={i} len={len}");
        }
    }
}

#[test]
fn l5_gfp2_pow_vartime() {
    let mut rng = DetRng::new(b"l5_gfp2_pow_vartime");

    let zero_exp: [u64; 0] = [];
    let exp_empty_word: [u64; 1] = [0];
    for _ in 0..8 {
        let a = rng.random_fp2_level5();
        assert!(eq(&a.pow_vartime(&zero_exp), &Fp2::<Level5>::one()));
        assert!(eq(&a.pow_vartime(&exp_empty_word), &Fp2::<Level5>::one()));
    }

    for _ in 0..8 {
        let a = rng.random_fp2_level5();
        assert!(eq(&a.pow_vartime(&[1]), &a));
        assert!(eq(&a.pow_vartime(&[2]), &a.sqr()));
        assert!(eq(&a.pow_vartime(&[3]), &a.sqr().mul(&a)));
        let a2 = a.sqr();
        let a4 = a2.sqr();
        let a5 = a4.mul(&a);
        assert!(eq(&a.pow_vartime(&[5]), &a5));
    }

    for _ in 0..4 {
        let a = rng.random_fp2_level5();
        let mut sq = a.clone();
        for _ in 0..64 {
            sq = sq.sqr();
        }
        let exp = [0u64, 1u64];
        assert!(eq(&a.pow_vartime(&exp), &sq));
    }

    for _ in 0..8 {
        let a = rng.random_fp2_level5();
        let m: u64 = (rng.random_u32() as u64) & 0x3FF;
        let n: u64 = (rng.random_u32() as u64) & 0x3FF;
        let am = a.pow_vartime(&[m]);
        let an = a.pow_vartime(&[n]);
        let amn = a.pow_vartime(&[m + n]);
        assert!(eq(&amn, &am.mul(&an)), "m={m} n={n}");
    }
}
