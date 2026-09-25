//! Properties of the fixed-precision integers against `num-bigint`, at
//! each container size.

mod mp_common;

use mp_common::*;
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use sqisign_rs::mp::Ibz;

fn check<const N: usize>() {
    let mut rng = DetRng::new(&format!("props{N}"));
    for _ in 0..200 {
        let sa = SIZES[rng.u32() as usize % SIZES.len()];
        let sb = SIZES[rng.u32() as usize % SIZES.len()];
        let a: Ibz<N> = rng.ibz(sa, true);
        let b: Ibz<N> = rng.ibz(sb, true);
        let (ba, bb) = (to_big(&a), to_big(&b));
        assert_eq!(to_big(&a.add(&b)), &ba + &bb);
        assert_eq!(to_big(&a.sub(&b)), &ba - &bb);
        assert_eq!(to_big(&a.mul(&b)), &ba * &bb);
        assert_eq!(to_big(&a.neg()), -&ba);
        assert_eq!(to_big(&a.abs()), ba.abs());
        assert_eq!(a.cmp(&b), ba.cmp(&bb));
        assert_eq!(a.bitsize(), ba.bits() as i32);
        assert_eq!(a.bitsize_ct(), ba.bits() as i32);
        if !bb.is_zero() {
            let (q, r) = a.div(&b);
            assert_eq!(to_big(&q), &ba / &bb, "trunc quotient");
            assert_eq!(to_big(&r), &ba % &bb, "trunc remainder");
            let m = a.modulo(&b);
            assert_eq!(to_big(&m), ba.mod_floor(&bb.abs()), "non-negative modulus");
        }
        let sh = rng.u32() % 300;
        // The reference's div_2exp is a logical shift of the active limbs with
        // the bound reduced by the shift. That is the floor division as long
        // as the zero fill from the top of the active limbs stays above the
        // sign bit of the result's top limb; outside that domain (negative
        // values shifted across limbs) the reference's result is not the
        // floor, and the port reproduces it (see mp_crossvalidate).
        let nwords = ((a.get_bound() + 63) / 64) as i64;
        let new_bound = if (sh as i32) < a.get_bound() {
            a.get_bound() - sh as i32
        } else {
            1
        } as i64;
        let new_words = (new_bound + 63) / 64;
        if ba >= BigInt::zero() || 64 * new_words - 1 < 64 * nwords - sh as i64 {
            assert_eq!(
                to_big(&a.div_2exp(sh)),
                ba.clone() >> sh,
                "floor shift of {a:?} by {sh}"
            );
        }
        assert_eq!(to_big(&a.mul_2exp(sh)), ba.clone() << sh);
        let g = a.gcd(&b);
        assert_eq!(to_big(&g), ba.gcd(&bb));
        let (g2, u, v) = a.xgcd(&b);
        assert_eq!(to_big(&g2), ba.gcd(&bb));
        assert_eq!(&to_big(&u) * &ba + &to_big(&v) * &bb, ba.gcd(&bb), "bezout");
        let m = b.abs();
        if !m.is_zero() {
            let bm = to_big(&m);
            match a.invmod(&m) {
                Some(inv) => {
                    assert!(ba.gcd(&bm).is_one());
                    let bi = to_big(&inv);
                    assert!(bi >= BigInt::zero() && bi < bm);
                    assert!((&bi * &ba).mod_floor(&bm).is_one() || bm.is_one());
                }
                None => assert!(!ba.gcd(&bm).is_one()),
            }
        }
        if !a.is_zero() {
            assert_eq!(a.two_adic(), ba.trailing_zeros().unwrap() as i32);
        }
        let mod2 = a.mod2exp(sh);
        assert_eq!(to_big(&mod2), ba.mod_floor(&(BigInt::one() << sh)));
        // sqrt floor of the magnitude
        let aa = a.abs();
        let s = aa.sqrt_floor();
        assert_eq!(to_big(&s), ba.abs().sqrt());
        // constant-time helpers
        assert_eq!(a.ct_lt_mask(&b) == u64::MAX, ba < bb);
        assert_eq!(a.ct_nonzero_mask() == u64::MAX, !ba.is_zero());
        let off = (rng.u32() % 700) as i32;
        let window = (ba.abs() >> off) & BigInt::from(u64::MAX);
        assert_eq!(BigInt::from(a.extract_u64(off)), window);
        let shift = (rng.u32() % 601) as i32 - 300;
        let out_bits = (rng.u32() % 900) as i32 + 64;
        let shifted = a.ct_shift(shift, 300, out_bits);
        let expect = if shift >= 0 {
            ba.clone() << shift
        } else {
            ba.clone() >> (-shift)
        };
        assert_eq!(to_big(&shifted), wrap(&expect, out_bits), "ct_shift");
        // the wide product buffer has na + nb + 1 limbs and the output must fit it
        let wide_bits = 64 * (((a.get_bound() + 63) / 64) + ((b.get_bound() + 63) / 64) + 1);
        let ob = out_bits.min(wide_bits);
        let p = (rng.u32() % 500) as i32 + 64;
        let hm = a.ct_highmul_p(&b, p, ob);
        assert_eq!(to_big(&hm), wrap(&((&ba * &bb) >> p), ob), "highmul_p");
        let hs = a.ct_highmul_s(&b, shift, 300, ob);
        let prod = &ba * &bb;
        let expect = if shift >= 0 {
            prod << shift
        } else {
            prod >> (-shift)
        };
        assert_eq!(
            to_big(&hs),
            wrap(&wrap(&expect, wide_bits), ob),
            "highmul_s"
        );
        let k = (rng.u32() % 2001) as i32 - 1000;
        let bound = a.get_bound() + 12;
        assert_eq!(
            to_big(&a.mul_by_int_and_set_bound(k, bound)),
            &ba * BigInt::from(k)
        );
    }
    // modular arithmetic on odd moduli
    for _ in 0..20 {
        let bits = 64 + (rng.u32() % 300) as i32;
        let mut m: Ibz<N> = rng.ibz(bits, false);
        m.limbs_set_odd();
        let bm = to_big(&m);
        let x: Ibz<N> = rng.ibz(bits, false).modulo(&m);
        let e: Ibz<N> = rng.ibz(bits, false).modulo(&m);
        let r = x.pow_mod(&e, &m);
        assert_eq!(to_big(&r), to_big(&x).modpow(&to_big(&e), &bm));
        // crt
        let mut m2: Ibz<N> = rng.ibz(bits, false);
        m2.limbs_set_odd();
        if m.gcd(&m2).is_one() {
            let a1 = rng.ibz::<N>(bits, true).modulo(&m);
            let a2 = rng.ibz::<N>(bits, true).modulo(&m2);
            let x = Ibz::crt(&a1, &a2, &m, &m2);
            let bx = to_big(&x);
            assert_eq!(bx.mod_floor(&bm), to_big(&a1));
            assert_eq!(bx.mod_floor(&to_big(&m2)), to_big(&a2));
        }
    }
    // primes: sieve-and-test loop, sqrt mod p, sqrt of -1
    let mut d = default_domain(&format!("props-prime{N}"));
    for _ in 0..4 {
        let bits = 40 + (rng.u32() % 200) as i32;
        let p = loop {
            let mut c: Ibz<N> = rng.ibz(bits, false);
            c.limbs_set_odd();
            if c.probab_prime(64, &mut d) {
                break c;
            }
        };
        let bp = to_big(&p);
        assert!(bp.bits() as i32 == bits);
        // trial division check against small factors as a sanity check
        for q in [3u32, 5, 7, 11, 13, 17, 19, 23, 29, 31] {
            assert!(!(&bp % q).is_zero());
        }
        assert!(!p.mul(&p).probab_prime(64, &mut d));
        let a: Ibz<N> = rng.ibz(bits, false).modulo(&p);
        let sq = a.mul(&a).modulo(&p);
        let r = sq.sqrt_mod_p(&p, &mut d).expect("square root of a square");
        assert_eq!(r.mul(&r).modulo(&p), sq);
        if (p.get() & 3) == 1 {
            let r = Ibz::sqrt_m1_mod_verified(&p).expect("root of -1");
            assert_eq!(r.mul(&r).add(&Ibz::one()).modulo(&p), Ibz::zero());
            let r2 = Ibz::sqrt_m1_mod(&p, &mut d).expect("root of -1");
            assert_eq!(r2.mul(&r2).add(&Ibz::one()).modulo(&p), Ibz::zero());
        }
        assert_eq!(a.legendre(&p), jacobi(&to_big(&a), &bp));
    }
    // strings round-trip
    let x: Ibz<N> = rng.ibz(500, true);
    let mut buf = [0u8; 64 * 64 + 2];
    let dec = x.to_decimal(&mut buf).to_string();
    assert_eq!(dec, to_big(&x).to_string());
    assert_eq!(Ibz::<N>::from_str_radix(&dec, 10).unwrap(), x);
}

trait SetOdd {
    fn limbs_set_odd(&mut self);
}

impl<const N: usize> SetOdd for Ibz<N> {
    fn limbs_set_odd(&mut self) {
        if self.is_even() {
            *self = self.add(&Ibz::one());
        }
    }
}

/// Two's complement truncation to `bits` bits as a signed value.
fn wrap(v: &BigInt, bits: i32) -> BigInt {
    let m = BigInt::one() << bits;
    let r = v.mod_floor(&m);
    if r >= (BigInt::one() << (bits - 1)) {
        r - m
    } else {
        r
    }
}

fn jacobi(a: &BigInt, n: &BigInt) -> i32 {
    let mut a = a.mod_floor(n);
    let mut n = n.clone();
    let mut sign = 1;
    while !a.is_zero() {
        while a.is_even() {
            a >>= 1;
            let r = (&n % 8u32).to_u32_digits().1.first().copied().unwrap_or(0);
            if r == 3 || r == 5 {
                sign = -sign;
            }
        }
        core::mem::swap(&mut a, &mut n);
        let a3 = (&a % 4u32) == BigInt::from(3);
        let n3 = (&n % 4u32) == BigInt::from(3);
        if a3 && n3 {
            sign = -sign;
        }
        a = a.mod_floor(&n);
    }
    if n.is_one() {
        sign
    } else {
        0
    }
}

#[test]
fn props_30() {
    check::<30>();
}

#[test]
fn props_40() {
    check::<40>();
}

#[test]
fn props_55() {
    check::<55>();
}

#[test]
fn constants() {
    let z = Ibz::<30>::zero();
    assert!(z.is_zero() && z.is_positive() && z.get_bound() == 1);
    assert_eq!(to_big(&Ibz::<30>::set(-7, 32)), BigInt::from(-7));
    assert_eq!(
        to_big(&Ibz::<30>::from_le_bytes(&[1, 2, 3])),
        BigInt::from(0x030201)
    );
    let _ = Sign::Plus;
}
