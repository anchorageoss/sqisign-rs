//! Algebraic-property tests for `Fp` and `Fp2`, instantiated for every
//! prime the crate knows. Each test checks one identity over many
//! deterministic random inputs. Byte-exact agreement with the C reference
//! is checked separately in `c_crossvalidate.rs`.

mod common;

macro_rules! fp_props {
    ($modname:ident, $prime:ty $(, $cfg:meta)?) => {
        $(#[cfg($cfg)])?
        mod $modname {
            use crate::common::{eq, eq2, DetRng, ITER};
            use hybrid_array::typenum::Unsigned;
            use sqisign_verify::fp::{Fp, Fp2};
            use sqisign_verify::params::Prime;
            use subtle::Choice;

            type P = $prime;
            type F = Fp<P>;
            type F2 = Fp2<P>;

            fn label(s: &str) -> Vec<u8> {
                format!("{}/{}", stringify!($modname), s).into_bytes()
            }

            #[test]
            fn constants_are_consistent() {
                let one = F::one();
                let zero = F::zero();
                assert!(bool::from(zero.ct_is_zero()));
                assert!(!bool::from(one.ct_is_zero()));
                assert!(eq(&F::from_small(1), &one));
                assert!(eq(&one.add(&one), &F::from_small(2)));
                // p encodes as itself, and p is c * 2^e - 1
                let p = P::prime_le_bytes();
                assert_eq!(p.len(), <P as Prime>::FpEncodedBytes::USIZE);
                assert_eq!(p[0] & 1, 1, "p is odd");
                assert!(F::decode(p).is_none(), "p itself is not a canonical element");
                let mut pm1 = p.to_vec();
                pm1[0] -= 1;
                let x = F::decode(&pm1).expect("p - 1 is canonical");
                assert!(eq(&x.add(&one), &zero), "p - 1 + 1 = 0");
                assert!(eq(&x, &one.neg()), "p - 1 = -1");
            }

            #[test]
            fn add_sub_neg() {
                let mut rng = DetRng::new(&label("add_sub_neg"));
                let zero = F::zero();
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let b = rng.random_fp::<P>();
                    let c = rng.random_fp::<P>();
                    assert!(eq(&a.add(&b), &b.add(&a)));
                    assert!(eq(&a.add(&b).add(&c), &a.add(&b.add(&c))));
                    assert!(eq(&a.add(&zero), &a));
                    assert!(eq(&a.sub(&a), &zero));
                    assert!(eq(&a.sub(&b), &a.add(&b.neg())));
                    assert!(eq(&a.neg().neg(), &a));
                    assert!(eq(&a.add(&a.neg()), &zero));
                    assert!(eq(&a.sub(&b).add(&b), &a));
                }
            }

            #[test]
            fn mul_sqr() {
                let mut rng = DetRng::new(&label("mul_sqr"));
                let one = F::one();
                let zero = F::zero();
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let b = rng.random_fp::<P>();
                    let c = rng.random_fp::<P>();
                    assert!(eq(&a.mul(&b), &b.mul(&a)));
                    assert!(eq(&a.mul(&b).mul(&c), &a.mul(&b.mul(&c))));
                    assert!(eq(&a.mul(&b.add(&c)), &a.mul(&b).add(&a.mul(&c))));
                    assert!(eq(&a.mul(&one), &a));
                    assert!(eq(&a.mul(&zero), &zero));
                    assert!(eq(&a.sqr(), &a.mul(&a)));
                    assert!(eq(&a.neg().sqr(), &a.sqr()));
                }
            }

            #[test]
            fn small_multipliers() {
                let mut rng = DetRng::new(&label("small"));
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let k = rng.random_u32();
                    assert!(eq(&a.mul_small(k), &a.mul(&F::from_small(k as u64))));
                    assert!(eq(&a.half().add(&a.half()), &a));
                    assert!(eq(&a.div3().mul_small(3), &a));
                    assert!(eq(&a.mul_small(0), &F::zero()));
                    assert!(eq(&a.mul_small(1), &a));
                }
            }

            #[test]
            fn inverse() {
                let mut rng = DetRng::new(&label("inverse"));
                let one = F::one();
                assert!(bool::from(F::zero().inv().ct_is_zero()), "1/0 is 0");
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    if bool::from(a.ct_is_zero()) {
                        continue;
                    }
                    let ai = a.inv();
                    assert!(eq(&a.mul(&ai), &one));
                    assert!(eq(&ai.inv(), &a));
                    // 1/a = a * (a^((p-3)/4))^4
                    let pro = a.exp3div4();
                    let pro4 = pro.sqr().sqr();
                    assert!(eq(&a.mul(&pro4), &ai));
                }
            }

            #[test]
            fn square_roots() {
                let mut rng = DetRng::new(&label("sqrt"));
                let mut squares = 0;
                let mut nonsquares = 0;
                assert!(bool::from(F::zero().is_square()));
                assert!(bool::from(F::one().is_square()));
                // -1 is a non-square since p = 3 mod 4
                assert!(!bool::from(F::one().neg().is_square()));
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let sq = a.sqr();
                    assert!(bool::from(sq.is_square()));
                    let r = sq.sqrt();
                    assert!(eq(&r, &a) || eq(&r, &a.neg()), "sqrt(a^2) = +-a");
                    // sqrt is x * x^((p-3)/4)
                    assert!(eq(&r, &sq.mul(&sq.exp3div4())));
                    if bool::from(a.is_square()) {
                        squares += 1;
                        assert!(eq(&a.sqrt().sqr(), &a));
                    } else {
                        nonsquares += 1;
                        assert!(!eq(&a.sqrt().sqr(), &a));
                    }
                }
                assert!(squares > 0 && nonsquares > 0, "sampler hit both classes");
            }

            #[test]
            fn encode_decode() {
                let mut rng = DetRng::new(&label("encode"));
                let n = <P as Prime>::FpEncodedBytes::USIZE;
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let bytes = a.encode();
                    assert_eq!(bytes.len(), n);
                    let back = F::decode(&bytes).expect("canonical bytes decode");
                    assert!(eq(&a, &back));
                    assert!(F::decode(&bytes[..n - 1]).is_none(), "short input rejected");
                    // decode_reduce agrees with decode on canonical input
                    assert!(eq(&F::decode_reduce(&bytes), &a));
                    // ... and reduces p + a to a
                    let mut wide = vec![0u8; n + 1];
                    let mut carry = 0u16;
                    for i in 0..n {
                        let s = bytes[i] as u16 + P::prime_le_bytes()[i] as u16 + carry;
                        wide[i] = s as u8;
                        carry = s >> 8;
                    }
                    wide[n] = carry as u8;
                    assert!(eq(&F::decode_reduce(&wide), &a));
                }
                // values >= p are rejected
                assert!(F::decode(P::prime_le_bytes()).is_none());
                let all_ff = vec![0xffu8; n];
                assert!(F::decode(&all_ff).is_none());
                // decode_reduce of an empty string is zero, and of long
                // strings is the integer mod p (checked against Horner form)
                assert!(bool::from(F::decode_reduce(&[]).ct_is_zero()));
                let mut long = vec![0u8; 3 * n + 7];
                rng.fill(&mut long);
                let mut acc = F::zero();
                let two56 = F::from_small(1 << 8);
                for &b in long.iter().rev() {
                    acc = acc.mul(&two56).add(&F::from_small(b as u64));
                }
                assert!(eq(&F::decode_reduce(&long), &acc));
            }

            #[test]
            fn constant_time_helpers() {
                let mut rng = DetRng::new(&label("ct"));
                for _ in 0..ITER {
                    let a = rng.random_fp::<P>();
                    let b = rng.random_fp::<P>();
                    assert!(eq(&F::select(&a, &b, Choice::from(0)), &a));
                    assert!(eq(&F::select(&a, &b, Choice::from(1)), &b));
                    let (mut x, mut y) = (a.clone(), b.clone());
                    x.cswap(&mut y, Choice::from(0));
                    assert!(eq(&x, &a) && eq(&y, &b));
                    x.cswap(&mut y, Choice::from(1));
                    assert!(eq(&x, &b) && eq(&y, &a));
                    assert!(eq(&a, &a));
                    assert_eq!(eq(&a, &b), eq(&b, &a));
                }
            }

            #[test]
            fn fp2_arithmetic() {
                let mut rng = DetRng::new(&label("fp2"));
                let one = F2::one();
                let i = F2::i_element();
                assert!(eq2(&i.sqr(), &one.neg()), "i^2 = -1");
                for _ in 0..ITER {
                    let a = rng.random_fp2::<P>();
                    let b = rng.random_fp2::<P>();
                    let c = rng.random_fp2::<P>();
                    assert!(eq2(&a.mul(&b), &b.mul(&a)));
                    assert!(eq2(&a.mul(&b).mul(&c), &a.mul(&b.mul(&c))));
                    assert!(eq2(&a.mul(&b.add(&c)), &a.mul(&b).add(&a.mul(&c))));
                    assert!(eq2(&a.sqr(), &a.mul(&a)));
                    assert!(eq2(&a.add_one(), &a.add(&one)));
                    assert!(eq2(&a.half().add(&a.half()), &a));
                    assert!(eq2(&a.mul_small(7), &a.mul(&F2::from_small(7))));
                    assert!(eq2(&a.mul(&a.conjugate()), &F2 { re: a.re.sqr().add(&a.im.sqr()), im: F::zero() }));
                    if !bool::from(a.ct_is_zero()) {
                        assert!(eq2(&a.mul(&a.inv()), &one));
                    }
                }
            }

            #[test]
            fn fp2_square_roots() {
                let mut rng = DetRng::new(&label("fp2 sqrt"));
                let mut squares = 0;
                for _ in 0..ITER {
                    let a = rng.random_fp2::<P>();
                    let sq = a.sqr();
                    assert!(bool::from(sq.is_square()));
                    let r = sq.sqrt();
                    assert!(eq2(&r.sqr(), &sq));
                    // the round-2 convention is the same root up to sign, with
                    // even real part (or zero real part and even imaginary part)
                    let c = sq.sqrt_canonical_even();
                    assert!(eq2(&c, &r) || eq2(&c, &r.neg()));
                    if bool::from(c.re.ct_is_zero()) {
                        assert_eq!(c.im.encode()[0] & 1, 0);
                    } else {
                        assert_eq!(c.re.encode()[0] & 1, 0);
                    }
                    let mut v = sq.clone();
                    assert!(bool::from(v.sqrt_verify()));
                    assert!(eq2(&v, &r));
                    if bool::from(a.is_square()) {
                        squares += 1;
                    } else {
                        let mut w = a.clone();
                        assert!(!bool::from(w.sqrt_verify()));
                    }
                }
                assert!(squares > 0);
                let x = F2 { re: F::one().neg(), im: F::zero() };
                assert!(bool::from(x.is_square()), "-1 is a square in Fp2");
                assert!(eq2(&x.sqrt().sqr(), &x));
            }

            #[test]
            fn fp2_helpers() {
                let mut rng = DetRng::new(&label("fp2 helpers"));
                let i = F2::i_element();
                for _ in 0..ITER {
                    let a = rng.random_fp2::<P>();
                    let b = rng.random_fp2::<P>();
                    assert!(eq2(&a.frob(), &a.conjugate()));
                    assert!(eq2(&a.frob().frob(), &a));
                    assert!(eq2(&a.mul_by_i(Choice::from(1)), &a.mul(&i)));
                    assert!(eq2(&a.mul_by_i(Choice::from(0)), &a.mul(&i.neg())));
                    // less_than: strict total order on encodings, im-major
                    let lt = bool::from(a.less_than(&b));
                    let gt = bool::from(b.less_than(&a));
                    assert!(!(lt && gt));
                    assert!(!bool::from(a.less_than(&a)));
                    let ea = a.encode();
                    let eb = b.encode();
                    let expect = ea.iter().rev().cmp(eb.iter().rev()) == core::cmp::Ordering::Less;
                    assert_eq!(lt, expect);
                    assert_eq!(lt || gt || eq2(&a, &b), true);
                }
                // im dominates re
                let lo = F2 { re: F::one().neg(), im: F::zero() };
                let hi = F2 { re: F::zero(), im: F::one() };
                assert!(bool::from(lo.less_than(&hi)));
                assert!(!bool::from(hi.less_than(&lo)));
            }

            #[test]
            fn fp2_batched_inverse_and_pow() {
                let mut rng = DetRng::new(&label("fp2 batch"));
                let one = F2::one();
                let mut xs: Vec<F2> = (0..9).map(|_| rng.random_fp2::<P>()).collect();
                let orig = xs.clone();
                let mut t1 = vec![F2::zero(); xs.len()];
                let mut t2 = vec![F2::zero(); xs.len()];
                F2::batched_inv(&mut xs, &mut t1, &mut t2);
                for (x, o) in xs.iter().zip(orig.iter()) {
                    assert!(eq2(&x.mul(o), &one));
                }
                let a = rng.random_fp2::<P>();
                assert!(eq2(&a.pow_vartime(&[0]), &one));
                assert!(eq2(&a.pow_vartime(&[1]), &a));
                assert!(eq2(&a.pow_vartime(&[5]), &a.sqr().sqr().mul(&a)));
                assert!(eq2(&a.pow_vartime(&[0, 1]), &{
                    let mut r = a.clone();
                    for _ in 0..64 {
                        r = r.sqr();
                    }
                    r
                }));
            }

            #[test]
            fn fp2_encode_decode() {
                let mut rng = DetRng::new(&label("fp2 encode"));
                let n = <P as Prime>::FpEncodedBytes::USIZE;
                for _ in 0..ITER {
                    let a = rng.random_fp2::<P>();
                    let bytes = a.encode();
                    assert_eq!(bytes.len(), 2 * n);
                    let back = F2::decode(&bytes).expect("canonical");
                    assert!(eq2(&a, &back));
                    assert!(F2::decode(&bytes[..2 * n - 1]).is_none());
                    let mut bad = bytes.to_vec();
                    bad[..n].copy_from_slice(P::prime_le_bytes());
                    assert!(F2::decode(&bad).is_none(), "re = p rejected");
                }
            }
        }
    };
}

fp_props!(p324_3, sqisign_verify::params::P324_3);
fp_props!(p500_27, sqisign_verify::params::P500_27);
fp_props!(p664_17, sqisign_verify::params::P664_17);
