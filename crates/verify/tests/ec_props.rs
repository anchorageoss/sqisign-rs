//! Property tests for the curve layer, at every prime: ladders against each
//! other, deterministic bases and hints, pairings (bilinearity, order,
//! discrete logarithms), isogeny chains landing on canonical models, and the
//! `E[p + 1]` basis. Curves other than `E0` are reached through 4-isogeny
//! chains from `E0`, so those are exercised too.

mod common;

macro_rules! ec_props {
    ($modname:ident, $prime:ty $(, $cfg:meta)?) => {
        $(#[cfg($cfg)])?
        mod $modname {
            use crate::common::DetRng;
            use sqisign_verify::ec::basis::*;
            use sqisign_verify::ec::isogeny::*;
            use sqisign_verify::ec::normalize::*;
            use sqisign_verify::ec::pairing::*;
            use sqisign_verify::ec::point::*;
            use sqisign_verify::ec::{EcBasis, EcCurve, EcIsogEven, EcPoint, MAX_ORDER_WORDS};
            use sqisign_verify::fp::Fp2;
            use sqisign_verify::params::Prime;
            use subtle::Choice;

            type P = $prime;
            const F: u32 = <P as Prime>::TWO_ADIC_EXPONENT;
            const W: usize = <P as Prime>::ORDER_WORDS;

            fn label(s: &str) -> Vec<u8> {
                format!("ec/{}/{}", stringify!($modname), s).into_bytes()
            }

            fn scalar(rng: &mut DetRng, bits: usize) -> [u64; MAX_ORDER_WORDS] {
                let mut s = [0u64; MAX_ORDER_WORDS];
                let mut buf = [0u8; 8 * MAX_ORDER_WORDS];
                rng.fill(&mut buf);
                for i in 0..W {
                    s[i] = u64::from_le_bytes(buf[8 * i..8 * i + 8].try_into().unwrap());
                }
                for i in 0..W {
                    let lo = i * 64;
                    if lo >= bits {
                        s[i] = 0;
                    } else if lo + 64 > bits {
                        s[i] &= (1u64 << (bits - lo)) - 1;
                    }
                }
                s
            }

            fn eq(a: &EcPoint<P>, b: &EcPoint<P>) -> bool {
                bool::from(a.ct_equal(b))
            }

            fn e0() -> EcCurve<P> {
                let mut e = EcCurve::<P>::init();
                e.normalize_curve_and_a24();
                e
            }

            /// A curve reached from E0 by a chain of even length, with its
            /// full-torsion basis and hint.
            fn random_curve(rng: &mut DetRng) -> (EcCurve<P>, EcBasis<P>, u8) {
                let mut e = e0();
                let b0 = ec_curve_to_basis_2f_from_hint(&mut e, F, 0).expect("basis");
                let length = F - (F % 2) - 2;
                let m = scalar(rng, F as usize);
                let k = ec_ladder3pt(&m[..W], &b0.p, &b0.q, &b0.pmq, &e).expect("ladder");
                let k = ec_dbl_iter(&k, (F - length) as usize, &mut e);
                let phi = EcIsogEven {
                    curve: e.clone(),
                    kernel: k,
                    length,
                };
                let mut e1 = iso_isogeny_2chain(&phi).expect("kernel has the right order");
                let (b1, hint) = ec_curve_to_basis_2f_to_hint(&mut e1, F, 0).expect("basis");
                (e1, b1, hint)
            }

            #[test]
            fn e0_precomputed_basis_is_a_basis() {
                let mut e = e0();
                let b = ec_curve_to_basis_2f_from_hint(&mut e, F, 0).expect("basis");
                assert!(test_basis_order_twof(&b, &e, F as usize));
                let (b2, hint) = ec_curve_to_basis_2f_to_hint(&mut e, F, 0).expect("basis");
                assert_eq!(hint, 0);
                assert!(eq(&b.p, &b2.p) && eq(&b.q, &b2.q) && eq(&b.pmq, &b2.pmq));
                // the difference point is P - Q or P + Q: either way P + Q via xadd has full order
                let ppq = xadd(&b.p, &b.q, &b.pmq);
                assert!(test_point_order_twof(&ppq, &e, F as usize));
                // Q is above (0, 0)
                let q2 = ec_dbl_iter(&b.q, (F - 1) as usize, &mut e);
                assert!(bool::from(q2.x.ct_is_zero()));
                // Weil pairing of the basis has full order
                let w = weil(F, &b.p, &b.q, &b.pmq, &mut e);
                let mut t = w.clone();
                for _ in 0..F - 1 {
                    t = t.sqr();
                }
                assert!(bool::from(t.ct_equal(&Fp2::<P>::one().neg())), "e(P,Q)^(2^(f-1)) = -1");
            }

            #[test]
            fn ladders_agree() {
                let mut rng = DetRng::new(&label("ladders"));
                let mut e = e0();
                let b = ec_curve_to_basis_2f_from_hint(&mut e, F, 0).expect("basis");
                for _ in 0..3 {
                    let a = scalar(&mut rng, F as usize);
                    let c = scalar(&mut rng, F as usize);
                    // [a]P + [c]Q two ways
                    let r1 = ec_biscalar_mul(&a[..W], &c[..W], F as usize, &b, &e).expect("biscalar");
                    let r2 = ec_biscalar_mul_verif(&a[..W], &c[..W], F as usize, &b, &mut e.clone()).expect("verif");
                    assert!(eq(&r1, &r2), "biscalar constant-time vs variable-time");
                    // [a]P via ladder, and ([a]P) + [c]Q via ladder3pt needs x([a]P - Q): use biscalar with c = 0
                    let ap = ec_mul(&b.p, &a[..W], F as usize, &mut e);
                    let ap2 = ec_biscalar_mul(&a[..W], &[0u64; MAX_ORDER_WORDS][..W], F as usize, &b, &e).expect("biscalar");
                    assert!(eq(&ap, &ap2), "ladder vs biscalar");
                    // P + [c]Q via ladder3pt vs biscalar(1, c)
                    let one = {
                        let mut s = [0u64; MAX_ORDER_WORDS];
                        s[0] = 1;
                        s
                    };
                    let l3 = ec_ladder3pt(&c[..W], &b.p, &b.q, &b.pmq, &e).expect("ladder3pt");
                    let l3b = ec_biscalar_mul(&one[..W], &c[..W], F as usize, &b, &e).expect("biscalar");
                    assert!(eq(&l3, &l3b), "ladder3pt vs biscalar");
                    // doublings agree
                    let d1 = ec_dbl(&b.p, &e);
                    let d2 = xdbl_e0(&b.p);
                    let d3 = xdbl_a24(&b.p, &e.a24, true);
                    assert!(eq(&d1, &d2) && eq(&d1, &d3));
                }
            }

            #[test]
            fn bary_coordinates() {
                let mut rng = DetRng::new(&label("bary"));
                let (e, b, _) = random_curve(&mut rng);
                let ppq = xadd(&b.p, &b.q, &b.pmq);
                let uvw = ec_points_to_bary_coordinates(&b.p, &b.q, &b.pmq);
                let sum = EcPoint::new(uvw.u.sub(&uvw.v), uvw.w.clone());
                let diff = EcPoint::new(uvw.u.add(&uvw.v), uvw.w.clone());
                assert!(eq(&sum, &ppq), "u - v : w is P + Q");
                assert!(eq(&diff, &b.pmq), "u + v : w is P - Q");
                let _ = e;
            }

            #[test]
            fn chain_lands_on_canonical_model_with_basis_and_hint() {
                let mut rng = DetRng::new(&label("chain"));
                for _ in 0..2 {
                    let (mut e1, b1, hint) = random_curve(&mut rng);
                    assert!(bool::from(e1.c.ct_is_one()));
                    assert!(test_basis_order_twof(&b1, &e1, F as usize));
                    let q2 = ec_dbl_iter(&b1.q, (F - 1) as usize, &mut e1);
                    assert!(bool::from(q2.x.ct_is_zero()), "Q above (0, 0)");
                    // hint round trip
                    let b2 = ec_curve_to_basis_2f_from_hint(&mut e1, F, hint).expect("basis");
                    assert!(eq(&b1.p, &b2.p) && eq(&b1.q, &b2.q) && eq(&b1.pmq, &b2.pmq));
                    // e-normalised variant agrees with the e-torsion recomputation
                    let e_small = F / 2;
                    let (b3, hint3) = ec_curve_to_basis_2f_to_hint(&mut e1, F, e_small).expect("basis");
                    assert_eq!(hint3, hint);
                    let b3s = ec_dbl_iter_basis(&b3, (F - e_small) as usize, &mut e1);
                    let b4 = ec_curve_to_basis_2f_from_hint(&mut e1, e_small, hint).expect("basis");
                    assert!(eq(&b3s.p, &b4.p) && eq(&b3s.q, &b4.q) && eq(&b3s.pmq, &b4.pmq));
                    // canonical model: the A coefficient is the maximum of the six
                    // any order-4 point not above (0, 0) is a theta null point of
                    // the curve and reproduces the canonical coefficient
                    let four = ec_dbl_iter(&b1.p, (F - 2) as usize, &mut e1);
                    let (e_can, _) = ec_theta_to_montgomery(&four).expect("valid theta point");
                    assert!(bool::from(e_can.a.ct_equal(&e1.a)), "already canonical");
                    assert!(bool::from(e_can.j_inv().ct_equal(&e1.j_inv())));
                    // an order-4 point above (0, 0) is rejected or gives another model
                    let four_q = ec_dbl_iter(&b1.q, (F - 2) as usize, &mut e1);
                    if let Some((e_bad, _)) = ec_theta_to_montgomery(&four_q) {
                        assert!(!bool::from(e_bad.j_inv().ct_equal(&e1.j_inv())) || bool::from(e_bad.a.ct_equal(&e1.a)));
                    }
                    let _ = ec_apply_isomorphism::<P>;
                }
            }

            #[test]
            fn chain_rejects_bad_kernels() {
                let mut rng = DetRng::new(&label("bad kernel"));
                let mut e = e0();
                let b0 = ec_curve_to_basis_2f_from_hint(&mut e, F, 0).expect("basis");
                let length = F - (F % 2) - 2;
                let m = scalar(&mut rng, F as usize);
                let k = ec_ladder3pt(&m[..W], &b0.p, &b0.q, &b0.pmq, &e).expect("ladder");
                // wrong order: one doubling too many
                let k_short = ec_dbl_iter(&k, (F - length + 1) as usize, &mut e);
                assert!(iso_isogeny_2chain(&EcIsogEven { curve: e.clone(), kernel: k_short, length }).is_none());
                // kernel above (0, 0): Q itself
                let kq = ec_dbl_iter(&b0.q, (F - length) as usize, &mut e);
                assert!(iso_isogeny_2chain(&EcIsogEven { curve: e.clone(), kernel: kq, length }).is_none());
                // odd length
                let k_ok = ec_dbl_iter(&k, (F - length) as usize, &mut e);
                assert!(iso_isogeny_2chain(&EcIsogEven { curve: e.clone(), kernel: k_ok, length: length - 1 }).is_none());
            }

            #[test]
            fn pairings_and_dlog() {
                let mut rng = DetRng::new(&label("pairing"));
                let (mut e1, b1, _) = random_curve(&mut rng);
                // bilinearity of Weil: e([a]P, Q) = e(P, Q)^a
                let w = weil(F, &b1.p, &b1.q, &b1.pmq, &mut e1);
                let a = scalar(&mut rng, 64);
                let ap = ec_mul(&b1.p, &a[..W], 64, &mut e1);
                // need x([a]P - Q): [a]P - Q = [a]P + [-1]Q; use biscalar with (a, 2^F - 1)
                let mut minus_one = [u64::MAX; MAX_ORDER_WORDS];
                for i in 0..W {
                    let lo = i * 64;
                    if lo >= F as usize {
                        minus_one[i] = 0;
                    } else if lo + 64 > F as usize {
                        minus_one[i] &= (1u64 << (F as usize - lo)) - 1;
                    }
                }
                let apmq = ec_biscalar_mul(&a[..W], &minus_one[..W], F as usize, &b1, &e1).expect("biscalar");
                let w_a = weil(F, &ap, &b1.q, &apmq, &mut e1);
                assert!(bool::from(w_a.ct_equal(&w.pow_vartime(&a[..1]))), "e([a]P, Q) = e(P, Q)^a");
                // reduced Tate of the basis has full order 2^F
                let t = reduced_tate(F, &b1.p, &b1.q, &b1.pmq, &mut e1);
                let mut t2 = t.clone();
                for _ in 0..F - 1 {
                    t2 = t2.sqr();
                }
                assert!(bool::from(t2.ct_equal(&Fp2::<P>::one().neg())));
                // dlog: R = [r1]P + [r2]Q, S = [s1]P + [s2]Q, full and partial torsion.
                // x-only points fix (R, S) only up to a common sign, so the
                // dlog is defined up to negating all four scalars modulo 2^e;
                // which sign comes out depends on the canonical square root
                // in the shared-difference computation (the toy prime picks
                // the other one, P19 B1). Both answers are accepted.
                for &e_dl in &[F, F - 7] {
                    let r1 = scalar(&mut rng, e_dl as usize);
                    let r2 = scalar(&mut rng, e_dl as usize);
                    let s1 = scalar(&mut rng, e_dl as usize);
                    let s2 = scalar(&mut rng, e_dl as usize);
                    let mk = |x: &[u64], y: &[u64], e1: &mut EcCurve<P>| {
                        let pt = ec_biscalar_mul(x, y, e_dl as usize, &b1, e1).expect("in span");
                        ec_dbl_iter(&pt, (F - e_dl) as usize, e1)
                    };
                    let r = mk(&r1[..W], &r2[..W], &mut e1);
                    let s = mk(&s1[..W], &s2[..W], &mut e1);
                    // x(R - S) = [2^(F-e)] ([r1 - s1] P + [r2 - s2] Q)
                    let mut d1 = [0u64; MAX_ORDER_WORDS];
                    let mut d2 = [0u64; MAX_ORDER_WORDS];
                    let mut borrow = 0u64;
                    for i in 0..W {
                        let (x, o1) = r1[i].overflowing_sub(s1[i]);
                        let (x, o2) = x.overflowing_sub(borrow);
                        d1[i] = x;
                        borrow = (o1 | o2) as u64;
                    }
                    borrow = 0;
                    for i in 0..W {
                        let (x, o1) = r2[i].overflowing_sub(s2[i]);
                        let (x, o2) = x.overflowing_sub(borrow);
                        d2[i] = x;
                        borrow = (o1 | o2) as u64;
                    }
                    let rms = mk(&d1[..W], &d2[..W], &mut e1);
                    let rs = EcBasis::new(r, s, rms);
                    if !test_basis_order_twof(&rs, &e1, e_dl as usize) {
                        continue; // not a basis of E[2^e]; skip this draw
                    }
                    let out = ec_dlog_2_tate(&b1, &rs, &mut e1, e_dl).expect("in span");
                    let neg = |x: &[u64]| -> Vec<u64> {
                        // 2^e - x mod 2^e over W words
                        let mut v = vec![0u64; W];
                        let mut borrow = 0u64;
                        for i in 0..W {
                            let (d1, b1_) = 0u64.overflowing_sub(x[i]);
                            let (d2, b2_) = d1.overflowing_sub(borrow);
                            v[i] = d2;
                            borrow = (b1_ | b2_) as u64;
                        }
                        for (i, w) in v.iter_mut().enumerate() {
                            let lo = 64 * i as u32;
                            if lo >= e_dl {
                                *w = 0;
                            } else if lo + 64 > e_dl {
                                *w &= (1u64 << (e_dl - lo)) - 1;
                            }
                        }
                        v
                    };
                    let got: Vec<Vec<u64>> = (0..4).map(|k| out[k][..W].to_vec()).collect();
                    let want: Vec<Vec<u64>> = [&r1[..W], &r2[..W], &s1[..W], &s2[..W]].iter().map(|x| x.to_vec()).collect();
                    let want_neg: Vec<Vec<u64>> = want.iter().map(|x| neg(x)).collect();
                    assert!(
                        got == want || got == want_neg,
                        "dlog is (r1, r2, s1, s2) up to a common sign: got {got:?}, want {want:?}"
                    );
                }
            }

            #[test]
            fn full_rational_basis() {
                let mut rng = DetRng::new(&label("p+1"));
                for curve in [e0(), random_curve(&mut rng).0] {
                    let mut e = curve;
                    let b = ec_basis_p_plus_one(&mut e).expect("basis found");
                    // [p + 1] R = 0 and R, T of full order: check via the 2^f and cofactor parts
                    let cof = [<P as Prime>::COFACTOR];
                    for pt in [&b.p, &b.q] {
                        let t = ec_mul(pt, &cof, <P as Prime>::COFACTOR_BITLENGTH, &mut e);
                        assert!(test_point_order_twof(&t, &e, F as usize), "2-part of full order");
                        let u = ec_dbl_iter(pt, F as usize, &mut e);
                        let v = ec_mul(&u, &cof, <P as Prime>::COFACTOR_BITLENGTH, &mut e);
                        assert!(bool::from(v.is_zero()), "[p + 1] R = 0");
                        assert!(!bool::from(u.is_zero()), "odd part non-trivial");
                    }
                    // the difference point is consistent: R + T via xadd has full 2-power order too
                    let sum = xadd(&b.p, &b.q, &b.pmq);
                    let t = ec_mul(&sum, &cof, <P as Prime>::COFACTOR_BITLENGTH, &mut e);
                    assert!(test_point_order_twof(&t, &e, F as usize));
                    let _ = Choice::from(1);
                }
            }
        }
    };
}

ec_props!(p324_3, sqisign_verify::params::P324_3);
ec_props!(p500_27, sqisign_verify::params::P500_27);
ec_props!(p664_17, sqisign_verify::params::P664_17);
