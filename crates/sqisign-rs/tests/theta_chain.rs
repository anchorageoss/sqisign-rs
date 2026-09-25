//! Property tests for the `(2^n, 2^n)`-isogeny chains, with kernels built
//! from endomorphisms of `E0` (see `theta_common`).
//!
//! For a kernel `{([u] P, theta(P))}` the isogeny `F: E0 x E0 -> E3 x E4`
//! restricted to `E0 x 0` has components of degrees `u` and `v = 2^n - u`,
//! so Weil pairings of pushed points are known powers of the pairing on
//! `E0`. The dual chain, with kernel `F(E0[2^(n + 2)] x 0)`, must land back
//! on `E0 x E0` and again satisfy pairing identities. Malformed kernels
//! must be rejected.

mod theta_common;

use num_bigint::BigUint;
use sqisign_verify::ec::pairing::weil;
use sqisign_verify::ec::point::test_point_order_twof;
use sqisign_verify::ec::{EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::theta::chain::{theta_chain_compute_and_eval, ChainOutput, MAX_POINTS};
use sqisign_verify::theta::{
    ChainMode, ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints,
};
use theta_common::*;

/// Per-prime constants the tests need.
pub struct Params {
    pub f: u32,
    pub words: usize,
    pub p_bits: u32,
    pub gens: Gens,
}

fn e0_couple<L: FpBackend>() -> ThetaCoupleCurve<L> {
    ThetaCoupleCurve {
        e1: e0::<L>(),
        e2: e0::<L>(),
    }
}

fn j1728<L: FpBackend>() -> Fp2<L> {
    Fp2::from_small(1728)
}

fn eq2<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

fn pow<L: FpBackend>(base: &Fp2<L>, exp: &BigUint, f: u32, words: usize) -> Fp2<L> {
    let e = exp % pow2(f);
    base.pow_vartime(&big_to_limbs(&e, words))
}

/// Weil pairing of order `2^f` of a triple `(A, B, A - B)` on `curve`.
fn pairing<L: FpBackend>(
    f: u32,
    a: &EcPoint<L>,
    b: &EcPoint<L>,
    amb: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> Fp2<L> {
    let mut e = curve.clone();
    weil(f, a, b, amb, &mut e)
}

/// The images of `(P0, 0), (Q0, 0), (P0 - Q0, 0)` (first factor) or
/// `(0, P0), ...` (second factor) as a triple on `E3` or `E4`.
fn component<L: FpBackend>(imgs: &[ThetaCouplePoint<L>], first: bool) -> [EcPoint<L>; 3] {
    let pick = |c: &ThetaCouplePoint<L>| if first { c.p1.clone() } else { c.p2.clone() };
    [pick(&imgs[0]), pick(&imgs[1]), pick(&imgs[2])]
}

fn run_chain<L: FpBackend>(
    n: u32,
    e12: &ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    pts: &[ThetaCouplePoint<L>],
    mode: ChainMode,
) -> Option<ChainOutput<L>> {
    let mut dom = e12.clone();
    theta_chain_compute_and_eval(n as u16, &mut dom, ker, pts, mode)
}

/// Full check of one Kani kernel of length `n`.
fn check_length<L: FpBackend + PrimePrecomp>(prm: &Params, n: u32, seed: u32) {
    let p = prime::<L>();
    let (u, endo) = find_endomorphism(n, &p, &prm.gens, prm.f, seed);
    let v = pow2(n) - &u;
    let kk = kani_kernel::<L>(n, &u, &endo, prm.f, prm.words, &prm.gens);
    let e12 = e0_couple::<L>();
    let b0 = e0_basis::<L>();

    // both kernel generators have order 2^(n + 2)
    for t in [&kk.ker.t1, &kk.ker.t2] {
        assert!(test_point_order_twof(&t.p1, &e12.e1, (n + 2) as usize));
        assert!(test_point_order_twof(&t.p2, &e12.e2, (n + 2) as usize));
    }

    // points to push: basis of E0[2^f] on each factor, basis of E0[2^(n+2)] on the first
    let pts = [
        ThetaCouplePoint::on_e1(b0.p.clone()),
        ThetaCouplePoint::on_e1(b0.q.clone()),
        ThetaCouplePoint::on_e1(b0.pmq.clone()),
        ThetaCouplePoint::on_e2(b0.p.clone()),
        ThetaCouplePoint::on_e2(b0.q.clone()),
        ThetaCouplePoint::on_e2(b0.pmq.clone()),
        ThetaCouplePoint::on_e1(kk.basis_n2.p.clone()),
        ThetaCouplePoint::on_e1(kk.basis_n2.q.clone()),
    ];
    assert_eq!(pts.len(), MAX_POINTS);

    let out = run_chain(n, &e12, &kk.ker, &pts, ChainMode::Both).expect("chain of length {n}");
    let imgs: Vec<ThetaCouplePoint<L>> = out
        .images
        .iter()
        .map(|i| i.clone().expect("image"))
        .collect();
    let e3 = &out.codomain.e1;
    let e4 = &out.codomain.e2;
    assert!(e3.is_a24_computed_and_normalized || bool::from(e3.c.ct_is_one()));

    let e_0 = pairing(prm.f, &b0.p, &b0.q, &b0.pmq, &e12.e1);
    let e_u = pow(&e_0, &u, prm.f, prm.words);
    let e_v = pow(&e_0, &v, prm.f, prm.words);

    // images of the E0[2^f] bases have full order and pairings e0^u, e0^v
    for (first, curve) in [(true, e3), (false, e4)] {
        for offset in [0usize, 3] {
            let tri = component(&imgs[offset..offset + 3], first);
            for pt in &tri {
                assert!(
                    test_point_order_twof(pt, curve, prm.f as usize),
                    "image order, n = {n}"
                );
            }
            let e = pairing(prm.f, &tri[0], &tri[1], &tri[2], curve);
            assert!(eq2(&e, &e_u) || eq2(&e, &e_v), "image pairing, n = {n}");
        }
    }
    // the two factors of F restricted to E0 x 0 have degrees u and v (one each)
    {
        let a = pairing(prm.f, &imgs[0].p1, &imgs[1].p1, &imgs[2].p1, e3);
        let b = pairing(prm.f, &imgs[0].p2, &imgs[1].p2, &imgs[2].p2, e4);
        assert!((eq2(&a, &e_u) && eq2(&b, &e_v)) || (eq2(&a, &e_v) && eq2(&b, &e_u)));
    }

    // the other modes agree with mode Both on what they normalise
    {
        let ver = run_chain(n, &e12, &kk.ker, &pts[..6], ChainMode::Verify).expect("verify mode");
        assert!(curves_equal(&ver.codomain.e1, e3));
        for (a, b) in ver.images.iter().flatten().zip(imgs.iter()) {
            assert!(points_equal(&a.p1, &b.p1));
        }
        let m1 = run_chain(n, &e12, &kk.ker, &pts[..3], ChainMode::E1).expect("E1 mode");
        assert!(curves_equal(&m1.codomain.e1, e3));
        for (a, b) in m1.images.iter().flatten().zip(imgs.iter()) {
            assert!(points_equal(&a.p1, &b.p1));
        }
        let m2 = run_chain(n, &e12, &kk.ker, &pts[3..6], ChainMode::E2).expect("E2 mode");
        assert!(curves_equal(&m2.codomain.e2, e4));
        for (a, b) in m2.images.iter().flatten().zip(imgs[3..].iter()) {
            assert!(points_equal(&a.p2, &b.p2));
        }
    }

    // the dual chain: kernel F(E0[2^(n+2)] x 0), codomain E0 x E0
    {
        let dual_ker = ThetaKernelCouplePoints {
            t1: imgs[6].clone(),
            t2: imgs[7].clone(),
        };
        let dual_pts = [
            ThetaCouplePoint::on_e1(imgs[0].p1.clone()),
            ThetaCouplePoint::on_e1(imgs[1].p1.clone()),
            ThetaCouplePoint::on_e1(imgs[2].p1.clone()),
            ThetaCouplePoint::on_e2(imgs[0].p2.clone()),
            ThetaCouplePoint::on_e2(imgs[1].p2.clone()),
            ThetaCouplePoint::on_e2(imgs[2].p2.clone()),
        ];
        let back =
            run_chain(n, &out.codomain, &dual_ker, &dual_pts, ChainMode::Both).expect("dual chain");
        assert!(
            eq2(&back.codomain.e1.j_inv(), &j1728::<L>()),
            "dual codomain E1, n = {n}"
        );
        assert!(
            eq2(&back.codomain.e2.j_inv(), &j1728::<L>()),
            "dual codomain E2, n = {n}"
        );
        let bimgs: Vec<ThetaCouplePoint<L>> = back.images.iter().flatten().cloned().collect();
        assert_eq!(bimgs.len(), 6);
        let e_uu = pow(&e_0, &(&u * &u), prm.f, prm.words);
        let e_vv = pow(&e_0, &(&v * &v), prm.f, prm.words);
        let e_uv = pow(&e_0, &(&u * &v), prm.f, prm.words);
        for offset in [0usize, 3] {
            // F^(F(R, 0)) = ([deg alpha] R, gamma^ alpha R) up to the order of the
            // factors: one component pairs to e0^(deg^2), the other to e0^(uv)
            let c1 = component(&bimgs[offset..offset + 3], true);
            let c2 = component(&bimgs[offset..offset + 3], false);
            let a = pairing(prm.f, &c1[0], &c1[1], &c1[2], &back.codomain.e1);
            let b = pairing(prm.f, &c2[0], &c2[1], &c2[2], &back.codomain.e2);
            let is_sq = |x: &Fp2<L>| eq2(x, &e_uu) || eq2(x, &e_vv);
            assert!(
                (is_sq(&a) && eq2(&b, &e_uv)) || (eq2(&a, &e_uv) && is_sq(&b)),
                "dual pairing, n = {n}"
            );
        }
    }

    // negative: a kernel of the right size but wrong shape is not accepted
    {
        // [2] ker with length n - 1: a valid (2^(n-1), 2^(n-1))-isogeny whose
        // codomain is a Jacobian, not a product
        let half = ThetaKernelCouplePoints {
            t1: kk.ker.t1.double(&e12),
            t2: kk.ker.t2.double(&e12),
        };
        for mode in [ChainMode::Both, ChainMode::Verify] {
            assert!(
                run_chain(n - 1, &e12, &half, &[], mode).is_none(),
                "truncated kernel accepted"
            );
        }
        // generators of order 2^(n+1) for a chain of length n
        assert!(
            run_chain(n, &e12, &half, &[], ChainMode::Verify).is_none(),
            "short kernel accepted"
        );
    }
    // negative: a pushed point with two non-zero components
    let bad = [ThetaCouplePoint::new(b0.p.clone(), b0.q.clone())];
    assert!(run_chain(n, &e12, &kk.ker, &bad, ChainMode::Both).is_none());
}

/// Rejections that need no valid kernel.
fn check_rejections<L: FpBackend + PrimePrecomp>(prm: &Params) {
    let e12 = e0_couple::<L>();
    let b0 = e0_basis::<L>();
    let n = prm.f - 2;
    let p = prime::<L>();
    let (u, endo) = find_endomorphism(n, &p, &prm.gens, prm.f, 3);
    let kk = kani_kernel::<L>(n, &u, &endo, prm.f, prm.words, &prm.gens);
    // lengths 0 and 1 are not supported
    assert!(run_chain(0, &e12, &kk.ker, &[], ChainMode::Both).is_none());
    let (u1, endo1) = find_endomorphism(n_min(prm.p_bits), &p, &prm.gens, prm.f, 0);
    let k1 = kani_kernel::<L>(n_min(prm.p_bits), &u1, &endo1, prm.f, prm.words, &prm.gens);
    let short = ThetaKernelCouplePoints {
        t1: k1.ker.t1.double_iter((n_min(prm.p_bits) - 1) as u16, &e12),
        t2: k1.ker.t2.double_iter((n_min(prm.p_bits) - 1) as u16, &e12),
    };
    assert!(run_chain(1, &e12, &short, &[], ChainMode::Both).is_none());
    // too many points
    let many: Vec<ThetaCouplePoint<L>> = (0..MAX_POINTS + 1)
        .map(|_| ThetaCouplePoint::on_e1(b0.p.clone()))
        .collect();
    assert!(run_chain(n, &e12, &kk.ker, &many, ChainMode::Both).is_none());
    // a product kernel <(P', 0), (0, Q')> is not a gluing kernel
    let prod = ThetaKernelCouplePoints {
        t1: ThetaCouplePoint::on_e1(kk.basis_n2.p.clone()),
        t2: ThetaCouplePoint::on_e2(kk.basis_n2.q.clone()),
    };
    for mode in [ChainMode::Both, ChainMode::Verify] {
        assert!(
            run_chain(n, &e12, &prod, &[], mode).is_none(),
            "product kernel accepted"
        );
    }
    // a non-isotropic kernel: (P', 0) and (Q', Q'), both order 2^(n+2)
    let noniso = ThetaKernelCouplePoints {
        t1: ThetaCouplePoint::on_e1(kk.basis_n2.p.clone()),
        t2: ThetaCouplePoint::new(kk.basis_n2.q.clone(), kk.basis_n2.q.clone()),
    };
    assert!(
        run_chain(n, &e12, &noniso, &[], ChainMode::Verify).is_none(),
        "non-isotropic kernel accepted"
    );
    // identity components
    let zero = ThetaKernelCouplePoints {
        t1: ThetaCouplePoint::new(EcPoint::identity(), EcPoint::identity()),
        t2: kk.ker.t2.clone(),
    };
    assert!(run_chain(n, &e12, &zero, &[], ChainMode::Verify).is_none());
}

/// Smallest length reachable by the endomorphism construction.
fn n_min(p_bits: u32) -> u32 {
    p_bits / 2 + 2
}

macro_rules! chain_tests {
    ($modname:ident, $lp:path, $vp:path, $sp:path $(, $cfg:meta)?) => {
        $(#[cfg($cfg)])?
        mod $modname {
            use super::*;
            use $lp as L;
            use $sp as s;
            use $vp as v;

            fn params() -> Params {
                Params {
                    f: v::TWO_ADIC_EXPONENT,
                    words: v::NWORDS_ORDER,
                    p_bits: s::P_BITS,
                    gens: Gens {
                        g2: mat_from_bytes(&s::ACTION_GEN2),
                        g3: mat_from_bytes(&s::ACTION_GEN3),
                        g4: mat_from_bytes(&s::ACTION_GEN4),
                    },
                }
            }

            #[test]
            fn length_min() {
                let p = params();
                check_length::<L>(&p, n_min(p.p_bits), 0);
            }

            #[test]
            fn length_mid() {
                let p = params();
                let n = (n_min(p.p_bits) + p.f - 2) / 2;
                check_length::<L>(&p, n, 1);
            }

            #[test]
            fn length_max() {
                let p = params();
                check_length::<L>(&p, p.f - 2, 2);
            }

            #[test]
            fn rejections() {
                check_rejections::<L>(&params());
            }
        }
    };
}

chain_tests!(
    p324_3,
    sqisign_verify::params::P324_3,
    sqisign_verify::precomp::p324_3,
    sqisign_rs::precomp::p324_3
);
chain_tests!(
    p500_27,
    sqisign_verify::params::P500_27,
    sqisign_verify::precomp::p500_27,
    sqisign_rs::precomp::p500_27
);
chain_tests!(
    p664_17,
    sqisign_verify::params::P664_17,
    sqisign_verify::precomp::p664_17,
    sqisign_rs::precomp::p664_17
);
