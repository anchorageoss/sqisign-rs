//! The `E0`-specific derivations: the deterministic basis of `E0[2^f]`
//! (spec Appendix B, Algorithm B.1), the two-dimensional discrete logarithm
//! in `E0[2^f]`, point halving into `E0(Fp4)`, and the action matrices of
//! `i`, `j`, `k`, `(i + j)/2` and `(1 + k)/2`.

use crate::curve::{Affine, Curve, Jac};
use crate::field::{fp2_eq, fp2_lex_less, fp2_small, prime, FieldOps, Fp2Ops, Fp4, Fp4Ctx};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::params::Prime;

type P2<L> = Option<Jac<Fp2<L>>>;

pub struct E0Basis<L: FpBackend> {
    pub p: Affine<Fp2<L>>,
    pub q: Affine<Fp2<L>>,
    pub f: u64,
}

/// A 2x2 matrix over `Z/2^f` in the reference's layout:
/// `[[a(P0), a(Q0)], [b(P0), b(Q0)]]` with `alpha(R) = a(R) P0 + b(R) Q0`.
pub type Mat = [[BigUint; 2]; 2];

/// Algorithm B.1 / `torsion_basis.py`: the canonical basis of `E0[2^f]`.
pub fn e0_basis<L: FpBackend>() -> E0Basis<L> {
    let ops = Fp2Ops::<L>::new();
    let c = Curve::new(&ops);
    let f = L::TWO_ADIC_EXPONENT as u64;
    let cof = L::COFACTOR;
    let two_f = BigUint::one() << f;
    let two_fm1 = BigUint::one() << (f - 1);
    let zero_zero: P2<L> = c.lift(&Some((Fp2::zero(), Fp2::zero())));

    // Points of exact order 2^f with x = 1 + m i, in order of increasing m.
    let mut candidates = (1u64..).filter_map(|m| {
        let x = fp2_small::<L>(1, m);
        let rhs = x.sqr().mul(&x).add(&x);
        if !bool::from(rhs.is_square()) {
            return None;
        }
        let y = rhs.sqrt();
        let y_neg = y.neg();
        let y = if fp2_lex_less(&y_neg, &y) { y_neg } else { y };
        debug_assert!(c.on_curve(&x, &y));
        let r = c.lift(&Some((x, y)));
        // 2^f divides the order iff [c 2^(f-1)] R != O
        c.mul(&r, &(BigUint::from(cof) << (f - 1)))?;
        // odd part d of the order: the smallest divisor d of c with [d 2^f] R = O
        let d = (1..=cof)
            .filter(|d| cof % d == 0)
            .find(|d| c.mul(&r, &(BigUint::from(*d) << f)).is_none())
            .expect("order divides c 2^f");
        Some(c.mul(&r, &BigUint::from(d)))
    });

    let p = candidates.next().expect("first basis point");
    let p2 = c.mul(&p, &two_fm1);
    assert!(p2.is_some());
    let q = candidates
        .find(|q| {
            // (P, Q) is a basis iff their 2^(f-1) multiples are distinct
            // non-zero points of E0[2], i.e. e(P, Q) has full order 2^f.
            let q2 = c.mul(q, &two_fm1);
            q2.is_some() && !c.eq(&p2, &q2)
        })
        .expect("second basis point");

    // Arrange for [2^(f-1)] Q = (0, 0).
    let q2 = c.mul(&q, &two_fm1);
    let (p, q) = if c.eq(&q2, &zero_zero) {
        (p, q)
    } else if c.eq(&p2, &zero_zero) {
        (q, p)
    } else {
        (p.clone(), c.add(&q, &p))
    };
    debug_assert!(c.eq(&c.mul(&q, &two_fm1), &zero_zero));
    debug_assert!(c.mul(&p, &two_f).is_none() && c.mul(&q, &two_f).is_none());

    E0Basis {
        p: c.to_affine(&p),
        q: c.to_affine(&q),
        f,
    }
}

/// x-coordinate of `P - Q`.
pub fn x_of_difference<L: FpBackend>(b: &E0Basis<L>) -> Fp2<L> {
    let ops = Fp2Ops::<L>::new();
    let c = Curve::new(&ops);
    let p = c.lift(&b.p);
    let q = c.lift(&b.q);
    c.to_affine(&c.sub(&p, &q)).expect("P != Q").0
}

/// Two-dimensional Pohlig-Hellman: `(a, b)` with `S = a P + b Q` in
/// `E0[2^f]`, bit by bit from the least significant, comparing
/// `[2^(f-1-k)] (S - a_k P - b_k Q)` against the three non-zero points of
/// `E0[2]`.
pub fn dlog2<L: FpBackend>(b: &E0Basis<L>, s: &Affine<Fp2<L>>) -> (BigUint, BigUint) {
    let ops = Fp2Ops::<L>::new();
    let c = Curve::new(&ops);
    let f = b.f;
    let p = c.lift(&b.p);
    let q = c.lift(&b.q);
    let mut s = c.lift(s);
    // [2^k] P and [2^k] Q for k < f
    let mut pk = Vec::with_capacity(f as usize);
    let mut qk = Vec::with_capacity(f as usize);
    let (mut tp, mut tq) = (p.clone(), q.clone());
    for _ in 0..f {
        pk.push(tp.clone());
        qk.push(tq.clone());
        tp = c.double(&tp);
        tq = c.double(&tq);
    }
    let p2 = &pk[f as usize - 1];
    let q2 = &qk[f as usize - 1];
    let pq2 = c.add(p2, q2);
    let mut a = BigUint::zero();
    let mut bb = BigUint::zero();
    for k in 0..f {
        let u = c.mul_pow2(&s, f - 1 - k);
        let (e1, e2) = if u.is_none() {
            (false, false)
        } else if c.eq(&u, p2) {
            (true, false)
        } else if c.eq(&u, q2) {
            (false, true)
        } else if c.eq(&u, &pq2) {
            (true, true)
        } else {
            panic!("point not in the span of the basis");
        };
        if e1 {
            a.set_bit(k, true);
            s = c.sub(&s, &pk[k as usize]);
        }
        if e2 {
            bb.set_bit(k, true);
            s = c.sub(&s, &qk[k as usize]);
        }
    }
    assert!(s.is_none(), "residual after Pohlig-Hellman");
    (a, bb)
}

/// The endomorphism `i: (x, y) -> (-x, i y)` on an `Fp2` point.
pub fn endo_i<L: FpBackend>(pt: &Affine<Fp2<L>>) -> Affine<Fp2<L>> {
    pt.as_ref()
        .map(|(x, y)| (x.neg(), y.mul(&Fp2::<L>::i_element())))
}

/// The Frobenius endomorphism `j: (x, y) -> (x^p, y^p)` on an `Fp2` point.
pub fn endo_j<L: FpBackend>(pt: &Affine<Fp2<L>>) -> Affine<Fp2<L>> {
    pt.as_ref().map(|(x, y)| (x.conjugate(), y.conjugate()))
}

/// A half `H` of an `Fp2`-rational point of `E0`, in `E0(Fp4)`. With
/// `x^3 + x = x (x - i)(x + i)` and `r_1^2 = x`, `r_2^2 = x - i`,
/// `r_3^2 = x + i`, `r_1 r_2 r_3 = y`, one half has
/// `x_H = x + r_1 r_2 + r_2 r_3 + r_3 r_1`. Which of the four halves comes
/// out does not matter to the callers.
pub fn halve<L: FpBackend>(ctx: &Fp4Ctx<L>, pt: &Affine<Fp2<L>>) -> Affine<Fp4<L>> {
    let (x, y) = pt.as_ref().expect("finite point");
    let c = Curve::new(ctx);
    let i = Fp2::<L>::i_element();
    let r1 = ctx.sqrt(&ctx.lift(x)).expect("x is a square in Fp4");
    let r2 = ctx
        .sqrt(&ctx.lift(&x.sub(&i)))
        .expect("x - i is a square in Fp4");
    let mut r3 = ctx
        .sqrt(&ctx.lift(&x.add(&i)))
        .expect("x + i is a square in Fp4");
    let prod = ctx.mul(&ctx.mul(&r1, &r2), &r3);
    if !ctx.eq(&prod, &ctx.lift(y)) {
        r3 = ctx.neg(&r3);
        let prod = ctx.mul(&ctx.mul(&r1, &r2), &r3);
        assert!(ctx.eq(&prod, &ctx.lift(y)), "r1 r2 r3 = y");
    }
    let xh = ctx.add(
        &ctx.add(
            &ctx.add(&ctx.lift(x), &ctx.mul(&r1, &r2)),
            &ctx.mul(&r2, &r3),
        ),
        &ctx.mul(&r3, &r1),
    );
    let rhs = ctx.add(&ctx.mul(&ctx.sqr(&xh), &xh), &xh);
    let mut yh = ctx.sqrt(&rhs).expect("x_H lifts to E0(Fp4)");
    let target = c.lift(&Some((ctx.lift(x), ctx.lift(y))));
    let twice = c.double(&c.lift(&Some((xh.clone(), yh.clone()))));
    if !c.eq(&twice, &target) {
        yh = ctx.neg(&yh);
        let twice = c.double(&c.lift(&Some((xh.clone(), yh.clone()))));
        assert!(c.eq(&twice, &target), "2H = P");
    }
    Some((xh, yh))
}

/// `(i + j)(H)` and `(1 + k)(H)` for a point `H` in `E0(Fp4)`, brought back
/// to `Fp2`. Both endomorphisms kill `E0[2]`, so for `H` a half of `R` the
/// results are `((i + j)/2)(R)` and `((1 + k)/2)(R)`.
pub fn halved_generators<L: FpBackend>(
    ctx: &Fp4Ctx<L>,
    h: &Affine<Fp4<L>>,
) -> (Affine<Fp2<L>>, Affine<Fp2<L>>) {
    let c = Curve::new(ctx);
    let (xh, yh) = h.as_ref().expect("finite");
    let i4 = ctx.lift(&Fp2::<L>::i_element());
    let ih = Some((ctx.neg(xh), ctx.mul(yh, &i4)));
    let jh = Some((ctx.frobenius(xh), ctx.frobenius(yh)));
    let (xj, yj) = jh.as_ref().expect("finite");
    let kh = Some((ctx.neg(xj), ctx.mul(yj, &i4)));
    let sum_ij = c.to_affine(&c.add(&c.lift(&ih), &c.lift(&jh)));
    let sum_1k = c.to_affine(&c.add(&c.lift(h), &c.lift(&kh)));
    let descend = |pt: Affine<Fp4<L>>| -> Affine<Fp2<L>> {
        pt.map(|(x, y)| {
            (
                ctx.descend(&x).expect("result is Fp2-rational"),
                ctx.descend(&y).expect("result is Fp2-rational"),
            )
        })
    };
    (descend(sum_ij), descend(sum_1k))
}

/// All six action matrices for the basis, in the reference's layout.
pub struct Actions {
    pub i: Mat,
    pub j: Mat,
    pub k: Mat,
    pub gen2: Mat,
    pub gen3: Mat,
    pub gen4: Mat,
}

fn matrix<L: FpBackend>(b: &E0Basis<L>, on_p: &Affine<Fp2<L>>, on_q: &Affine<Fp2<L>>) -> Mat {
    let (ap, bp) = dlog2(b, on_p);
    let (aq, bq) = dlog2(b, on_q);
    [[ap, aq], [bp, bq]]
}

pub fn actions<L: FpBackend>(b: &E0Basis<L>) -> Actions {
    let ctx = Fp4Ctx::<L>::new();
    let i = matrix(b, &endo_i::<L>(&b.p), &endo_i::<L>(&b.q));
    let j = matrix(b, &endo_j::<L>(&b.p), &endo_j::<L>(&b.q));
    let kp = endo_i::<L>(&endo_j::<L>(&b.p));
    let kq = endo_i::<L>(&endo_j::<L>(&b.q));
    let k = matrix(b, &kp, &kq);
    let hp = halve(&ctx, &b.p);
    let hq = halve(&ctx, &b.q);
    let (ij_p, k1_p) = halved_generators(&ctx, &hp);
    let (ij_q, k1_q) = halved_generators(&ctx, &hq);
    let gen3 = matrix(b, &ij_p, &ij_q);
    let gen4 = matrix(b, &k1_p, &k1_q);
    Actions {
        gen2: i.clone(),
        i,
        j,
        k,
        gen3,
        gen4,
    }
}

/// Sanity checks on the computed matrices: `M_i^2 = -1`, `M_j^2 = -p`,
/// `M_k = M_j M_i` in the stored (transposed) layout, and the halved
/// generators double to the integral ones.
pub fn check_actions<L: Prime>(f: u64, a: &Actions) {
    let m = BigUint::one() << f;
    let p = prime::<L>() % &m;
    let neg = |x: &BigUint| (&m - (x % &m)) % &m;
    let mul = |x: &Mat, y: &Mat| -> Mat {
        let e = |r: usize, c: usize| (&x[r][0] * &y[0][c] + &x[r][1] * &y[1][c]) % &m;
        [[e(0, 0), e(0, 1)], [e(1, 0), e(1, 1)]]
    };
    let scalar =
        |s: &BigUint| -> Mat { [[s.clone(), BigUint::zero()], [BigUint::zero(), s.clone()]] };
    let eq = |x: &Mat, y: &Mat| (0..2).all(|r| (0..2).all(|c| &x[r][c] % &m == &y[r][c] % &m));
    assert!(
        eq(&mul(&a.i, &a.i), &scalar(&neg(&BigUint::one()))),
        "i^2 = -1"
    );
    assert!(eq(&mul(&a.j, &a.j), &scalar(&neg(&p))), "j^2 = -p");
    // In the stored layout M_alpha holds coordinates as columns, so
    // composition alpha o beta is M_alpha M_beta; k = i o j.
    assert!(eq(&mul(&a.i, &a.j), &a.k), "k = i j");
    let dbl = |x: &Mat| -> Mat {
        [
            [(&x[0][0] * 2u32) % &m, (&x[0][1] * 2u32) % &m],
            [(&x[1][0] * 2u32) % &m, (&x[1][1] * 2u32) % &m],
        ]
    };
    let addm = |x: &Mat, y: &Mat| -> Mat {
        [
            [(&x[0][0] + &y[0][0]) % &m, (&x[0][1] + &y[0][1]) % &m],
            [(&x[1][0] + &y[1][0]) % &m, (&x[1][1] + &y[1][1]) % &m],
        ]
    };
    assert!(eq(&dbl(&a.gen3), &addm(&a.i, &a.j)), "2 gen3 = i + j");
    assert!(
        eq(&dbl(&a.gen4), &addm(&scalar(&BigUint::one()), &a.k)),
        "2 gen4 = 1 + k"
    );
    let _ = fp2_eq::<sqisign_verify::params::P324_3>;
}
