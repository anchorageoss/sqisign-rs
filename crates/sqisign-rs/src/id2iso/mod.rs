//! Ideal-to-isogeny translation (spec Section 4.6): the round-3 Qlapoty
//! method, a port of the SQIsign reference's `src/id2iso`.
//!
//! [`ideal_to_isogeny_qlapoty`] turns a left `O0`-ideal `I` into the
//! codomain of the isogeny `E0 -> E0 / I` together with the images of the
//! precomputed basis of `E0[2^f]`. It solves the Qlapoty norm equation
//! (`quat::qlapoty`) for two equivalent ideals whose norms sum to `2^e`,
//! then evaluates the isogeny through a `(2^e, 2^e)`-isogeny of `E0 x E0`
//! given by Kani's lemma (`sqisign_verify::theta`). The images of the basis
//! are a first-class output because PRISM's signatures are made of them.
//!
//! The helpers ([`endomorphism_application_even_basis`],
//! [`kernel_dlogs_to_ideal_even`], [`change_of_basis_matrix_tate`]) are the
//! reference's other `id2iso` entry points, used by key generation and
//! signing.
//!
//! Randomness enters only through Qlapoty (`rng`), so a seeded generator
//! gives the reference's result value for value. Nothing here is constant
//! time in the ideal: the quaternion arithmetic is variable time, the
//! scalar multiplications on the curves are constant time in their scalars.

use crate::mp::{Ibz, Rng};
use crate::quat::qlapoty::{qlapoty_with_stats, QlapotyStats};
use crate::quat::{Mat2x2, QuatAlg, QuatAlgElem, QuatIdeal, Vec2};
use sqisign_verify::ec::pairing::ec_dlog_2_tate;
use sqisign_verify::ec::point::ec_biscalar_mul;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint, MAX_ORDER_WORDS};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::params::Prime;
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::theta::chain::theta_chain_compute_and_eval;
use sqisign_verify::theta::{
    ChainMode, ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints,
};

pub mod params;

/// The action of the endomorphisms of `E0` on the precomputed basis of
/// `E0[2^f]`, as `2x2` matrices over `Z / 2^f` in the reference's layout
/// (`theta(P0) = m00 P0 + m10 Q0`).
#[derive(Clone, Debug)]
pub struct E0Actions<const N: usize> {
    /// Action of `i`.
    pub action_i: Mat2x2<N>,
    /// Action of `j`.
    pub action_j: Mat2x2<N>,
    /// Action of the second `O0` generator, `i`.
    pub action_gen2: Mat2x2<N>,
    /// Action of the third `O0` generator, `(i + j) / 2`.
    pub action_gen3: Mat2x2<N>,
    /// Action of the fourth `O0` generator, `(1 + k) / 2`.
    pub action_gen4: Mat2x2<N>,
}

fn mat_from_bytes<const N: usize, const B: usize>(m: &crate::precomp::Mat2x2<B>) -> Mat2x2<N> {
    Mat2x2([
        [Ibz::from_le_bytes(&m[0][0]), Ibz::from_le_bytes(&m[0][1])],
        [Ibz::from_le_bytes(&m[1][0]), Ibz::from_le_bytes(&m[1][1])],
    ])
}

impl<const N: usize> E0Actions<N> {
    /// Load the matrices from the generated little-endian constants.
    pub fn from_bytes<const B: usize>(
        i: &crate::precomp::Mat2x2<B>,
        j: &crate::precomp::Mat2x2<B>,
        gen2: &crate::precomp::Mat2x2<B>,
        gen3: &crate::precomp::Mat2x2<B>,
        gen4: &crate::precomp::Mat2x2<B>,
    ) -> Self {
        Self {
            action_i: mat_from_bytes(i),
            action_j: mat_from_bytes(j),
            action_gen2: mat_from_bytes(gen2),
            action_gen3: mat_from_bytes(gen3),
            action_gen4: mat_from_bytes(gen4),
        }
    }
}

/// Counters of an ideal-to-isogeny run, for measuring failures and retries
/// (spec Section 9.4).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Id2IsoStats {
    /// Qlapoty's loop counters.
    pub qlapoty: QlapotyStats,
    /// Norm equations Qlapoty could not solve (the reference aborts here).
    pub qlapoty_failures: u64,
    /// Kani chains that did not split into a product of elliptic curves.
    pub chain_failures: u64,
}

/// Result of an ideal-to-isogeny translation.
#[derive(Clone, Debug)]
pub struct Id2IsoOutput<L: Prime, const N: usize> {
    /// `beta1`, an element of the ideal of norm `d1 n(I)` (Qlapoty output).
    pub beta1: QuatAlgElem<N>,
    /// `d1 = n(beta1) / n(I)`, odd.
    pub d1: Ibz<N>,
    /// The codomain `E0 / I`, in the canonical form of the theta splitting.
    pub codomain: EcCurve<L>,
    /// Images of the precomputed basis `(P0, Q0, P0 - Q0)` of `E0[2^f]`.
    pub basis: EcBasis<L>,
}

/// The curve `E0: y^2 = x^3 + x`.
#[inline]
pub fn e0_curve<L: FpBackend>() -> EcCurve<L> {
    EcCurve::default()
}

/// The precomputed basis `(P0, Q0, P0 - Q0)` of `E0[2^f]`
/// (`BASIS_EVEN` in the reference).
pub fn e0_basis_even<L: FpBackend + PrimePrecomp>() -> EcBasis<L> {
    let dec = |b: &[u8]| EcPoint::from_x(Fp2::<L>::decode(b).expect("precomputed x-coordinate"));
    EcBasis::new(
        dec(L::e0_basis_px()),
        dec(L::e0_basis_qx()),
        dec(L::e0_basis_pmqx()),
    )
}

/// A non-negative integer below `2^(64 MAX_ORDER_WORDS)` as scalar words.
fn to_scalar<const N: usize>(x: &Ibz<N>) -> [u64; MAX_ORDER_WORDS] {
    let mut s = [0u64; MAX_ORDER_WORDS];
    x.to_digits(&mut s);
    s
}

/// Apply a `2x2` matrix to a basis of `E[2^f]`, in place: for
/// `mat = [[a, c], [b, d]]` the new basis is `([a]P + [b]Q, [c]P + [d]Q)`,
/// and the third point is recomputed when `last_point` is set. The matrix
/// is reduced modulo `2^f` in place. `None` if a scalar multiplication is
/// not defined (a basis point with a zero coordinate).
pub fn matrix_application_even_basis<L: FpBackend, const N: usize>(
    bas: &mut EcBasis<L>,
    e: &EcCurve<L>,
    mat: &mut Mat2x2<N>,
    f: u32,
    last_point: bool,
) -> Option<()> {
    let pow_two = Ibz::<N>::one().mul_2exp(f);
    let tmp_bas = bas.clone();
    for row in mat.0.iter_mut() {
        for m in row.iter_mut() {
            *m = m.modulo(&pow_two);
        }
    }
    let w = L::ORDER_WORDS;
    let f = f as usize;
    let s0 = to_scalar(&mat.0[0][0]);
    let s1 = to_scalar(&mat.0[1][0]);
    bas.p = ec_biscalar_mul(&s0[..w], &s1[..w], f, &tmp_bas, e)?;
    let s0 = to_scalar(&mat.0[0][1]);
    let s1 = to_scalar(&mat.0[1][1]);
    bas.q = ec_biscalar_mul(&s0[..w], &s1[..w], f, &tmp_bas, e)?;
    if last_point {
        let d0 = mat.0[0][0].sub(&mat.0[0][1]).modulo(&pow_two);
        let d1 = mat.0[1][0].sub(&mat.0[1][1]).modulo(&pow_two);
        let s0 = to_scalar(&d0);
        let s1 = to_scalar(&d1);
        bas.pmq = ec_biscalar_mul(&s0[..w], &s1[..w], f, &tmp_bas, e)?;
    }
    Some(())
}

/// Apply an endomorphism `theta` of `E0` (an element of `O0` up to an odd
/// content) to a basis of `E[2^f]`, in place, through its action matrix.
/// The matrix of `theta ∈ O0` on the precomputed basis of `E0[2^f]`,
/// modulo `2^f`: `theta` written as `content * primitive` in `O0`'s basis,
/// the generators' action matrices combined. `content` must be odd.
pub fn endomorphism_action_matrix<const N: usize>(
    theta: &QuatAlgElem<N>,
    f: u32,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
) -> Mat2x2<N> {
    let two_e = Ibz::<N>::one().mul_2exp(f);
    let (mut coeffs, content) = theta.make_primitive(&alg.o0);
    for c in coeffs.0.iter_mut() {
        *c = c.modulo(&two_e);
    }
    debug_assert!(content.is_odd());
    let content = content.modulo(&two_e);
    let zero = Ibz::<N>::zero();
    let mut mat = Mat2x2([[zero, zero], [zero, zero]]);
    let gens = [&act.action_gen2, &act.action_gen3, &act.action_gen4];
    for i in 0..2 {
        mat.0[i][i] = mat.0[i][i].add(&coeffs.0[0]);
        for j in 0..2 {
            for (g, c) in gens.iter().zip(coeffs.0[1..].iter()) {
                let tmp = g.0[i][j].mul(c);
                mat.0[i][j] = mat.0[i][j].add(&tmp).modulo(&two_e);
            }
            mat.0[i][j] = mat.0[i][j].mul(&content).modulo(&two_e);
        }
    }
    mat
}

/// Apply `theta ∈ O0` to `bas`, the image of the precomputed basis of
/// `E0[2^f]` on `e`: `bas <- ψ(theta(B0))` through the action matrix.
pub fn endomorphism_application_even_basis<L: FpBackend, const N: usize>(
    bas: &mut EcBasis<L>,
    e: &EcCurve<L>,
    theta: &QuatAlgElem<N>,
    f: u32,
    last_point: bool,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
) -> Option<()> {
    let mut mat = endomorphism_action_matrix(theta, f, alg, act);
    matrix_application_even_basis(bas, e, &mut mat, f, last_point)
}

/// The `O0`-ideal of norm `2^f` whose kernel is generated by
/// `vec2[0] P0 + vec2[1] Q0` on the precomputed basis of `E0[2^f]`
/// (KernelToIdeal): the inert part and the split part (`1 + i` or `i - 1`
/// when the generator was divisible by `1 + i`, `1` otherwise). `None` if
/// the point does not have full order.
pub fn kernel_dlogs_to_ideal_even<const N: usize>(
    vec2: &Vec2<N>,
    f: u32,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
) -> Option<(QuatIdeal<N>, QuatAlgElem<N>)> {
    let mut two_pow = Ibz::<N>::one().mul_2exp(f);
    // a P + b [j + (1 + k)/2] P = [i] P
    let mut mat = Mat2x2([[vec2.0[0], Ibz::zero()], [vec2.0[1], Ibz::zero()]]);
    let v = act.action_j.eval(vec2);
    mat.0[0][1] = v.0[0];
    mat.0[1][1] = v.0[1];
    let v = act.action_gen4.eval(vec2);
    mat.0[0][1] = mat.0[0][1].add(&v.0[0]).modulo(&two_pow);
    mat.0[1][1] = mat.0[1][1].add(&v.0[1]).modulo(&two_pow);
    let (inv, ok) = mat.inv_mod(&two_pow);
    if !ok {
        return None;
    }
    let v = act.action_i.eval(vec2);
    let mut v = inv.eval(&v);
    v.0[0] = v.0[0].modulo(&two_pow);
    v.0[1] = v.0[1].modulo(&two_pow);

    // gen = a - i + b (j + (1 + k)/2)
    let mut gen = QuatAlgElem::<N>::set(2, 0, -2, 0, 0);
    gen.coord.0[0] = v.0[0].add(&v.0[0]).add(&v.0[1]);
    gen.coord.0[2] = v.0[1].add(&v.0[1]);
    gen.coord.0[3] = v.0[1];
    for c in gen.coord.0.iter_mut() {
        c.set_bound(f as i32 + 10);
    }

    // divisible by (1 + i)? then remove that factor and record it
    let mut split_test = QuatAlgElem::<N>::set(1, 1, 1, 0, 0);
    split_test = gen.mul(&split_test, alg);
    let (_, n) = split_test.make_primitive(&alg.o0);
    let split;
    if !n.is_odd() {
        let (c, exact) = split_test.coord.scalar_div(&n);
        debug_assert!(exact);
        split_test.coord = c;
        gen = split_test;
        split = QuatAlgElem::set(1, 1, -1, 0, 0);
        two_pow = Ibz::<N>::one().mul_2exp(f - 1);
    } else {
        split = QuatAlgElem::set(1, 1, 0, 0, 0);
    }
    let ideal = QuatIdeal::create_o0_pow_two(&gen, &two_pow, alg)?;
    debug_assert!(ideal.norm == two_pow);
    Some((ideal, split))
}

fn mat_from_dlogs<L: Prime, const N: usize>(d: &[[u64; MAX_ORDER_WORDS]; 4]) -> Mat2x2<N> {
    let w = L::ORDER_WORDS;
    Mat2x2([
        [Ibz::from_digits(&d[0][..w]), Ibz::from_digits(&d[2][..w])],
        [Ibz::from_digits(&d[1][..w]), Ibz::from_digits(&d[3][..w])],
    ])
}

/// Change of basis: `mat` with `(mat v) . B2 = v . B1` for a basis `B1` of
/// `E[2^f]` and a basis `B2` of the full `E[2^e]`, by Tate pairings. `None`
/// if `B1` is not in the span of `B2`.
pub fn change_of_basis_matrix_tate<L: FpBackend + PrimePrecomp, const N: usize>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    e: &mut EcCurve<L>,
    f: u32,
) -> Option<Mat2x2<N>> {
    let d = ec_dlog_2_tate(b2, b1, e, f)?;
    Some(mat_from_dlogs::<L, N>(&d))
}

/// Change of basis the other way: `mat` with `(mat v) . B1 = [2^(e-f)] v . B2`
/// for a basis `B1` of the full `E[2^e]` and a basis `B2` of `E[2^f]`, by
/// inverting [`change_of_basis_matrix_tate`]'s outcome modulo `2^f`.
pub fn change_of_basis_matrix_tate_invert<L: FpBackend + PrimePrecomp, const N: usize>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    e: &mut EcCurve<L>,
    f: u32,
) -> Option<Mat2x2<N>> {
    let d = ec_dlog_2_tate(b1, b2, e, f)?;
    let mut m: Mat2x2<N> = mat_from_dlogs::<L, N>(&d);
    let [[mut r1, mut s1], [mut r2, mut s2]] = m.0;
    let ok = Ibz::invmat(&mut r1, &mut r2, &mut s1, &mut s2, f as i32);
    debug_assert!(ok);
    m.0 = [[r1, s1], [r2, s2]];
    Some(m)
}

/// IdealToIsogeny (spec Algorithm 4.6, Qlapoty): the codomain of the
/// isogeny `E0 -> E0 / I` and the images of the precomputed basis of
/// `E0[2^f]`, for a left `O0`-ideal `I` of odd norm. The norm equation
/// gives `beta1, beta2` in `I` with `n(beta1) + n(beta2) = n(I) 2^e`; the
/// isogeny is read off the `(2^e, 2^e)`-isogeny of `E0 x E0` with kernel
/// `{([d1^-1 theta] P, P)}`. `None` if the norm equation has no solution
/// or the chain does not split (both counted in `stats`); the caller
/// retries with fresh randomness as key generation does.
pub fn ideal_to_isogeny_qlapoty<L: FpBackend + PrimePrecomp, const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
    rng: &mut impl Rng,
    stats: &mut Id2IsoStats,
) -> Option<Id2IsoOutput<L, N>> {
    let f = L::TWO_ADIC_EXPONENT;
    let used_torsion = alg.qlapoty_used_power_of_two;
    debug_assert!(used_torsion + 2 <= f);

    let Some((beta1, d1, theta)) = qlapoty_with_stats(ideal, alg, rng, &mut stats.qlapoty) else {
        stats.qlapoty_failures += 1;
        return None;
    };
    debug_assert!(d1.is_odd());

    // theta / d1 on the basis of E0 (mod 2^f)
    let two_e = Ibz::<N>::one().mul_2exp(f);
    let d1_inv = d1.invmod(&two_e)?;
    let e0 = e0_curve::<L>();
    let basis_even = e0_basis_even::<L>();

    let apply = theta.scalar_mul(&d1_inv);
    let mut basis_temp = basis_even.clone();
    endomorphism_application_even_basis(&mut basis_temp, &e0, &apply, f, false, alg, act)?;
    let ker = ThetaKernelCouplePoints::from_bases(&basis_temp, &basis_even);
    let mut e0e0 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };

    // beta1 / d1 on the basis of E0: the points to push
    let apply = beta1.scalar_mul(&d1_inv);
    let mut basis_temp = basis_even;
    endomorphism_application_even_basis(&mut basis_temp, &e0, &apply, f, true, alg, act)?;
    let push = [
        ThetaCouplePoint::on_e1(basis_temp.p),
        ThetaCouplePoint::on_e1(basis_temp.q),
        ThetaCouplePoint::on_e1(basis_temp.pmq),
    ];

    let Some(out) =
        theta_chain_compute_and_eval(used_torsion as u16, &mut e0e0, &ker, &push, ChainMode::E2)
    else {
        stats.chain_failures += 1;
        return None;
    };
    // with the fixed theta structure the codomain of interest is always E2
    let img = |i: usize| out.images[i].as_ref().expect("pushed point").p2.clone();
    Some(Id2IsoOutput {
        beta1,
        d1,
        codomain: out.codomain.e2,
        basis: EcBasis::new(img(0), img(1), img(2)),
    })
}

/// The codomain and basis images only (the reference's
/// `dim2id2iso_arbitrary_isogeny_evaluation`).
pub fn arbitrary_isogeny_evaluation<L: FpBackend + PrimePrecomp, const N: usize>(
    ideal: &QuatIdeal<N>,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
    rng: &mut impl Rng,
) -> Option<(EcCurve<L>, EcBasis<L>)> {
    let out = ideal_to_isogeny_qlapoty::<L, N>(ideal, alg, act, rng, &mut Id2IsoStats::default())?;
    Some((out.codomain, out.basis))
}
