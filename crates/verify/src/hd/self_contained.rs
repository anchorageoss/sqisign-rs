//! The compact verifier, end to end and self-contained: from the public
//! key, the commitment curve, the hints, the challenge scalar and the
//! response scalars to accept or reject, with no data beyond the
//! signature's.
//!
//! 1. the challenge isogeny and the rescaled basis of `E_chall`
//!    ([`recover_challenge`]) and the response images
//!    ([`recover_response`]);
//! 2. the norm equation `2^e − q = a_1^2 + a_2^2` ([`norm_equation_2f_minus_q`])
//!    and `m = v_2(a_2)`;
//! 3. the canonical four-torsion, the product theta null point, the starting
//!    symplectic matrices ([`starting_two_symplectic_matrices`]) and the
//!    gluing changes of basis;
//! 4. the two gluing chains ([`KaniGluingChainHalf`]) and the optimal-strategy
//!    chains ([`run_strategy_chain`]) on each half;
//! 5. the projective middle-codomain match ([`middle_codomain_matches`]) and
//!    the paper's image condition on `(P_com, 0, 0, 0)`.
//!
//! Level I parameters come from [`HdLevel`]; the layer is the round-2
//! compact verifier made generic over the prime.

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::{EcCurve, JacPoint};
use crate::fp::{Fp2, FpBackend};

use super::params::HdLevel;
use super::uint::U4;

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::basis::hd_torsion_basis;
use crate::hd::canonical::make_canonical;
use crate::hd::chain::middle_codomain_matches;
use crate::hd::challenge::recover_challenge;
use crate::hd::gluing_chain::{
    jac_mul_u128, point_matrix_product_k, KaniGluingChainHalf, TuplePoint4,
};
use crate::hd::hd_verify::HdReject;
use crate::hd::isogeny::apply_plain_image;
use crate::hd::kani::{
    gluing_bc_dim4_f1, gluing_bc_dim4_f2, inverse_mod_pow2, norm_equation_2f_minus_q,
    starting_two_symplectic_matrices,
};
use crate::hd::point::ThetaPointDim4;
use crate::hd::product_theta::{product_theta_dim2, ThetaStructureDim1};
use crate::hd::response::{recover_response, ResponseScalars};
use crate::hd::strategy::{run_strategy_chain, StrategyChain};
use crate::hd::structure::ThetaStructureDim4;

/// Upper bound on `m = v_2(a_2)` accepted by [`build_setup`]: the gluing base
/// change reduces values mod `2^(m+2)` and multiplies them as `i128`, with
/// the largest products about `2^(3m+5)`; `m ≤ 40` keeps them in range, and
/// stays far below the `r − 3` cap the half-chain step counts need. A
/// legitimate response has a tiny `m` (`v_2` of a random even number), so
/// this rejects only a crafted `q`.
const MAX_GLUING_M: usize = 40;

/// Signature-derived inputs to [`hd_verify`] (everything comes from the wire
/// signature + recovered public-key/commitment curves; no oracle data).
/// A parsed compact signature with its public key: what [`hd_verify`]
/// consumes. `chal_limbs` is the recomputed challenge, `λ` bits.
pub struct HdSignature<'a, L: FpBackend> {
    /// Public-key curve Montgomery `A` (stage 1).
    pub a_pk: Fp2<L>,
    /// Commitment curve Montgomery `A` (stage 1).
    pub a_com: Fp2<L>,
    /// Public-key `2^f`-torsion basis hints.
    pub hint_pk_p: u32,
    pub hint_pk_q: u32,
    /// Commitment `2^f`-torsion basis hints.
    pub hint_com_p: u32,
    pub hint_com_q: u32,
    /// The signed message.
    /// The signature's challenge, little-endian limbs (≥ `hd_challenge_len` bytes).
    pub chal_limbs: &'a [u64],
    /// The signature's challenge as little-endian bytes (for the binding check).
    /// Response scalars `(a, b, c_or_d)`.
    pub resp_a: i128,
    pub resp_b: i128,
    pub resp_c_or_d: i128,
    /// Response degree `q` (the full ~136-bit value).
    /// The response degree, `q < 2^e`.
    pub q: U4,
}

// small helpers (mirroring the validated tests)

fn jac_dbl_n<L: HdLevel>(p: &JacPoint<L>, n: u32, c: &EcCurve<L>) -> JacPoint<L> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, c);
    }
    acc
}

/// `e_4(U,V)` (biextension `weil`, PARI's inverse).
fn weil4<L: HdLevel>(u: &JacPoint<L>, v: &JacPoint<L>, c: &mut EcCurve<L>) -> Fp2<L> {
    let uv = jac_add(u, &v.neg(), c);
    weil(2, &u.to_xz(), &v.to_xz(), &uv.to_xz(), c)
}

fn m0_from_canon(mt: &[[u8; 2]; 2], mu: &[[u8; 2]; 2]) -> [[u8; 4]; 4] {
    [
        [mt[0][0], 0, mt[1][0], 0],
        [0, mu[0][0], 0, mu[1][0]],
        [mt[0][1], 0, mt[1][1], 0],
        [0, mu[0][1], 0, mu[1][1]],
    ]
}

fn dim1_null<L: HdLevel>(u1: &JacPoint<L>) -> [Fp2<L>; 2] {
    let (x, _) = crate::hd::basis::jac_to_affine(u1);
    ThetaStructureDim1::<L>::from_torsion(&x, &Fp2::<L>::one())
        .null()
        .clone()
}

// Low 128 bits of a `U256`. Byte-oriented so it does not depend on the
// `crypto-bigint` word size (u64 on 64-bit, u32 on 32-bit/`no_std` targets).

/// The self-derived gluing-chain inputs (all completion-independent).
struct Setup<L: FpBackend> {
    e_com: EcCurve<L>,
    e_chal: EcCurve<L>,
    points_m: [JacPoint<L>; 4],
    /// The full `2^e`-torsion commitment basis point `P_com` (stage-6 input `T`).
    p_com: JacPoint<L>,
    r_com: JacPoint<L>,
    s_com: JacPoint<L>,
    phi_r: JacPoint<L>,
    phi_s: JacPoint<L>,
    zero12: [Fp2<L>; 4],
    m0: [[u8; 4]; 4],
    e4: Fp2<L>,
    a1: u128,
    a2: u128,
    q4: u128,
    m: usize,
}

/// `B_Kpp = kernel_basis(M, e1, R_com, S_com, φR, lamb·φS)` (modulus `2^(e1+2)=2^70`).
fn b_kpp<L: HdLevel>(s: &Setup<L>, m_full: &[[u128; 8]; 8]) -> [TuplePoint4<L>; 4] {
    let mask = (1u128 << L::R) - 1;
    let lamb = inverse_mod_pow2(s.q4, mask);
    let lamb_s2 = jac_mul_u128(&s.phi_s, lamb, &s.e_chal);
    point_matrix_product_k(
        m_full, &s.r_com, &s.s_com, &s.phi_r, &lamb_s2, mask, &s.e_com, &s.e_chal,
    )
}

/// One self-derived half-chain: the gluing chain plus the optimal-strategy plain
/// chain (kernels + codomains). Returns `None` if any step is uncomputable.
fn run_half<L: HdLevel>(
    s: &Setup<L>,
    m_full: &[[u128; 8]; 8],
    m_glue: &[[i64; 8]; 8],
    dual: bool,
) -> Option<(KaniGluingChainHalf<L>, StrategyChain<L>)> {
    let chain = KaniGluingChainHalf::new(
        &s.points_m,
        &s.zero12,
        &s.m0,
        &s.e4,
        s.a1,
        s.a2,
        s.q4,
        s.m,
        m_full,
        m_glue,
        dual,
        &s.e_com,
        &s.e_chal,
    )?;
    let basis = b_kpp(s, m_full);
    let post: [ThetaPointDim4<L>; 4] = core::array::from_fn(|i| chain.evaluate(&basis[i]));
    let start = ThetaStructureDim4::<L>::new(chain.codomain_null().clone());
    let n_plain = (L::E1 - s.m as u32 - 1) as usize;
    let sc = run_strategy_chain(&start, &post, n_plain)?;
    Some((chain, sc))
}

/// Projective equality of the affine `x`-coordinates of two Jacobian points
/// (`x(A) = x(B) ⟺ A = ±B`): `A.x·B.z² = B.x·A.z²`.
fn x_eq<L: HdLevel>(a: &JacPoint<L>, b: &JacPoint<L>) -> bool {
    let (az2, bz2) = (a.z.sqr(), b.z.sqr());
    bool::from(a.x.mul(&bz2).ct_equal(&b.x.mul(&az2)))
}

/// Evaluate a dim-4 theta point on `C1 = Hadamard(C2)` through `F2 =
/// F2_dual.dual()`: the dim-4 dual plain steps (reverse order, each
/// `precomp = inv(C_{i-1})`, `hadamard=true`) then the splitting (dual gluing).
fn eval_f2_dual<L: HdLevel>(
    chain2: &KaniGluingChainHalf<L>,
    codomains2: &[ThetaPointDim4<L>],
    y: &ThetaPointDim4<L>,
) -> Option<TuplePoint4<L>> {
    let n = codomains2.len() + 1; // total F2_dual dim-4 steps (1 gluing + n-1 plain)
    let mut coords: [Fp2<L>; 16] = y.coords().clone();
    for i in (1..n).rev() {
        // Dual of forward step i (C_{i-1} → C_i); precomp = inv(C_{i-1}).
        let k = i - 1;
        let cprev = if k == 0 {
            chain2.codomain_null().coords()
        } else {
            codomains2[k - 1].coords()
        };
        if cprev.iter().any(|x| bool::from(x.ct_is_zero())) {
            return None;
        }
        // One batched inversion of the codomain null (Montgomery's trick),
        // rather than 16 individual inversions per dual step.
        let mut inv: [Fp2<L>; 16] = cprev.clone();
        let mut t1: [Fp2<L>; 16] = core::array::from_fn(|_| Fp2::<L>::zero());
        let mut t2: [Fp2<L>; 16] = core::array::from_fn(|_| Fp2::<L>::zero());
        crate::hd::field::batched_inv(&mut inv, &mut t1, &mut t2);
        let hs = hadamard(&pointwise_square(&coords));
        let prod: [Fp2<L>; 16] = core::array::from_fn(|t| hs[t].mul(&inv[t]));
        coords = hadamard(&prod); // hadamard = true
    }
    chain2.splitting_eval(&ThetaPointDim4::new(coords))
}

/// Stage 6 - the HD-image check. Evaluate `T = (P_com, 0, 0, 0)` through
/// `F = F2 ∘ F1` and verify `F(T) = (±a₁·P_com, ±a₂·P_com, *, 0_{E_chal})`.
fn hd_image_check<L: HdLevel>(
    s: &Setup<L>,
    chain1: &KaniGluingChainHalf<L>,
    sc1: &StrategyChain<L>,
    chain2: &KaniGluingChainHalf<L>,
    sc2: &StrategyChain<L>,
) -> Option<bool> {
    let id = JacPoint::<L>::identity();
    let t4 = TuplePoint4::new(s.p_com.clone(), id.clone(), id.clone(), id);

    // F1: gluing chain, then the plain (2,2,2,2) steps. The strategy loop
    // already built each isogeny; reuse its stored image precomputation `1/O`
    // (`apply_plain_image`) instead of rebuilding the isogenies via from_kernel.
    let mut pt = chain1.evaluate(&t4);
    for inv in &sc1.image_precomp {
        pt = apply_plain_image(inv, &pt);
    }
    // F2 = F2_dual.dual(): pt is on C1 = Hadamard(C2).
    let ft = eval_f2_dual(chain2, &sc2.codomains, &pt)?;

    let a1p = jac_mul_u128(&s.p_com, s.a1, &s.e_com);
    let a2p = jac_mul_u128(&s.p_com, s.a2, &s.e_com);
    let ok0 = x_eq(&ft.c[0], &a1p);
    let ok1 = x_eq(&ft.c[1], &a2p);
    let ok3 = bool::from(ft.c[3].z.ct_is_zero()); // FT[3] = 0_{E_chal}
    Some(ok0 && ok1 && ok3)
}

/// Recover the response basis and assemble the completion-independent setup.
/// Returns `None` if challenge/response recovery fails or the norm equation has
/// no solution (a malformed/forged `q`).
fn build_setup<L: HdLevel>(sig: &HdSignature<L>) -> Result<Setup<L>, HdReject> {
    let chal = recover_challenge::<L>(&sig.a_pk, sig.hint_pk_p, sig.hint_pk_q, sig.chal_limbs)
        .ok_or(HdReject::ChallengeRecovery)?;
    let q4 = sig.q.low_u128();
    let rsp = recover_response::<L>(
        &chal,
        &sig.a_com,
        sig.hint_com_p,
        sig.hint_com_q,
        ResponseScalars {
            a: sig.resp_a,
            b: sig.resp_b,
            c_or_d: sig.resp_c_or_d,
            q: q4,
        },
    )
    .ok_or(HdReject::ResponseRecovery)?;

    // 5b.4: norm equation 2^e - q = a1² + a2² (a1 odd, a2 even); m = v₂(a2).
    let (a1u, a2u) = norm_equation_2f_minus_q(L::E_EMBED, &sig.q).ok_or(HdReject::NormEquation)?;
    let a1 = a1u.low_u128();
    let a2 = a2u.low_u128();
    let m = a2.trailing_zeros() as usize;
    // `m = v₂(a2)` must lie in the protocol's small range. Out-of-range values
    // would otherwise drive the dim-4 gluing base change's `i128` products past
    // `i128::MAX` (overflow / panic) and the half-chain step counts (`67 - m`)
    // below into a `u32` underflow. See [`MAX_GLUING_M`].
    if m > MAX_GLUING_M {
        return Err(HdReject::GluingBound);
    }
    let lamb4 = q4 & 3;

    let mut e_com = EcCurve::<L>::from_a(&sig.a_com).ok_or(HdReject::BadCurve)?;
    e_com.normalize_a24();
    let mut e_chal = chal.e_chal.clone();
    e_chal.normalize_a24();

    // The full 2^e-torsion commitment basis (same recovery as recover_response,
    // before the 2^(e-r) rescaling to R_com); P_com is the stage-6 input.
    let (p_com, _q_com) = hd_torsion_basis::<L>(&sig.a_com, sig.hint_com_p, sig.hint_com_q)
        .ok_or(HdReject::BadCurve)?;

    let to4 = L::R - 2;
    let p1_4 = jac_dbl_n::<L>(&rsp.r_com, to4, &e_com);
    let q1_4 = jac_dbl_n(&rsp.s_com, to4, &e_com);
    let (t1, t2, mt) = make_canonical(&p1_4, &q1_4, &mut e_com).ok_or(HdReject::CanonicalCom)?;
    let r2_4 = jac_dbl_n(&rsp.phi_rsp_r_com, to4, &e_chal);
    let mut s2_4 = jac_dbl_n(&rsp.phi_rsp_s_com, to4, &e_chal);
    if lamb4 == 3 {
        s2_4 = jac_add(&jac_dbl(&s2_4, &e_chal), &s2_4, &e_chal); // 3·S2_4
    }
    let (u1, _u2, mu) = make_canonical(&r2_4, &s2_4, &mut e_chal).ok_or(HdReject::CanonicalChal)?;

    let e4 = weil4(&t1, &t2, &mut e_com).inv();
    let zero12 = product_theta_dim2(&dim1_null(&t1), &dim1_null(&u1));
    let m0 = m0_from_canon(&mt, &mu);

    let k = L::R - 3 - m as u32;
    let points_m = [
        jac_dbl_n(&rsp.r_com, k, &e_com),
        jac_dbl_n(&rsp.s_com, k, &e_com),
        jac_dbl_n(&rsp.phi_rsp_r_com, k, &e_chal),
        jac_dbl_n(&rsp.phi_rsp_s_com, k, &e_chal),
    ];

    Ok(Setup {
        e_com,
        e_chal,
        points_m,
        p_com,
        r_com: rsp.r_com,
        s_com: rsp.s_com,
        phi_r: rsp.phi_rsp_r_com,
        phi_s: rsp.phi_rsp_s_com,
        zero12,
        m0,
        e4,
        a1,
        a2,
        q4,
        m,
    })
}

/// The output of one half-chain: its gluing chain plus the strategy-derived
/// plain chain.
type HalfResult<L> = (KaniGluingChainHalf<L>, StrategyChain<L>);

/// Run both independent half-chains (F1, F2_dual). Serial by default; with the
/// `parallel` feature, on two threads (see below). Either way the results are
/// bit-identical - the halves do not interact until the middle-codomain match.
fn run_both_halves<L: HdLevel>(
    s: &Setup<L>,
    m1: &[[u128; 8]; 8],
    mg1: &[[i64; 8]; 8],
    m2: &[[u128; 8]; 8],
    mg2: &[[i64; 8]; 8],
) -> Option<(HalfResult<L>, HalfResult<L>)> {
    let r1 = run_half(s, m1, mg1, false)?;
    let r2 = run_half(s, m2, mg2, true)?;
    Some((r1, r2))
}

/// Self-derive both half-chains, then run stages 5 (middle-codomain match) and
/// 6 (HD-image). Accept only if BOTH pass.
fn self_derived_check<L: HdLevel>(s: &Setup<L>) -> Result<(), HdReject> {
    let mask = (1u128 << L::R) - 1;
    let (m1, m2) =
        starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask).ok_or(HdReject::ChainFailed)?;
    let mg1 = gluing_bc_dim4_f1(s.a1, s.a2, s.q4, s.m, &m1);
    let mg2 = gluing_bc_dim4_f2(s.a1, s.a2, s.q4, s.m, &m2);

    let ((chain1, sc1), (chain2, sc2)) =
        run_both_halves(s, &m1, &mg1, &m2, &mg2).ok_or(HdReject::ChainFailed)?;

    // Stage 5: the projective middle-codomain match.
    let c1 = sc1.last_codomain().ok_or(HdReject::ChainFailed)?;
    let c2 = sc2.last_codomain().ok_or(HdReject::ChainFailed)?;
    if !middle_codomain_matches(c1, c2) {
        return Err(HdReject::MiddleCodomainMismatch);
    }

    // Stage 6: the HD-image condition F(T) = (±a₁P, ±a₂P, *, 0).
    if !hd_image_check(s, &chain1, &sc1, &chain2, &sc2).unwrap_or(false) {
        return Err(HdReject::HdImageMismatch);
    }
    Ok(())
}

/// Run the end-to-end self-contained Level-1 verification (all 6 FastVerify
/// stages).
///
/// Performs (0/1) the challenge binding, then self-derives stages 2-4 from the
/// signature, and checks (5) the dim-4 middle-codomain match and (6) the
/// HD-image condition. Returns `Ok(())` on accept or the rejection reason. No
/// oracle data is consulted.
pub fn hd_verify<L: HdLevel>(sig: &HdSignature<L>) -> Result<(), HdReject> {
    let setup = build_setup(sig)?;
    self_derived_check(&setup)
}

/// Convenience: `true` iff the signature verifies (self-contained).
#[inline]
pub fn hd_verify_bool<L: HdLevel>(sig: &HdSignature<L>) -> bool {
    hd_verify::<L>(sig).is_ok()
}

/// Diagnostic: the stage-6 HD-image `F(T)` for `T = (P_com, 0, 0, 0)`, together
/// with the self-derived `a₁·P_com` and `a₂·P_com` it is checked against.
/// Returns `None` if the chain is uncomputable. Exposed for oracle validation;
/// the verify path uses the internal boolean check.
#[allow(clippy::type_complexity)]
pub fn hd_image<L: HdLevel>(
    sig: &HdSignature<L>,
) -> Option<(TuplePoint4<L>, JacPoint<L>, JacPoint<L>)> {
    let s = build_setup(sig).ok()?;
    let mask = (1u128 << L::R) - 1;
    let (m1, m2) = starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask)?;
    let mg1 = gluing_bc_dim4_f1(s.a1, s.a2, s.q4, s.m, &m1);
    let mg2 = gluing_bc_dim4_f2(s.a1, s.a2, s.q4, s.m, &m2);
    let (chain1, sc1) = run_half(&s, &m1, &mg1, false)?;
    let (chain2, sc2) = run_half(&s, &m2, &mg2, true)?;

    let id = JacPoint::<L>::identity();
    let t4 = TuplePoint4::new(s.p_com.clone(), id.clone(), id.clone(), id);
    let mut pt = chain1.evaluate(&t4);
    for inv in &sc1.image_precomp {
        pt = apply_plain_image(inv, &pt);
    }
    let ft = eval_f2_dual(&chain2, &sc2.codomains, &pt)?;
    let a1p = jac_mul_u128(&s.p_com, s.a1, &s.e_com);
    let a2p = jac_mul_u128(&s.p_com, s.a2, &s.e_com);
    Some((ft, a1p, a2p))
}
