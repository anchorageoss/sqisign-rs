//! Phase 5b.6 - the **end-to-end self-contained** Level-1 dim-4 verifier.
//!
//! [`hd_verify_l1`] takes only signature-derived data and runs the whole
//! dimension-4 FastVerify pipeline with **no oracle input on the verify path**:
//!
//! 1. challenge binding - recompute `chal` from the curves + message and compare;
//! 2. stages 2-3 - recover the challenge isogeny ([`recover_challenge_l1`]) and
//!    the response basis ([`recover_response_l1`]);
//! 3. 5b.4 - the Kani norm equation `2^e - q = a1² + a2²`
//!    ([`norm_equation_2f_minus_q`]) gives `a1, a2` (and `m = v₂(a2)`);
//! 4. the canonical dim-1 bases, product theta null, the starting symplectic
//!    matrices ([`starting_two_symplectic_matrices`], self-consistent completion)
//!    and the dim-4 gluing change of basis;
//! 5. the two gluing chains ([`KaniGluingChainHalf`]) → gluing codomain + the
//!    post-gluing kernel basis;
//! 6. the optimal-strategy chain loop ([`run_strategy_chain`]) on each half;
//! 7. the projective middle-codomain match.
//!
//! Stage 6 (the HD-image / point-evaluation condition, Phase 5b.7) is out of
//! scope; the accept/reject decision here is the middle-codomain match, which is
//! the strongest check the half-chain construction supports.

use crypto_bigint::U256;

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::{EcCurve, JacPoint};
use crate::{Fp2, Level1};

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::basis::hd_torsion_basis_l1;
use crate::hd::canonical::make_canonical;
use crate::hd::chain::middle_codomain_matches;
use crate::hd::challenge::recover_challenge_l1;
use crate::hd::gluing_chain::{
    jac_mul_u128, point_matrix_product_k, KaniGluingChainHalf, TuplePoint4,
};
use crate::hd::hd_verify::{hd_challenge_from_curves, hd_challenge_len, HdReject, MAX_CHAL_BYTES};
use crate::hd::isogeny::apply_plain_image;
use crate::hd::kani::{
    gluing_bc_dim4_f1, gluing_bc_dim4_f2, inverse_mod_pow2, norm_equation_2f_minus_q,
    starting_two_symplectic_matrices,
};
use crate::hd::point::ThetaPointDim4;
use crate::hd::product_theta::{product_theta_dim2, ThetaStructureDim1};
use crate::hd::response::{recover_response_l1, ResponseScalars};
use crate::hd::strategy::{run_strategy_chain, StrategyChain};
use crate::hd::structure::ThetaStructureDim4;

/// Level-1 parameters of the half-chain (`KaniEndoHalf` naming): the embedding
/// dimension `e = 136`, the available `2^f`-torsion exponent `f = 70`, and the
/// half-chain length `e1 = e2 = ⌈e/2⌉ = 68`.
const E_EMBED: u32 = 136;
const F_TORSION: u32 = 70;
const E1: u32 = 68;
type F = Fp2<Level1>;

/// Signature-derived inputs to [`hd_verify_l1`] (everything comes from the wire
/// signature + recovered public-key/commitment curves; no oracle data).
pub struct HdSignatureL1<'a> {
    /// Public-key curve Montgomery `A` (stage 1).
    pub a_pk: F,
    /// Commitment curve Montgomery `A` (stage 1).
    pub a_com: F,
    /// Public-key `2^f`-torsion basis hints.
    pub hint_pk_p: u32,
    pub hint_pk_q: u32,
    /// Commitment `2^f`-torsion basis hints.
    pub hint_com_p: u32,
    pub hint_com_q: u32,
    /// The signed message.
    pub message: &'a [u8],
    /// The signature's challenge, little-endian limbs (≥ [`hd_challenge_len`] bytes).
    pub chal_limbs: &'a [u64],
    /// The signature's challenge as little-endian bytes (for the binding check).
    pub claimed_chal: &'a [u8],
    /// Response scalars `(a, b, c_or_d)`.
    pub resp_a: i128,
    pub resp_b: i128,
    pub resp_c_or_d: i128,
    /// Response degree `q` (the full ~136-bit value).
    pub q: U256,
}

// small helpers (mirroring the validated tests)

fn jac_dbl_n(p: &JacPoint<Level1>, n: u32, c: &EcCurve<Level1>) -> JacPoint<Level1> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, c);
    }
    acc
}

/// `e_4(U,V)` (biextension `weil`, PARI's inverse).
fn weil4(u: &JacPoint<Level1>, v: &JacPoint<Level1>, c: &mut EcCurve<Level1>) -> F {
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

fn dim1_null(u1: &JacPoint<Level1>) -> [F; 2] {
    let (x, _) = crate::hd::basis::jac_to_affine(u1);
    ThetaStructureDim1::<Level1>::from_torsion(&x, &F::one())
        .null()
        .clone()
}

// Low 128 bits of a `U256`. Byte-oriented so it does not depend on the
// `crypto-bigint` word size (u64 on 64-bit, u32 on 32-bit/`no_std` targets).
fn low_u128(x: &U256) -> u128 {
    let b = x.to_le_bytes();
    u128::from_le_bytes(b[..16].try_into().unwrap())
}

fn u256_to_u128(x: &U256) -> u128 {
    low_u128(x)
}

/// The self-derived gluing-chain inputs (all completion-independent).
struct Setup {
    e_com: EcCurve<Level1>,
    e_chal: EcCurve<Level1>,
    points_m: [JacPoint<Level1>; 4],
    /// The full `2^e`-torsion commitment basis point `P_com` (stage-6 input `T`).
    p_com: JacPoint<Level1>,
    r_com: JacPoint<Level1>,
    s_com: JacPoint<Level1>,
    phi_r: JacPoint<Level1>,
    phi_s: JacPoint<Level1>,
    zero12: [F; 4],
    m0: [[u8; 4]; 4],
    e4: F,
    a1: u128,
    a2: u128,
    q4: u128,
    m: usize,
}

/// `B_Kpp = kernel_basis(M, e1, R_com, S_com, φR, lamb·φS)` (modulus `2^(e1+2)=2^70`).
fn b_kpp(s: &Setup, m_full: &[[u128; 8]; 8]) -> [TuplePoint4<Level1>; 4] {
    let mask = (1u128 << F_TORSION) - 1;
    let lamb = inverse_mod_pow2(s.q4, mask);
    let lamb_s2 = jac_mul_u128(&s.phi_s, lamb, &s.e_chal);
    point_matrix_product_k(
        m_full, &s.r_com, &s.s_com, &s.phi_r, &lamb_s2, mask, &s.e_com, &s.e_chal,
    )
}

/// One self-derived half-chain: the gluing chain plus the optimal-strategy plain
/// chain (kernels + codomains). Returns `None` if any step is uncomputable.
fn run_half(
    s: &Setup,
    m_full: &[[u128; 8]; 8],
    m_glue: &[[i64; 8]; 8],
    dual: bool,
) -> Option<(KaniGluingChainHalf<Level1>, StrategyChain<Level1>)> {
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
    let post: [ThetaPointDim4<Level1>; 4] = core::array::from_fn(|i| chain.evaluate(&basis[i]));
    let start = ThetaStructureDim4::new(chain.codomain_null().clone());
    let n_plain = (E1 - s.m as u32 - 1) as usize;
    let sc = run_strategy_chain(&start, &post, n_plain)?;
    Some((chain, sc))
}

/// Projective equality of the affine `x`-coordinates of two Jacobian points
/// (`x(A) = x(B) ⟺ A = ±B`): `A.x·B.z² = B.x·A.z²`.
fn x_eq(a: &JacPoint<Level1>, b: &JacPoint<Level1>) -> bool {
    let (az2, bz2) = (a.z.sqr(), b.z.sqr());
    bool::from(a.x.mul(&bz2).ct_equal(&b.x.mul(&az2)))
}

/// Evaluate a dim-4 theta point on `C1 = Hadamard(C2)` through `F2 =
/// F2_dual.dual()`: the dim-4 dual plain steps (reverse order, each
/// `precomp = inv(C_{i-1})`, `hadamard=true`) then the splitting (dual gluing).
fn eval_f2_dual(
    chain2: &KaniGluingChainHalf<Level1>,
    codomains2: &[ThetaPointDim4<Level1>],
    y: &ThetaPointDim4<Level1>,
) -> Option<TuplePoint4<Level1>> {
    let n = codomains2.len() + 1; // total F2_dual dim-4 steps (1 gluing + n-1 plain)
    let mut coords: [F; 16] = y.coords().clone();
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
        let mut inv: [F; 16] = cprev.clone();
        let mut t1: [F; 16] = core::array::from_fn(|_| F::zero());
        let mut t2: [F; 16] = core::array::from_fn(|_| F::zero());
        crate::hd::field::batched_inv(&mut inv, &mut t1, &mut t2);
        let hs = hadamard(&pointwise_square(&coords));
        let prod: [F; 16] = core::array::from_fn(|t| hs[t].mul(&inv[t]));
        coords = hadamard(&prod); // hadamard = true
    }
    chain2.splitting_eval(&ThetaPointDim4::new(coords))
}

/// Stage 6 - the HD-image check. Evaluate `T = (P_com, 0, 0, 0)` through
/// `F = F2 ∘ F1` and verify `F(T) = (±a₁·P_com, ±a₂·P_com, *, 0_{E_chal})`.
fn hd_image_check(
    s: &Setup,
    chain1: &KaniGluingChainHalf<Level1>,
    sc1: &StrategyChain<Level1>,
    chain2: &KaniGluingChainHalf<Level1>,
    sc2: &StrategyChain<Level1>,
) -> Option<bool> {
    let id = JacPoint::<Level1>::identity();
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
fn build_setup(sig: &HdSignatureL1) -> Option<Setup> {
    let chal = recover_challenge_l1(&sig.a_pk, sig.hint_pk_p, sig.hint_pk_q, sig.chal_limbs)?;
    let q4 = low_u128(&sig.q);
    let rsp = recover_response_l1(
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
    )?;

    // 5b.4: norm equation 2^e - q = a1² + a2² (a1 odd, a2 even); m = v₂(a2).
    let (a1u, a2u) = norm_equation_2f_minus_q(E_EMBED, &sig.q)?;
    let a1 = u256_to_u128(&a1u);
    let a2 = u256_to_u128(&a2u);
    let m = a2.trailing_zeros() as usize;
    // A valid response degree gives a nonzero even `a2`, so `m = v₂(a2) ≤ 67` and
    // the half-chain step counts below (`67 - m`) stay in range. Reject an
    // out-of-range `m` (e.g. `a2 = 0` when `2^e - q` is a perfect square) so the
    // counts cannot underflow.
    if m > (F_TORSION - 3) as usize {
        return None;
    }
    let lamb4 = q4 & 3;

    let mut e_com = EcCurve::<Level1>::from_a(&sig.a_com)?;
    e_com.normalize_a24();
    let mut e_chal = chal.e_chal.clone();
    e_chal.normalize_a24();

    // The full 2^e-torsion commitment basis (same recovery as recover_response,
    // before the 2^(e-r) rescaling to R_com); P_com is the stage-6 input.
    let (p_com, _q_com) = hd_torsion_basis_l1(&sig.a_com, sig.hint_com_p, sig.hint_com_q)?;

    let to4 = F_TORSION - 2;
    let p1_4 = jac_dbl_n(&rsp.r_com, to4, &e_com);
    let q1_4 = jac_dbl_n(&rsp.s_com, to4, &e_com);
    let (t1, t2, mt) = make_canonical(&p1_4, &q1_4, &mut e_com)?;
    let r2_4 = jac_dbl_n(&rsp.phi_rsp_r_com, to4, &e_chal);
    let mut s2_4 = jac_dbl_n(&rsp.phi_rsp_s_com, to4, &e_chal);
    if lamb4 == 3 {
        s2_4 = jac_add(&jac_dbl(&s2_4, &e_chal), &s2_4, &e_chal); // 3·S2_4
    }
    let (u1, _u2, mu) = make_canonical(&r2_4, &s2_4, &mut e_chal)?;

    let e4 = weil4(&t1, &t2, &mut e_com).inv();
    let zero12 = product_theta_dim2(&dim1_null(&t1), &dim1_null(&u1));
    let m0 = m0_from_canon(&mt, &mu);

    let k = F_TORSION - 3 - m as u32; // 67 - m
    let points_m = [
        jac_dbl_n(&rsp.r_com, k, &e_com),
        jac_dbl_n(&rsp.s_com, k, &e_com),
        jac_dbl_n(&rsp.phi_rsp_r_com, k, &e_chal),
        jac_dbl_n(&rsp.phi_rsp_s_com, k, &e_chal),
    ];

    Some(Setup {
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
type HalfResult = (KaniGluingChainHalf<Level1>, StrategyChain<Level1>);

/// Run both independent half-chains (F1, F2_dual). Serial by default; with the
/// `parallel` feature, on two threads (see below). Either way the results are
/// bit-identical - the halves do not interact until the middle-codomain match.
#[cfg(not(feature = "parallel"))]
fn run_both_halves(
    s: &Setup,
    m1: &[[u128; 8]; 8],
    mg1: &[[i64; 8]; 8],
    m2: &[[u128; 8]; 8],
    mg2: &[[i64; 8]; 8],
) -> Option<(HalfResult, HalfResult)> {
    let r1 = run_half(s, m1, mg1, false)?;
    let r2 = run_half(s, m2, mg2, true)?;
    Some((r1, r2))
}

/// `parallel`: F1 and F2_dual are independent until the middle-codomain
/// comparison, so run them on two scoped threads (borrowing `s`) and join
/// before stages 5-6. This only changes wall-clock latency; the total work and
/// every result are identical to the serial path.
#[cfg(feature = "parallel")]
fn run_both_halves(
    s: &Setup,
    m1: &[[u128; 8]; 8],
    mg1: &[[i64; 8]; 8],
    m2: &[[u128; 8]; 8],
    mg2: &[[i64; 8]; 8],
) -> Option<(HalfResult, HalfResult)> {
    std::thread::scope(|scope| {
        let h1 = scope.spawn(|| run_half(s, m1, mg1, false));
        // Run F2_dual on the current thread while F1 runs on the spawned one.
        let r2 = run_half(s, m2, mg2, true);
        let r1 = h1.join().expect("F1 half-chain thread panicked");
        Some((r1?, r2?))
    })
}

/// Self-derive both half-chains, then run stages 5 (middle-codomain match) and
/// 6 (HD-image). Accept only if BOTH pass.
fn self_derived_check(s: &Setup) -> Result<(), HdReject> {
    let mask = (1u128 << F_TORSION) - 1;
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
pub fn hd_verify_l1(sig: &HdSignatureL1) -> Result<(), HdReject> {
    // (0/1) challenge binding: recompute chal from the curves + message.
    let n = hd_challenge_len::<Level1>();
    let mut chal = [0u8; MAX_CHAL_BYTES];
    if !hd_challenge_from_curves(&sig.a_com, &sig.a_pk, sig.message, &mut chal[..n]) {
        return Err(HdReject::BadCurve);
    }
    if sig.claimed_chal.len() != n || chal[..n] != *sig.claimed_chal {
        return Err(HdReject::ChallengeMismatch);
    }

    // (2-6) self-derive the chain and check both the middle codomain and HD-image.
    let setup = build_setup(sig).ok_or(HdReject::ChainFailed)?;
    self_derived_check(&setup)
}

/// Convenience: `true` iff the signature verifies (self-contained).
#[inline]
pub fn hd_verify_l1_bool(sig: &HdSignatureL1) -> bool {
    hd_verify_l1(sig).is_ok()
}

/// Diagnostic: the stage-6 HD-image `F(T)` for `T = (P_com, 0, 0, 0)`, together
/// with the self-derived `a₁·P_com` and `a₂·P_com` it is checked against.
/// Returns `None` if the chain is uncomputable. Exposed for oracle validation;
/// the verify path uses the internal boolean check.
#[allow(clippy::type_complexity)]
pub fn hd_image_l1(
    sig: &HdSignatureL1,
) -> Option<(TuplePoint4<Level1>, JacPoint<Level1>, JacPoint<Level1>)> {
    let s = build_setup(sig)?;
    let mask = (1u128 << F_TORSION) - 1;
    let (m1, m2) = starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask)?;
    let mg1 = gluing_bc_dim4_f1(s.a1, s.a2, s.q4, s.m, &m1);
    let mg2 = gluing_bc_dim4_f2(s.a1, s.a2, s.q4, s.m, &m2);
    let (chain1, sc1) = run_half(&s, &m1, &mg1, false)?;
    let (chain2, sc2) = run_half(&s, &m2, &mg2, true)?;

    let id = JacPoint::<Level1>::identity();
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
