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

use alloc::vec::Vec;

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::{EcCurve, JacPoint};
use crate::fp::{Fp2, FpBackend};

use super::params::HdLevel;
use super::uint::U4;

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::canonical::make_canonical;
use crate::hd::chain::middle_codomain_matches;
use crate::hd::challenge::recover_challenge_prof;
use crate::hd::gluing_chain::{
    jac_mul_u128, point_matrix_product_k_fixed, FixedBase, KaniGluingChainHalf, TuplePoint4,
};
use crate::hd::hd_verify::HdReject;
use crate::hd::isogeny::apply_plain_image;
use crate::hd::kani::{
    gluing_bc_dim4_f1, gluing_bc_dim4_f2, inverse_mod_pow2, norm_equation_2f_minus_q,
    starting_two_symplectic_matrices,
};
use crate::hd::point::ThetaPointDim4;
use crate::hd::product_theta::{product_theta_dim2, ThetaStructureDim1};
use crate::hd::profile::{
    CycleProfiler, NoProfiler, Profiler, SLOTS, SLOT_BASIS, SLOT_COMPARE, SLOT_GLUING, SLOT_IMAGE,
    SLOT_SETUP, SUB_SLOTS,
};
use crate::hd::response::{recover_response_prof, ResponseScalars};
use crate::hd::strategy::{run_strategy_chain_prof, StrategyChain};
use crate::hd::structure::ThetaStructureDim4;
use crate::hd::two_point::{
    run_strategy_chain_two_point_checked_prof, run_strategy_chain_two_point_prof,
};

/// Upper bound on `m = v_2(a_2)` accepted by [`build_setup`]: the gluing base
/// change reduces values mod `2^(m+2)` and multiplies them as `i128`, with
/// the largest products about `2^(3m+5)`; `m ≤ 40` keeps them in range, and
/// stays far below the `r − 3` cap the half-chain step counts need. A
/// legitimate response has a tiny `m` (`v_2` of a random even number), so
/// this rejects only a crafted `q`.
const MAX_GLUING_M: usize = 40;

/// How many kernel points the plain chain of each half-chain carries from
/// the gluing codomain on (SECENG-2421).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum KernelPoints {
    /// The four generators `g0` to `g3`, as before the two-point variant
    /// existed: the default, and the path [`hd_verify`] runs.
    #[default]
    Four,
    /// `g0` and `g1` only, with `g2` and `g3` formed by the permutation
    /// [`PI`](crate::hd::two_point::PI) at every step, after a check at the
    /// gluing: the gluing images of `g2` and `g3` are compared projectively
    /// with `PI` of those of `g0` and `g1`, and a half-chain whose images do
    /// not match runs with four points (reported by
    /// [`HdVerifyReport::fell_back`]).
    TwoWithImageCheck,
    /// `g0` and `g1` only, with `g2` and `g3` formed by `PI` at every step,
    /// without the check at the gluing.
    TwoWithoutImageCheck,
}

/// The options of [`hd_verify_with`]. The default is the four-point
/// verifier.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HdVerifyOptions {
    /// The kernel points carried after the gluing.
    pub kernel_points: KernelPoints,
}

/// What [`hd_verify_with`] reports beside the verdict.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HdVerifyReport {
    /// Per half-chain (`F1`, then `F2_dual`): `true` when the half ran with
    /// four points because the check of [`KernelPoints::TwoWithImageCheck`]
    /// failed. Both `false` in the other modes, and for a signature rejected
    /// before its half-chains.
    pub fell_back: [bool; 2],
}

/// The codomain null points of both half-chains, from [`hd_codomains`].
#[derive(Clone, Debug)]
pub struct HdCodomains<L: FpBackend> {
    /// The gluing codomain of `F1` and of `F2_dual`.
    pub gluing: [ThetaPointDim4<L>; 2],
    /// The plain steps' codomains of `F1` and of `F2_dual`, in chain order.
    pub plain: [Vec<ThetaPointDim4<L>>; 2],
    /// The gluing images of the four full-order generators `g0` to `g3` of
    /// `F1` and of `F2_dual`: the points the plain chains start from.
    pub gluing_images: [[ThetaPointDim4<L>; 4]; 2],
    /// As in [`HdVerifyReport`].
    pub fell_back: [bool; 2],
}

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
    /// `lambda = q^-1 mod 2^r`.
    lambda: u128,
    /// The doubling chains of `R_com`, `S_com`, `phi R` and `phi S`, the
    /// base points of the kernel generators, computed once.
    tables: [FixedBase<L>; 4],
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
    point_matrix_product_k_fixed(m_full, &s.tables, s.lambda, mask, &s.e_com, &s.e_chal)
}

/// One self-derived half-chain: the gluing chain plus the optimal-strategy plain
/// chain (kernels + codomains). Returns `None` if any step is uncomputable.
/// The third value is `true` when `kernel_points` was
/// [`KernelPoints::TwoWithImageCheck`] and the half ran with four points.
fn run_half<L: HdLevel, P: Profiler>(
    s: &Setup<L>,
    m_full: &[[u128; 8]; 8],
    m_glue: &[[i64; 8]; 8],
    dual: bool,
    kernel_points: KernelPoints,
    prof: &mut P,
) -> Option<(KaniGluingChainHalf<L>, StrategyChain<L>, bool)> {
    let chain = KaniGluingChainHalf::new_prof(
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
        prof,
    )?;
    let t0 = prof.now();
    let basis = b_kpp(s, m_full);
    let tb = prof.now();
    prof.sub_add(14, tb - t0);
    let start = ThetaStructureDim4::<L>::new(chain.codomain_null().clone());
    let n_plain = (L::E1 - s.m as u32 - 1) as usize;
    let t1 = prof.now();
    prof.sub_add(15, t1 - tb);
    prof.add(SLOT_SETUP, t1 - t0);
    match kernel_points {
        KernelPoints::Four => {
            let post: [ThetaPointDim4<L>; 4] = core::array::from_fn(|i| chain.evaluate(&basis[i]));
            prof.add(SLOT_GLUING, prof.now() - t1);
            let sc = run_strategy_chain_prof(&start, &post, n_plain, prof)?;
            Some((chain, sc, false))
        }
        KernelPoints::TwoWithImageCheck => {
            let post: [ThetaPointDim4<L>; 4] = core::array::from_fn(|i| chain.evaluate(&basis[i]));
            prof.add(SLOT_GLUING, prof.now() - t1);
            let (sc, fell_back) =
                run_strategy_chain_two_point_checked_prof(&start, &post, n_plain, prof)?;
            Some((chain, sc, fell_back))
        }
        KernelPoints::TwoWithoutImageCheck => {
            let pair = [chain.evaluate(&basis[0]), chain.evaluate(&basis[1])];
            prof.add(SLOT_GLUING, prof.now() - t1);
            let sc = run_strategy_chain_two_point_prof(&start, &pair, n_plain, prof)?;
            Some((chain, sc, false))
        }
    }
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
        // The dual-step weights as products of the other coordinates of the
        // codomain null (1/C up to a common scalar), no inversion.
        let inv: [Fp2<L>; 16] = crate::hd::field::products_of_others(cprev);
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
    build_setup_prof(sig, &mut NoProfiler)
}

/// [`build_setup`] with a [`Profiler`]: the challenge and response recovery
/// and the commitment basis to `SLOT_BASIS`, the norm equation, the
/// canonical torsion, the product null point and the chain points to
/// `SLOT_SETUP`.
fn build_setup_prof<L: HdLevel, P: Profiler>(
    sig: &HdSignature<L>,
    prof: &mut P,
) -> Result<Setup<L>, HdReject> {
    let t0 = prof.now();
    let chal = recover_challenge_prof::<L, P>(
        &sig.a_pk,
        sig.hint_pk_p,
        sig.hint_pk_q,
        sig.chal_limbs,
        prof,
    )
    .ok_or(HdReject::ChallengeRecovery)?;
    let q4 = sig.q.low_u128();
    let rsp = recover_response_prof::<L, P>(
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
        prof,
    )
    .ok_or(HdReject::ResponseRecovery)?;
    let t1 = prof.now();
    prof.add(SLOT_BASIS, t1 - t0);

    // 5b.4: norm equation 2^e - q = a1² + a2² (a1 odd, a2 even); m = v₂(a2).
    let ts = prof.now();
    let (a1u, a2u) = norm_equation_2f_minus_q(L::E_EMBED, &sig.q).ok_or(HdReject::NormEquation)?;
    prof.sub_add(11, prof.now() - ts);
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

    // P_com, the stage-6 input, is the first point of the basis recover_response
    // derived (the same hinted basis was computed a second time here before
    // SECENG-2428).
    let p_com = rsp.p_com.clone();
    let t3 = prof.now();
    prof.add(SLOT_SETUP, t3 - t1);

    // The doubling chains of the four base points, once (SECENG-2429): the
    // 4-torsion points, points_m and every kernel generator read them.
    let mask_r = (1u128 << L::R) - 1;
    let lambda = inverse_mod_pow2(q4, mask_r);
    let bits = L::R as usize;
    let tables = [
        FixedBase::new(&rsp.r_com, bits, &e_com),
        FixedBase::new(&rsp.s_com, bits, &e_com),
        FixedBase::new(&rsp.phi_rsp_r_com, bits, &e_chal),
        FixedBase::new(&rsp.phi_rsp_s_com, bits, &e_chal),
    ];
    let to4 = (L::R - 2) as usize;
    let p1_4 = tables[0].pts[to4].clone();
    let q1_4 = tables[1].pts[to4].clone();
    let (t1, t2, mt) = make_canonical(&p1_4, &q1_4, &mut e_com).ok_or(HdReject::CanonicalCom)?;
    let r2_4 = tables[2].pts[to4].clone();
    let mut s2_4 = tables[3].pts[to4].clone();
    if lamb4 == 3 {
        s2_4 = jac_add(&jac_dbl(&s2_4, &e_chal), &s2_4, &e_chal); // 3·S2_4
    }
    let (u1, _u2, mu) = make_canonical(&r2_4, &s2_4, &mut e_chal).ok_or(HdReject::CanonicalChal)?;

    let e4 = weil4(&t1, &t2, &mut e_com).inv();
    let zero12 = product_theta_dim2(&dim1_null(&t1), &dim1_null(&u1));
    let m0 = m0_from_canon(&mt, &mu);

    let k = (L::R - 3 - m as u32) as usize;
    let points_m: [JacPoint<L>; 4] = core::array::from_fn(|i| tables[i].pts[k].clone());

    let t4 = prof.now();
    prof.add(SLOT_SETUP, t4 - t3);
    prof.sub_add(12, t4 - t3);
    Ok(Setup {
        e_com,
        e_chal,
        points_m,
        p_com,
        lambda,
        tables,
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
/// plain chain, and whether it fell back to four points.
type HalfResult<L> = (KaniGluingChainHalf<L>, StrategyChain<L>, bool);

/// Run both independent half-chains (F1, F2_dual). Serial by default; with the
/// `parallel` feature, on two threads (see below). Either way the results are
/// bit-identical - the halves do not interact until the middle-codomain match.
fn run_both_halves<L: HdLevel, P: Profiler>(
    s: &Setup<L>,
    m1: &[[u128; 8]; 8],
    mg1: &[[i64; 8]; 8],
    m2: &[[u128; 8]; 8],
    mg2: &[[i64; 8]; 8],
    kernel_points: KernelPoints,
    prof: &mut P,
) -> Option<(HalfResult<L>, HalfResult<L>)> {
    let r1 = run_half(s, m1, mg1, false, kernel_points, prof)?;
    let r2 = run_half(s, m2, mg2, true, kernel_points, prof)?;
    Some((r1, r2))
}

/// Self-derive both half-chains, then run stages 5 (middle-codomain match) and
/// 6 (HD-image). Accept only if BOTH pass.
fn self_derived_check<L: HdLevel, P: Profiler>(
    s: &Setup<L>,
    opts: &HdVerifyOptions,
    prof: &mut P,
) -> (Result<(), HdReject>, HdVerifyReport) {
    let none = HdVerifyReport::default();
    let t0 = prof.now();
    let mask = (1u128 << L::R) - 1;
    let (m1, m2) = match starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask) {
        Some(m) => m,
        None => return (Err(HdReject::ChainFailed), none),
    };
    let mg1 = gluing_bc_dim4_f1(s.a1, s.a2, s.q4, s.m, &m1);
    let mg2 = gluing_bc_dim4_f2(s.a1, s.a2, s.q4, s.m, &m2);
    let t1 = prof.now();
    prof.add(SLOT_SETUP, t1 - t0);
    prof.sub_add(13, t1 - t0);

    let ((chain1, sc1, fb1), (chain2, sc2, fb2)) =
        match run_both_halves(s, &m1, &mg1, &m2, &mg2, opts.kernel_points, prof) {
            Some(r) => r,
            None => return (Err(HdReject::ChainFailed), none),
        };
    // From here on the report carries the halves' fall-back flags, whatever
    // the verdict.
    let report = HdVerifyReport {
        fell_back: [fb1, fb2],
    };

    // Stage 5: the projective middle-codomain match.
    let t1 = prof.now();
    let (c1, c2) = match (sc1.last_codomain(), sc2.last_codomain()) {
        (Some(c1), Some(c2)) => (c1, c2),
        _ => return (Err(HdReject::ChainFailed), report),
    };
    let matches = middle_codomain_matches(c1, c2);
    let t2 = prof.now();
    prof.add(SLOT_COMPARE, t2 - t1);
    if !matches {
        return (Err(HdReject::MiddleCodomainMismatch), report);
    }

    // Stage 6: the HD-image condition F(T) = (±a₁P, ±a₂P, *, 0).
    let image_ok = hd_image_check(s, &chain1, &sc1, &chain2, &sc2).unwrap_or(false);
    prof.add(SLOT_IMAGE, prof.now() - t2);
    if !image_ok {
        return (Err(HdReject::HdImageMismatch), report);
    }
    (Ok(()), report)
}

/// Run the end-to-end self-contained Level-1 verification (all 6 FastVerify
/// stages).
///
/// Performs (0/1) the challenge binding, then self-derives stages 2-4 from the
/// signature, and checks (5) the dim-4 middle-codomain match and (6) the
/// HD-image condition. Returns `Ok(())` on accept or the rejection reason. No
/// oracle data is consulted. This is [`hd_verify_with`] at the default
/// options, the four-point verifier.
pub fn hd_verify<L: HdLevel>(sig: &HdSignature<L>) -> Result<(), HdReject> {
    hd_verify_with::<L>(sig, &HdVerifyOptions::default()).0
}

/// [`hd_verify`] with options: the verdict, and the report of the run. The
/// report's fall-back flags are those of the half-chains whatever the
/// verdict once both half-chains have run; a signature rejected before or
/// while they run reports no fall-back.
pub fn hd_verify_with<L: HdLevel>(
    sig: &HdSignature<L>,
    opts: &HdVerifyOptions,
) -> (Result<(), HdReject>, HdVerifyReport) {
    hd_verify_prof(sig, opts, &mut NoProfiler)
}

fn hd_verify_prof<L: HdLevel, P: Profiler>(
    sig: &HdSignature<L>,
    opts: &HdVerifyOptions,
    prof: &mut P,
) -> (Result<(), HdReject>, HdVerifyReport) {
    let setup = match build_setup_prof(sig, prof) {
        Ok(s) => s,
        Err(e) => return (Err(e), HdVerifyReport::default()),
    };
    self_derived_check(&setup, opts, prof)
}

/// [`hd_verify_with`] timed per phase with the caller's cycle counter
/// `clock` (SECENG-2424): the verdict, the report, and the cycles of each
/// slot of [`crate::hd::profile`] (both half-chains added), in the order of
/// `SLOT_NAMES`. The timing calls sit on the phase boundaries of the
/// shipped code; the computation is the same. Decoding the wire format is
/// the caller's and is not in the slots.
pub fn hd_verify_profiled<L: HdLevel>(
    sig: &HdSignature<L>,
    opts: &HdVerifyOptions,
    clock: &dyn Fn() -> u64,
) -> (Result<(), HdReject>, HdVerifyReport, [u64; SLOTS]) {
    let (verdict, report, slots, _) = hd_verify_profiled_full::<L>(sig, opts, clock);
    (verdict, report, slots)
}

/// [`hd_verify_profiled`] with the sub-slots of [`crate::hd::profile::SUB_NAMES`]
/// as a fourth value (SECENG-2428); the slots are the same.
#[allow(clippy::type_complexity)] // the verdict, the report, the slots and the sub-slots.
pub fn hd_verify_profiled_full<L: HdLevel>(
    sig: &HdSignature<L>,
    opts: &HdVerifyOptions,
    clock: &dyn Fn() -> u64,
) -> (
    Result<(), HdReject>,
    HdVerifyReport,
    [u64; SLOTS],
    [u64; SUB_SLOTS],
) {
    let mut prof = CycleProfiler::new(clock);
    let (verdict, report) = hd_verify_prof(sig, opts, &mut prof);
    (verdict, report, prof.slots, prof.subs)
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
    let (chain1, sc1, _) = run_half(&s, &m1, &mg1, false, KernelPoints::Four, &mut NoProfiler)?;
    let (chain2, sc2, _) = run_half(&s, &m2, &mg2, true, KernelPoints::Four, &mut NoProfiler)?;

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

/// Diagnostic: the codomain null points of both half-chains under `opts`,
/// for comparing the two-point and four-point chains step by step. Returns
/// `None` if the setup or a half-chain is not computable. The verdict is not
/// evaluated here; [`hd_verify_with`] does that.
pub fn hd_codomains<L: HdLevel>(
    sig: &HdSignature<L>,
    opts: &HdVerifyOptions,
) -> Option<HdCodomains<L>> {
    let s = build_setup(sig).ok()?;
    let mask = (1u128 << L::R) - 1;
    let (m1, m2) = starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask)?;
    let mg1 = gluing_bc_dim4_f1(s.a1, s.a2, s.q4, s.m, &m1);
    let mg2 = gluing_bc_dim4_f2(s.a1, s.a2, s.q4, s.m, &m2);
    let ((chain1, sc1, fb1), (chain2, sc2, fb2)) = run_both_halves(
        &s,
        &m1,
        &mg1,
        &m2,
        &mg2,
        opts.kernel_points,
        &mut NoProfiler,
    )?;
    let basis1 = b_kpp(&s, &m1);
    let basis2 = b_kpp(&s, &m2);
    Some(HdCodomains {
        gluing: [
            chain1.codomain_null().clone(),
            chain2.codomain_null().clone(),
        ],
        plain: [sc1.codomains, sc2.codomains],
        gluing_images: [
            core::array::from_fn(|i| chain1.evaluate(&basis1[i])),
            core::array::from_fn(|i| chain2.evaluate(&basis2[i])),
        ],
        fell_back: [fb1, fb2],
    })
}

/// The setup's derived integers and points (SECENG-2429): what the kernel
/// generators are built from, and the generators themselves. For the record
/// and the tests; not on the verification path.
pub struct HdKernelSetup<L: FpBackend> {
    /// The norm equation's `a1` (odd) and `a2`, and `m = v2(a2)`.
    pub a1: u128,
    pub a2: u128,
    pub m: usize,
    /// `lambda = q^-1 mod 2^r`.
    pub lambda: u128,
    /// `[P1_m, Q1_m, R2_m, S2_m]` of order `2^(m+3)`.
    pub points_m: [JacPoint<L>; 4],
    /// The product null point, the dimension-2 base-change matrix over
    /// `Z/4` and `e_4(T1, T2)^-1`.
    pub zero12: [Fp2<L>; 4],
    pub m0: [[u8; 4]; 4],
    pub e4: Fp2<L>,
    /// The four full-order generators of `F1` and of `F2_dual` (`b_kpp` of
    /// `M1` and of `M2`).
    pub generators: [[TuplePoint4<L>; 4]; 2],
}

/// [`HdKernelSetup`] of a signature; `None` where the verifier rejects
/// before the kernel points.
pub fn hd_kernel_setup<L: HdLevel>(sig: &HdSignature<L>) -> Option<HdKernelSetup<L>> {
    let s = build_setup(sig).ok()?;
    let mask = (1u128 << L::R) - 1;
    let (m1, m2) = starting_two_symplectic_matrices(s.a1, s.a2, s.q4, mask)?;
    let generators = [b_kpp(&s, &m1), b_kpp(&s, &m2)];
    Some(HdKernelSetup {
        a1: s.a1,
        a2: s.a2,
        m: s.m,
        lambda: s.lambda,
        points_m: s.points_m.clone(),
        zero12: s.zero12.clone(),
        m0: s.m0,
        e4: s.e4.clone(),
        generators,
    })
}
