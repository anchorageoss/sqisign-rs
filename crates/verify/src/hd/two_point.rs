//! The two-point kernel stack of the dimension-4 plain chain (SECENG-2421).
//!
//! From the gluing codomain on, the four kernel generators of a half-chain
//! come in two pairs related by one fixed permutation [`PI`] of the 16 theta
//! coordinates, the exchange of index bits 0 and 2 and of bits 1 and 3 under
//! `k = i0 + 2 i1 + 4 i2 + 8 i3`: coordinate `a` of `g2` is coordinate
//! `PI[a]` of `g0`, projectively, and likewise `g3` from `g1`. The loop here
//! doubles and pushes `g0` and `g1` only and forms `g2`, `g3` by [`permute_pi`]
//! wherever a step reads them. It is the walk of
//! [`run_strategy_chain`](crate::hd::strategy::run_strategy_chain) with a stack
//! of pairs; the strategy, the step ([`IsogenyDim4::from_kernel`]), the images
//! and the codomains are those of the four-point loop. The four-point loop
//! stays the default; see `KernelPoints` in `self_contained`.
//!
//! In debug builds every codomain null point from the gluing codomain on is
//! checked to be fixed by [`PI`] ([`null_fixed_by_pi`], six comparisons);
//! release builds compile the check out.

use alloc::vec;
use alloc::vec::Vec;

use crate::{Fp2, FpBackend};

use crate::hd::isogeny::IsogenyDim4;
use crate::hd::point::{ThetaPointDim4, THETA_DIM4_N};
use crate::hd::profile::{NoProfiler, Profiler, SLOT_PLAIN_CODOMAIN, SLOT_PLAIN_PUSHES};
use crate::hd::strategy::{
    optimised_strategy, run_strategy_chain_prof, StrategyChain, STRATEGY_MUL_C,
};
use crate::hd::structure::ThetaStructureDim4;

/// `PI[a]`: the index `a` with bits 0 and 2 exchanged and bits 1 and 3
/// exchanged. An involution; its fixed points are 0, 5, 10 and 15.
pub const PI: [usize; THETA_DIM4_N] = [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15];

/// The six two-cycles of [`PI`], each as `(a, PI[a])` with `a < PI[a]`.
pub const PI_PAIRS: [(usize, usize); 6] = [(1, 4), (2, 8), (3, 12), (6, 9), (7, 13), (11, 14)];

/// The point whose coordinate `a` is coordinate `PI[a]` of `p`. Since `PI`
/// is an involution, this is also the point whose coordinate `PI[a]` is
/// coordinate `a` of `p`.
#[inline]
pub fn permute_pi<L: FpBackend>(p: &ThetaPointDim4<L>) -> ThetaPointDim4<L> {
    let c = p.coords();
    ThetaPointDim4::new(core::array::from_fn(|a| c[PI[a]].clone()))
}

/// The fixed points of [`PI`].
pub const PI_FIXED: [usize; 4] = [0, 5, 10, 15];

/// Whether `o` is fixed by [`PI`] projectively, that is `o[PI[a]] = t o[a]`
/// for all `a` and one `t`, in six comparisons: with a reference pair
/// `(k0, PI[k0])`, `o[a] * o[PI[k0]] == o[PI[a]] * o[k0]` for each of the
/// six pairs `(a, PI[a])` of [`PI_PAIRS`]. The reference `k0` is the lowest
/// fixed point of `PI` with a non-zero coordinate, so that `t = 1` and the
/// six comparisons read `o[a] == o[PI[a]]`, which is the whole condition.
/// When the four fixed coordinates are all zero, `k0` is the pivot of `o`
/// and one more comparison, `o[k0]^2 == o[PI[k0]]^2`, pins `t^2 = 1`, which
/// with the six makes the condition complete in that case too.
pub fn null_fixed_by_pi<L: FpBackend>(o: &ThetaPointDim4<L>) -> bool {
    let c = o.coords();
    let fixed = PI_FIXED
        .iter()
        .copied()
        .find(|&f| !bool::from(c[f].ct_is_zero()));
    let k0 = fixed.unwrap_or_else(|| o.pivot());
    let r = &c[k0];
    let rp = &c[PI[k0]];
    let pairs = PI_PAIRS
        .iter()
        .all(|&(a, b)| bool::from(c[a].mul(rp).ct_equal(&c[b].mul(r))));
    let unit = fixed.is_some() || bool::from(r.sqr().ct_equal(&rp.sqr()));
    pairs && unit
}

/// Whether the gluing images `post` of `g0` to `g3` are paired by [`PI`]:
/// `post[2]` equals `PI` of `post[0]` and `post[3]` equals `PI` of `post[1]`,
/// projectively. This is the per-signature check of the two-point verifier.
pub fn gluing_images_paired<L: FpBackend>(post: &[ThetaPointDim4<L>; 4]) -> bool {
    post[2].projective_eq(&permute_pi(&post[0])) && post[3].projective_eq(&permute_pi(&post[1]))
}

/// The two-point chain from the gluing images `post` of `g0` to `g3` when
/// [`gluing_images_paired`] holds, and otherwise the four-point chain
/// ([`run_strategy_chain`]) from the same images. The flag is `true` in the
/// second case. Both chains are `l` plain steps from `start`.
pub fn run_strategy_chain_two_point_checked<L: FpBackend>(
    start: &ThetaStructureDim4<L>,
    post: &[ThetaPointDim4<L>; 4],
    l: usize,
) -> Option<(StrategyChain<L>, bool)> {
    run_strategy_chain_two_point_checked_prof(start, post, l, &mut NoProfiler)
}

/// [`run_strategy_chain_two_point_checked`] with a [`Profiler`].
pub fn run_strategy_chain_two_point_checked_prof<L: FpBackend, P: Profiler>(
    start: &ThetaStructureDim4<L>,
    post: &[ThetaPointDim4<L>; 4],
    l: usize,
    prof: &mut P,
) -> Option<(StrategyChain<L>, bool)> {
    if gluing_images_paired(post) {
        let pair = [post[0].clone(), post[1].clone()];
        Some((
            run_strategy_chain_two_point_prof(start, &pair, l, prof)?,
            false,
        ))
    } else {
        Some((run_strategy_chain_prof(start, post, l, prof)?, true))
    }
}

/// The step kernel `[K0, K1, PI(K0), PI(K1)]` from the pair `[K0, K1]`.
#[inline]
fn kernel_from_pair<L: FpBackend>(pair: &[ThetaPointDim4<L>; 2]) -> [ThetaPointDim4<L>; 4] {
    [
        pair[0].clone(),
        pair[1].clone(),
        permute_pi(&pair[0]),
        permute_pi(&pair[1]),
    ]
}

/// Derive an `l`-step plain `(2,2,2,2)`-chain from the gluing images of `g0`
/// and `g1` by walking the optimal strategy with a stack of pairs.
///
/// * `start` is the gluing codomain, the structure `basis` lives on.
/// * `basis` are the images of `g0` and `g1`, of order `2^(l+2)` on `start`.
/// * `l` is the number of plain steps.
///
/// Each step's kernel is `[K0, K1, PI(K0), PI(K1)]` from the pair on top of
/// the stack, doubled down to order 8; the saved pairs are pushed through the
/// step's image. The output is the [`StrategyChain`] of the four-point loop,
/// its `kernels` holding the four points handed to each step. Returns `None`
/// where the four-point loop does: a structure that cannot be doubled on, or
/// a step whose codomain or image is not computable.
pub fn run_strategy_chain_two_point<L: FpBackend>(
    start: &ThetaStructureDim4<L>,
    basis: &[ThetaPointDim4<L>; 2],
    l: usize,
) -> Option<StrategyChain<L>> {
    run_strategy_chain_two_point_prof(start, basis, l, &mut NoProfiler)
}

/// [`run_strategy_chain_two_point`] with a [`Profiler`], the slots as in
/// [`run_strategy_chain_prof`](crate::hd::strategy::run_strategy_chain_prof):
/// the doublings and the images of the saved pairs to `SLOT_PLAIN_PUSHES`,
/// `IsogenyDim4::from_kernel` and the next structure's `precompute` to
/// `SLOT_PLAIN_CODOMAIN`.
pub fn run_strategy_chain_two_point_prof<L: FpBackend, P: Profiler>(
    start: &ThetaStructureDim4<L>,
    basis: &[ThetaPointDim4<L>; 2],
    l: usize,
    prof: &mut P,
) -> Option<StrategyChain<L>> {
    debug_assert!(
        null_fixed_by_pi(start.null_point()),
        "the gluing codomain null point is not fixed by PI"
    );
    if l == 0 {
        return Some(StrategyChain {
            kernels: Vec::new(),
            codomains: Vec::new(),
            image_precomp: Vec::new(),
        });
    }
    let strategy = optimised_strategy(l, STRATEGY_MUL_C);

    let mut strat_idx = 0usize;
    let mut level: Vec<u32> = vec![0];
    let mut kernel_elements: Vec<[ThetaPointDim4<L>; 2]> = vec![basis.clone()];
    let mut cur = start.clone();
    let mut cur_ready = cur.precompute();

    let mut kernels: Vec<[ThetaPointDim4<L>; 4]> = Vec::with_capacity(l);
    let mut codomains: Vec<ThetaPointDim4<L>> = Vec::with_capacity(l);
    let mut image_precomp: Vec<[Fp2<L>; THETA_DIM4_N]> = Vec::with_capacity(l);

    for k in 0..l {
        let mut prev: u32 = level.iter().sum();
        let target = (l - 1 - k) as u32;
        let t0 = prof.now();
        while prev != target {
            if !cur_ready {
                return None;
            }
            let s = strategy[strat_idx];
            level.push(s);
            prev += s;
            let src = kernel_elements.last().expect("non-empty kernel stack");
            let next: [ThetaPointDim4<L>; 2] =
                core::array::from_fn(|i| cur.double_iter(&src[i], s));
            kernel_elements.push(next);
            strat_idx += 1;
        }
        let t1 = prof.now();
        prof.add(SLOT_PLAIN_PUSHES, t1 - t0);

        // The descended top-of-stack pair gives the order-8 kernel of this step.
        let ker = kernel_from_pair(kernel_elements.last().expect("non-empty kernel stack"));
        let iso = IsogenyDim4::from_kernel(&ker)?;
        let t2 = prof.now();
        prof.add(SLOT_PLAIN_CODOMAIN, t2 - t1);
        debug_assert!(
            null_fixed_by_pi(iso.codomain_null()),
            "a codomain null point of the two-point chain is not fixed by PI"
        );
        kernels.push(ker);
        codomains.push(iso.codomain_null().clone());
        image_precomp.push(iso.inv_dual_null().clone());

        kernel_elements.pop();
        level.pop();

        // Push the remaining saved pairs through the isogeny.
        for elt in kernel_elements.iter_mut() {
            for slot in elt.iter_mut() {
                *slot = iso.image(slot)?;
            }
        }
        let t3 = prof.now();
        prof.add(SLOT_PLAIN_PUSHES, t3 - t2);
        cur = iso.codomain().clone();
        cur_ready = cur.precompute();
        prof.add(SLOT_PLAIN_CODOMAIN, prof.now() - t3);
    }

    Some(StrategyChain {
        kernels,
        codomains,
        image_precomp,
    })
}
