//! Phase 5b.6 - the dimension-4 optimal-strategy chain loop.
//!
//! A SQIsignHD half-chain `F1` is one gluing step (Phase 4) followed by a run of
//! `L` plain `(2,2,2,2)`-steps. The gluing exposes, on its codomain, four
//! full-order (`2^{L+2}`) kernel generators; the rest of the chain is *derived*
//! from them by an optimal-strategy walk that, at each step, doubles the saved
//! kernel bases down to the order-8 kernel of the current step, computes the
//! plain isogeny ([`IsogenyDim4::from_kernel`]), and pushes the remaining bases
//! through it ([`IsogenyDim4::image`]). This is the dimension-4 analogue of the
//! dim-2 `IsogenyChainDim2` loop and of `isogeny_chain_dim4.py`'s plain part -
//! the bookkeeping that turns a single kernel basis into all `L` per-step
//! kernels (the ~938 doublings / ~1882 pushforwards of a Level-1 half-chain).
//!
//! # Strategy-independence (why a self-consistent completion is fine here)
//!
//! The codomain of step `k` and its order-8 kernel are determined by the chain
//! and the position `k`, **not** by the strategy (which only chooses the
//! doubling/pushforward *path*). So any valid optimal strategy reproduces the
//! same per-step kernels and codomains. Likewise, the per-step kernels are the
//! same subgroup generators regardless of the symplectic completion used to
//! build the starting basis (5b.4 decision) - different completions give the
//! same chain up to a global symplectic change of basis, so the middle-codomain
//! match (the accept/reject invariant) is preserved.
//!
//! # `alloc`
//!
//! The optimal-strategy walk keeps a stack of intermediate kernel bases (depth
//! and size data-dependent, bounded by the chain length). It therefore uses
//! heap `Vec`s - the crate's only `alloc` use, off the constant-time path.

use alloc::vec;
use alloc::vec::Vec;

use crate::{Fp2, FpBackend};

use crate::hd::isogeny::IsogenyDim4;
use crate::hd::point::{ThetaPointDim4, THETA_DIM4_N};
use crate::hd::structure::ThetaStructureDim4;

/// Doubling-to-image cost ratio fed to [`optimised_strategy`]. Measured
/// post-Phase-7 at Level 1 (`tests/hd_primitive_bench.rs`): one dim-4 doubling
/// ≈ 5.7 µs, one isogeny image ≈ 3.6 µs, so a doubling costs ≈ 1.6 image
/// evaluations. (Phase 7's batch inversion made `from_kernel`/`image` cheaper
/// relative to doubling; the previous `1.0` assumed they were equal and so
/// over-spent on doublings.) `from_kernel` is a fixed per-step cost and does not
/// enter the doubling/image trade-off.
const STRATEGY_MUL_C: f64 = 1.6;

/// Optimal `(2,2)`-chain strategy of `n` leaves (Algorithm 60, SIKE spec;
/// `utilities/strategy.py::optimised_strategy`). Returns the `n - 1` internal
/// nodes (each the number of leaves to the right of that node), the sequence
/// the chain loop walks depth-first. `mul_c` is the doubling cost relative to a
/// unit isogeny evaluation.
///
/// Any valid optimal strategy yields the same isogeny chain; the precise cost
/// model only affects efficiency, so a plain `mul_c` is used (the reference's
/// gluing-aware `precompute_strategy_with_first_eval` is for the *whole* chain
/// including the costly gluing first-step, which this loop does not recompute).
pub fn optimised_strategy(n: usize, mul_c: f64) -> Vec<u32> {
    if n <= 1 {
        return Vec::new();
    }
    let eval_c = 1.0_f64;
    // s[i] = strategy (internal-node list) for a subtree of i leaves; c[i] cost.
    let mut s: Vec<Vec<u32>> = Vec::with_capacity(n + 1);
    let mut c: Vec<f64> = Vec::with_capacity(n + 1);
    for _ in 0..=n {
        s.push(Vec::new());
        c.push(0.0);
    }
    for i in 2..=n {
        // b = argmin over 1..i of c[i-b] + c[b] + b*mul_c + (i-b)*eval_c
        // (ties: smallest b, matching Python's min over ascending b).
        let mut best_b = 1usize;
        let mut best_cost = f64::INFINITY;
        for b in 1..i {
            let cost = c[i - b] + c[b] + (b as f64) * mul_c + ((i - b) as f64) * eval_c;
            if cost < best_cost {
                best_cost = cost;
                best_b = b;
            }
        }
        let mut si = Vec::with_capacity(i - 1);
        si.push(best_b as u32);
        si.extend_from_slice(&s[i - best_b]);
        si.extend_from_slice(&s[best_b]);
        s[i] = si;
        c[i] = best_cost;
    }
    core::mem::take(&mut s[n])
}

/// The result of [`run_strategy_chain`]: the per-step order-8 kernels derived by
/// the strategy walk (one `[_; 4]` per plain step, in chain order) and the
/// corresponding per-step codomain theta null points.
pub struct StrategyChain<L: FpBackend> {
    pub kernels: Vec<[ThetaPointDim4<L>; 4]>,
    pub codomains: Vec<ThetaPointDim4<L>>,
    /// Per-step image precomputation `1/O` (in chain order). Lets a caller push
    /// a point through the forward chain (`apply_plain_image`) without rebuilding
    /// the isogenies - used by the stage-6 HD-image evaluation.
    pub image_precomp: Vec<[Fp2<L>; THETA_DIM4_N]>,
}

impl<L: FpBackend> StrategyChain<L> {
    /// The last codomain (the half-chain's middle theta null point), if any.
    #[inline]
    pub fn last_codomain(&self) -> Option<&ThetaPointDim4<L>> {
        self.codomains.last()
    }
}

/// Derive an `L`-step plain `(2,2,2,2)`-chain from a single full-order kernel
/// basis by walking the optimal strategy.
///
/// * `start` is the theta structure the `basis` lives on (the gluing codomain).
/// * `basis` are four kernel generators of order `2^{L+2}` on `start`.
/// * `l` is the number of plain steps.
///
/// Returns the per-step kernels and codomains, or `None` if any step's codomain
/// is not computable (the non-suitable-null case Phase 4 would handle; not
/// expected at Level 1).
pub fn run_strategy_chain<L: FpBackend>(
    start: &ThetaStructureDim4<L>,
    basis: &[ThetaPointDim4<L>; 4],
    l: usize,
) -> Option<StrategyChain<L>> {
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
    let mut kernel_elements: Vec<[ThetaPointDim4<L>; 4]> = vec![basis.clone()];
    // The doubling formula needs the structure's precomputed inverses; a
    // structure with a non-suitable null cannot be doubled on (the base-change
    // fallback, out of scope - return `None`).
    let mut cur = start.clone();
    let mut cur_ready = cur.precompute();

    let mut kernels: Vec<[ThetaPointDim4<L>; 4]> = Vec::with_capacity(l);
    let mut codomains: Vec<ThetaPointDim4<L>> = Vec::with_capacity(l);
    let mut image_precomp: Vec<[Fp2<L>; THETA_DIM4_N]> = Vec::with_capacity(l);

    for k in 0..l {
        let mut prev: u32 = level.iter().sum();
        let target = (l - 1 - k) as u32;
        while prev != target {
            if !cur_ready {
                return None; // non-suitable null met a required doubling
            }
            let s = strategy[strat_idx];
            level.push(s);
            prev += s;
            let src = kernel_elements.last().expect("non-empty kernel stack");
            let next: [ThetaPointDim4<L>; 4] =
                core::array::from_fn(|i| cur.double_iter(&src[i], s));
            kernel_elements.push(next);
            strat_idx += 1;
        }

        // The descended top-of-stack is the order-8 kernel of this step.
        let ker = kernel_elements
            .last()
            .expect("non-empty kernel stack")
            .clone();
        let iso = IsogenyDim4::from_kernel(&ker)?;
        kernels.push(ker);
        codomains.push(iso.codomain_null().clone());
        image_precomp.push(iso.inv_dual_null().clone());

        kernel_elements.pop();
        level.pop();

        // Push the remaining saved bases through the isogeny.
        let mut pushed: Vec<[ThetaPointDim4<L>; 4]> = Vec::with_capacity(kernel_elements.len());
        for elt in &kernel_elements {
            let mut img: [ThetaPointDim4<L>; 4] = core::array::from_fn(|_| elt[0].clone());
            for (i, slot) in img.iter_mut().enumerate() {
                *slot = iso.image(&elt[i])?;
            }
            pushed.push(img);
        }
        kernel_elements = pushed;
        cur = iso.codomain().clone();
        cur_ready = cur.precompute();
    }

    Some(StrategyChain {
        kernels,
        codomains,
        image_precomp,
    })
}
