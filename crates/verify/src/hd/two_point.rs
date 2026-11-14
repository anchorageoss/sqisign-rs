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
//! of pairs and the same strategy. The codomain null points along the chain
//! and their Hadamard transforms are fixed by [`PI`], so each step computes
//! its dual null point, the inverse and the next doubling constants on the
//! ten orbits of [`PI`] (`isogeny_from_pair`, `hadamard_pi_assign`,
//! `ThetaStructureDim4::precompute_pi`), and the doubling constants only
//! when a doubling follows; a step that meets a vanishing coordinate runs
//! the full computation of the four-point loop
//! ([`IsogenyDim4::from_kernel`]). The codomains and the images are those
//! of the four-point loop, which stays the default; see `KernelPoints` in
//! `self_contained`.
//!
//! The checked entry tests, at the gluing, that the images of `g2`, `g3`
//! are the permutations of those of `g0`, `g1` and that the gluing
//! codomain's null point is fixed by [`PI`], and runs the four-point loop
//! otherwise. In debug builds every later codomain null point is checked to
//! be fixed too; release builds compile that out.

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

/// One index of each of the ten orbits of [`PI`]: the four fixed points,
/// then the lower index of each two-cycle.
pub const ORBIT_REPS: [usize; 10] = [0, 5, 10, 15, 1, 2, 3, 6, 7, 11];

/// The orbit (its position in [`ORBIT_REPS`]) of each index.
pub const ORBIT_OF: [usize; THETA_DIM4_N] = [0, 4, 5, 6, 4, 1, 7, 8, 5, 7, 2, 9, 6, 8, 9, 3];

/// The sixteen coordinates of a vector fixed by [`PI`] from its values on
/// the ten orbits.
#[inline]
pub fn expand_orbits<L: FpBackend>(v: &[Fp2<L>; 10]) -> [Fp2<L>; THETA_DIM4_N] {
    core::array::from_fn(|a| v[ORBIT_OF[a]].clone())
}

/// The Hadamard transform of a vector fixed by [`PI`] coordinate by
/// coordinate, in place; the result is fixed by [`PI`] too. The transform
/// is done on the pairs of index bits that [`PI`] exchanges, (0, 2) then
/// (1, 3). In each stage the two-by-two blocks whose other two bits are
/// equal are symmetric (three distinct outputs: four additions and one
/// doubling), and the remaining two blocks are transposes of each other
/// (one of them computed, eight additions): 32 additions and 4 doublings
/// in all, against 64 additions for [`hadamard_assign`](crate::hd::arith::hadamard_assign).
pub fn hadamard_pi_assign<L: FpBackend>(a: &mut [Fp2<L>; THETA_DIM4_N]) {
    // one stage on bits (lo, hi) with strides (s, t): `sym` are the bases
    // of the two symmetric blocks, `gen` the base of the computed general
    // block and `mirror` the base of its transpose
    let mut stage = |s: usize, t: usize, sym: [usize; 2], gen: usize, mirror: usize| {
        for c in sym {
            let (p, q, r) = (a[c].clone(), a[c + s].clone(), a[c + s + t].clone());
            let sum = p.add(&r);
            let q2 = q.add(&q);
            a[c] = sum.add(&q2);
            a[c + s + t] = sum.sub(&q2);
            let d = p.sub(&r);
            a[c + t] = d.clone();
            a[c + s] = d;
        }
        let (p, q, r, w) = (
            a[gen].clone(),
            a[gen + s].clone(),
            a[gen + t].clone(),
            a[gen + s + t].clone(),
        );
        let (pq, rw, pmq, rmw) = (p.add(&q), r.add(&w), p.sub(&q), r.sub(&w));
        a[gen] = pq.add(&rw);
        a[gen + s] = pmq.add(&rmw);
        a[gen + t] = pq.sub(&rw);
        a[gen + s + t] = pmq.sub(&rmw);
        // the transposed block: its (i, j) entry is the (j, i) entry
        a[mirror] = a[gen].clone();
        a[mirror + s] = a[gen + t].clone();
        a[mirror + t] = a[gen + s].clone();
        a[mirror + s + t] = a[gen + s + t].clone();
    };
    // bits 0 and 2: blocks indexed by bits (1, 3)
    stage(1, 4, [0, 10], 2, 8);
    // bits 1 and 3: blocks indexed by bits (0, 2)
    stage(2, 8, [0, 5], 1, 4);
}

/// The edges of the tree through the ten orbits along which the two-point
/// step reads its codomain: `(k, j)` is the flip of index bit `k` at index
/// `j`, measured by kernel point `k` (`g0` for bit 0, `g1` for bit 1; a
/// flip of bit 2 or 3 is a flip of bit 0 or 1 at the permuted index).
/// The first eight edges are a path over nine orbits,
/// `{0} {1,4} {5} {7,13} {15} {11,14} {10} {2,8} {3,12}`; the ninth hangs
/// the last orbit `{6,9}` on the eighth vertex `{2,8}`. No path visits all
/// ten: every single-bit flip joins one of the four orbits `{1,4}`,
/// `{2,8}`, `{7,13}`, `{11,14}` to one of the six others.
const TREE_EDGES: [(usize, usize); 9] = [
    (0, 0),
    (0, 4),
    (1, 5),
    (1, 13),
    (0, 15),
    (0, 11),
    (1, 10),
    (0, 2),
    (0, 8),
];

/// The orbits of the tree's vertices: the nine of the path in order, then
/// the leaf.
const TREE_ORBITS: [usize; 10] = [0, 4, 1, 8, 3, 9, 2, 5, 6, 7];

/// A vector with ratios `x[v'] / x[v] = num / den` along each edge
/// `v → v'` of the tree, up to a common scalar, on the ten orbits (in the
/// order of [`ORBIT_REPS`]): prefix products of the numerators and suffix
/// products of the denominators along the path, the leaf's denominator
/// folded into the suffixes. 25 multiplications.
fn along_tree<L: FpBackend>(num: &[Fp2<L>; 9], den: &[Fp2<L>; 9]) -> [Fp2<L>; 10] {
    // suf[i] = den[8] * prod_{i <= t < 8} den[t], for the path vertices
    let mut suf: [Fp2<L>; 9] = core::array::from_fn(|_| den[8].clone());
    for i in (0..8).rev() {
        suf[i] = den[i].mul(&suf[i + 1]);
    }
    // pre = prod_{t < i} num[t]
    let mut x: [Fp2<L>; 10] = core::array::from_fn(|_| Fp2::zero());
    x[TREE_ORBITS[0]] = suf[0].clone();
    let mut pre = num[0].clone();
    for i in 1..9 {
        if i == 7 {
            // the leaf, from the eighth vertex: its own suffix is den[7]
            x[TREE_ORBITS[9]] = pre.mul(&den[7]).mul(&num[8]);
        }
        x[TREE_ORBITS[i]] = pre.mul(&suf[i]);
        if i < 8 {
            pre.mul_assign(&num[i]);
        }
    }
    x
}

/// The two-point step from the order-8 pair `[K0, K1]`: the codomain null
/// point and the image weights computed on the ten orbits of [`PI`], from
/// the squared-dual coordinates of `K0` and `K1` alone (those of
/// `PI(K0)`, `PI(K1)` are their permutations). `None` when a coordinate
/// the tree reads vanishes; the caller then runs the full computation.
fn isogeny_from_pair<L: FpBackend>(pair: &[ThetaPointDim4<L>; 2]) -> Option<IsogenyDim4<L>> {
    let hsk: [[Fp2<L>; THETA_DIM4_N]; 2] = core::array::from_fn(|k| {
        let mut a = pair[k].coords().clone();
        crate::hd::arith::to_squared_theta_assign(&mut a);
        a
    });
    // O[j ^ 2^k] / O[j] = hsk[k][j ^ 2^k] / hsk[k][j] for the dual null O
    let from: [Fp2<L>; 9] = core::array::from_fn(|e| {
        let (k, j) = TREE_EDGES[e];
        hsk[k][j].clone()
    });
    let to: [Fp2<L>; 9] = core::array::from_fn(|e| {
        let (k, j) = TREE_EDGES[e];
        hsk[k][j ^ (1 << k)].clone()
    });
    if from
        .iter()
        .chain(to.iter())
        .any(|c| bool::from(c.ct_is_zero()))
    {
        return None;
    }
    let dual_null = along_tree(&to, &from);
    let inv_dual_null = along_tree(&from, &to);
    let mut null = expand_orbits(&dual_null);
    hadamard_pi_assign(&mut null);
    Some(IsogenyDim4::from_parts(
        ThetaPointDim4::new(null),
        expand_orbits(&inv_dual_null),
    ))
}

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
/// [`gluing_images_paired`] holds and the null point of `start`, the gluing
/// codomain, is fixed by [`PI`] ([`null_fixed_by_pi`]): four evaluations
/// through the gluing and six comparisons per half. Otherwise the
/// four-point chain ([`run_strategy_chain`]) from the same images; the flag
/// is `true` in that case. Both chains are `l` plain steps from `start`.
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
    if gluing_images_paired(post) && null_fixed_by_pi(start.null_point()) {
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
    // The doubling constants of a codomain are computed only when a
    // doubling follows. `cur_pi`: the null point of `cur` came from the
    // orbit computation, so its constants can be computed on the orbits.
    let mut cur_pi = false;

    let mut kernels: Vec<[ThetaPointDim4<L>; 4]> = Vec::with_capacity(l);
    let mut codomains: Vec<ThetaPointDim4<L>> = Vec::with_capacity(l);
    let mut image_precomp: Vec<[Fp2<L>; THETA_DIM4_N]> = Vec::with_capacity(l);

    for k in 0..l {
        let mut prev: u32 = level.iter().sum();
        let target = (l - 1 - k) as u32;
        let t0 = prof.now();
        while prev != target {
            let ready = if cur_pi {
                cur.precompute_pi()
            } else {
                cur.precompute()
            };
            if !ready {
                return None;
            }
            let s = strategy[strat_idx];
            level.push(s);
            prev += s;
            // The doublings run in place on a copy of the top pair.
            let mut next = kernel_elements
                .last()
                .expect("non-empty kernel stack")
                .clone();
            for pt in next.iter_mut() {
                cur.double_iter_assign(pt, s);
            }
            kernel_elements.push(next);
            strat_idx += 1;
        }
        let t1 = prof.now();
        prof.add(SLOT_PLAIN_PUSHES, t1 - t0);

        // The descended top-of-stack pair gives the order-8 kernel of this
        // step: the codomain on the orbits, or the full computation from
        // the four points when a coordinate it reads vanishes.
        let top = kernel_elements.last().expect("non-empty kernel stack");
        let ker = kernel_from_pair(top);
        let (iso, on_orbits) = match isogeny_from_pair(top) {
            Some(iso) => (iso, true),
            None => (IsogenyDim4::from_kernel(&ker)?, false),
        };
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

        // Push the remaining saved pairs through the isogeny, in place.
        for elt in kernel_elements.iter_mut() {
            for slot in elt.iter_mut() {
                if !iso.image_assign(slot) {
                    return None;
                }
            }
        }
        let t3 = prof.now();
        prof.add(SLOT_PLAIN_PUSHES, t3 - t2);
        cur = iso.into_codomain();
        cur_pi = on_orbits;
        prof.add(SLOT_PLAIN_CODOMAIN, prof.now() - t3);
    }

    Some(StrategyChain {
        kernels,
        codomains,
        image_precomp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hd::arith::hadamard;
    use crate::params::P324_3;

    type F = Fp2<P324_3>;

    fn value(seed: u64) -> F {
        // a spread of field elements from small integers
        let a = F::from_small(seed.wrapping_mul(0x9e37_79b9).wrapping_add(17) % 1_000_003);
        let b = F::from_small(seed.wrapping_mul(0x85eb_ca6b).wrapping_add(5) % 999_983);
        a.sqr().add(&b.mul(&F::i_element())).add(&b.sqr().sqr())
    }

    #[test]
    fn orbit_tables_are_consistent() {
        for a in 0..THETA_DIM4_N {
            assert_eq!(ORBIT_OF[a], ORBIT_OF[PI[a]]);
            assert!(ORBIT_REPS[ORBIT_OF[a]] == a || ORBIT_REPS[ORBIT_OF[a]] == PI[a]);
        }
        // every tree edge is a single-bit flip between the orbits listed
        let mut seen = [false; 10];
        seen[TREE_ORBITS[0]] = true;
        for (e, &(k, j)) in TREE_EDGES.iter().enumerate() {
            let from = if e < 8 {
                TREE_ORBITS[e]
            } else {
                TREE_ORBITS[7]
            };
            let to = if e < 8 {
                TREE_ORBITS[e + 1]
            } else {
                TREE_ORBITS[9]
            };
            assert_eq!(ORBIT_OF[j], from, "edge {e}");
            assert_eq!(ORBIT_OF[j ^ (1 << k)], to, "edge {e}");
            assert!(!seen[to]);
            seen[to] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn hadamard_of_a_fixed_vector() {
        for seed in 0..50u64 {
            let v10: [F; 10] = core::array::from_fn(|i| value(31 * seed + i as u64));
            let v = expand_orbits(&v10);
            let want = hadamard(&v);
            let mut got = v.clone();
            hadamard_pi_assign(&mut got);
            for a in 0..THETA_DIM4_N {
                assert!(
                    bool::from(got[a].ct_equal(&want[a])),
                    "seed {seed} index {a}"
                );
                assert!(bool::from(got[a].ct_equal(&got[PI[a]])));
            }
        }
    }

    #[test]
    fn along_tree_has_the_edge_ratios() {
        let num: [F; 9] = core::array::from_fn(|i| value(1000 + i as u64));
        let den: [F; 9] = core::array::from_fn(|i| value(2000 + i as u64));
        let x = along_tree(&num, &den);
        for e in 0..9 {
            let from = if e < 8 {
                TREE_ORBITS[e]
            } else {
                TREE_ORBITS[7]
            };
            let to = if e < 8 {
                TREE_ORBITS[e + 1]
            } else {
                TREE_ORBITS[9]
            };
            // x[to] * den = x[from] * num
            assert!(bool::from(
                x[to].mul(&den[e]).ct_equal(&x[from].mul(&num[e]))
            ));
        }
    }
}
