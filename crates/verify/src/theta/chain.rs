//! The (2^n, 2^n)-isogeny chain driver (spec Algorithm `Isogeny22Chain`):
//! gluing, `n - 2` generic steps along a balanced strategy, and splitting.

pub use super::gluing::GeneralCouplePoint;
use super::gluing::{start_chain, GluingChangeCoordMatrix, STACK};
use super::isogeny::{theta_isogeny_compute, ThetaIsogeny};
use super::splitting::splitting_to_elliptic_product;
pub use super::splitting::MAX_POINTS;
use super::structure::{
    hadamard, hadamard_assign, invert_point, pointwise_product, pointwise_product_assign,
    to_squared_theta, to_squared_theta_assign,
};
use super::{
    ChainMode, ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints, ThetaPoint,
    ThetaStructure,
};
use crate::fp::{Fp2, FpBackend};
use subtle::Choice;

/// Upper bound on the surfaces a chain can record for its return trip:
/// `e' - 2` for the largest torsion exponent in use (`e' = 321` at level
/// V), rounded up.
pub const MAX_RECORDED_STEPS: usize = 320;

/// What the return trip through a chain needs from the way out (P24): for
/// each intermediate surface `A_j` (`j = 1..=len`, `A_1` the gluing's
/// codomain), the coordinate-wise inverse of its regular theta null point,
/// and the gluing's change of coordinates.
///
/// The dual of a `(2, 2)`-step `f: A -> B` at a point `Q` of `B` with dual
/// coordinates `q` is `f_hat(Q) = H(q^2) / theta_A(0)` in regular
/// coordinates of `A`: `f_hat o f = [2]`, and the duplication formula on
/// `A` reads `theta_A(2P) theta_A(0) = H(theta_B^dual(f(P))^2)` (spec
/// Algorithm 4.114 factored through the step). A step back is therefore
/// four squarings, two Hadamard transforms and four multiplications per
/// point, with no kernel and no codomain arithmetic. The last step back is
/// the dual of the gluing, `A_1 -> E_1 x E_2`, by the same formula in the
/// coordinates the gluing adapted to its kernel: `theta_ad(P) = M (x_1 x_2 :
/// x_1 z_2 : z_1 x_2 : z_1 z_2)`, whose null point is `M (1 : 0 : 0 : 0)`,
/// the first column of `M`.
///
/// Both formulas divide by a null point, so the return trip is defined
/// only when every recorded null point has four non-zero coordinates
/// ([`Self::all_nonzero`]); a zero would mean an intermediate product of
/// elliptic curves, which the forward checks do not exclude.
pub struct DualCache<L: FpBackend> {
    inv_null: [ThetaPoint<L>; MAX_RECORDED_STEPS],
    len: usize,
    matrix: GluingChangeCoordMatrix<L>,
    all_nonzero: Choice,
}

impl<L: FpBackend> Default for DualCache<L> {
    fn default() -> Self {
        Self::new()
    }
}

impl<L: FpBackend> DualCache<L> {
    /// An empty cache.
    pub fn new() -> Self {
        Self {
            inv_null: core::array::from_fn(|_| ThetaPoint::zero()),
            len: 0,
            matrix: GluingChangeCoordMatrix {
                m: core::array::from_fn(|_| core::array::from_fn(|_| Fp2::zero())),
            },
            all_nonzero: Choice::from(1),
        }
    }

    /// Number of recorded surfaces.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether nothing was recorded.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Whether every recorded null point, the gluing's adapted null point
    /// included, has four non-zero coordinates, so that the return trip is
    /// defined.
    pub fn all_nonzero(&self) -> bool {
        bool::from(self.all_nonzero & !Self::has_zero(&self.adapted_null()))
    }

    /// The gluing's change of coordinates.
    pub fn matrix(&self) -> &GluingChangeCoordMatrix<L> {
        &self.matrix
    }

    fn has_zero(p: &ThetaPoint<L>) -> Choice {
        p.x.ct_is_zero() | p.y.ct_is_zero() | p.z.ct_is_zero() | p.t.ct_is_zero()
    }

    /// Record a surface by its regular theta null point. `false` if the
    /// cache is full.
    fn push_null(&mut self, null: &ThetaPoint<L>) -> bool {
        if self.len >= MAX_RECORDED_STEPS {
            return false;
        }
        self.all_nonzero &= !Self::has_zero(null);
        self.inv_null[self.len] = invert_point(null);
        self.len += 1;
        true
    }

    /// The null point of the product `E_1 x E_2` in the gluing's adapted
    /// coordinates, `M (1 : 0 : 0 : 0)`.
    fn adapted_null(&self) -> ThetaPoint<L> {
        ThetaPoint {
            x: self.matrix.m[0][0].clone(),
            y: self.matrix.m[1][0].clone(),
            z: self.matrix.m[2][0].clone(),
            t: self.matrix.m[3][0].clone(),
        }
    }

    /// Pull `pts`, given in dual coordinates on the codomain of the step
    /// after the last recorded surface, back through the duals of the
    /// recorded generic steps. On return `pts` are in dual coordinates on
    /// `A_1`, the gluing's codomain. Requires [`Self::all_nonzero`].
    pub fn pull_back(&self, pts: &mut [ThetaPoint<L>]) {
        for inv in self.inv_null[..self.len].iter().rev() {
            for p in pts.iter_mut() {
                // regular coordinates on the previous surface, then dual
                to_squared_theta_assign(p);
                pointwise_product_assign(p, inv);
                hadamard_assign(p);
            }
        }
    }

    /// Whether the dual of the gluing carries `q` (dual coordinates on
    /// `A_1`) to the couple point `expected` on `E_1 x E_2`: `theta_ad(f_hat
    /// Q) theta_ad(0) = H(q^2)`, compared projectively in the adapted
    /// coordinates, where `theta_ad(expected) = M (x_1 x_2 : x_1 z_2 : z_1
    /// x_2 : z_1 z_2)` (the Segre embedding is injective, so this is
    /// equality on both Kummer lines). `false` if either side is the zero
    /// vector. Requires [`Self::all_nonzero`].
    pub fn gluing_dual_equals(&self, q: &ThetaPoint<L>, expected: &ThetaCouplePoint<L>) -> bool {
        let lhs = to_squared_theta(q);
        let rhs = pointwise_product(
            &self.adapted_null(),
            &self.matrix.couple_point_to_theta(expected),
        );
        let lhs_zero =
            lhs.x.ct_is_zero() & lhs.y.ct_is_zero() & lhs.z.ct_is_zero() & lhs.t.ct_is_zero();
        let rhs_zero =
            rhs.x.ct_is_zero() & rhs.y.ct_is_zero() & rhs.z.ct_is_zero() & rhs.t.ct_is_zero();
        bool::from(!(lhs_zero | rhs_zero) & proj_equal_choice(&lhs, &rhs))
    }
}

/// Projective equality of two theta points, every pair of coordinates
/// cross-multiplied (correct also when a first coordinate is zero); the
/// zero vector is projectively equal to everything.
pub fn proj_equal_choice<L: FpBackend>(a: &ThetaPoint<L>, b: &ThetaPoint<L>) -> Choice {
    let ca = [&a.x, &a.y, &a.z, &a.t];
    let cb = [&b.x, &b.y, &b.z, &b.t];
    let mut eq = Choice::from(1);
    for i in 0..4 {
        for j in (i + 1)..4 {
            eq &= ca[i].mul(cb[j]).ct_equal(&ca[j].mul(cb[i]));
        }
    }
    eq
}

/// Result of a chain: the codomain product and the images of the pushed
/// points, in order.
pub struct ChainOutput<L: FpBackend> {
    /// `E3 x E4`, with the factors selected by the mode in canonical form.
    /// A factor the mode does not select is left at the default curve, and
    /// the corresponding image components are meaningless.
    pub codomain: ThetaCoupleCurve<L>,
    /// Images of the pushed points, `Some` for the first `points.len()`.
    pub images: [Option<ThetaCouplePoint<L>>; MAX_POINTS],
}

/// A theta structure with the images of the pushed points, as they are
/// after some number of steps: the structure's dual coordinates.
pub struct ChainStage<L: FpBackend> {
    /// The structure after the step.
    pub structure: ThetaStructure<L>,
    /// The pushed points, in the structure's dual coordinates.
    pub points: [ThetaPoint<L>; MAX_POINTS],
}

/// The generic steps of a chain along a balanced strategy: `steps` more
/// `(2, 2)`-isogenies from `structure`, whose kernel stack (`theta_q1`,
/// `theta_q2` with the remaining doublings in `todo`, top at `current`)
/// and pushed points `pts` are given. The first step reads regular
/// coordinates, later steps dual coordinates. With `verify`, every step
/// checks its kernel except (`skip_last_check`) the last one, whose check
/// the splitting performs. `snapshot_after` copies the state after that
/// many steps (counted from 1) into `snapshot`. `record`, if given, receives
/// the regular null points of the codomains of the first `until` steps.
#[allow(clippy::too_many_arguments)]
fn generic_steps<L: FpBackend>(
    steps: u16,
    mut step: ThetaIsogeny<L>,
    theta_q1: &mut [ThetaPoint<L>; STACK],
    theta_q2: &mut [ThetaPoint<L>; STACK],
    todo: &mut [u16; STACK],
    mut current: isize,
    pts: &mut [ThetaPoint<L>],
    verify: bool,
    skip_last_check: bool,
    snapshot_after: u16,
    snapshot: Option<&mut ChainStage<L>>,
    first_reads_dual: bool,
    progress: &mut u16,
    record: Option<(&mut DualCache<L>, u16)>,
) -> Option<ThetaStructure<L>> {
    let mut snapshot = snapshot;
    let mut record = record;
    let mut i = 1u16;
    *progress = 0;
    while current >= 0 && todo[current as usize] != 0 {
        let is_not_step2 = i != 1 || first_reads_dual;
        while todo[current as usize] != 1 {
            debug_assert!(todo[current as usize] >= 2);
            current += 1;
            if current as usize >= STACK {
                return None;
            }
            let c = current as usize;
            let num_dbls = todo[c - 1] / 2;
            theta_q1[c] = step.codomain.double_iter(&theta_q1[c - 1], num_dbls);
            theta_q2[c] = step.codomain.double_iter(&theta_q2[c - 1], num_dbls);
            todo[c] = todo[c - 1] - num_dbls;
        }
        let c = current as usize;
        let verify_step = verify && !(skip_last_check && i == steps);
        step = theta_isogeny_compute(&theta_q1[c], &theta_q2[c], is_not_step2, verify_step)?;
        *progress = i;
        for p in pts.iter_mut() {
            step.eval_assign(p, is_not_step2);
        }
        for j in 0..c {
            step.eval_assign(&mut theta_q1[j], is_not_step2);
            step.eval_assign(&mut theta_q2[j], is_not_step2);
            debug_assert!(todo[j] > 0);
            todo[j] -= 1;
        }
        current -= 1;
        if i == snapshot_after {
            if let Some(snap) = snapshot.as_deref_mut() {
                snap.structure = step.codomain.clone();
                for (dst, src) in snap.points.iter_mut().zip(pts.iter()) {
                    *dst = src.clone();
                }
            }
        }
        // the return trip's data for this step's codomain: its regular
        // null point, the Hadamard transform of the dual one
        if let Some((cache, until)) = record.as_mut() {
            if i <= *until && !cache.push_null(&hadamard(&step.codomain.dual_null_point())) {
                return None;
            }
        }
        i += 1;
    }
    if i != steps + 1 {
        return None;
    }
    debug_assert_eq!(current, -1);
    Some(step.codomain)
}

/// Run the chain up to (not including) the splitting step. `snapshot`, if
/// given, receives the state after `n - 1` steps (the penultimate
/// structure with the pushed points in its dual coordinates). `record`, if
/// given, receives the return trip's data for `A_1 .. A_(n-2)`.
#[allow(clippy::too_many_arguments)]
fn chain_core<L: FpBackend>(
    n: u16,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
    general: &[GeneralCouplePoint<L>],
    mode: ChainMode,
    check_last: bool,
    snapshot: Option<&mut ChainStage<L>>,
    record: Option<&mut DualCache<L>>,
) -> Option<(ThetaStructure<L>, [ThetaPoint<L>; MAX_POINTS], usize)> {
    if n < 2 || points.len() + general.len() > MAX_POINTS {
        return None;
    }
    let mut record = record;
    e12.e1.normalize_a24();
    e12.e2.normalize_a24();
    let verify = mode == ChainMode::Verify;

    let mut space = 1usize;
    let mut i = 1u16;
    while i < n {
        space += 1;
        i *= 2;
    }
    if space > STACK {
        return None;
    }
    let mut todo = [0u16; STACK];
    todo[0] = n;
    let mut pts: [ThetaPoint<L>; MAX_POINTS] = core::array::from_fn(|_| ThetaPoint::zero());
    let npts = points.len() + general.len();

    // gluing
    let stage = start_chain(
        &mut todo,
        e12,
        ker,
        points,
        general,
        &mut pts[..npts],
        verify,
    )?;
    let mut theta_q1 = stage.theta_q1;
    let mut theta_q2 = stage.theta_q2;
    // the entry consumed by the gluing is dropped; for n = 1 nothing is left
    let current: isize = stage.current as isize - 1;

    // the return trip's data for the gluing: its codomain's regular null
    // point (`A_1`) and its change of coordinates
    if let Some(cache) = record.as_mut() {
        cache.matrix = stage.first_step.matrix.clone();
        if !cache.push_null(&stage.first_step.codomain) {
            return None;
        }
    }

    // the generic structure after the gluing, in regular coordinates
    let step = ThetaIsogeny {
        codomain: ThetaStructure {
            inv_dual_null_point: invert_point(&stage.first_step.codomain),
            dbl_data: ThetaPoint::zero(),
            inv_sqr_null_point: invert_point(&to_squared_theta(&stage.first_step.codomain)),
            precomputation: true,
        },
    };
    // the gluing was step 1; the snapshot after n - 1 steps is after n - 2
    // generic steps
    let codomain = generic_steps(
        n - 1,
        step,
        &mut theta_q1,
        &mut theta_q2,
        &mut todo,
        current,
        &mut pts[..npts],
        verify,
        !check_last,
        n - 2,
        snapshot,
        false,
        &mut 0u16,
        // the surfaces `A_2 .. A_(n-2)`: the return trip starts on `A_(n-1)`
        record.map(|c| (c, n.saturating_sub(3))),
    )?;
    Some((codomain, pts, npts))
}

/// A `(2^n, 2^n)`-isogeny chain that does not end on a product: the final
/// structure and the images of the pushed points (dual coordinates), plus
/// the state after `n - 1` steps. Every step's kernel is checked. Used by
/// the ECLIPSE verifier, whose `Phi_1` ends on
/// a Jacobian and whose `Phi_2` is the same chain one step short. `None` if
/// `n < 3`, if the kernel is malformed or if a step's check fails.
pub fn theta_chain_two_stage<L: FpBackend>(
    n: u16,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
) -> Option<(ChainStage<L>, ChainStage<L>)> {
    if n < 3 {
        return None;
    }
    let mut mid = ChainStage {
        structure: ThetaStructure::from_null_point(&ThetaPoint::zero()),
        points: core::array::from_fn(|_| ThetaPoint::zero()),
    };
    let (structure, points, _) = chain_core(
        n,
        e12,
        ker,
        points,
        &[],
        ChainMode::Verify,
        true,
        Some(&mut mid),
        None,
    )?;
    Some((mid, ChainStage { structure, points }))
}

/// [`theta_chain_two_stage`] with, in addition, the return trip's data
/// recorded into `cache` (P24): the regular null points of `A_1 ..
/// A_(n-2)` and the gluing's change of coordinates, with which
/// [`DualCache::pull_back`] evaluates the dual of the first `n - 1` steps
/// on points of `A_(n-1)` and [`DualCache::gluing_dual_equals`] the dual
/// of the gluing. `None` also when `n - 2 > MAX_RECORDED_STEPS`.
pub fn theta_chain_two_stage_recorded<L: FpBackend>(
    n: u16,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
    cache: &mut DualCache<L>,
) -> Option<(ChainStage<L>, ChainStage<L>)> {
    if n < 3 || (n as usize) - 2 > MAX_RECORDED_STEPS {
        return None;
    }
    let mut mid = ChainStage {
        structure: ThetaStructure::from_null_point(&ThetaPoint::zero()),
        points: core::array::from_fn(|_| ThetaPoint::zero()),
    };
    let (structure, points, _) = chain_core(
        n,
        e12,
        ker,
        points,
        &[],
        ChainMode::Verify,
        true,
        Some(&mut mid),
        Some(cache),
    )?;
    debug_assert_eq!(cache.len(), n as usize - 2);
    Some((mid, ChainStage { structure, points }))
}

/// [`theta_chain_two_stage`] with, in addition, general couple points (both
/// components non-zero) pushed through the chain; their images follow the
/// special points' in both stages' `points`. Used by the literal
/// Algorithm 3 test, whose step 5 needs `Phi(P_vk, -P_sig)`.
pub fn theta_chain_two_stage_general<L: FpBackend>(
    n: u16,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
    general: &[GeneralCouplePoint<L>],
) -> Option<(ChainStage<L>, ChainStage<L>)> {
    if n < 3 {
        return None;
    }
    let mut mid = ChainStage {
        structure: ThetaStructure::from_null_point(&ThetaPoint::zero()),
        points: core::array::from_fn(|_| ThetaPoint::zero()),
    };
    let (structure, points, _) = chain_core(
        n,
        e12,
        ker,
        points,
        general,
        ChainMode::Verify,
        true,
        Some(&mut mid),
        None,
    )?;
    Some((mid, ChainStage { structure, points }))
}

/// A `(2^n, 2^n)`-isogeny chain starting from a theta structure rather than
/// a product: the kernel generators `k1`, `k2` (order `2^(n+2)`, dual
/// coordinates of `structure`) and the pushed points are theta points; the
/// first step reads its inputs as regular coordinates (`first_reads_dual`
/// false, the convention of the step after a gluing) or as dual ones (the
/// convention of every later step); the chain then proceeds as usual and ends with the splitting into a product
/// of elliptic curves, both factors in canonical form. Every step's kernel
/// is checked. `None` if `n < 2`, on a malformed kernel or if the codomain
/// is not a product; `progress` reports the last step completed (`n + 1`
/// once all steps are done, `n + 2` after the splitting).
pub fn theta_chain_from_structure<L: FpBackend>(
    n: u16,
    structure: &ThetaStructure<L>,
    k1: &ThetaPoint<L>,
    k2: &ThetaPoint<L>,
    points: &[ThetaPoint<L>],
    first_reads_dual: bool,
    progress: &mut u16,
) -> Option<ChainOutput<L>> {
    *progress = 0;
    if n < 2 || points.len() > MAX_POINTS {
        return None;
    }
    let mut space = 1usize;
    let mut i = 1u16;
    while i < n {
        space += 1;
        i *= 2;
    }
    if space > STACK {
        return None;
    }
    let mut todo = [0u16; STACK];
    todo[0] = n;
    let mut theta_q1: [ThetaPoint<L>; STACK] = core::array::from_fn(|_| ThetaPoint::zero());
    let mut theta_q2: [ThetaPoint<L>; STACK] = core::array::from_fn(|_| ThetaPoint::zero());
    theta_q1[0] = k1.clone();
    theta_q2[0] = k2.clone();
    let mut pts: [ThetaPoint<L>; MAX_POINTS] = core::array::from_fn(|_| ThetaPoint::zero());
    let npts = points.len();
    for (dst, src) in pts.iter_mut().zip(points.iter()) {
        *dst = src.clone();
    }
    let step = ThetaIsogeny {
        codomain: structure.clone(),
    };
    let codomain = generic_steps(
        n,
        step,
        &mut theta_q1,
        &mut theta_q2,
        &mut todo,
        0,
        &mut pts[..npts],
        true,
        true,
        0,
        None,
        first_reads_dual,
        progress,
        None,
    )?;
    *progress = n + 1;
    let dual = codomain.dual_null_point();
    let (e34, images) = splitting_to_elliptic_product(&dual, &pts[..npts], ChainMode::Both)?;
    *progress = n + 2;
    Some(ChainOutput {
        codomain: e34,
        images,
    })
}

/// Compute the `(2^n, 2^n)`-isogeny from `E1 x E2` with kernel
/// `[4] <ker.t1, ker.t2>` (the generators have order `2^(n+2)`) and push
/// `points` through it. Each pushed point must have a zero component. The
/// input curves get their `A24` normalised. `None` if `n < 2` (a single
/// gluing step has no splitting basis, as in the reference), if the kernel
/// is malformed, or if the codomain is not a product of elliptic curves.
pub fn theta_chain_compute_and_eval<L: FpBackend>(
    n: u16,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
    mode: ChainMode,
) -> Option<ChainOutput<L>> {
    let (codomain, pts, npts) = chain_core(n, e12, ker, points, &[], mode, false, None, None)?;
    let dual = codomain.dual_null_point();
    let (e34, images) = splitting_to_elliptic_product(&dual, &pts[..npts], mode)?;
    Some(ChainOutput {
        codomain: e34,
        images,
    })
}
