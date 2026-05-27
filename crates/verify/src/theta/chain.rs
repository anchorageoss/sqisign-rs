//! splitting to compute a chain of 2-isogenies between elliptic products.

use crate::ec::EcBasis;
use crate::fp::FpBackend;
use rand_core::RngCore;

use super::basis_change::apply_isomorphism;
use super::couple::double_couple_jac_point_iter;
use super::gluing::{gluing_compute, gluing_eval_basis, gluing_eval_point_special_case};
use super::isogeny::{
    theta_isogeny_compute, theta_isogeny_compute_2, theta_isogeny_compute_4, theta_isogeny_eval,
};
use super::splitting::{
    splitting_compute, theta_point_to_montgomery_point, theta_product_structure_to_elliptic_product,
};
use super::theta_structure::{double_iter, theta_precomputation};
use super::{
    ThetaCoupleCurve, ThetaCoupleJacPoint, ThetaCouplePoint, ThetaIsogeny, ThetaKernelCouplePoints,
    ThetaPoint, ThetaStructure, HD_EXTRA_TORSION,
};

/// Maximum stack depth for the doubling strategy tree.
/// `ceil(log2(n)) + 2` for any n up to 512 covers all security levels.
const MAX_SPACE: usize = 16;

/// Maximum number of image points that can be evaluated through the chain.
const MAX_PTS: usize = 4;

#[inline]
#[allow(clippy::too_many_arguments)]
fn theta_chain_compute_impl<L: FpBackend + crate::precomp::LevelPrecomp>(
    n: u32,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    extra_torsion: bool,
    pts_in: &mut [ThetaCouplePoint<L>],
    verify: bool,
    randomize: bool,
    rng: Option<&mut dyn RngCore>,
) -> Option<ThetaCoupleCurve<L>> {
    let num_p = pts_in.len();
    debug_assert!(num_p <= MAX_PTS);

    // Lift the kernel basis to Jacobian coordinates
    let mut bas1 = EcBasis {
        p: ker.t1.p1.clone(),
        q: ker.t2.p1.clone(),
        pmq: ker.t1m2.p1.clone(),
    };
    let mut bas2 = EcBasis {
        p: ker.t1.p2.clone(),
        q: ker.t2.p2.clone(),
        pmq: ker.t1m2.p2.clone(),
    };

    let (jac_t1_p1, jac_t2_p1, _) = crate::ec::basis::lift_basis(&mut bas1, &mut e12.e1);
    let (jac_t1_p2, jac_t2_p2, _) = crate::ec::basis::lift_basis(&mut bas2, &mut e12.e2);

    let xy_t1 = ThetaCoupleJacPoint {
        p1: jac_t1_p1,
        p2: jac_t1_p2,
    };
    let xy_t2 = ThetaCoupleJacPoint {
        p1: jac_t2_p1,
        p2: jac_t2_p2,
    };

    let extra: u32 = if extra_torsion { HD_EXTRA_TORSION } else { 0 };

    // Image points through the chain
    let mut pts: [ThetaPoint<L>; MAX_PTS] = core::array::from_fn(|_| ThetaPoint::default());

    // Compute strategy tree space: ceil(log2(n)) + 1
    let mut space: usize = 1;
    {
        let mut i: u32 = 1;
        while i < n {
            space += 1;
            i = i.saturating_mul(2);
        }
    }
    debug_assert!(space <= MAX_SPACE);

    let mut todo = [0u16; MAX_SPACE];
    todo[0] = (n - 2 + extra) as u16;

    let mut current: i32 = 0;

    // Jacobian kernel points for the doubling strategy before gluing
    let mut jac_q1: [ThetaCoupleJacPoint<L>; MAX_SPACE] =
        core::array::from_fn(|_| ThetaCoupleJacPoint::default());
    let mut jac_q2: [ThetaCoupleJacPoint<L>; MAX_SPACE] =
        core::array::from_fn(|_| ThetaCoupleJacPoint::default());
    jac_q1[0] = xy_t1;
    jac_q2[0] = xy_t2;

    // Doubling strategy: push kernel points down to order 2^1
    while todo[current as usize] != 1 {
        current += 1;
        let prev = (current - 1) as usize;
        let cur = current as usize;
        let t = todo[prev];
        // Gluing is expensive, so use asymmetric splitting near the end
        let num_dbls: u16 = if t >= 16 { t / 2 } else { t - 1 };
        jac_q1[cur] = double_couple_jac_point_iter(&jac_q1[prev], num_dbls as u32, e12);
        jac_q2[cur] = double_couple_jac_point_iter(&jac_q2[prev], num_dbls as u32, e12);
        todo[cur] = t - num_dbls;
    }

    // Theta kernel points (populated after gluing)
    let mut theta_q1: [ThetaPoint<L>; MAX_SPACE] = core::array::from_fn(|_| ThetaPoint::default());
    let mut theta_q2: [ThetaPoint<L>; MAX_SPACE] = core::array::from_fn(|_| ThetaPoint::default());

    // Gluing step
    let first_step_opt = gluing_compute(
        e12,
        &jac_q1[current as usize],
        &jac_q2[current as usize],
        verify,
    );
    let first_step = first_step_opt?;

    for j in 0..num_p {
        pts[j] = gluing_eval_point_special_case(&pts_in[j], &first_step)?;
    }
    for j in 0..current as usize {
        let (tq1, tq2) = gluing_eval_basis(&jac_q1[j], &jac_q2[j], &first_step);
        theta_q1[j] = tq1;
        theta_q2[j] = tq2;
        todo[j] -= 1;
    }

    current -= 1;

    // Set up theta structure for the first codomain
    let mut theta = ThetaStructure {
        null_point: first_step.codomain,
        precomputation: false,
        ..ThetaStructure::default()
    };
    theta_precomputation(&mut theta);

    // Remaining isogeny steps
    let mut step: Option<ThetaIsogeny<L>> = None;
    let mut i: u32 = 1;
    while current >= 0 && todo[current as usize] != 0 {
        // Doubling strategy within theta model
        while todo[current as usize] != 1 {
            current += 1;
            let prev = (current - 1) as usize;
            let cur = current as usize;
            let t = todo[prev];
            let num_dbls = t / 2;
            theta_q1[cur] = double_iter(&theta_q1[prev], &mut theta, num_dbls as u32);
            theta_q2[cur] = double_iter(&theta_q2[prev], &mut theta, num_dbls as u32);
            todo[cur] = t - num_dbls;
        }

        let cur = current as usize;
        let this_step_opt = if i == n - 2 {
            theta_isogeny_compute(&theta, &theta_q1[cur], &theta_q2[cur], false, false, verify)
        } else if i == n - 1 {
            theta_isogeny_compute(&theta, &theta_q1[cur], &theta_q2[cur], true, false, false)
        } else {
            theta_isogeny_compute(&theta, &theta_q1[cur], &theta_q2[cur], false, true, verify)
        };
        let this_step = this_step_opt?;

        for pt in pts[..num_p].iter_mut() {
            *pt = theta_isogeny_eval(&this_step, pt);
        }

        theta = this_step.codomain.clone();

        for j in 0..current as usize {
            theta_q1[j] = theta_isogeny_eval(&this_step, &theta_q1[j]);
            theta_q2[j] = theta_isogeny_eval(&this_step, &theta_q2[j]);
            todo[j] -= 1;
        }

        step = Some(this_step);
        current -= 1;
        i += 1;
    }

    if !extra_torsion {
        if n >= 3 {
            // SAFETY: when n >= 3, todo[0] = n-2 >= 1, so the while loop above
            // executes at least once, setting step = Some(...).
            let last_step = step
                .as_ref()
                .expect("invariant: step must be set when n >= 3");
            theta_q1[0] = theta_isogeny_eval(last_step, &theta_q1[0]);
            theta_q2[0] = theta_isogeny_eval(last_step, &theta_q2[0]);
        }

        // Penultimate step: from 4-torsion
        let step_4 = theta_isogeny_compute_4(&theta, &theta_q1[0], &theta_q2[0], false, false);
        for pt in pts[..num_p].iter_mut() {
            *pt = theta_isogeny_eval(&step_4, pt);
        }
        theta = step_4.codomain.clone();
        theta_q1[0] = theta_isogeny_eval(&step_4, &theta_q1[0]);
        theta_q2[0] = theta_isogeny_eval(&step_4, &theta_q2[0]);

        // Ultimate step: from 2-torsion
        let step_2 = theta_isogeny_compute_2(&theta, &theta_q1[0], &theta_q2[0], true, false);
        for pt in pts[..num_p].iter_mut() {
            *pt = theta_isogeny_eval(&step_2, pt);
        }
        theta = step_2.codomain.clone();
    }

    // Splitting step
    let last_step = splitting_compute(&theta, if extra_torsion { 8 } else { -1 }, randomize, rng)?;

    let split_b = last_step.b.clone();
    let (e1, e2) = theta_product_structure_to_elliptic_product(&split_b)?;

    for j in 0..num_p {
        let p_iso = apply_isomorphism(&last_step.basis_change, &pts[j]);
        let couple = theta_point_to_montgomery_point(&p_iso, &last_step.b)?;
        pts_in[j] = couple;
    }

    Some(ThetaCoupleCurve { e1, e2 })
}

/// Compute a theta chain and evaluate image points.
#[inline]
pub fn theta_chain_compute_and_eval<L: FpBackend + crate::precomp::LevelPrecomp>(
    n: u32,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    extra_torsion: bool,
    pts: &mut [ThetaCouplePoint<L>],
) -> Option<ThetaCoupleCurve<L>> {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, pts, false, false, None)
}

/// Compute a theta chain with verification checks.
///
/// Used in signature verification. Runs additional validity checks at each
/// isogeny step to detect malformed inputs.
#[inline]
pub fn theta_chain_compute_and_eval_verify<L: FpBackend + crate::precomp::LevelPrecomp>(
    n: u32,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    extra_torsion: bool,
    pts: &mut [ThetaCouplePoint<L>],
) -> Option<ThetaCoupleCurve<L>> {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, pts, true, false, None)
}

/// Compute a theta chain with randomized splitting (signing only).
///
/// Applies a random normalization transform to the final splitting step
/// to hide the splitting choice.
#[inline]
pub fn theta_chain_compute_and_eval_randomized<L: FpBackend + crate::precomp::LevelPrecomp>(
    n: u32,
    e12: &mut ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    extra_torsion: bool,
    pts: &mut [ThetaCouplePoint<L>],
    rng: &mut dyn RngCore,
) -> Option<ThetaCoupleCurve<L>> {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, pts, false, true, Some(rng))
}
