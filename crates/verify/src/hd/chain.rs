//! Driving a sequence of plain dimension-4 `(2,2,2,2)`-isogeny steps into a
//! half-chain, and the middle-codomain match that is the heart of SQIsignHD
//! verification.
//!
//! A full SQIsignHD chain `F = F2 ∘ F1` is built as two half-chains `F1` and
//! `F2_dual`. Each half-chain is one gluing step (Phase 4) followed by a run of
//! plain `IsogenyDim4` steps. Because a plain step's codomain depends only on
//! its kernel (see [`crate::hd::isogeny`]), this driver folds
//! [`IsogenyDim4::from_kernel`] over the per-step kernels: step `k`'s codomain
//! is the theta structure on which step `k+1`'s kernel lives, so applying the
//! steps in order reproduces the whole sequence of codomains. The kernels
//! themselves are produced by the reference's optimal-strategy point
//! doublings and pushforwards (not part of this phase).
//!
//! The gluing step is **not** implemented here (Phase 4); this driver covers
//! the plain-step portion of each half-chain.

use crate::FpBackend;

use crate::hd::arith::hadamard;
use crate::hd::isogeny::IsogenyDim4;
use crate::hd::point::ThetaPointDim4;

/// Apply a sequence of plain isogeny steps and return the final codomain theta
/// null point. Returns `None` if any step's codomain is not computable (the
/// non-generic case handled by Phase 4). No heap.
#[inline]
pub fn run_half_chain<L: FpBackend>(
    kernels: &[[ThetaPointDim4<L>; 4]],
) -> Option<ThetaPointDim4<L>> {
    let mut last: Option<ThetaPointDim4<L>> = None;
    for k8 in kernels {
        let iso = IsogenyDim4::from_kernel(k8)?;
        last = Some(iso.codomain_null().clone());
    }
    last
}

/// Like [`run_half_chain`] but writes each step's codomain theta null point into
/// `out` (no heap; the caller supplies a buffer at least `kernels.len()` long).
/// Returns the number of codomains written, or `None` if a step is not
/// computable. Intended for validation/inspection; the hot path uses
/// [`run_half_chain`].
#[inline]
pub fn run_half_chain_collect<L: FpBackend>(
    kernels: &[[ThetaPointDim4<L>; 4]],
    out: &mut [ThetaPointDim4<L>],
) -> Option<usize> {
    assert!(out.len() >= kernels.len(), "output buffer too small");
    for (i, k8) in kernels.iter().enumerate() {
        let iso = IsogenyDim4::from_kernel(k8)?;
        out[i] = iso.codomain_null().clone();
    }
    Some(kernels.len())
}

/// The middle-codomain check (`verify_middle_codomain`): `F1`'s last codomain
/// theta null point must equal the Hadamard of `F2_dual`'s last codomain theta
/// null point, **projectively**.
#[inline]
pub fn middle_codomain_matches<L: FpBackend>(
    f1_last: &ThetaPointDim4<L>,
    f2_dual_last: &ThetaPointDim4<L>,
) -> bool {
    let hc2 = ThetaPointDim4::new(hadamard(f2_dual_last.coords()));
    f1_last.projective_eq(&hc2)
}
