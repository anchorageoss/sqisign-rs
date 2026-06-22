//! Phase 5b.5 - dimension-1 theta structures and the product-theta wrappers
//! that form the dimension-4 gluing's domain.
//!
//! The HD gluing maps from a **product** of lower-dimensional theta structures:
//! a dim-1 theta structure on each of `E_com` and `E_chal`, combined into a
//! dim-2 product (the domain of the `m` dim-2 `(2,2)`-isogeny steps), and then
//! two copies of the resulting dim-2 codomain combined into a dim-4 product
//! (the gluing's domain). This module ports those structures
//! (`theta_structures/Theta_dim1.py`, `montgomery_theta.py`, and the
//! `product_theta_point*` helpers of `Theta_dim2.py` / `theta_helpers_dim4.py`).
//!
//! # What is reused vs new
//!
//! The dim-2 `(2,2)`-theta module in `sqisign-verify` (`crate::theta`) is the
//! *same mathematical object* as the gluing's dim-2 sub-chain, but it is built
//! over **elliptic products** with the gluing's own theta-null convention - not
//! the HD dim-1 convention `(X+Z : X-Z)` keyed to a canonical 4-torsion basis.
//! So the dim-1 null and the product (tensor) maps are **new thin wrappers**;
//! the only reuse is the `Fp2` field arithmetic. The heavy reuse target - the
//! dim-2 `(2,2)`-isogeny *chain* itself - is the strategy loop (Phase 5b.6),
//! not this sub-phase.
//!
//! # The maps (all theta points are projective, defined up to a global scalar)
//!
//! * **dim-1 null** from a 4-torsion point `P = (X : Z)`:
//!   `(X + Z, X - Z)` (`torsion_to_theta_null_point`).
//! * **dim-1 point** `(X : Z)` with null `(a, b)`:
//!   `(a·(X - Z), b·(X + Z))`, and the identity `(0 : 0)` maps to the null
//!   (`montgomery_point_to_theta_point`).
//! * **dim-2 product** of dim-1 theta points `t, u`:
//!   `P[k] = t[k mod 2]·u[k div 2]` (the Kronecker product `t ⊗ u`).
//! * **dim-4 product** of dim-2 theta points `s₁, s₂`:
//!   `P[k] = s₁[k mod 4]·s₂[k div 4]` (the Kronecker product `s₁ ⊗ s₂`).

use crate::{Fp2, FpBackend};

/// A dimension-1 theta structure: an elliptic curve point group in the theta
/// model, characterised by its 2-coordinate theta null point. Built from a
/// canonical `4`-torsion point `P` (with the partner `Q` satisfying `Q[0] = -1`
/// in the SQIsignHD convention).
#[derive(Clone, Debug)]
pub struct ThetaStructureDim1<L: FpBackend> {
    null: [Fp2<L>; 2],
}

impl<L: FpBackend> ThetaStructureDim1<L> {
    /// The theta null point `(X + Z, X - Z)` of the dim-1 structure induced by
    /// the 4-torsion point `P = (X : Z)` (`torsion_to_theta_null_point`).
    #[inline]
    pub fn from_torsion(x: &Fp2<L>, z: &Fp2<L>) -> Self {
        Self {
            null: [x.add(z), x.sub(z)],
        }
    }

    /// Build directly from a known null point `(a, b)`.
    #[inline]
    pub fn from_null(null: [Fp2<L>; 2]) -> Self {
        Self { null }
    }

    /// The theta null point `(a, b)`.
    #[inline]
    pub fn null(&self) -> &[Fp2<L>; 2] {
        &self.null
    }

    /// Map a Montgomery point `(X : Z)` to its dim-1 theta point
    /// `(a·(X - Z), b·(X + Z))` (`montgomery_point_to_theta_point`). The
    /// identity `(0 : 0)` maps to the theta null point.
    #[inline]
    pub fn montgomery_to_theta(&self, x: &Fp2<L>, z: &Fp2<L>) -> [Fp2<L>; 2] {
        if bool::from(x.ct_is_zero()) && bool::from(z.ct_is_zero()) {
            return self.null.clone();
        }
        let (a, b) = (&self.null[0], &self.null[1]);
        [a.mul(&x.sub(z)), b.mul(&x.add(z))]
    }
}

/// The dim-2 product theta point `t ⊗ u`: `P[k] = t[k mod 2]·u[k div 2]`
/// (`Theta_dim2.product_theta_point`). Applied to the two dim-1 **null** points
/// it gives the `ProductThetaStructureDim2` null.
#[inline]
pub fn product_theta_dim2<L: FpBackend>(t: &[Fp2<L>; 2], u: &[Fp2<L>; 2]) -> [Fp2<L>; 4] {
    core::array::from_fn(|k| t[k & 1].mul(&u[(k >> 1) & 1]))
}

/// The dim-4 product theta point `s₁ ⊗ s₂`: `P[k] = s₁[k mod 4]·s₂[k div 4]`
/// (`theta_helpers_dim4.product_theta_point_dim2_dim4`). Applied to the two
/// dim-2 **null** points it gives the `ProductThetaStructureDim2To4` null.
#[inline]
pub fn product_theta_dim2to4<L: FpBackend>(s1: &[Fp2<L>; 4], s2: &[Fp2<L>; 4]) -> [Fp2<L>; 16] {
    core::array::from_fn(|k| s1[k & 3].mul(&s2[(k >> 2) & 3]))
}

/// Convenience: the dim-2 product **null** point of two dim-1 structures
/// (`ProductThetaStructureDim2`).
#[inline]
pub fn product_null_dim2<L: FpBackend>(
    t1: &ThetaStructureDim1<L>,
    t2: &ThetaStructureDim1<L>,
) -> [Fp2<L>; 4] {
    product_theta_dim2(t1.null(), t2.null())
}
