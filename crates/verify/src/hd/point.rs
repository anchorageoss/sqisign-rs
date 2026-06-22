//! [`ThetaPointDim4`]: a point in the dimension-4 level-2 theta model, given by
//! 16 projective 𝔽p² coordinates.

use crate::{Fp2, FpBackend};

/// Number of theta coordinates in dimension 4 (`2⁴`).
pub const THETA_DIM4_N: usize = 16;

/// A point in the dimension-4 theta model: 16 projective 𝔽p² coordinates,
/// indexed by `k = i₀ + 2·i₁ + 4·i₂ + 8·i₃`.
///
/// The coordinates are **projective**: scaling all 16 by a common non-zero
/// `Fp2` value yields the same point (see [`Self::projective_eq`]).
#[derive(Clone, Debug)]
pub struct ThetaPointDim4<L: FpBackend> {
    coords: [Fp2<L>; THETA_DIM4_N],
}

impl<L: FpBackend> ThetaPointDim4<L> {
    /// Construct from 16 coordinates.
    #[inline]
    pub fn new(coords: [Fp2<L>; THETA_DIM4_N]) -> Self {
        Self { coords }
    }

    /// The 16 coordinates.
    #[inline]
    pub fn coords(&self) -> &[Fp2<L>; THETA_DIM4_N] {
        &self.coords
    }

    /// The pivot: index of the first non-zero coordinate (capped at 15, to
    /// match the sage `ThetaPointDim4.__eq__` convention). For any genuine
    /// (non-zero) theta point this is the smallest `k` with `P[k] ≠ 0`.
    ///
    /// Not constant-time. The dimension-4 verification path operates entirely
    /// on public data (signature, public key), so leaking the pivot index is
    /// not a concern; constant-time hardening is deferred to a later phase.
    #[inline]
    pub fn pivot(&self) -> usize {
        let mut k0 = 0;
        while k0 < THETA_DIM4_N - 1 && bool::from(self.coords[k0].ct_is_zero()) {
            k0 += 1;
        }
        k0
    }

    /// Projective equality.
    ///
    /// Returns `true` iff `self` and `other` represent the same projective
    /// point, i.e. `self[l]·other[k₀] == other[l]·self[k₀]` for every `l`,
    /// where `k₀` is the pivot of `self`. This mirrors
    /// `ThetaPointDim4.__eq__` in the sage reference exactly. **Never** compare
    /// theta points coordinate-by-coordinate.
    #[inline]
    pub fn projective_eq(&self, other: &Self) -> bool {
        let k0 = self.pivot();
        let p = &self.coords;
        let q = &other.coords;
        for l in 0..THETA_DIM4_N {
            let lhs = p[l].mul(&q[k0]);
            let rhs = q[l].mul(&p[k0]);
            if !bool::from(lhs.ct_equal(&rhs)) {
                return false;
            }
        }
        true
    }

    /// Return a new point with every coordinate multiplied by `lambda` (a
    /// change of projective representative; the point itself is unchanged).
    #[inline]
    pub fn scale(&self, lambda: &Fp2<L>) -> Self {
        Self {
            coords: core::array::from_fn(|k| self.coords[k].mul(lambda)),
        }
    }

    /// Pivot-normalise to the canonical representative whose pivot coordinate
    /// is `1`. Returns `(pivot, normalised_point)`.
    ///
    /// This matches the `normalized` / `pivot` fields emitted by the Phase 0
    /// oracle: divide every coordinate by the first non-zero one. Two equal
    /// projective points have the same pivot and identical normalised
    /// coordinates.
    #[inline]
    pub fn normalize(&self) -> (usize, Self) {
        let k0 = self.pivot();
        let inv = crate::hd::field::inv(&self.coords[k0]);
        (k0, self.scale(&inv))
    }

    /// `true` iff every coordinate is zero (a degenerate, non-projective
    /// point). Genuine theta points are never zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|c| bool::from(c.ct_is_zero()))
    }
}
