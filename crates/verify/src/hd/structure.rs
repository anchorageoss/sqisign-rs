//! [`ThetaStructureDim4`]: an abelian fourfold in the level-2 theta model,
//! characterised by its theta null point, with the precomputed constants and
//! the doubling formula.

use crate::{Fp2, FpBackend};

use crate::hd::arith::{
    hadamard, hadamard_assign, pointwise_mul_assign, pointwise_square, pointwise_square_assign,
    to_squared_theta_assign,
};
use crate::hd::point::{ThetaPointDim4, THETA_DIM4_N as N};

/// A dimension-4 theta structure: its theta null point (the theta point at the
/// origin) plus the inverses precomputed for the doubling formula.
///
/// Construction does not precompute; call [`Self::precompute`] (or use
/// [`Self::double`], which asserts it was done) before doubling.
#[derive(Clone, Debug)]
pub struct ThetaStructureDim4<L: FpBackend> {
    null_point: ThetaPointDim4<L>,
    precomputed: bool,
    /// `1 / O[k]` up to a common scalar (the products of the other
    /// coordinates of the null point).
    inv_null: [Fp2<L>; N],
    /// `1 / H(S(O))[k]` up to a common scalar, likewise.
    inv_dual_sq: [Fp2<L>; N],
}

impl<L: FpBackend> ThetaStructureDim4<L> {
    /// Build a theta structure from its theta null point.
    #[inline]
    pub fn new(null_point: ThetaPointDim4<L>) -> Self {
        Self {
            null_point,
            precomputed: false,
            inv_null: core::array::from_fn(|_| Fp2::zero()),
            inv_dual_sq: core::array::from_fn(|_| Fp2::zero()),
        }
    }

    /// The theta null point.
    #[inline]
    pub fn null_point(&self) -> &ThetaPointDim4<L> {
        &self.null_point
    }

    /// The theta null point (the identity of the group law). Alias of
    /// [`Self::null_point`], matching the sage `.zero()` accessor.
    #[inline]
    pub fn zero(&self) -> &ThetaPointDim4<L> {
        &self.null_point
    }

    /// The Hadamard-transformed theta structure, whose null point is `H(O)`.
    ///
    /// The middle-codomain check in verification compares
    /// `C1.zero()` against `C2.hadamard().zero()` projectively; this provides
    /// the `C2.hadamard()` side. Because `H` is an involution up to the scalar
    /// `16`, for a codomain `C2` (whose standard null point is `H(O₂)`) this
    /// returns a structure whose null point is projectively `O₂`.
    #[inline]
    pub fn hadamard(&self) -> Self {
        Self::new(ThetaPointDim4::new(hadamard(self.null_point.coords())))
    }

    /// `true` iff this structure supports the projective doubling formula:
    /// every null-point coordinate and every squared-dual coordinate `H(S(O))`
    /// is non-zero. When this fails the sage reference falls back to a random
    /// symplectic base change; that fallback is out of scope for this phase.
    #[inline]
    pub fn has_suitable_doubling(&self) -> bool {
        let o = self.null_point.coords();
        let hso = hadamard(&pointwise_square(o));
        o.iter().all(|c| !bool::from(c.ct_is_zero()))
            && hso.iter().all(|c| !bool::from(c.ct_is_zero()))
    }

    /// Precompute the inverses used by [`Self::double`]:
    /// `inv_null = 1/O` and `inv_dual_sq = 1/H(S(O))`.
    ///
    /// Returns `false` (and leaves the structure un-precomputed) if the
    /// structure does not support doubling (see [`Self::has_suitable_doubling`]).
    pub fn precompute(&mut self) -> bool {
        if self.precomputed {
            return true;
        }
        let o = self.null_point.coords();
        let so = pointwise_square(o);
        let hso = hadamard(&so);

        let suitable = o.iter().all(|c| !bool::from(c.ct_is_zero()))
            && hso.iter().all(|c| !bool::from(c.ct_is_zero()));
        if !suitable {
            return false;
        }

        // The doubling weights as products of the other coordinates, no
        // inversion: the doubled point is the same up to a common scalar.
        self.inv_null = crate::hd::field::products_of_others(o);
        self.inv_dual_sq = crate::hd::field::products_of_others(&hso);
        self.precomputed = true;
        true
    }

    /// Double a theta point on this structure: `P ↦ 2·P`.
    ///
    /// Ports `ThetaPointDim4.double` (suitable-null-point path) from the sage
    /// reference:
    /// `2P = O⁻¹ ⊙ H( H(S(O))⁻¹ ⊙ S(H(S(P))) )`,
    /// where `S` is coordinate-wise squaring, `H` the Hadamard transform, and
    /// `⊙` coordinate-wise multiplication. The result is projective.
    ///
    /// Requires [`Self::precompute`] to have been called.
    #[inline]
    pub fn double(&self, p: &ThetaPointDim4<L>) -> ThetaPointDim4<L> {
        let mut q = p.clone();
        self.double_assign(&mut q);
        q
    }

    /// [`Self::double`] in place: `p <- 2·p`.
    #[inline]
    pub fn double_assign(&self, p: &mut ThetaPointDim4<L>) {
        debug_assert!(
            self.precomputed,
            "ThetaStructureDim4::double called before precompute()"
        );
        // U_χ(P) = S(H(S(P))), then ⊙ 1/H(S(O)), H, ⊙ 1/O, all in place.
        let a = p.coords_mut();
        to_squared_theta_assign(a);
        pointwise_square_assign(a);
        pointwise_mul_assign(a, &self.inv_dual_sq);
        hadamard_assign(a);
        pointwise_mul_assign(a, &self.inv_null);
    }

    /// Iterated doubling: `P ↦ 2ⁿ·P`. Returns `P` unchanged for `n == 0`.
    #[inline]
    pub fn double_iter(&self, p: &ThetaPointDim4<L>, n: u32) -> ThetaPointDim4<L> {
        let mut out = p.clone();
        self.double_iter_assign(&mut out, n);
        out
    }

    /// [`Self::double_iter`] in place: `p <- 2ⁿ·p`.
    #[inline]
    pub fn double_iter_assign(&self, p: &mut ThetaPointDim4<L>, n: u32) {
        for _ in 0..n {
            self.double_assign(p);
        }
    }
}
