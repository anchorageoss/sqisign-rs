//! The dimension-4 **gluing** isogeny that starts each half-chain (the Kani
//! embedding), ported from the sage `GluingIsogenyDim4`
//! (`isogenies/gluing_isogeny_dim4.py`).
//!
//! A gluing step differs from a plain [`crate::hd::IsogenyDim4`] step in two ways,
//! both stemming from its codomain having **zero dual theta-null coordinates**:
//!
//! 1. **Codomain.** The kernel is given by five 8-torsion points whose
//!    translation directions are `[1,2,4,8,3]` - four single-bit generators
//!    plus one **diagonal** generator (`3 = bit0 ⊕ bit1`, the sum of the first
//!    two). The diagonal edge lets the spanning-tree reconstruction
//!    (`codomain_dual_null_and_inverse`) cover all 16 indices even though
//!    the zero dual coordinates make some single-bit edges unusable
//!    (zero denominators). The resulting dual null `O` legitimately contains
//!    zeros; the codomain null point is `H(O)` as usual.
//!
//! 2. **Image evaluation.** The plain `image(P) = H(H(S(P)) ⊙ O⁻¹)` formula
//!    divides by `O`, which is impossible where `O` is zero. `special_image`
//!    instead recovers those coordinates from **translates** `P + Tₖ` of `P` by
//!    4-torsion points above the kernel: a coordinate `i` with `O[i] = 0` is
//!    read off from a translate at index `i ⊕ indₖ` where `O[i ⊕ indₖ] ≠ 0`,
//!    correcting for the projective factor between `P` and `P + Tₖ`.
//!
//! The product → dim-4 change of theta coordinates that sets up the gluing
//! domain and kernel (the `N_dim4` Kani base change, derived from `a₁,a₂,q,m`)
//! is part of the chain orchestration; here the gluing domain and kernel are
//! taken as given (the oracle supplies the base-changed inputs at Level 1).

use crate::{Fp2, FpBackend};

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::isogeny::codomain_dual_null_and_inverse;
use crate::hd::point::{ThetaPointDim4, THETA_DIM4_N as N};
use crate::hd::structure::ThetaStructureDim4;

/// The canonical kernel directions of a dim-4 gluing isogeny: four single-bit
/// generators and one diagonal (`3 = bit0 ⊕ bit1`).
pub const GLUING_KERNEL_DIRS: [usize; 5] = [1, 2, 4, 8, 3];

/// Maximum number of translates supported by [`GluingIsogenyDim4::special_image`]
/// (the gluing uses two).
const MAX_TRANS: usize = 4;

/// A computed dimension-4 gluing isogeny: its codomain structure, the codomain
/// **dual** theta null point (with zeros), and the per-index inverse dual
/// constants used by [`Self::special_image`].
#[derive(Clone, Debug)]
pub struct GluingIsogenyDim4<L: FpBackend> {
    codomain: ThetaStructureDim4<L>,
    /// Codomain dual theta null point `O` (some coordinates are zero).
    dual_null: [Fp2<L>; N],
    /// `1 / O[i]` where `O[i] ≠ 0`, else `0`.
    inv_dual_null: [Fp2<L>; N],
    /// `true` where `O[i] = 0`.
    dual_null_zero: [bool; N],
}

impl<L: FpBackend> GluingIsogenyDim4<L> {
    /// Compute the gluing isogeny from its five kernel 8-torsion theta points
    /// and their translation directions (use [`GLUING_KERNEL_DIRS`]).
    ///
    /// Returns `None` if the dual null point cannot be reconstructed (no root
    /// spans the 4-cube).
    pub fn from_kernel(k8: &[ThetaPointDim4<L>; 5], dirs: &[usize; 5]) -> Option<Self> {
        let hsk: [[Fp2<L>; N]; 5] =
            core::array::from_fn(|k| hadamard(&pointwise_square(k8[k].coords())));

        // One batched inversion yields both the dual null `O` (with its zero
        // coordinates) and its inverse `1/O` for the image precomputation.
        let (dual_null, inv_dual_null, dual_null_zero, _has_zero) =
            codomain_dual_null_and_inverse(&hsk, dirs)?;

        let codomain_null = ThetaPointDim4::new(hadamard(&dual_null));

        Some(Self {
            codomain: ThetaStructureDim4::new(codomain_null),
            dual_null,
            inv_dual_null,
            dual_null_zero,
        })
    }

    /// The codomain theta structure.
    #[inline]
    pub fn codomain(&self) -> &ThetaStructureDim4<L> {
        &self.codomain
    }

    /// The codomain theta null point.
    #[inline]
    pub fn codomain_null(&self) -> &ThetaPointDim4<L> {
        self.codomain.null_point()
    }

    /// Number of zero dual theta-null coordinates of the codomain (6 for the
    /// Level-1 gluing steps).
    #[inline]
    pub fn dual_zero_count(&self) -> usize {
        self.dual_null_zero.iter().filter(|&&z| z).count()
    }

    /// Evaluate the gluing isogeny on a domain theta point `p`, given its
    /// translates `p + Tₖ` (`l_trans`) by 4-torsion points above the kernel and
    /// their indices `l_trans_ind` (e.g. `[1, 2]`). Ports `special_image`.
    pub fn special_image(
        &self,
        p: &ThetaPointDim4<L>,
        l_trans: &[ThetaPointDim4<L>],
        l_trans_ind: &[usize],
    ) -> ThetaPointDim4<L> {
        let nt = l_trans.len();
        assert!(nt <= MAX_TRANS && nt == l_trans_ind.len());

        let hs_p = hadamard(&pointwise_square(p.coords()));
        let mut hsl: [[Fp2<L>; N]; MAX_TRANS] =
            core::array::from_fn(|_| core::array::from_fn(|_| Fp2::zero()));
        for (k, q) in l_trans.iter().enumerate() {
            hsl[k] = hadamard(&pointwise_square(q.coords()));
        }
        let o = &self.dual_null;

        // λₖ⁻¹: the projective factor relating P + Tₖ to the reference
        // representative. Pick any (j, j⊕indₖ) with both HSL[k][j] ≠ 0 and
        // O[j⊕indₖ] ≠ 0, and set λₖ⁻¹ = (HS_P[j⊕indₖ]·O[j]) / (HSL[k][j]·O[j⊕indₖ]).
        let mut lambda_inv: [Fp2<L>; MAX_TRANS] = core::array::from_fn(|_| Fp2::one());
        for (k, ind) in l_trans_ind.iter().enumerate() {
            for j in 0..N {
                let jpk = j ^ ind;
                if !bool::from(hsl[k][j].ct_is_zero()) && !bool::from(o[jpk].ct_is_zero()) {
                    let num = hs_p[jpk].mul(&o[j]);
                    let den = hsl[k][j].mul(&o[jpk]);
                    lambda_inv[k] = num.mul(&den.inv());
                    break;
                }
            }
        }

        // U_χ(f(P)): direct where O[i] ≠ 0, else via a translate.
        let mut u: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        for (i, ui) in u.iter_mut().enumerate() {
            if !self.dual_null_zero[i] {
                *ui = hs_p[i].mul(&self.inv_dual_null[i]);
            } else {
                for (k, ind) in l_trans_ind.iter().enumerate() {
                    let ipk = i ^ ind;
                    if !self.dual_null_zero[ipk] {
                        *ui = self.inv_dual_null[ipk].mul(&hsl[k][ipk]).mul(&lambda_inv[k]);
                        break;
                    }
                }
            }
        }

        ThetaPointDim4::new(hadamard(&u))
    }
}
