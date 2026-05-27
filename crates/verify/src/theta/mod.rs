//!
//! Implements the SQIsign2D-West approach for fast verification by working
//! with 2-dimensional isogenies in theta coordinates on abelian surfaces.

use crate::ec::{EcBasis, EcCurve, EcPoint, JacPoint};
use crate::fp::{Fp2, FpBackend};

pub mod basis_change;
pub mod chain;
pub mod couple;
pub mod gluing;
pub mod isogeny;
pub mod splitting;
pub mod theta_structure;

/// Additional bits of 2-power torsion consumed by the (2,2)-isogeny chain
/// beyond the target isogeny degree. Computing a degree-2ⁿ isogeny via
/// the theta model requires 2ⁿ⁺ᴴᴰ_ᴱˣᵀᴿᴬ_ᵀᴼᴿˢᴵᴼᴺ torsion points.
pub const HD_EXTRA_TORSION: u32 = 2;

/// A point on an elliptic product E1 x E2 in (X : Z) coordinates.
#[derive(Clone, Debug)]
pub struct ThetaCouplePoint<L: FpBackend> {
    pub p1: EcPoint<L>,
    pub p2: EcPoint<L>,
}

/// A triple (T1, T2, T1-T2) of couple points forming a kernel.
#[derive(Clone, Debug)]
pub struct ThetaKernelCouplePoints<L: FpBackend> {
    pub t1: ThetaCouplePoint<L>,
    pub t2: ThetaCouplePoint<L>,
    pub t1m2: ThetaCouplePoint<L>,
}

/// A point on an elliptic product E1 x E2 in (X : Y : Z) Jacobian coordinates.
#[derive(Clone, Debug)]
pub struct ThetaCoupleJacPoint<L: FpBackend> {
    pub p1: JacPoint<L>,
    pub p2: JacPoint<L>,
}

/// An elliptic product E1 x E2.
#[derive(Clone, Debug)]
pub struct ThetaCoupleCurve<L: FpBackend> {
    pub e1: EcCurve<L>,
    pub e2: EcCurve<L>,
}

/// An elliptic product E1 x E2 with torsion bases B1, B2.
#[derive(Clone, Debug)]
pub struct ThetaCoupleCurveWithBasis<L: FpBackend> {
    pub e1: EcCurve<L>,
    pub e2: EcCurve<L>,
    pub b1: EcBasis<L>,
    pub b2: EcBasis<L>,
}

/// A point in the theta model with projective coordinates (x : y : z : t).
#[derive(Clone, Debug)]
pub struct ThetaPoint<L: FpBackend> {
    pub x: Fp2<L>,
    pub y: Fp2<L>,
    pub z: Fp2<L>,
    pub t: Fp2<L>,
}

/// A compact theta point with two coordinates, used when components repeat.
#[derive(Clone, Debug)]
pub struct ThetaPointCompact<L: FpBackend> {
    pub x: Fp2<L>,
    pub y: Fp2<L>,
}

/// A theta structure: null point plus 8 precomputed 𝔽p² values for
/// efficient doubling and (2,2)-isogeny computation.
#[derive(Clone, Debug)]
pub struct ThetaStructure<L: FpBackend> {
    pub null_point: ThetaPoint<L>,
    pub precomputation: bool,

    // Precomputed from to_squared_theta(null_point) = (XX, YY, ZZ, TT)
    pub cap_xyz0: Fp2<L>, // XX * YY * ZZ
    pub cap_yzt0: Fp2<L>, // YY * ZZ * TT
    pub cap_xzt0: Fp2<L>, // XX * ZZ * TT
    pub cap_xyt0: Fp2<L>, // XX * YY * TT

    // Precomputed from null_point = (x, y, z, t) directly
    pub xyz0: Fp2<L>, // x * y * z
    pub yzt0: Fp2<L>, // y * z * t
    pub xzt0: Fp2<L>, // x * z * t
    pub xyt0: Fp2<L>, // x * y * t
}

/// A 2×2 𝔽p² matrix for the action-by-translation in gluing.
#[derive(Clone, Debug)]
pub struct TranslationMatrix<L: FpBackend> {
    pub g00: Fp2<L>,
    pub g01: Fp2<L>,
    pub g10: Fp2<L>,
    pub g11: Fp2<L>,
}

/// A 4×4 𝔽p² matrix for theta basis changes.
#[derive(Clone, Debug)]
pub struct BasisChangeMatrix<L: FpBackend> {
    pub m: [[Fp2<L>; 4]; 4],
}

/// Precomputed basis change matrix: 4×4 of u8 indices into FP2_CONSTANTS.
#[derive(Clone, Debug)]
pub struct PrecompBasisChangeMatrix {
    pub m: [[u8; 4]; 4],
}

/// A gluing (2,2) theta isogeny from an elliptic product.
#[derive(Clone, Debug)]
pub struct ThetaGluing<L: FpBackend> {
    pub domain: ThetaCoupleCurve<L>,
    pub xy_k1_8: ThetaCoupleJacPoint<L>,
    pub image_k1_8: ThetaPointCompact<L>,
    pub basis_change: BasisChangeMatrix<L>,
    pub precomputation: ThetaPoint<L>,
    pub codomain: ThetaPoint<L>,
}

/// A standard (2,2) theta isogeny between theta structures.
#[derive(Clone, Debug)]
pub struct ThetaIsogeny<L: FpBackend> {
    pub t1_8: ThetaPoint<L>,
    pub t2_8: ThetaPoint<L>,
    pub hadamard_bool_1: bool,
    pub hadamard_bool_2: bool,
    pub domain: ThetaStructure<L>,
    pub precomputation: ThetaPoint<L>,
    pub codomain: ThetaStructure<L>,
}

/// A splitting isomorphism from a theta structure back to an elliptic product.
#[derive(Clone, Debug)]
pub struct ThetaSplitting<L: FpBackend> {
    pub basis_change: BasisChangeMatrix<L>,
    pub b: ThetaStructure<L>,
}

impl<L: FpBackend> Default for ThetaPoint<L> {
    #[inline]
    fn default() -> Self {
        Self {
            x: Fp2::zero(),
            y: Fp2::zero(),
            z: Fp2::zero(),
            t: Fp2::zero(),
        }
    }
}

impl<L: FpBackend> Default for ThetaPointCompact<L> {
    #[inline]
    fn default() -> Self {
        Self {
            x: Fp2::zero(),
            y: Fp2::zero(),
        }
    }
}

impl<L: FpBackend> Default for ThetaStructure<L> {
    #[inline]
    fn default() -> Self {
        Self {
            null_point: ThetaPoint::default(),
            precomputation: false,
            cap_xyz0: Fp2::zero(),
            cap_yzt0: Fp2::zero(),
            cap_xzt0: Fp2::zero(),
            cap_xyt0: Fp2::zero(),
            xyz0: Fp2::zero(),
            yzt0: Fp2::zero(),
            xzt0: Fp2::zero(),
            xyt0: Fp2::zero(),
        }
    }
}

impl<L: FpBackend> Default for BasisChangeMatrix<L> {
    #[inline]
    fn default() -> Self {
        Self {
            m: core::array::from_fn(|_| core::array::from_fn(|_| Fp2::zero())),
        }
    }
}

impl<L: FpBackend> Default for TranslationMatrix<L> {
    #[inline]
    fn default() -> Self {
        Self {
            g00: Fp2::zero(),
            g01: Fp2::zero(),
            g10: Fp2::zero(),
            g11: Fp2::zero(),
        }
    }
}

impl<L: FpBackend> Default for ThetaCouplePoint<L> {
    #[inline]
    fn default() -> Self {
        Self {
            p1: EcPoint::identity(),
            p2: EcPoint::identity(),
        }
    }
}

impl<L: FpBackend> Default for ThetaCoupleJacPoint<L> {
    #[inline]
    fn default() -> Self {
        Self {
            p1: JacPoint::identity(),
            p2: JacPoint::identity(),
        }
    }
}

impl<L: FpBackend> Default for ThetaCoupleCurve<L> {
    #[inline]
    fn default() -> Self {
        Self {
            e1: EcCurve::default(),
            e2: EcCurve::default(),
        }
    }
}
