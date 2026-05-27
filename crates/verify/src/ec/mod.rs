//!
//! Provides Montgomery curve types and arithmetic in projective (X:Z)
//! coordinates, Jacobian (X:Y:Z) coordinates, isogeny types, and
//! torsion-basis types. All types are generic over the security level.

use crate::fp::{Fp2, FpBackend};
use crate::params::SecurityLevel;
use zeroize::Zeroize;

pub mod basis;
pub mod curve;
pub mod isogeny;
pub mod jacobian;
pub mod pairing;
pub mod point;

pub use isogeny::{EcKps2, EcKps4};

/// Projective point on the Kummer line in Montgomery (X:Z) coordinates.
#[derive(Clone, Debug)]
pub struct EcPoint<L: SecurityLevel> {
    pub x: Fp2<L>,
    pub z: Fp2<L>,
}

/// Jacobian point (X:Y:Z) on a Montgomery curve, representing the
/// affine point `(X/Z^2, Y/Z^3)`.
#[derive(Clone, Debug)]
pub struct JacPoint<L: SecurityLevel> {
    pub x: Fp2<L>,
    pub y: Fp2<L>,
    pub z: Fp2<L>,
}

/// Three components (u, v, w) encoding the (X:Z) coordinates of both
/// the addition and subtraction of two distinct points:
/// `P+Q = (u-v : w)` and `P-Q = (u+v : w)`.
#[derive(Clone, Debug)]
pub struct AddComponents<L: SecurityLevel> {
    pub u: Fp2<L>,
    pub v: Fp2<L>,
    pub w: Fp2<L>,
}

/// A basis `{P, Q, P-Q}` of a torsion subgroup, in (X:Z) coordinates.
#[derive(Clone, Debug)]
pub struct EcBasis<L: SecurityLevel> {
    pub p: EcPoint<L>,
    pub q: EcPoint<L>,
    pub pmq: EcPoint<L>,
}

/// Montgomery curve `y^2 = x^3 + (A/C)x^2 + x` in projective form.
///
/// Stores `(A : C)` and a cached copy of `(A+2C : 4C)` for fast
/// doubling. The `is_a24_computed_and_normalized` flag tracks whether
/// `a24` holds `((A+2C)/(4C) : 1)`.
#[derive(Clone, Debug)]
pub struct EcCurve<L: SecurityLevel> {
    pub a: Fp2<L>,
    pub c: Fp2<L>,
    pub a24: EcPoint<L>,
    pub is_a24_computed_and_normalized: bool,
}

/// Isomorphism of Montgomery curves: `(X:Z) -> (Nx*X + Nz*Z : D*Z)`.
#[derive(Clone, Debug)]
pub struct EcIsomorphism<L: SecurityLevel> {
    pub nx: Fp2<L>,
    pub nz: Fp2<L>,
    pub d: Fp2<L>,
}

/// Even-degree isogeny (degree 2^length).
#[derive(Clone, Debug)]
pub struct EcIsogEven<L: SecurityLevel> {
    pub curve: EcCurve<L>,
    pub kernel: EcPoint<L>,
    pub length: u32,
}

impl<L: FpBackend> Default for EcPoint<L> {
    #[inline]
    fn default() -> Self {
        EcPoint::identity()
    }
}

impl<L: FpBackend> EcPoint<L> {
    /// The point at infinity `(1 : 0)`.
    #[inline]
    pub fn identity() -> Self {
        Self {
            x: Fp2::one(),
            z: Fp2::zero(),
        }
    }

    #[inline]
    pub fn new(x: Fp2<L>, z: Fp2<L>) -> Self {
        Self { x, z }
    }
}

impl<L: FpBackend> Default for JacPoint<L> {
    #[inline]
    fn default() -> Self {
        JacPoint::identity()
    }
}

impl<L: FpBackend> JacPoint<L> {
    /// The identity `(0 : 1 : 0)` in Jacobian coordinates.
    #[inline]
    pub fn identity() -> Self {
        Self {
            x: Fp2::zero(),
            y: Fp2::one(),
            z: Fp2::zero(),
        }
    }

    #[inline]
    pub fn new(x: Fp2<L>, y: Fp2<L>, z: Fp2<L>) -> Self {
        Self { x, y, z }
    }
}

impl<L: FpBackend> Default for EcCurve<L> {
    #[inline]
    fn default() -> Self {
        Self {
            a: Fp2::zero(),
            c: Fp2::one(),
            a24: EcPoint::identity(),
            is_a24_computed_and_normalized: false,
        }
    }
}

impl<L: FpBackend> EcBasis<L> {
    #[inline]
    pub fn new(p: EcPoint<L>, q: EcPoint<L>, pmq: EcPoint<L>) -> Self {
        Self { p, q, pmq }
    }
}

impl<L: FpBackend> AddComponents<L> {
    #[inline]
    pub fn new(u: Fp2<L>, v: Fp2<L>, w: Fp2<L>) -> Self {
        Self { u, v, w }
    }
}

impl<L: SecurityLevel> Zeroize for EcPoint<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.x.zeroize();
        self.z.zeroize();
    }
}

impl<L: SecurityLevel> Zeroize for EcBasis<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.p.zeroize();
        self.q.zeroize();
        self.pmq.zeroize();
    }
}

impl<L: SecurityLevel> Zeroize for EcCurve<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.a.zeroize();
        self.c.zeroize();
        self.a24.zeroize();
        self.is_a24_computed_and_normalized = false;
    }
}
