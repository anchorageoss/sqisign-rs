//! Montgomery curves over `Fp2` in x-only projective coordinates, following
//! the SQIsign round-3 reference (`src/ec`) and spec Section 4.4.
//!
//! Types are generic over the prime `L`. Everything is x-only: a point is
//! `(X : Z)` on the Kummer line, a torsion basis is `(P, Q, P - Q)`, and a
//! curve is `(A : C)` with a cached, optionally normalised `A24 = (A + 2C :
//! 4C)`.

use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;
use zeroize::Zeroize;

pub mod basis;
pub mod curve;
pub mod isogeny;
#[cfg(feature = "compact")]
pub mod jacobian;
pub(crate) mod mp;
pub mod normalize;
pub mod pairing;
pub mod point;

/// Upper bound on the 64-bit words of an order-sized scalar across all
/// parameter sets (`ceil(669 / 64)`).
pub const MAX_ORDER_WORDS: usize = 11;

/// Upper bound on the bit length of scalars fed to the ladders.
pub const MAX_ORDER_BITS: usize = 64 * MAX_ORDER_WORDS;

/// A point on the Kummer line in projective `(X : Z)` coordinates.
#[derive(Clone, Debug)]
pub struct EcPoint<L: Prime> {
    /// Projective `X`.
    pub x: Fp2<L>,
    /// Projective `Z`; zero for the point at infinity.
    pub z: Fp2<L>,
}

/// Jacobian point `(X : Y : Z)` on a Montgomery curve, the affine point
/// `(X/Z^2, Y/Z^3)` (the dimension-4 layer needs full points).
#[cfg(feature = "compact")]
#[derive(Clone, Debug)]
pub struct JacPoint<L: Prime> {
    /// `X`.
    pub x: Fp2<L>,
    /// `Y`.
    pub y: Fp2<L>,
    /// `Z`.
    pub z: Fp2<L>,
}

/// Three components `(u, v, w)` encoding the `(X : Z)` coordinates of both
/// the sum and the difference of two distinct points: `P + Q = (u - v : w)`,
/// `P - Q = (u + v : w)`.
#[cfg(feature = "compact")]
#[derive(Clone, Debug)]
pub struct AddComponents<L: Prime> {
    /// `u`.
    pub u: Fp2<L>,
    /// `v`.
    pub v: Fp2<L>,
    /// `w`.
    pub w: Fp2<L>,
}

/// A basis `(P, Q)` of a torsion subgroup, carried as `x(P), x(Q), x(P - Q)`.
#[derive(Clone, Debug)]
pub struct EcBasis<L: Prime> {
    /// `x(P)`.
    pub p: EcPoint<L>,
    /// `x(Q)`.
    pub q: EcPoint<L>,
    /// `x(P - Q)`.
    pub pmq: EcPoint<L>,
}

/// Montgomery curve `C y^2 = C x^3 + A x^2 + C x`, stored as `(A : C)` with a
/// cache of `A24 = (A + 2C : 4C)` that may be normalised to `((A+2)/4 : 1)`.
#[derive(Clone, Debug)]
pub struct EcCurve<L: Prime> {
    /// Projective `A`.
    pub a: Fp2<L>,
    /// Projective `C`, never zero.
    pub c: Fp2<L>,
    /// `(A + 2C : 4C)`, or `((A + 2C)/4C : 1)` when normalised.
    pub a24: EcPoint<L>,
    /// Whether `a24` holds the normalised form.
    pub is_a24_computed_and_normalized: bool,
}

/// An isogeny of degree `2^length` given by its domain and a kernel generator.
#[derive(Clone, Debug)]
pub struct EcIsogEven<L: Prime> {
    /// Domain curve.
    pub curve: EcCurve<L>,
    /// A generator of the kernel, of order `2^length`.
    pub kernel: EcPoint<L>,
    /// Length of the walk as a 2-isogeny chain; must be even.
    pub length: u32,
}

/// A change of coordinates `(X : Z) -> (a X + b Z : c X + d Z)` between
/// Montgomery models of the same curve.
#[derive(Clone, Debug)]
pub struct EcChangeCoordMatrix<L: Prime> {
    /// Top-left entry.
    pub a: Fp2<L>,
    /// Top-right entry.
    pub b: Fp2<L>,
    /// Bottom-left entry.
    pub c: Fp2<L>,
    /// Bottom-right entry.
    pub d: Fp2<L>,
}

/// Barycentric coordinates `(u, v, w)` of a pair of points: `P + Q = (u - v :
/// w)` and `P - Q = (u + v : w)`.
#[derive(Clone, Debug)]
pub struct EcBaryCoordinates<L: Prime> {
    /// `u`.
    pub u: Fp2<L>,
    /// `v`.
    pub v: Fp2<L>,
    /// `w`.
    pub w: Fp2<L>,
}

impl<L: FpBackend> Default for EcPoint<L> {
    #[inline]
    fn default() -> Self {
        Self::identity()
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

    /// `(x : z)`.
    #[inline]
    pub fn new(x: Fp2<L>, z: Fp2<L>) -> Self {
        Self { x, z }
    }

    /// `(x : 1)`.
    #[inline]
    pub fn from_x(x: Fp2<L>) -> Self {
        Self { x, z: Fp2::one() }
    }
}

impl<L: FpBackend> Default for EcCurve<L> {
    /// The curve `E0: y^2 = x^3 + x`, i.e. `(A : C) = (0 : 1)`, with `A24`
    /// not yet computed.
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
    /// Bundle three points into a basis.
    #[inline]
    pub fn new(p: EcPoint<L>, q: EcPoint<L>, pmq: EcPoint<L>) -> Self {
        Self { p, q, pmq }
    }
}

impl<L: Prime> Zeroize for EcPoint<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.x.zeroize();
        self.z.zeroize();
    }
}

impl<L: Prime> Zeroize for EcBasis<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.p.zeroize();
        self.q.zeroize();
        self.pmq.zeroize();
    }
}

impl<L: Prime> Zeroize for EcCurve<L> {
    #[inline]
    fn zeroize(&mut self) {
        self.a.zeroize();
        self.c.zeroize();
        self.a24.zeroize();
        self.is_a24_computed_and_normalized = false;
    }
}
