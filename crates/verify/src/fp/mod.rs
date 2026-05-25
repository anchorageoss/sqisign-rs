//! SQIsign.
//!
//! All arithmetic uses Montgomery form. SQIsign primes satisfy
//! `p = 3 mod 4`, so the extension is `Fp2 = Fp[i] / (i^2 + 1)`.
//!
//! Per-level dispatch is via the [`FpBackend`] trait: each security level
//! implements the primitive arithmetic on the raw limb arrays, and
//! [`Fp`] / [`Fp2`] methods forward to those implementations. Generic
//! code can write `where L: FpBackend` to operate on any level.
//!
//! The Level 1 backend in [`level1`] uses a 5-limb unsaturated radix-`2^51`
//! representation for the 251-bit prime `p = 5 * 2^248 - 1`. Limb layouts
//! and internal helpers are an implementation detail of each backend; the
//! [`Fp`] / [`Fp2`] surface is the only stable API.

use crate::params::SecurityLevel;
use hybrid_array::Array;
use subtle::Choice;
use zeroize::Zeroize;

#[allow(clippy::module_inception)]
pub mod fp;
pub mod fp2;
pub mod level1;
pub mod level3;
pub mod level5;

/// A prime-field element.
///
/// The internal representation is opaque per-level: each [`FpBackend`]
/// chooses a limb layout suited to its prime. User code interacts with
/// `Fp` only through its methods; the limbs are not part of the public
/// API.
#[derive(Clone, Debug)]
pub struct Fp<L: SecurityLevel> {
    pub(crate) limbs: Array<u64, L::FpLimbs>,
}

/// A quadratic-extension field element `a = re + im * i`, where `i^2 = -1`.
///
/// `Fp2` is laid out as two independent `Fp` values. Karatsuba
/// multiplication uses three `Fp` multiplications and five `Fp`
/// adds/subs.
#[derive(Clone, Debug)]
pub struct Fp2<L: SecurityLevel> {
    pub re: Fp<L>,
    pub im: Fp<L>,
}

/// Per-level field-arithmetic backend.
///
/// Each `SecurityLevel` implements this trait with the prime-specific
/// limb routines. The `Fp` and `Fp2` user-facing methods dispatch through
/// this trait, so generic code can write `where L: FpBackend` to operate
/// on any level.
///
/// All methods operate on the raw limb storage. Inputs and outputs are
/// in the backend's internal Montgomery form unless noted otherwise.
/// Arithmetic operations (`add`, `sub`, `neg`, `mul`, `sqr`, `mul_small`)
/// leave results lazily reduced in `[0, 2p)`; the final reduction to
/// `[0, p)` happens during encoding and comparison.
pub trait FpBackend: SecurityLevel {
    /// `out <- 0`.
    fn set_zero(out: &mut Array<u64, Self::FpLimbs>);

    /// `out <- 1` in Montgomery form.
    fn set_one(out: &mut Array<u64, Self::FpLimbs>);

    /// `out <- val` in Montgomery form, treating `val` as an unsigned
    /// integer that fits in the field.
    fn set_small(out: &mut Array<u64, Self::FpLimbs>, val: u64);

    /// Constant-time equality test. Returns `Choice(1)` if `a == b`
    /// (after full reduction), `Choice(0)` otherwise.
    fn is_equal(a: &Array<u64, Self::FpLimbs>, b: &Array<u64, Self::FpLimbs>) -> Choice;

    /// Constant-time zero test. Returns `Choice(1)` if `a == 0`
    /// (after full reduction), `Choice(0)` otherwise.
    fn is_zero(a: &Array<u64, Self::FpLimbs>) -> Choice;

    /// `out <- a` (copy).
    fn copy(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// `out <- a + b mod 2p`.
    fn add(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    );

    /// `out <- a - b mod 2p`.
    fn sub(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    );

    /// `out <- -a mod 2p`.
    fn neg(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// Montgomery multiplication: `out <- a * b * R^{-1} mod 2p`.
    fn mul(
        out: &mut Array<u64, Self::FpLimbs>,
        a: &Array<u64, Self::FpLimbs>,
        b: &Array<u64, Self::FpLimbs>,
    );

    /// Specialized Montgomery squaring: `out <- a^2 * R^{-1} mod 2p`.
    fn sqr(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// `out <- 1 / a mod p`. If `a == 0` the output is `0` (no panic).
    fn inv(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// `out <- sqrt(a) mod p`. The caller is responsible for ensuring
    /// `a` is a quadratic residue; on non-QR inputs the output is
    /// well-defined but is not a square root of `a`. The result is
    /// determined only up to sign.
    fn sqrt(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// Returns `Choice(1)` if `a` is a quadratic residue (or zero) in
    /// Fp, `Choice(0)` otherwise.
    fn is_square(a: &Array<u64, Self::FpLimbs>) -> Choice;

    /// `out <- a / 2 mod p`.
    fn half(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// `out <- a / 3 mod p`.
    fn div3(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// Square root progenitor: `out <- a^((p-3)/4) mod p`. Combined
    /// with one extra multiplication this yields `sqrt(a)` when `p = 3
    /// mod 4`.
    fn exp3div4(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>);

    /// `out <- a * val mod 2p` for a small (32-bit) integer multiplier.
    fn mul_small(out: &mut Array<u64, Self::FpLimbs>, a: &Array<u64, Self::FpLimbs>, val: u32);

    /// Serialize `a` to its canonical little-endian byte form. Writes
    /// exactly `Self::FpEncodedBytes::USIZE` bytes.
    fn encode(out: &mut [u8], a: &Array<u64, Self::FpLimbs>);

    /// Deserialize an `Fp` element from `Self::FpEncodedBytes::USIZE`
    /// canonical little-endian bytes. Returns `Choice(1)` if the input
    /// represented an integer in `[0, p)`, `Choice(0)` otherwise. On
    /// out-of-range input the output is zeroed.
    fn decode(out: &mut Array<u64, Self::FpLimbs>, bytes: &[u8]) -> Choice;

    /// Decode a possibly-longer little-endian byte string with full
    /// modular reduction. Used to map a hash output uniformly into Fp.
    fn decode_reduce(out: &mut Array<u64, Self::FpLimbs>, bytes: &[u8]);

    /// Constant-time conditional swap: if `ctl` is set, swap `a`
    /// and `b`; otherwise leave them unchanged.
    fn cswap(a: &mut Array<u64, Self::FpLimbs>, b: &mut Array<u64, Self::FpLimbs>, ctl: Choice);

    /// Constant-time conditional select: if `ctl` is clear, set
    /// `out <- a0`; if `ctl` is set, set `out <- a1`.
    fn select(
        out: &mut Array<u64, Self::FpLimbs>,
        a0: &Array<u64, Self::FpLimbs>,
        a1: &Array<u64, Self::FpLimbs>,
        ctl: Choice,
    );
}

impl<L: SecurityLevel> Zeroize for Fp<L> {
    fn zeroize(&mut self) {
        self.limbs.as_mut_slice().zeroize();
    }
}

impl<L: SecurityLevel> Zeroize for Fp2<L> {
    fn zeroize(&mut self) {
        self.re.zeroize();
        self.im.zeroize();
    }
}
