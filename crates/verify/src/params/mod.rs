//! Parameter sets.
//!
//! Every prime used by SQIsign and PRISM has the form `p = c * 2^e - 1`.
//! The [`Prime`] trait carries the facts about `p` that the field layer
//! needs; protocol-level parameters (challenge lengths, encodings, ...)
//! will be layered on top in later crates and tickets, so that a prime
//! shared between two parameter sets (round-3 level III and round-2
//! level V both use `27 * 2^500 - 1`) has exactly one field backend.
//!
//! Marker structs are named after the SQIsign reference's build variants
//! (`p324_3` is `3 * 2^324 - 1`) and the round-3 NIST levels are aliases.

use hybrid_array::ArraySize;

/// Facts about a prime `p = c * 2^e - 1` needed by the field layer.
///
/// Only independent, per-prime data belongs here. Anything derivable is
/// computed where it is used.
pub trait Prime: Default + Clone + core::fmt::Debug + 'static {
    /// Number of 64-bit limbs in an `Fp` element (the backend's layout:
    /// unsaturated in the portable backends, saturated in the assembly ones).
    type FpLimbs: ArraySize;
    /// Byte length of a canonical `Fp` encoding, `ceil(log2(p) / 8)`.
    type FpEncodedBytes: ArraySize;
    /// Byte length of a canonical `Fp2` encoding, `2 * FpEncodedBytes`.
    type Fp2EncodedBytes: ArraySize;

    /// The odd cofactor `c`.
    const COFACTOR: u64;
    /// The two-adic exponent `e`: the full `2^e`-torsion of a supersingular
    /// curve over `Fp2` is rational.
    const TWO_ADIC_EXPONENT: u32;
    /// Bit length of the cofactor `c`.
    const COFACTOR_BITLENGTH: usize;
    /// 64-bit words needed for an integer of `log2 p` bits (the reference's
    /// `NWORDS_ORDER`): the width of scalars in the curve layer.
    const ORDER_WORDS: usize;

    /// `p` as canonical little-endian bytes of length `FpEncodedBytes`.
    fn prime_le_bytes() -> &'static [u8];
}

macro_rules! prime_set {
    ($(#[$doc:meta])* $name:ident, $module:ident, limbs = $limbs:ty, bytes = $bytes:ty, bytes2 = $bytes2:ty, c = $c:expr, e = $e:expr, words = $words:expr) => {
        $(#[$doc])*
        #[allow(non_camel_case_types)] // named after the reference's build variant
        #[derive(Default, Clone, Debug)]
        pub struct $name;

        impl Prime for $name {
            type FpLimbs = $limbs;
            type FpEncodedBytes = $bytes;
            type Fp2EncodedBytes = $bytes2;
            const COFACTOR: u64 = $c;
            const TWO_ADIC_EXPONENT: u32 = $e;
            const COFACTOR_BITLENGTH: usize = (64 - ($c as u64).leading_zeros()) as usize;
            const ORDER_WORDS: usize = $words;
            fn prime_le_bytes() -> &'static [u8] {
                &crate::fp::$module::PRIME_LE_BYTES
            }
        }
    };
}

prime_set!(
    /// `p = 3 * 2^324 - 1` (326 bits). SQIsign round 3, NIST level I.
    P324_3, p324_3,
    limbs = hybrid_array::sizes::U6, bytes = hybrid_array::sizes::U41, bytes2 = hybrid_array::sizes::U82,
    c = 3, e = 324, words = 6
);
/// Limbs of a `P500_27` element: nine unsaturated 57-bit limbs in the
/// generated radix backend (other architectures), eight saturated ones in
/// the x86-64 backend.
#[cfg(not(target_arch = "x86_64"))]
pub type P500_27Limbs = hybrid_array::sizes::U9;
/// Limbs of a `P500_27` element on x86-64.
#[cfg(target_arch = "x86_64")]
pub type P500_27Limbs = hybrid_array::sizes::U8;

prime_set!(
    /// `p = 27 * 2^500 - 1` (505 bits). SQIsign round 3, NIST level III
    /// (and round 2, NIST level V).
    P500_27, p500_27,
    limbs = P500_27Limbs, bytes = hybrid_array::sizes::U64, bytes2 = hybrid_array::sizes::U128,
    c = 27, e = 500, words = 8
);
prime_set!(
    /// `p = 17 * 2^664 - 1` (669 bits). SQIsign round 3, NIST level V.
    P664_17, p664_17,
    limbs = hybrid_array::sizes::U11, bytes = hybrid_array::sizes::U84, bytes2 = hybrid_array::sizes::U168,
    c = 17, e = 664, words = 11
);

/// NIST level I (round-3 parameters).
pub type Level1 = P324_3;
/// NIST level III (round-3 parameters).
pub type Level3 = P500_27;
/// NIST level V (round-3 parameters).
pub type Level5 = P664_17;

pub mod sqisign_v3;
