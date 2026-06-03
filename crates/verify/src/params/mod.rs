//!
//! Defines the [`SecurityLevel`] trait and the marker structs [`Level1`],
//! [`Level3`], [`Level5`] for NIST security levels I, III, V. Downstream
//! crates are generic over `L: SecurityLevel` and the compiler
//! monomorphizes one specialized copy per level from a single source.
//!
//! All three levels ([`Level1`], [`Level3`], [`Level5`]) have full
//! `SecurityLevel` implementations.

use hybrid_array::ArraySize;

pub mod level1;
pub mod level3;
pub mod level5;

/// Marker trait for SQIsign security levels (NIST I, III, V).
///
/// Associated types encode array/buffer sizes via `hybrid-array` typenum.
/// Associated constants encode scalar parameters used in arithmetic and
/// loop bounds.
///
/// # Design rules
///
/// - Only INDEPENDENT, per-level parameters belong here.
/// - Derived sizes are computed in downstream crate impl blocks.
/// - Do not add trait bounds beyond `Default + Clone + Debug + 'static`.
pub trait SecurityLevel: Default + Clone + core::fmt::Debug + 'static {
    /// Number of 64-bit limbs in a prime-field element `Fp`.
    ///
    /// For Level 1 this is `U5`: 5 unsaturated radix-2^51 limbs
    /// (51 * 5 = 255 bits of storage for a 251-bit prime).
    type FpLimbs: ArraySize;

    /// Number of 64-bit limbs for `mp`-layer scalar/order intermediates.
    /// This is `NWORDS_ORDER`, not the field width.
    type MpLimbs: ArraySize;

    /// Byte length of a serialized canonical `Fp` element.
    type FpEncodedBytes: ArraySize;

    /// Byte length of a serialized canonical `Fp2` element (`2 *
    /// FpEncodedBytes`).
    type Fp2EncodedBytes: ArraySize;

    /// Byte length of a serialized public key.
    type PkLen: ArraySize;

    /// Byte length of a serialized standard signature.
    type SigLen: ArraySize;

    /// Byte length of a serialized expanded signature.
    type ExpandedSigLen: ArraySize;

    /// Byte length of a serialized compressed signature.
    type CompressedSigLen: ArraySize;

    /// Byte length of a serialized secret key.
    type SkLen: ArraySize;

    /// The prime `p` as a static byte slice (little-endian canonical
    /// encoding, length `FP_ENCODED_BYTES`).
    fn prime_le_bytes() -> &'static [u8];

    /// Security parameter `lambda` in bits (128, 192, or 256).
    const LAMBDA: u32;

    /// Exponent `f` of the 2-power torsion available on the starting
    /// curve E0. The full 2ᶠ-torsion E0\[2ᶠ\] ≅ (ℤ/2ᶠ)² is
    /// rational over `Fp2`. Equal to [`TORSION_EVEN_POWER`](Self::TORSION_EVEN_POWER).
    const F_CHR: u32;

    /// Bit-length of the response isogeny degree `2^E_RSP`. The
    /// response isogeny in the sigma protocol has degree exactly
    /// `2^E_RSP`. Must satisfy `E_RSP < F_CHR`.
    const E_RSP: u32;

    /// Bit-length of the challenge space. The verifier's challenge
    /// scalar is drawn from `{0, ..., 2^E_CHL - 1}`. Equal to
    /// `LAMBDA` for standard security.
    const E_CHL: u32;

    /// Maximum number of hash-to-challenge iterations before aborting.
    /// Each iteration squeezes from SHAKE256 and checks whether the
    /// result falls in the valid challenge range.
    const HASH_ITERATIONS: u32;

    /// Number of 64-bit limbs for order/scalar arrays. Determines the
    /// width of multi-precision scalar arithmetic in the EC layer.
    /// Level 1 = 4 (256 bits), Level 3 = 6 (384 bits), Level 5 = 8
    /// (512 bits).
    const NWORDS_ORDER: usize;

    /// The 2-adic valuation of `p + 1`, i.e. the largest `f` such that
    /// 2ᶠ | (p + 1). This is the exponent of the available even
    /// torsion on the supersingular curve E0 over `Fp2`.
    const TORSION_EVEN_POWER: u32;

    /// Bit-length of the odd cofactor `(p + 1) / 2^TORSION_EVEN_POWER`.
    /// For Level 1: `(p+1)/2^248 = 5`, which has bit-length 3.
    const P_COFACTOR_FOR_2F_BITLENGTH: usize;

    /// Length (in bits) of the response isogeny in the sigma protocol.
    /// Equal to [`E_RSP`](Self::E_RSP).
    const SQISIGN_RESPONSE_LENGTH: u32;
}

/// NIST Security Level I (128-bit post-quantum security).
///
/// Prime: `p = 5 * 2^248 - 1`, encoded in 32 bytes.
#[derive(Default, Clone, Debug)]
pub struct Level1;

/// NIST Security Level III (192-bit post-quantum security).
///
/// Prime: `p = 65 * 2^376 - 1`, encoded in 48 bytes.
#[derive(Default, Clone, Debug)]
pub struct Level3;

/// NIST Security Level V (256-bit post-quantum security).
///
/// Prime: `p = 27 * 2^500 - 1`, encoded in 64 bytes.
#[derive(Default, Clone, Debug)]
pub struct Level5;
