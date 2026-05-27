//!
//! Prime: `p = 27 * 2^500 - 1` (505 bits).
//! Field uses 9 limbs of 57-bit radix (unsaturated Montgomery form).

use super::{Level5, SecurityLevel};
use hybrid_array::sizes::{U128, U129, U257, U292, U420, U576, U64, U8, U9};

/// The Level 5 prime `p = 27 * 2^500 - 1` as 64 little-endian bytes.
///
/// In hex: `0x01af...ff` (505 bits). The bottom 62 bytes are `0xff`,
/// byte 62 is `0xAF`, and byte 63 is `0x01`.
///
/// Derivation: `p + 1 = 27 * 2^500 = 0x1B * 2^500`. Bit 500 lies at
/// bit-position 4 within byte 62 (since 500 = 62 * 8 + 4). Subtracting
/// 1 sets all 500 low bits to 1 and decrements the upper portion from
/// `0x1B` to `0x1A`, which straddles bytes 62 and 63 as `0xAF` and
/// `0x01` respectively.
pub const PRIME_LE_BYTES: [u8; 64] = {
    let mut bytes = [0xffu8; 64];
    bytes[62] = 0xAF;
    bytes[63] = 0x01;
    bytes
};

impl SecurityLevel for Level5 {
    /// 9 limbs × 57-bit radix = 513 bits of storage for the 505-bit prime.
    type FpLimbs = U9;
    /// 8 limbs × 64 bits = 512-bit scalars for order arithmetic.
    type MpLimbs = U8;
    /// `p` fits in 64 bytes (505 bits).
    type FpEncodedBytes = U64;
    /// Two `Fp` elements = 128 bytes.
    type Fp2EncodedBytes = U128;
    /// Public key: 1-byte header + 2 × 64 bytes for the `Fp2` j-invariant.
    type PkLen = U129;
    /// Signature: compressed response isogeny encoding (292 bytes).
    type SigLen = U292;
    /// Expanded signature (420 bytes).
    type ExpandedSigLen = U420;
    /// Compressed signature (257 bytes).
    type CompressedSigLen = U257;
    /// Secret key: ideal norm + generator coords + basis-change matrix.
    /// Actual content is 572 bytes; U576 is the next upstream hybrid-array
    /// size. The 4 trailing bytes are zero-padded.
    type SkLen = U576;

    fn prime_le_bytes() -> &'static [u8] {
        &PRIME_LE_BYTES
    }

    /// 256-bit post-quantum security.
    const LAMBDA: u32 = 256;

    /// `p + 1 = 27 * 2^500`, so the full `2^500`-torsion is available.
    const F_CHR: u32 = 500;
    /// Response isogeny has degree `2^253`.
    const E_RSP: u32 = 253;
    /// Challenge scalar is 256 bits (matching `LAMBDA`).
    const E_CHL: u32 = 256;
    /// Up to 512 SHAKE256 squeeze attempts to find a valid challenge.
    const HASH_ITERATIONS: u32 = 512;
    /// 8 limbs × 64 = 512-bit scalar width.
    const NWORDS_ORDER: usize = 8;
    /// `v_2(p + 1) = 500`.
    const TORSION_EVEN_POWER: u32 = 500;
    /// `(p + 1) / 2^500 = 27 = 0b11011`, which is 5 bits.
    const P_COFACTOR_FOR_2F_BITLENGTH: usize = 5;
    /// Response isogeny length = 253 bits (same as `E_RSP`).
    const SQISIGN_RESPONSE_LENGTH: u32 = 253;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level5_prime_is_correct() {
        let bytes = Level5::prime_le_bytes();
        assert_eq!(bytes.len(), 64);
        for &b in &bytes[..62] {
            assert_eq!(b, 0xFF, "low 62 bytes of p must all be 0xFF");
        }
        assert_eq!(bytes[62], 0xAF, "byte 62 of p must be 0xAF");
        assert_eq!(bytes[63], 0x01, "top byte of p must be 0x01");
    }

    #[test]
    fn level5_prime_is_3_mod_4() {
        let bytes = Level5::prime_le_bytes();
        assert_eq!(bytes[0] & 0b11, 3, "p mod 4 must be 3");
    }

    /// Verify the bit-length of p. The most significant byte is 0x01
    /// (bit 504 set), giving 505 bits total. This matches `27 * 2^500 - 1`.
    #[test]
    fn level5_prime_bitlength() {
        let bytes = Level5::prime_le_bytes();
        // Byte 63 = 0x01 means bit 504 is set, bits 505-511 are zero.
        assert_eq!(bytes[63], 0x01);
        // Byte 62 = 0xAF = 0b1010_1111, so bit 503 (bit 7 of byte 62) is set.
        // The topmost set bit is bit 504 (in byte 63), giving 505-bit prime.
        assert_eq!(bytes[62] & 0x80, 0x80, "bit 503 must be set");
    }

    const _: () = assert!(Level5::F_CHR > Level5::LAMBDA);
    const _: () = assert!(Level5::E_RSP > 0);

    #[test]
    fn level5_protocol_exponents_in_range() {
        assert_eq!(Level5::LAMBDA, 256);
        assert_eq!(Level5::F_CHR, 500);
        assert_eq!(Level5::E_RSP, 253);
    }
}
