//!
//! Prime: `p = 65 * 2^376 - 1` (383 bits).
//! Field uses 7 limbs of 55-bit radix (unsaturated Montgomery form).

use super::{Level3, SecurityLevel};
use hybrid_array::sizes::{U224, U432, U48, U6, U7, U96, U97};

/// The Level 3 prime `p = 65 * 2^376 - 1` as 48 little-endian bytes.
///
/// In hex: `0x40ff...ff` (383 bits). All bytes are `0xff` except the
/// top byte which is `0x40`.
pub const PRIME_LE_BYTES: [u8; 48] = {
    let mut bytes = [0xffu8; 48];
    bytes[47] = 0x40;
    bytes
};

impl SecurityLevel for Level3 {
    /// 7 limbs × 55-bit radix = 385 bits of storage for the 383-bit prime.
    type FpLimbs = U7;
    /// 6 limbs × 64 bits = 384-bit scalars for order arithmetic.
    type MpLimbs = U6;
    /// `p` fits in 48 bytes (383 bits).
    type FpEncodedBytes = U48;
    /// Two `Fp` elements = 96 bytes.
    type Fp2EncodedBytes = U96;
    /// Public key: 1-byte header + 2 × 48 bytes for the `Fp2` j-invariant.
    type PkLen = U97;
    /// Signature: compressed response isogeny encoding (224 bytes).
    type SigLen = U224;
    /// Secret key: ideal norm + generator coords + basis-change matrix (432 bytes).
    type SkLen = U432;

    fn prime_le_bytes() -> &'static [u8] {
        &PRIME_LE_BYTES
    }

    /// 192-bit post-quantum security.
    const LAMBDA: u32 = 192;

    /// `p + 1 = 65 × 2^376`, so the full `2^376`-torsion is available.
    const F_CHR: u32 = 376;
    /// Response isogeny has degree `2^192`.
    const E_RSP: u32 = 192;
    /// Challenge scalar is 192 bits (matching `LAMBDA`).
    const E_CHL: u32 = 192;
    /// Up to 256 SHAKE256 squeeze attempts to find a valid challenge.
    const HASH_ITERATIONS: u32 = 256;
    /// 6 limbs × 64 = 384-bit scalar width.
    const NWORDS_ORDER: usize = 6;
    /// `v_2(p + 1) = 376`.
    const TORSION_EVEN_POWER: u32 = 376;
    /// `(p + 1) / 2^376 = 65`, which is 7 bits.
    const P_COFACTOR_FOR_2F_BITLENGTH: usize = 7;
    /// Response isogeny length = 192 bits (same as `E_RSP`).
    const SQISIGN_RESPONSE_LENGTH: u32 = 192;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level3_prime_is_correct() {
        let bytes = Level3::prime_le_bytes();
        assert_eq!(bytes.len(), 48);
        for &b in &bytes[..47] {
            assert_eq!(b, 0xFF, "low 47 bytes of p must all be 0xFF");
        }
        assert_eq!(bytes[47], 0x40, "top byte of p must be 0x40");
    }

    #[test]
    fn level3_prime_is_3_mod_4() {
        let bytes = Level3::prime_le_bytes();
        assert_eq!(bytes[0] & 0b11, 3, "p mod 4 must be 3");
    }

    const _: () = assert!(Level3::F_CHR > Level3::LAMBDA);
    const _: () = assert!(Level3::E_RSP > 0);

    #[test]
    fn level3_protocol_exponents_in_range() {
        assert_eq!(Level3::LAMBDA, 192);
        assert_eq!(Level3::F_CHR, 376);
        assert_eq!(Level3::E_RSP, 192);
    }
}
