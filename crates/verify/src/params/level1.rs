//!
//! Prime: `p = 5 * 2^248 - 1` (251 bits).
//! Field uses 5 limbs of 51-bit radix (unsaturated Montgomery form).

use super::{Level1, SecurityLevel};
use hybrid_array::sizes::{U148, U288, U32, U4, U5, U64, U65};

/// The Level 1 prime `p = 5 * 2^248 - 1` as 32 little-endian bytes.
///
/// In hex this is `0x04ff..ff` (top byte `0x04`, then 31 bytes of `0xff`),
/// matching `p` in `sqisign_parameters.txt`.
pub const PRIME_LE_BYTES: [u8; 32] = {
    let mut bytes = [0xffu8; 32];
    bytes[31] = 0x04;
    bytes
};

impl SecurityLevel for Level1 {
    /// 5 limbs × 51-bit radix = 255 bits of storage for the 251-bit prime.
    type FpLimbs = U5;
    /// 4 limbs × 64 bits = 256-bit scalars for order arithmetic.
    type MpLimbs = U4;
    /// `p` fits in 32 bytes (251 bits).
    type FpEncodedBytes = U32;
    /// Two `Fp` elements = 64 bytes.
    type Fp2EncodedBytes = U64;
    /// Public key: 1-byte header + 2 × 32 bytes for the `Fp2` j-invariant.
    type PkLen = U65;
    /// Signature: compressed response isogeny encoding (148 bytes).
    type SigLen = U148;
    /// Secret key: ideal norm + generator coords + basis-change matrix (288 bytes).
    type SkLen = U288;

    fn prime_le_bytes() -> &'static [u8] {
        &PRIME_LE_BYTES
    }

    /// 128-bit post-quantum security.
    const LAMBDA: u32 = 128;

    /// `p + 1 = 5 × 2^248`, so the full `2^248`-torsion is available.
    const F_CHR: u32 = 248;
    /// Response isogeny has degree `2^126`.
    const E_RSP: u32 = 126;
    /// Challenge scalar is 128 bits (matching `LAMBDA`).
    const E_CHL: u32 = 128;
    /// Up to 64 SHAKE256 squeeze attempts to find a valid challenge.
    const HASH_ITERATIONS: u32 = 64;
    /// 4 limbs × 64 = 256-bit scalar width.
    const NWORDS_ORDER: usize = 4;
    /// `v_2(p + 1) = 248`.
    const TORSION_EVEN_POWER: u32 = 248;
    /// `(p + 1) / 2^248 = 5`, which is 3 bits.
    const P_COFACTOR_FOR_2F_BITLENGTH: usize = 3;
    /// Response isogeny length = 126 bits (same as `E_RSP`).
    const SQISIGN_RESPONSE_LENGTH: u32 = 126;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reconstruct `p` from `PRIME_LE_BYTES` and verify it equals
    /// `5 * 2^248 - 1`.
    ///
    /// We do the comparison limb-by-limb: little-endian, bytes 0..31 all
    /// `0xFF` except byte 31 which is `0x04`. The integer is therefore
    /// `4 * 2^248 + (2^248 - 1) = 5 * 2^248 - 1`.
    #[test]
    fn level1_prime_is_correct() {
        let bytes = Level1::prime_le_bytes();
        assert_eq!(bytes.len(), 32);
        for &b in &bytes[..31] {
            assert_eq!(b, 0xFF, "low 31 bytes of p must all be 0xFF");
        }
        assert_eq!(bytes[31], 0x04, "top byte of p must be 0x04");
    }

    /// `p = 5 * 2^248 - 1`, so `p mod 4 = (5*2^248 mod 4) - 1 mod 4 =
    /// 0 - 1 mod 4 = 3`. Required for the `Fp2 = Fp[i]/(i^2 + 1)`
    /// construction and the Fermat-style `Fp` square root.
    #[test]
    fn level1_prime_is_3_mod_4() {
        let bytes = Level1::prime_le_bytes();
        // Bottom byte determines the value mod 4 (since 256 = 0 mod 4).
        assert_eq!(bytes[0] & 0b11, 3, "p mod 4 must be 3");
    }

    /// Sanity: protocol exponents are in expected ranges. The two
    /// `const _: () = assert!(...)` blocks fire at compile time if the
    /// invariants are violated by future edits; the runtime asserts
    /// just make the constants observable in the test output.
    const _: () = assert!(Level1::F_CHR > Level1::LAMBDA);
    const _: () = assert!(Level1::E_RSP > 0);

    #[test]
    fn level1_protocol_exponents_in_range() {
        assert_eq!(Level1::LAMBDA, 128);
        assert_eq!(Level1::F_CHR, 248);
        assert_eq!(Level1::E_RSP, 126);
    }
}
