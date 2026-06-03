//!
//! Contains the base curve E0, torsion point bases, and other per-level
//! constant data needed by the EC and isogeny layers.
//!
//! Enable the `signing` feature to include quaternion data, endomorphism
//! action matrices, and torsion degree constants needed by the signing path.

use crate::params::SecurityLevel;

pub mod level1;
pub mod level3;
pub mod level5;

/// Level-specific precomputed constants needed by the verification and
/// EC layers. Implemented for each security level marker type.
pub trait LevelPrecomp: SecurityLevel {
    /// Canonical 𝔽p²-encoded bytes for the x-coordinate of the first
    /// generator of the 2ᶠ-torsion basis on E0.
    fn basis_e0_px_bytes() -> &'static [u8];

    /// Canonical 𝔽p²-encoded bytes for the x-coordinate of the second
    /// generator of the 2ᶠ-torsion basis on E0.
    fn basis_e0_qx_bytes() -> &'static [u8];

    /// The odd cofactor (p+1) / 2ᶠ as 64-bit limbs (little-endian).
    fn p_cofactor_for_2f() -> &'static [u64];

    /// Bit-length of the odd cofactor.
    fn p_cofactor_for_2f_bitlength() -> u32;

    /// Exponent f such that the torsion subgroup is ℤ/2ᶠ × ℤ/2ᶠ.
    fn torsion_even_power() -> u32;

    /// 10 precomputed 4×4 basis change matrices for splitting transforms.
    fn splitting_transforms() -> &'static [[[u8; 4]; 4]; 10];

    /// 6 precomputed 4×4 normalization matrices for splitting.
    fn normalization_transforms() -> &'static [[[u8; 4]; 4]; 6];

    /// Character evaluation table for splitting.
    fn chi_eval() -> &'static [[i32; 4]; 4];

    /// Pairs of indices for the 10 possible zero-positions in splitting.
    fn even_index() -> &'static [[i32; 2]; 10];
}

impl LevelPrecomp for crate::params::Level1 {
    fn basis_e0_px_bytes() -> &'static [u8] {
        &level1::BASIS_E0_PX_BYTES
    }
    fn basis_e0_qx_bytes() -> &'static [u8] {
        &level1::BASIS_E0_QX_BYTES
    }
    fn p_cofactor_for_2f() -> &'static [u64] {
        level1::P_COFACTOR_FOR_2F
    }
    fn p_cofactor_for_2f_bitlength() -> u32 {
        level1::P_COFACTOR_FOR_2F_BITLENGTH
    }
    fn torsion_even_power() -> u32 {
        level1::TORSION_EVEN_POWER
    }
    fn splitting_transforms() -> &'static [[[u8; 4]; 4]; 10] {
        &level1::SPLITTING_TRANSFORMS
    }
    fn normalization_transforms() -> &'static [[[u8; 4]; 4]; 6] {
        &level1::NORMALIZATION_TRANSFORMS
    }
    fn chi_eval() -> &'static [[i32; 4]; 4] {
        &level1::CHI_EVAL
    }
    fn even_index() -> &'static [[i32; 2]; 10] {
        &level1::EVEN_INDEX
    }
}

impl LevelPrecomp for crate::params::Level3 {
    fn basis_e0_px_bytes() -> &'static [u8] {
        &level3::BASIS_E0_PX_BYTES
    }
    fn basis_e0_qx_bytes() -> &'static [u8] {
        &level3::BASIS_E0_QX_BYTES
    }
    fn p_cofactor_for_2f() -> &'static [u64] {
        level3::P_COFACTOR_FOR_2F
    }
    fn p_cofactor_for_2f_bitlength() -> u32 {
        level3::P_COFACTOR_FOR_2F_BITLENGTH
    }
    fn torsion_even_power() -> u32 {
        level3::TORSION_EVEN_POWER
    }
    fn splitting_transforms() -> &'static [[[u8; 4]; 4]; 10] {
        &level3::SPLITTING_TRANSFORMS
    }
    fn normalization_transforms() -> &'static [[[u8; 4]; 4]; 6] {
        &level3::NORMALIZATION_TRANSFORMS
    }
    fn chi_eval() -> &'static [[i32; 4]; 4] {
        &level3::CHI_EVAL
    }
    fn even_index() -> &'static [[i32; 2]; 10] {
        &level3::EVEN_INDEX
    }
}

impl LevelPrecomp for crate::params::Level5 {
    fn basis_e0_px_bytes() -> &'static [u8] {
        &level5::BASIS_E0_PX_BYTES
    }
    fn basis_e0_qx_bytes() -> &'static [u8] {
        &level5::BASIS_E0_QX_BYTES
    }
    fn p_cofactor_for_2f() -> &'static [u64] {
        level5::P_COFACTOR_FOR_2F
    }
    fn p_cofactor_for_2f_bitlength() -> u32 {
        level5::P_COFACTOR_FOR_2F_BITLENGTH
    }
    fn torsion_even_power() -> u32 {
        level5::TORSION_EVEN_POWER
    }
    fn splitting_transforms() -> &'static [[[u8; 4]; 4]; 10] {
        &level5::SPLITTING_TRANSFORMS
    }
    fn normalization_transforms() -> &'static [[[u8; 4]; 4]; 6] {
        &level5::NORMALIZATION_TRANSFORMS
    }
    fn chi_eval() -> &'static [[i32; 4]; 4] {
        &level5::CHI_EVAL
    }
    fn even_index() -> &'static [[i32; 2]; 10] {
        &level5::EVEN_INDEX
    }
}
