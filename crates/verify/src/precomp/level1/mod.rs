//!
//! Contains EC layer constants (basis points, splitting transforms) and
//! signing-layer constants (quaternion data, endomorphism actions, torsion degrees).

pub mod e0_basis;
pub mod ec_params;

/// Exponent f such that the torsion subgroup is ℤ/2ᶠ × ℤ/2ᶠ.
pub const TORSION_EVEN_POWER: u32 = 248;

/// Bit-length of the odd cofactor (p+1) / 2ᶠ.
pub const P_COFACTOR_FOR_2F_BITLENGTH: u32 = 3;

/// The odd cofactor (p+1) / 2ᶠ as a single 64-bit limb.
/// For Level 1: `(p+1)/2^248 = 5`.
pub const P_COFACTOR_FOR_2F: &[u64] = &[5];

/// Canonical 𝔽p²-encoded bytes for `BASIS_E0_PX`, the x-coordinate of the
/// first generator of the 2ᶠ-torsion basis on E0.
pub const BASIS_E0_PX_BYTES: [u8; 64] = [
    0x78, 0x00, 0xb4, 0xae, 0x5e, 0xd9, 0x19, 0x21, 0x8b, 0xa7, 0xbf, 0x59, 0x1a, 0x99, 0xbe, 0x44,
    0xc4, 0x16, 0x62, 0xa6, 0xc3, 0x04, 0xcc, 0x83, 0x24, 0xb1, 0x82, 0xca, 0x7f, 0x87, 0x9b, 0x01,
    0x75, 0xd2, 0xf9, 0xc3, 0x3d, 0x13, 0x04, 0x8e, 0x74, 0x92, 0x42, 0x51, 0xae, 0xdd, 0xcf, 0xb2,
    0x2f, 0xe9, 0x67, 0x98, 0xaa, 0x0a, 0x15, 0x52, 0x42, 0xe0, 0xea, 0x49, 0xdb, 0x2a, 0x44, 0x04,
];

/// Canonical 𝔽p²-encoded bytes for `BASIS_E0_QX`, the x-coordinate of the
/// second generator of the 2ᶠ-torsion basis on E0.
pub const BASIS_E0_QX_BYTES: [u8; 64] = [
    0x1f, 0xeb, 0x93, 0x55, 0x2a, 0x25, 0x16, 0x7c, 0xf3, 0xe1, 0x4b, 0xa5, 0xf7, 0x78, 0x86, 0x87,
    0x1d, 0x04, 0x0d, 0x05, 0x17, 0x27, 0xdf, 0x9f, 0x71, 0x0b, 0x5c, 0x7d, 0x47, 0xfd, 0x5f, 0x04,
    0xee, 0xaa, 0xcd, 0xa0, 0xb7, 0x28, 0xe2, 0xfd, 0xae, 0xa5, 0x8e, 0x6f, 0x4b, 0x05, 0xff, 0x39,
    0x9a, 0xb3, 0x76, 0x36, 0xfb, 0xa8, 0x65, 0x44, 0xdc, 0x73, 0x18, 0xdf, 0xe9, 0xd4, 0x87, 0x04,
];

/// Indices encoding the set `{0, 1, i, -1, -i}`.
/// Used by `SPLITTING_TRANSFORMS` and `NORMALIZATION_TRANSFORMS` to
/// encode basis-change matrices over `Fp2` compactly as u8 indices.
pub const FP2_ZERO: u8 = 0;
pub const FP2_ONE: u8 = 1;
pub const FP2_I: u8 = 2;
pub const FP2_MINUS_ONE: u8 = 3;
pub const FP2_MINUS_I: u8 = 4;

/// Pairs of indices for the 10 possible zero-positions in splitting.
pub const EVEN_INDEX: [[i32; 2]; 10] = [
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3],
];

/// Character evaluation table for splitting.
pub const CHI_EVAL: [[i32; 4]; 4] = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]];

/// 10 precomputed 4×4 basis change matrices for splitting transforms.
/// Each entry is a u8 index into `{0, 1, i, -1, -i}`.
pub const SPLITTING_TRANSFORMS: [[[u8; 4]; 4]; 10] = [
    [
        [FP2_ONE, FP2_I, FP2_ONE, FP2_I],
        [FP2_ONE, FP2_MINUS_I, FP2_MINUS_ONE, FP2_I],
        [FP2_ONE, FP2_I, FP2_MINUS_ONE, FP2_MINUS_I],
        [FP2_MINUS_ONE, FP2_I, FP2_MINUS_ONE, FP2_I],
    ],
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        [FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO, FP2_ZERO],
    ],
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        [FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO],
    ],
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE],
    ],
    [
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
        [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
        [FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE],
    ],
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
    ],
    [
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
        [FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE],
    ],
    [
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
        [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
    ],
    [
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
        [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
        [FP2_MINUS_ONE, FP2_ONE, FP2_ONE, FP2_MINUS_ONE],
    ],
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
    ],
];

/// 6 precomputed 4×4 normalization matrices for splitting.
pub const NORMALIZATION_TRANSFORMS: [[[u8; 4]; 4]; 6] = [
    [
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
    ],
    [
        [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
        [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
    ],
    [
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
        [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
        [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
    ],
    [
        [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
        [FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE],
        [FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE],
        [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
    ],
    [
        [FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE],
        [FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I],
        [FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I],
        [FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE],
    ],
    [
        [FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE],
        [FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I],
        [FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I],
        [FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE],
    ],
];
