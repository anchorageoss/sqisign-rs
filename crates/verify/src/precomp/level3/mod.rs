//!
//! Contains EC layer constants (basis points, splitting transforms) and
//! signing-layer constants (quaternion data, endomorphism actions, torsion degrees).

pub mod e0_basis;
pub mod ec_params;

/// Exponent f such that the torsion subgroup is `Z/2^f x Z/2^f`.
pub const TORSION_EVEN_POWER: u32 = 376;

/// Bit-length of the odd cofactor `(p+1) / 2^f`.
pub const P_COFACTOR_FOR_2F_BITLENGTH: u32 = 7;

/// The odd cofactor `(p+1) / 2^f` as a single 64-bit limb.
/// For Level 3: `(p+1)/2^376 = 65 = 0x41`.
pub const P_COFACTOR_FOR_2F: &[u64] = &[0x41];

/// Canonical Fp2-encoded bytes for `BASIS_E0_PX`, the x-coordinate of the
/// first generator of the `2^f`-torsion basis on E0.
pub const BASIS_E0_PX_BYTES: [u8; 96] = [
    0x17, 0x43, 0xd2, 0xd1, 0x69, 0xed, 0xa2, 0x3e, 0xc4, 0xdb, 0x76, 0x0c, 0xc2, 0x11, 0xd3, 0xdb,
    0xc5, 0xc2, 0x7b, 0xd7, 0x87, 0x3a, 0xb7, 0x3f, 0x85, 0x2c, 0x5d, 0x5b, 0x3a, 0xa7, 0x59, 0xd0,
    0x56, 0xc4, 0x6e, 0xd2, 0x71, 0xf7, 0xe2, 0x48, 0xff, 0xdb, 0xb6, 0x7f, 0xc2, 0xa1, 0x98, 0x17,
    0xf0, 0x11, 0xdd, 0x9d, 0x1b, 0x46, 0xfc, 0xaa, 0x1f, 0x52, 0x23, 0x56, 0x1c, 0x20, 0xb5, 0xf4,
    0xdd, 0x2f, 0x51, 0xd7, 0x3a, 0x77, 0xa3, 0xaf, 0xa1, 0xf1, 0xcd, 0x1f, 0x6c, 0xcf, 0xd0, 0xae,
    0xc2, 0x65, 0x46, 0x2c, 0x00, 0x0e, 0xcd, 0x69, 0x8a, 0x34, 0x7d, 0x82, 0x5d, 0x9c, 0xb1, 0x2c,
];

/// Canonical Fp2-encoded bytes for `BASIS_E0_QX`, the x-coordinate of the
/// second generator of the `2^f`-torsion basis on E0.
pub const BASIS_E0_QX_BYTES: [u8; 96] = [
    0x23, 0xa7, 0x10, 0xc1, 0x17, 0x60, 0xc5, 0x5b, 0x68, 0xb8, 0x7d, 0xe5, 0x66, 0xdb, 0x7a, 0xd6,
    0x6c, 0xff, 0x70, 0x95, 0xbc, 0xd0, 0xd8, 0x24, 0x4a, 0x16, 0x12, 0xf4, 0xec, 0xb9, 0xe5, 0x4b,
    0xf3, 0xaf, 0x19, 0x68, 0x06, 0xad, 0x24, 0x4a, 0xc9, 0xd1, 0x31, 0x6e, 0xad, 0x13, 0x92, 0x12,
    0x35, 0x93, 0x7b, 0x97, 0x7e, 0xcb, 0x5f, 0x99, 0xcd, 0x77, 0x26, 0xfb, 0x1b, 0x4b, 0xf2, 0x22,
    0x5e, 0xe3, 0x58, 0x8d, 0x06, 0x9a, 0xb7, 0x4c, 0xa7, 0x82, 0xf7, 0x97, 0xba, 0x7a, 0x81, 0x31,
    0x04, 0xdc, 0x57, 0xea, 0x05, 0x4f, 0xf4, 0x35, 0x2a, 0xd4, 0x0f, 0xb1, 0x5c, 0x59, 0x2a, 0x03,
];

// These are level-independent but kept per-level for import consistency.

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

/// 10 precomputed 4x4 basis change matrices for splitting transforms.
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

/// 6 precomputed 4x4 normalization matrices for splitting.
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
