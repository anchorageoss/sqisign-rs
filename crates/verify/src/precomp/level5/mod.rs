//!
//! Contains EC layer constants and signing-layer constants
//! (quaternion data, endomorphism actions, torsion degrees).

pub mod e0_basis;
pub mod ec_params;

pub use ec_params::{P_COFACTOR_FOR_2F, P_COFACTOR_FOR_2F_BITLENGTH, TORSION_EVEN_POWER};

/// Canonical 𝔽p²-encoded bytes for `BASIS_E0_PX`, the x-coordinate of the
/// first generator of the 2ᶠ-torsion basis on E0.
pub const BASIS_E0_PX_BYTES: [u8; 128] = [
    0xc0, 0xb4, 0x87, 0xb1, 0xde, 0x5c, 0x02, 0x80, 0x81, 0xdf, 0x6f, 0x2b, 0x61, 0x2a, 0x84, 0xe7,
    0x04, 0x36, 0x08, 0xc4, 0xc4, 0x49, 0xf1, 0xc8, 0x9e, 0x45, 0x2d, 0x2e, 0x92, 0x09, 0xa0, 0x6b,
    0x68, 0xbc, 0xf9, 0x37, 0x67, 0xf0, 0x87, 0x35, 0x4a, 0x1f, 0xb0, 0x71, 0x38, 0x23, 0xbe, 0x6a,
    0x02, 0xc8, 0x10, 0xf0, 0x87, 0xe4, 0xd5, 0x13, 0x1f, 0xcb, 0x5f, 0x08, 0xe5, 0xaf, 0x9f, 0x00,
    0x3d, 0x06, 0x67, 0xe8, 0xff, 0xfb, 0xdd, 0x8c, 0xf4, 0x97, 0xd1, 0xcd, 0x82, 0x57, 0x8f, 0x12,
    0x81, 0x92, 0xed, 0x48, 0x8e, 0x19, 0xb7, 0x28, 0xac, 0xaa, 0xa9, 0xeb, 0xb2, 0x40, 0x19, 0xd0,
    0xd8, 0x7a, 0xd4, 0x43, 0x4a, 0x20, 0x4d, 0xff, 0x50, 0xf1, 0x14, 0xce, 0xdd, 0x85, 0x67, 0xe4,
    0xa1, 0xab, 0x8f, 0xd8, 0xa2, 0xe7, 0xe2, 0xd3, 0x80, 0xcf, 0xf3, 0x6e, 0x51, 0x2a, 0xc4, 0x00,
];

/// Canonical 𝔽p²-encoded bytes for `BASIS_E0_QX`, the x-coordinate of the
/// second generator of the 2ᶠ-torsion basis on E0.
pub const BASIS_E0_QX_BYTES: [u8; 128] = [
    0x09, 0xaa, 0xc8, 0xbf, 0x99, 0x40, 0x65, 0x6e, 0x5e, 0xef, 0x70, 0x2c, 0x57, 0xa3, 0x8a, 0xf6,
    0x19, 0x82, 0x60, 0xb4, 0xeb, 0xfc, 0x2a, 0x87, 0x1e, 0xa9, 0x09, 0x8b, 0xb5, 0x99, 0xcb, 0xf0,
    0x25, 0xb7, 0xc4, 0x4b, 0xe7, 0x71, 0x63, 0xe9, 0x66, 0x7e, 0x43, 0x64, 0x27, 0x6f, 0x3c, 0xd4,
    0xd1, 0x57, 0x86, 0x40, 0xf8, 0xd3, 0xdc, 0x3d, 0xcd, 0x59, 0x18, 0xe6, 0x1b, 0xe9, 0xbc, 0x00,
    0xe8, 0xfb, 0xff, 0x86, 0x70, 0x8b, 0xe5, 0x2b, 0xf5, 0xc8, 0x2e, 0xd7, 0x45, 0x51, 0xd7, 0x3b,
    0x99, 0x49, 0x30, 0x52, 0xa7, 0x43, 0x94, 0xfd, 0x39, 0xbd, 0x0e, 0xd7, 0x58, 0x0e, 0xed, 0x1b,
    0xb8, 0x4e, 0xb4, 0xc6, 0xd7, 0x1d, 0xe2, 0x90, 0x5f, 0x87, 0x66, 0x44, 0xd6, 0x31, 0x4f, 0x13,
    0x1e, 0x56, 0x5a, 0x50, 0x13, 0x08, 0xbf, 0xad, 0x10, 0x38, 0xde, 0x87, 0x99, 0x8c, 0xa4, 0x00,
];

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
