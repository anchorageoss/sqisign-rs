//!
//! `P` and `Q` generate `E0[2^f]` where `f = TORSION_EVEN_POWER`.
//! Only x-coordinates are stored (Montgomery model uses x-only arithmetic).
//! Fp elements are stored as Montgomery-form limb arrays (5 × 51-bit radix).
//!
//! Generated from SageMath precompute scripts. DO NOT EDIT.

/// Number of 64-bit limbs per `Fp` element (5 for Level 1).
pub const NWORDS_FIELD: usize = 5;

/// Real part of `P.x` where `P` is the first generator of `E0[2^248]`.
pub const BASIS_E0_P_X_RE: [u64; 5] = [
    0x5bcab12000c08,
    0x452654b56d052,
    0x26f81b5190a0a,
    0x36cfd66a361eb,
    0x12726610d11b,
];
/// Imaginary part of `P.x` where `P` is the first generator of `E0[2^248]`.
pub const BASIS_E0_P_X_IM: [u64; 5] = [
    0x6b96065c83efc,
    0x29da1d4a82cd9,
    0x190797ab98bdf,
    0x6841aa6eeee05,
    0x1377c5431166,
];
/// Real part of `Q.x` where `Q` is the second generator of `E0[2^248]`.
pub const BASIS_E0_Q_X_RE: [u64; 5] = [
    0x21dd55b97832f,
    0x210f2d30b26ad,
    0x680bcfcf6396,
    0x27b318ec126a7,
    0x4ffba5956012,
];
/// Imaginary part of `Q.x` where `Q` is the second generator of `E0[2^248]`.
pub const BASIS_E0_Q_X_IM: [u64; 5] = [
    0x74590149117e3,
    0x4982edefcc606,
    0x2ae3db0cc6884,
    0x7d0384872f5ec,
    0x4fbb0fcb5a52,
];
