//!
//! Fp elements are stored as Montgomery-form limb arrays.
//! Generated from SageMath precompute scripts. DO NOT EDIT.

pub const NWORDS_FIELD: usize = 7;

pub const BASIS_E0_P_X_RE: [u64; 7] = [
    0x94635b7b34b8c,
    0x431475975ec8c7,
    0x380f3b6b0f3d6c,
    0x2e90ddd88ba021,
    0x5eb0a59679b654,
    0x347706dc01cb41,
    0xb7765ed4a44a5,
];
pub const BASIS_E0_P_X_IM: [u64; 7] = [
    0x412c2c3df0cc54,
    0x2338803450b7d0,
    0x206883ec0e5d2f,
    0x407a7e72205c5d,
    0x187f5a00661d99,
    0x5905b6352b7e4d,
    0x3032c0ad99418,
];
pub const BASIS_E0_Q_X_RE: [u64; 7] = [
    0x7333ee4f4818b7,
    0x29c73aefc7681b,
    0x3db742e2128546,
    0x3f8774b65cc12a,
    0x332cf22a3425e2,
    0x4a219e343591d2,
    0x6d1dfdb6ea8ff,
];
pub const BASIS_E0_Q_X_IM: [u64; 7] = [
    0x1fdc82b11838c5,
    0x681f359137f9af,
    0x3eb05affc54924,
    0x509e310ef21e09,
    0x5a97b9d957fd56,
    0x6e7c043e0db389,
    0x4fbc3aab7429d,
];
