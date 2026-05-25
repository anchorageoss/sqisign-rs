//!
//! Fp elements are stored as Montgomery-form limb arrays.
//! Generated from SageMath precompute scripts. DO NOT EDIT.

pub const NWORDS_FIELD: usize = 9;

pub const BASIS_E0_P_X_RE: [u64; 9] = [
    0xa6f4f854fc265d,
    0x4c8b84aa5fb427,
    0x1309a22f7f8bedd,
    0x12326c230bb7339,
    0x1177007f8e443b7,
    0x3227e897204471,
    0x173c12694021af7,
    0xd9af272f428697,
    0x523eda847ef5,
];
pub const BASIS_E0_P_X_IM: [u64; 9] = [
    0xded07bbc792c63,
    0x11e1b26dec9cebc,
    0xc046644c2b6cd7,
    0x10fa781cd249b7d,
    0x1c100f6a2ab7eb,
    0x3268453a15b6a9,
    0x54d2827aa042c2,
    0x1976f2e8b7c96ec,
    0x16e01b2e8125f,
];
pub const BASIS_E0_Q_X_RE: [u64; 9] = [
    0xf73af643285709,
    0xf149be2f088d45,
    0xcd261395ea3c0a,
    0x3a51f18f48bd2c,
    0x20878d18902069,
    0x1dde2b7d4cfad79,
    0x1cfc83af281db52,
    0xcb86b4138f7754,
    0xb4deb1e3f8a7,
];
pub const BASIS_E0_Q_X_IM: [u64; 9] = [
    0x924c7a12f1ab1e,
    0x37608c2f01a03,
    0x15ab8f95ccf5c3e,
    0x99e325091f7251,
    0xc375ef1b0a8b52,
    0x1e7185439fe829e,
    0x1393f18a069901e,
    0x1171a261ad16dd5,
    0x6573978c1c85,
];
