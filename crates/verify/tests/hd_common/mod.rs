//! Shared helpers for the dim-4 (`hd`) integration tests: JSON loading and
//! parsing of the `[re_hex, im_hex]` 𝔽p² serialization used by the oracle vectors.

// Each integration-test binary that does `mod common;` pulls in the whole
// module but uses only a subset of the helpers; silence the resulting
// per-binary dead-code warnings.
#![allow(dead_code)]

use serde_json::Value;
use sqisign_verify::hd::{ThetaPointDim4, THETA_DIM4_N};
use sqisign_verify::{Fp2, Level1};

pub type F = Fp2<Level1>;
pub type Pt = ThetaPointDim4<Level1>;

/// The Phase 0 ground-truth vectors (per-step codomain theta null points).
pub const PHASE0_VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/test_vectors_l1.json"
);

/// Read and parse a JSON file, panicking with context on failure.
pub fn load(path: &str) -> Value {
    let s = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    serde_json::from_str(&s).expect("valid json")
}

/// Big-endian hex (minimal width, optional `0x`) -> little-endian 32-byte array.
pub fn le32(hexstr: &str) -> [u8; 32] {
    let s = hexstr.trim_start_matches("0x");
    let s = if s.len() % 2 == 1 {
        format!("0{s}")
    } else {
        s.to_string()
    };
    let be = hex::decode(&s).expect("valid hex");
    assert!(be.len() <= 32, "coordinate exceeds 32 bytes");
    let mut le = [0u8; 32];
    for (i, b) in be.iter().rev().enumerate() {
        le[i] = *b;
    }
    le
}

/// Parse a `[re_hex, im_hex]` pair into an `Fp2<Level1>`.
pub fn parse_fp2(pair: &Value) -> F {
    let re = pair[0].as_str().expect("re hex");
    let im = pair[1].as_str().expect("im hex");
    let mut buf = [0u8; 64];
    buf[..32].copy_from_slice(&le32(re));
    buf[32..].copy_from_slice(&le32(im));
    F::decode(&buf).expect("coordinate in [0,p)")
}

/// Parse a 16-element coords array into a theta point.
pub fn parse_coords(coords: &Value) -> Pt {
    let arr = coords.as_array().expect("coords array");
    assert_eq!(arr.len(), THETA_DIM4_N, "theta point must have 16 coords");
    ThetaPointDim4::new(core::array::from_fn(|k| parse_fp2(&arr[k])))
}

/// Parse a `{ "coords": [...] }` node into a theta point.
pub fn parse_node(node: &Value) -> Pt {
    parse_coords(&node["coords"])
}

/// Exact (non-projective) `Fp2` equality.
pub fn fp2_eq(a: &F, b: &F) -> bool {
    bool::from(a.ct_equal(b))
}

/// An all-zero theta point, used to pre-fill output buffers.
pub fn zero_point() -> Pt {
    ThetaPointDim4::new(core::array::from_fn(|_| F::zero()))
}

/// A deterministic, definitely-non-zero `Fp2` scalar for re-scaling tests.
pub fn lambda() -> F {
    let re = F::from_small(0x9e37_79b9_7f4a_7c15);
    let im = F::i_element().mul(&F::from_small(0xc2b2_ae3d_27d4_eb4f));
    re.add(&im)
}
