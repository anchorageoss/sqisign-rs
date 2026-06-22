//! Dimension-4 theta-model arithmetic for SQIsignHD verification.
//!
//! This is the foundational arithmetic layer of the dimension-4 SQIsignHD
//! verifier: theta points (16 𝔽p² coordinates), theta structures (an abelian
//! fourfold characterised by its theta null point), the Hadamard transform,
//! and theta doubling. It deliberately contains **no isogenies** - those are
//! built on top of this layer in later phases.
//!
//! It mirrors the structure of the dimension-2 theta module in
//! `sqisign-verify` (`crate::theta`), but with 16 coordinates and the
//! dimension-4 formulas, and it reuses the existing `Fp2` field arithmetic.
//!
//! # Projective coordinates (read this first)
//!
//! A theta point - and in particular a theta **null** point - is defined only
//! up to a common non-zero scalar: `(λ·P₀ : … : λ·P₁₅)` and `(P₀ : … : P₁₅)`
//! denote the same point. Equality therefore **must** cross-multiply against a
//! pivot ([`ThetaPointDim4::projective_eq`]); comparing coordinates directly is
//! a bug. This is the single most important correctness fact ported from the
//! Phase 0 oracle.
//!
//! # Index convention
//!
//! The 16 coordinates are indexed by a 4-bit value
//! `k = i₀ + 2·i₁ + 4·i₂ + 8·i₃`, matching the sage `Theta_dim4` reference
//! (`multindex_to_index`).
//!
//! # Reference
//!
//! Ported from the `Theta_dim4` sage package by Pierrick Dartois
//! (<https://github.com/Pierrick-Dartois/Theta_dim4>), specifically
//! `theta_structures/Theta_dim4.py` and `theta_structures/theta_helpers_dim4.py`.
//!
//! # Module organization
//!
//! This module was previously the standalone `theta4` crate; it now lives
//! inside `sqisign-verify` as the dimension-4 counterpart to the dim-2
//! [`crate::theta`] module. It reuses the crate's `Fp`/`Fp2`/`ec` primitives
//! directly (`crate::…`). The optimal-strategy chain loop ([`strategy`]) keeps a
//! small, data-dependent stack of intermediate kernel bases (bounded by the
//! chain length, ~67 at Level 1) in heap `Vec`s - the only `alloc` use, and off
//! the constant-time path (verification orchestration, not field arithmetic).

pub mod arith;
pub mod basis;
pub mod canonical;
pub mod challenge;
pub mod chain;
pub mod dim2;
pub mod dim4;
pub mod field;
pub mod gluing;
pub mod gluing_chain;
pub mod hd_verify;
pub mod isogeny;
pub mod kani;
pub mod point;
pub mod product_theta;
pub mod response;
pub mod self_contained;
pub mod strategy;
pub mod structure;
pub mod wire;
mod nqr_tables_l1;
mod nqr_tables_l3;
mod nqr_tables_l5;

pub use arith::{act_point, hadamard, pointwise_square, to_squared_theta};
pub use basis::{
    canonical_hints, canonical_hints_l1, hd_torsion_basis, hd_torsion_basis_l1, jac_to_affine,
    torsion_basis_2f_from_hint, HdNqr,
};
pub use canonical::make_canonical;
pub use challenge::{recover_challenge_l1, ChallengeRecovery};
pub use dim2::{
    apply_mat4, base_change_theta_dim2, hadamard2, squared_theta2, GluingThetaIsogenyDim2,
    IsogenyChainDim2, ThetaIsogenyDim2, ThetaStructureDim2, TuplePoint,
};
pub use chain::{middle_codomain_matches, run_half_chain, run_half_chain_collect};
pub use dim4::{apply_base_change_theta_dim4, base_change_theta_dim4};
pub use gluing::{GluingIsogenyDim4, GLUING_KERNEL_DIRS};
pub use gluing_chain::{
    jac_mul_u128, point_matrix_product_k, KaniGluingChainHalf, TuplePoint4,
};
pub use hd_verify::{
    hd_challenge, hd_challenge_from_curves, hd_challenge_len, hd_verify, hd_verify_checked,
    recover_response_cd, HdReject, HdVerifyInputs,
};
pub use isogeny::IsogenyDim4;
pub use kani::{
    complete_symplectic_dim4, gluing_bc_dim4_f1, gluing_bc_dim4_f2, gluing_dim2_f1, gluing_dim2_f2,
    inverse_mod_pow2, is_symplectic_dim4, kernel_matrix_f1, kernel_matrix_f2_dual, matrix_f,
    matrix_f_dual, norm_equation_2f_minus_q, starting_two_symplectic_matrices, sum_of_two_squares,
    F_MATRIX_L1,
};
pub use point::{ThetaPointDim4, THETA_DIM4_N};
pub use product_theta::{
    product_null_dim2, product_theta_dim2, product_theta_dim2to4, ThetaStructureDim1,
};
pub use response::{recover_response_l1, ResponseRecovery, ResponseScalars};
pub use self_contained::{hd_image_l1, hd_verify_l1, hd_verify_l1_bool, HdSignatureL1};
pub use wire::{
    encode_public_key, encode_signature, hd_verify_bytes_l1, hd_verify_bytes_l1_bool,
    hd_verify_l1_parsed, parse_public_key, parse_signature, ParsedPublicKey, ParsedSignature,
    PK_WIRE_BYTES, SIG_WIRE_BYTES,
};
pub use strategy::{optimised_strategy, run_strategy_chain, StrategyChain};
pub use structure::ThetaStructureDim4;
