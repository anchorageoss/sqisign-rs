//!
//! Provides big integer operations, compound types (vectors, matrices,
//! quaternion elements, lattices, ideals), and the algorithms needed
//! for the signing path of SQIsign v2.0.
//!
//! # Big Integer Status
//!
//! This crate currently uses `num-bigint` for arbitrary-precision integers.
//! This is a temporary choice for correctness-first development. The
//! long-term plan is to replace `num-bigint` with either:
//! - `crypto-bigint` (constant-time, no_std)
//! - A custom fixed-precision implementation based on Won Kim et al.
//!   (ePrint 2025/1649) worst-case bounds
//!
//! **This crate requires `std` due to `num-bigint`.** It is intentionally
//! excluded from the verification path, which remains `no_std`.

// SECURITY: num-bigint is NOT constant-time. Signing path only; not used in verification.

pub mod algebra;
pub mod dim2;
pub mod dim4;
pub mod dpe;
pub mod fast_modpow;
pub mod hnf;
pub mod ideal;
pub mod intbig;
pub mod integers;
pub mod lat_ball;
pub mod lattice;
pub mod lll;
pub mod montgomery;
pub mod normeq;
pub mod rational;
pub mod types;
