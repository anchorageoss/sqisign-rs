//!
//! This crate provides the signing protocol. Key generation and the
//! `SecretKey` type live in `sqisign-keygen`.
//!
//! # Note on `std` requirement
//!
//! This crate currently requires `std` because `sqisign-quaternion`
//! depends on `num-bigint`. When the quaternion layer is upgraded to
//! a `no_std`-compatible big integer library, this crate will become
//! `no_std` + `alloc`.

pub use crate::keygen::SecretKey;

use crate::id2iso::sign_precomp::HasSigningPrecomp;
use sqisign_verify::types::{PublicKey, Signature};

#[allow(clippy::module_inception)]
pub mod sign;

/// Sign a message using the SQIsign protocol.
///
/// Precomputed constants are constructed automatically from the
/// security level.
pub fn sign<L: HasSigningPrecomp + sqisign_verify::precomp::LevelPrecomp>(
    sk: &SecretKey<L>,
    pk: &PublicKey<L>,
    msg: &[u8],
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Result<Signature<L>, sqisign_verify::Error> {
    sign::protocols_sign(pk, sk, msg, rng)
}
