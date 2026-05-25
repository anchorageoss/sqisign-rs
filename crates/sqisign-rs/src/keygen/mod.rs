//!
//! This crate provides keypair generation and the `SecretKey` type.
//! It depends on the full SQIsign stack including quaternion algebra.
//!
//! # Note on `std` requirement
//!
//! This crate currently requires `std` because `sqisign-quaternion`
//! depends on `num-bigint`. When the quaternion layer is upgraded to
//! a `no_std`-compatible big integer library, this crate will become
//! `no_std` + `alloc`.

use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::quaternion::types::{IbzMat2x2, QuatLeftIdeal};
use sqisign_verify::ec::{EcBasis, EcCurve};
use sqisign_verify::fp::FpBackend;
use sqisign_verify::params::Level1;
use sqisign_verify::types::PublicKey;
use zeroize::Zeroize;

#[allow(clippy::module_inception)]
pub mod keygen;

pub(crate) mod sk_encoding;

/// Generate a fresh SQIsign keypair.
///
/// Precomputed constants are constructed automatically for the chosen
/// security level.
pub fn keypair<L: HasSigningPrecomp + sqisign_verify::precomp::LevelPrecomp>(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> (PublicKey<L>, SecretKey<L>) {
    let precomp = L::signing_precomp();
    keygen::protocols_keygen(rng, &precomp)
}

/// A secret key for SQIsign signing.
///
/// Contains the secret endomorphism ideal, the basis-change matrix
/// relating the canonical torsion basis to the image of the E0 basis
/// under the secret isogeny, and the canonical basis itself.
pub struct SecretKey<L: FpBackend + sqisign_verify::precomp::LevelPrecomp = Level1> {
    /// The public curve (same A-coefficient as the public key).
    pub curve: EcCurve<L>,
    /// The secret left ideal of O0 encoding the secret isogeny.
    pub secret_ideal: QuatLeftIdeal,
    /// 2x2 change-of-basis matrix `M` encoding the coordinates of
    /// `BA_can` in the basis `BA0_two`: `(M * v) . BA0_two = v . BA_can`,
    /// where `BA_can` is the canonical `2^f`-torsion basis on the public
    /// curve and `BA0_two` is the image of the E0 torsion basis
    /// through the secret isogeny. Entries are mod `2^TORSION_EVEN_POWER`.
    pub mat_ba_can_to_ba0_two: IbzMat2x2,
    /// Canonical `2^f`-torsion basis on the public curve.
    pub canonical_basis: EcBasis<L>,
}

impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> core::fmt::Debug for SecretKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SecretKey([REDACTED])")
    }
}

impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> core::fmt::Display for SecretKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SecretKey([REDACTED])")
    }
}

impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> Zeroize for SecretKey<L> {
    fn zeroize(&mut self) {
        self.curve.zeroize();
        self.secret_ideal.zeroize();
        self.mat_ba_can_to_ba0_two.zeroize();
        self.canonical_basis.zeroize();
    }
}

impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> Drop for SecretKey<L> {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> zeroize::ZeroizeOnDrop for SecretKey<L> {}
