//! SQIsign round 3 in pure Rust: key generation, signing, verification.
//!
//! SQIsign as submitted to the third round of NIST's additional-signatures
//! process (`the-sqisign` at tag `nist-v3`): the primes `3 * 2^324 - 1`,
//! `27 * 2^500 - 1` and `17 * 2^664 - 1`, 83 / 129 / 169-byte public keys,
//! 200 / 306 / 406-byte signatures. The 300 known-answer vectors of the
//! reference are reproduced byte for byte (`crates/kat`). Round 2 is not
//! supported: round-2 keys and signatures are not accepted and there is no
//! conversion (0.5.0 removed the round-2 implementation).
//!
//! Verification lives in [`sqisign_verify`], which this crate re-exports;
//! a verifier depends on that crate alone. This crate adds the signing
//! side: fixed-precision integers ([`mp`], no `num-bigint`), the quaternion
//! layer in inert representation with constant-time lattice reduction
//! ([`quat`]), the ideal-to-isogeny translation ([`id2iso`]), the signing
//! constants ([`precomp`]) and the protocol ([`sqisign`]), with
//! [`SigningKey`] as the typed surface.
//!
//! ```
//! use sqisign_rs::{generate, Level1, PublicKey, SigningKey, Verifier};
//!
//! # fn main() -> Result<(), sqisign_rs::Error> {
//! let mut rng = rand_core::OsRng;
//! let (pk, sk): (PublicKey<Level1>, SigningKey<Level1>) = generate(&mut rng);
//! let sig = sk.sign(b"hello world", &mut rng)?;
//! pk.verify(b"hello world", &sig)?;
//!
//! // the wire: 83-byte key, 200-byte signature at level I
//! let pk2 = PublicKey::<Level1>::from_bytes(&pk.to_bytes())?;
//! pk2.verify(b"hello world", &sig)?;
//! # Ok(())
//! # }
//! ```
//!
//! Levels III and V are the type parameter (`generate::<Level3>`). The
//! library is `no_std` with `alloc` (the secret-key encoding returns a
//! `Vec`); unit tests link `std`.

#![cfg_attr(not(test), no_std)]
#![warn(missing_docs)]

extern crate alloc;

#[cfg(feature = "compact")]
pub mod compact;
pub mod id2iso;
pub mod mp;
pub mod precomp;
pub mod quat;
pub mod secure_alloc;
pub mod sqisign;

#[cfg(feature = "compact")]
pub use compact::{generate_compact, CompactSignError, CompactSigningKey};
pub use secure_alloc::ZeroizingAllocator;
#[cfg(feature = "compact")]
pub use sqisign_verify::compact::{CompactLevel, CompactPublicKey, CompactSignature};
pub use sqisign_verify::{
    signature, CompressedSignature, Error, Level, Level1, Level3, Level5, PreparedPublicKey,
    PublicKey, Signature, SignatureEncoding, Verifier,
};

use alloc::vec::Vec;
use sqisign_verify::rng::Rng;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// A level with its signing-side parameters: the quaternion algebra and
/// the endomorphism data, and the secret key of the matching width.
pub trait SigningLevel: Level {
    /// The signing parameters ([`sqisign::Params`] at the level's width).
    type Params: 'static;
    /// The protocol-level secret key.
    type SecretKey: Clone + Zeroize;
    /// The parameters, built from the constants.
    fn params() -> Self::Params;
    /// `protocols_keygen`.
    fn keygen<R: Rng>(
        params: &Self::Params,
        rng: &mut R,
    ) -> Option<(sqisign_verify::sqisign::PublicKey<Self>, Self::SecretKey)>;
    /// `protocols_sign`.
    fn sign<R: Rng>(
        params: &Self::Params,
        pk: &sqisign_verify::sqisign::PublicKey<Self>,
        sk: &Self::SecretKey,
        msg: &[u8],
        rng: &mut R,
    ) -> Option<sqisign_verify::sqisign::Signature<Self>>;
    /// `secret_key_to_bytes`.
    fn secret_key_to_bytes(
        params: &Self::Params,
        sk: &Self::SecretKey,
        pk: &sqisign_verify::sqisign::PublicKey<Self>,
    ) -> Vec<u8>;
    /// `secret_key_from_bytes`, with the embedded public key.
    fn secret_key_from_bytes(
        params: &Self::Params,
        bytes: &[u8],
    ) -> Option<(Self::SecretKey, sqisign_verify::sqisign::PublicKey<Self>)>;
}

macro_rules! signing_level {
    ($prime:ty, $modname:ident, $level:ident) => {
        impl SigningLevel for $prime {
            type Params = sqisign::Params<{ crate::precomp::$modname::IBZ_NLIMBS }>;
            type SecretKey = sqisign::SecretKey<$prime, { crate::precomp::$modname::IBZ_NLIMBS }>;
            fn params() -> Self::Params {
                sqisign::$level()
            }
            fn keygen<R: Rng>(
                params: &Self::Params,
                rng: &mut R,
            ) -> Option<(sqisign_verify::sqisign::PublicKey<Self>, Self::SecretKey)> {
                sqisign::keygen::<$prime, { crate::precomp::$modname::IBZ_NLIMBS }>(params, rng)
            }
            fn sign<R: Rng>(
                params: &Self::Params,
                pk: &sqisign_verify::sqisign::PublicKey<Self>,
                sk: &Self::SecretKey,
                msg: &[u8],
                rng: &mut R,
            ) -> Option<sqisign_verify::sqisign::Signature<Self>> {
                sqisign::sign::<$prime, { crate::precomp::$modname::IBZ_NLIMBS }>(
                    params, pk, sk, msg, rng,
                )
            }
            fn secret_key_to_bytes(
                params: &Self::Params,
                sk: &Self::SecretKey,
                pk: &sqisign_verify::sqisign::PublicKey<Self>,
            ) -> Vec<u8> {
                sqisign::secret_key_to_bytes(params, sk, pk)
            }
            fn secret_key_from_bytes(
                params: &Self::Params,
                bytes: &[u8],
            ) -> Option<(Self::SecretKey, sqisign_verify::sqisign::PublicKey<Self>)> {
                sqisign::secret_key_from_bytes::<$prime, { crate::precomp::$modname::IBZ_NLIMBS }>(
                    params, bytes,
                )
            }
        }
    };
}

signing_level!(sqisign_verify::P324_3, p324_3, level1);
signing_level!(sqisign_verify::P500_27, p500_27, level3);
signing_level!(sqisign_verify::P664_17, p664_17, level5);

/// A `rand_core` generator as the protocol's entropy source.
struct Entropy<'a, R: rand_core::RngCore + rand_core::CryptoRng>(&'a mut R);

impl<R: rand_core::RngCore + rand_core::CryptoRng> Rng for Entropy<'_, R> {
    fn fill(&mut self, out: &mut [u8]) -> bool {
        self.0.try_fill_bytes(out).is_ok()
    }
}

/// A signing key: the secret key, its public key and the level's signing
/// parameters. Zeroised on drop.
pub struct SigningKey<L: SigningLevel> {
    sk: L::SecretKey,
    pk: PublicKey<L>,
    params: L::Params,
}

/// Generate a key pair. Loops on the negligible events the reference's key
/// generation retries on, and on a failure of `rng`.
pub fn generate<L: SigningLevel>(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> (PublicKey<L>, SigningKey<L>) {
    let params = L::params();
    let (pk, sk) = loop {
        if let Some(pair) = L::keygen(&params, &mut Entropy(rng)) {
            break pair;
        }
    };
    let pk = PublicKey::from_inner(pk);
    (pk.clone(), SigningKey { sk, pk, params })
}

impl<L: SigningLevel> SigningKey<L> {
    /// Sign `msg`. `Err(Error::InternalError)` on an entropy failure or on
    /// the negligible non-coprimality events the reference signs off on.
    pub fn sign(
        &self,
        msg: &[u8],
        rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
    ) -> Result<Signature<L>, Error> {
        L::sign(
            &self.params,
            self.pk.inner(),
            &self.sk,
            msg,
            &mut Entropy(rng),
        )
        .map(Signature::from_inner)
        .ok_or(Error::InternalError)
    }

    /// The public key.
    pub fn public_key(&self) -> &PublicKey<L> {
        &self.pk
    }

    /// The reference's secret-key encoding: the public key, the secret
    /// ideal in inert form, the basis-change matrix (270 / 417 / 549 bytes).
    pub fn to_bytes(&self) -> Vec<u8> {
        L::secret_key_to_bytes(&self.params, &self.sk, self.pk.inner())
    }

    /// Decode the reference's secret-key encoding.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::SK_BYTES {
            return Err(Error::InvalidLength);
        }
        let params = L::params();
        let (sk, pk) = L::secret_key_from_bytes(&params, bytes).ok_or(Error::MalformedInput)?;
        Ok(Self {
            sk,
            pk: PublicKey::from_inner(pk),
            params,
        })
    }
}

impl<L: SigningLevel> core::fmt::Debug for SigningKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SigningKey([REDACTED])")
    }
}

impl<L: SigningLevel> Zeroize for SigningKey<L> {
    fn zeroize(&mut self) {
        self.sk.zeroize();
    }
}

impl<L: SigningLevel> Drop for SigningKey<L> {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl<L: SigningLevel> ZeroizeOnDrop for SigningKey<L> {}

impl<L: SigningLevel> signature::RandomizedSigner<Signature<L>> for SigningKey<L> {
    fn try_sign_with_rng(
        &self,
        rng: &mut impl signature::rand_core::CryptoRngCore,
        msg: &[u8],
    ) -> Result<Signature<L>, signature::Error> {
        self.sign(msg, rng).map_err(|_| signature::Error::new())
    }
}
