//!
//! Pure Rust implementation of the SQIsign signature scheme (NIST PQC
//! Additional Signatures, Round 2/3 candidate).
//!
//! This crate provides key generation, signing, and re-exports everything
//! from `sqisign-verify` for verification. For verify-only usage (`no_std`),
//! depend on `sqisign-verify` directly.
//!
//! Two signature schemes are available, chosen at keygen time: the dimension-2
//! formats ([`generate`]) and the **compact** 108-byte format
//! ([`generate_compact`]). Verification autodetects the format from byte length
//! via [`AnySignature`]. At Level 1, dim-2 verification is ~4-7 ms and compact
//! verification ~33 ms (~20.5 ms with the `parallel` feature).
//!
//! ## Quick Start
//!
//! ```
//! use sqisign_rs::{generate, generate_compact, PublicKey, SigningKey, Verifier};
//!
//! # fn main() -> Result<(), sqisign_rs::Error> {
//! let mut rng = rand::rngs::OsRng;
//!
//! // Standard (dimension-2, 148 bytes at Level 1):
//! let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
//! let sig = sk.sign(b"hello world", &mut rng)?;
//! pk.verify(b"hello world", &sig)?;
//!
//! // Compact (smallest signature, 108 bytes at Level 1):
//! let (cpk, csk) = generate_compact(&mut rng);
//! let csig = csk.sign(b"hello world", &mut rng)?;
//! cpk.verify(b"hello world", &csig)?;
//! # Ok(())
//! # }
//! ```

pub mod alloc;
pub mod id2iso;
pub mod keygen;
pub mod precomp_signing;
pub mod quaternion;
pub mod sign;

// Re-export everything from sqisign-verify, including the unified
// `AnySignature` (which now autodetects dim-2 *and* dim-4 SQIsignHD by length)
// and `PublicKey` (which accepts both the 65-byte dim-2 and 64-byte HD key
// encodings). The dim-4 verifier itself lives in `sqisign_verify::hd`.
pub use sqisign_verify::*;

// Public API.
pub use keygen::SecretKey;

// The compact (108-byte-signature) scheme's signing-side entry points. The
// matching verification types `CompactPublicKey` / `CompactSignature` come from
// the `pub use sqisign_verify::*` re-export above.
pub use sign::{generate_compact, CompactSigningKey};

use hybrid_array::typenum::Unsigned;
use id2iso::sign_precomp::HasSigningPrecomp;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// A signing key that bundles everything needed to produce signatures.
///
/// Created by [`generate`]. Holds the secret key and public key.
pub struct SigningKey<L: sqisign_verify::fp::FpBackend + LevelPrecomp = Level1> {
    sk: SecretKey<L>,
    pk: PublicKey<L>,
}

impl<L: HasSigningPrecomp + LevelPrecomp> SigningKey<L> {
    /// Sign a message.
    #[inline]
    pub fn sign(
        &self,
        msg: &[u8],
        rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
    ) -> Result<Signature<L>, Error> {
        crate::sign::sign(&self.sk, &self.pk, msg, rng)
    }

    /// The public key corresponding to this signing key.
    #[inline]
    pub fn public_key(&self) -> &PublicKey<L> {
        &self.pk
    }

    /// Encode the signing key as `secret_key_bytes || public_key_bytes`.
    pub fn to_bytes(&self) -> Result<Vec<u8>, Error> {
        let sk_bytes = self.sk.to_bytes()?;
        let pk_bytes = self.pk.to_bytes();
        let mut out = Vec::with_capacity(L::SkLen::USIZE + L::PkLen::USIZE);
        out.extend_from_slice(&sk_bytes);
        out.extend_from_slice(&pk_bytes);
        Ok(out)
    }

    /// Decode a signing key from `secret_key_bytes || public_key_bytes`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let expected = L::SkLen::USIZE + L::PkLen::USIZE;
        if bytes.len() != expected {
            return Err(Error::InvalidLength);
        }
        let (sk_bytes, pk_bytes) = bytes.split_at(L::SkLen::USIZE);
        let mut sk = SecretKey::<L>::from_bytes(sk_bytes)?;
        let pk = PublicKey::<L>::from_bytes(pk_bytes)?;
        sk.populate_from_pk(&pk);
        Ok(SigningKey { sk, pk })
    }
}

impl<L: sqisign_verify::fp::FpBackend + LevelPrecomp> core::fmt::Debug for SigningKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SigningKey([REDACTED])")
    }
}

impl<L: sqisign_verify::fp::FpBackend + LevelPrecomp> core::fmt::Display for SigningKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SigningKey([REDACTED])")
    }
}

impl<L: id2iso::sign_precomp::HasSigningPrecomp + LevelPrecomp> SigningKey<L> {
    /// Return a wrapper that prints the raw signing key bytes as hex.
    ///
    /// # Security Warning
    ///
    /// The returned type implements [`Display`](core::fmt::Display) and
    /// will output **secret key material in plaintext**. Use only for
    /// debugging in secure, ephemeral environments. Never log this
    /// output in production, persist it to files, or transmit it over
    /// untrusted channels.
    #[inline]
    pub fn display_secret(&self) -> SigningKeyDisplay<'_, L> {
        SigningKeyDisplay(self)
    }
}

/// Wrapper returned by [`SigningKey::display_secret`] that prints raw
/// signing key bytes as lowercase hex.
///
/// # Security Warning
///
/// This will output secret key material in plaintext when formatted.
pub struct SigningKeyDisplay<'a, L: sqisign_verify::fp::FpBackend + LevelPrecomp>(
    &'a SigningKey<L>,
);

impl<L: id2iso::sign_precomp::HasSigningPrecomp + LevelPrecomp> core::fmt::Display
    for SigningKeyDisplay<'_, L>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self.0.to_bytes() {
            Ok(bytes) => {
                for &b in bytes.as_slice() {
                    write!(f, "{b:02x}")?;
                }
                Ok(())
            }
            Err(_) => f.write_str("<encoding error>"),
        }
    }
}

impl<L: sqisign_verify::fp::FpBackend + LevelPrecomp> Zeroize for SigningKey<L> {
    fn zeroize(&mut self) {
        self.sk.zeroize();
    }
}

impl<L: sqisign_verify::fp::FpBackend + LevelPrecomp> Drop for SigningKey<L> {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl<L: sqisign_verify::fp::FpBackend + LevelPrecomp> ZeroizeOnDrop for SigningKey<L> {}

impl<L: HasSigningPrecomp + LevelPrecomp> signature::RandomizedSigner<Signature<L>>
    for SigningKey<L>
{
    fn try_sign_with_rng(
        &self,
        rng: &mut impl signature::rand_core::CryptoRngCore,
        msg: &[u8],
    ) -> Result<Signature<L>, signature::Error> {
        self.sign(msg, rng).map_err(|_| signature::Error::new())
    }
}

impl<L: HasSigningPrecomp + LevelPrecomp> signature::Keypair for SigningKey<L> {
    type VerifyingKey = PublicKey<L>;

    fn verifying_key(&self) -> Self::VerifyingKey {
        self.pk.clone()
    }
}

/// Generate a fresh SQIsign keypair.
///
/// Returns the public key (for the verifier) and a signing key (for
/// the signer). Level 1 (128-bit security) is the default; specify
/// `generate::<Level3>` or `generate::<Level5>` for higher levels.
pub fn generate<L: HasSigningPrecomp + LevelPrecomp>(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> (PublicKey<L>, SigningKey<L>) {
    let precomp = L::signing_precomp();
    let (pk, sk) = keygen::keygen::protocols_keygen(rng, &precomp);
    let signing_key = SigningKey { sk, pk: pk.clone() };
    (pk, signing_key)
}
