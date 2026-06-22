//! The compact signing-side API: [`CompactSigningKey`] and [`generate_compact`].
//!
//! Mirrors the dim-2 [`crate::SigningKey`] / [`crate::generate`] pattern for the
//! compact (108-byte) scheme. The verification-side types [`CompactPublicKey`]
//! and [`CompactSignature`] live in `sqisign-verify` (re-exported from this
//! crate); this module adds the secret
//! key and the keygen/sign entry points.
//!
//! Compact keys are **distinct** from dim-2 keys: the schemes use different
//! torsion-basis conventions and are not interchangeable. A [`CompactSigningKey`]
//! signs into a [`CompactSignature`], verified by a [`CompactPublicKey`]; the
//! dim-2 [`crate::PublicKey`] does not verify compact signatures and vice versa.
//! The compact scheme is implemented at Level 1.

use zeroize::{Zeroize, ZeroizeOnDrop};

use sqisign_verify::hd::encode_public_key;
use sqisign_verify::{CompactPublicKey, CompactSignature, Error, Level1, SecurityLevel};

use crate::sign::dim4::{dim4_keygen, dim4_sign, Dim4PublicKey, Dim4SecretKey};

/// A compact signing key - produces 108-byte signatures.
///
/// Created by [`generate_compact`]. Distinct from the dim-2
/// [`crate::SigningKey`]; the two schemes' keys are not interchangeable.
pub struct CompactSigningKey<L: SecurityLevel = Level1> {
    sk: Dim4SecretKey,
    pk_dim4: Dim4PublicKey,
    pk: CompactPublicKey<L>,
}

impl CompactSigningKey<Level1> {
    /// Sign a message, producing a 108-byte [`CompactSignature`].
    ///
    /// The caller supplies a cryptographic RNG (use `OsRng` in production, or
    /// the NIST DRBG for KATs).
    #[inline]
    pub fn sign(
        &self,
        msg: &[u8],
        rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
    ) -> Result<CompactSignature<Level1>, Error> {
        let bytes = dim4_sign(&self.pk_dim4, &self.sk, msg, rng)?;
        CompactSignature::from_bytes(&bytes)
    }

    /// The compact public key corresponding to this signing key.
    #[inline]
    pub fn public_key(&self) -> &CompactPublicKey<Level1> {
        &self.pk
    }
}

impl signature::RandomizedSigner<CompactSignature<Level1>> for CompactSigningKey<Level1> {
    fn try_sign_with_rng(
        &self,
        rng: &mut impl signature::rand_core::CryptoRngCore,
        msg: &[u8],
    ) -> Result<CompactSignature<Level1>, signature::Error> {
        self.sign(msg, rng).map_err(|_| signature::Error::new())
    }
}

impl signature::Keypair for CompactSigningKey<Level1> {
    type VerifyingKey = CompactPublicKey<Level1>;

    fn verifying_key(&self) -> Self::VerifyingKey {
        self.pk.clone()
    }
}

impl<L: SecurityLevel> core::fmt::Debug for CompactSigningKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("CompactSigningKey([REDACTED])")
    }
}

impl<L: SecurityLevel> Zeroize for CompactSigningKey<L> {
    fn zeroize(&mut self) {
        // Scrub all secret-derived material in the key: the secret ideal, the
        // basis-change matrix, and the curve (`Dim4SecretKey::zeroize`). As on
        // the rest of the signing side, num-bigint heap copies are scrubbed only
        // with a zeroizing allocator; the logical values are cleared here.
        self.sk.zeroize();
    }
}

impl<L: SecurityLevel> Drop for CompactSigningKey<L> {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl<L: SecurityLevel> ZeroizeOnDrop for CompactSigningKey<L> {}

/// Generate a fresh compact keypair (the smallest signatures: 108 bytes at
/// Level 1).
///
/// Returns the compact public key (for verification) and the compact signing
/// key. The caller chooses the scheme at keygen time; compact keys verify only
/// compact signatures.
pub fn generate_compact(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> (CompactPublicKey<Level1>, CompactSigningKey<Level1>) {
    let (pk_dim4, sk) = dim4_keygen(rng);
    // Build the public type via its 64-byte encoding (avoids exposing an
    // internal constructor across the crate boundary).
    let pk_bytes = encode_public_key(&pk_dim4.a_pk, pk_dim4.hint_pk_p, pk_dim4.hint_pk_q)
        .expect("invariant: canonical HD hints fit the spare bits");
    let pk = CompactPublicKey::from_bytes(&pk_bytes)
        .expect("invariant: a just-encoded compact public key parses");
    let signing_key = CompactSigningKey {
        sk,
        pk_dim4,
        pk: pk.clone(),
    };
    (pk, signing_key)
}
