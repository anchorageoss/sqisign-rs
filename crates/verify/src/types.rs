//! The typed public surface: security levels, keys and signatures with
//! their byte encodings, and the RustCrypto [`signature`] traits.
//!
//! The protocol itself is [`crate::sqisign`], generic over the prime and
//! parameterised at run time by [`VerifyParams`]; this module fixes the
//! parameters per level at compile time so that a key or a signature
//! carries its level in its type and the encodings have fixed lengths.

use crate::compressed::{self, MAX_COMPRESSED_BYTES};
use crate::fp::FpBackend;
use crate::params::{Prime, P324_3, P500_27, P664_17};
use crate::precomp::PrimePrecomp;
use crate::sqisign::{
    self, PreparedPublicKey, VerifyParams, MAX_PUBLICKEY_BYTES, MAX_SIGNATURE_BYTES,
};
use crate::Error;
use hybrid_array::typenum::{Unsigned, U129, U169, U176, U200, U269, U306, U353, U406, U83};
use hybrid_array::{Array, ArraySize};

/// A NIST security level of SQIsign round 3: the prime with the level's
/// constants and encoded sizes fixed at compile time.
pub trait Level: FpBackend + PrimePrecomp + Prime + 'static {
    /// The verifier's constants (`sqisign::level1` and siblings).
    const PARAMS: VerifyParams;
    /// A short name for logs and benchmarks.
    const NAME: &'static str;
    /// Encoded public-key length.
    type PkLen: ArraySize;
    /// Encoded signature length.
    type SigLen: ArraySize;
    /// Encoded secret-key length (the signing crate's encoding).
    const SK_BYTES: usize;
    /// Encoded compressed-signature length.
    type CompressedLen: ArraySize;
}

impl Level for P324_3 {
    const PARAMS: VerifyParams = sqisign::level1();
    const NAME: &'static str = "level I (p324_3)";
    type PkLen = U83;
    type SigLen = U200;
    const SK_BYTES: usize = 270;
    type CompressedLen = U176;
}

impl Level for P500_27 {
    const PARAMS: VerifyParams = sqisign::level3();
    const NAME: &'static str = "level III (p500_27)";
    type PkLen = U129;
    type SigLen = U306;
    const SK_BYTES: usize = 417;
    type CompressedLen = U269;
}

impl Level for P664_17 {
    const PARAMS: VerifyParams = sqisign::level5();
    const NAME: &'static str = "level V (p664_17)";
    type PkLen = U169;
    type SigLen = U406;
    const SK_BYTES: usize = 549;
    type CompressedLen = U353;
}

/// NIST level I: `p = 3 * 2^324 - 1`, 83-byte keys, 200-byte signatures.
pub type Level1 = P324_3;
/// NIST level III: `p = 27 * 2^500 - 1`, 129-byte keys, 306-byte signatures.
pub type Level3 = P500_27;
/// NIST level V: `p = 17 * 2^664 - 1`, 169-byte keys, 406-byte signatures.
pub type Level5 = P664_17;

/// A public key: the curve `E_pk` and the hint of its canonical basis.
pub struct PublicKey<L: Level> {
    inner: sqisign::PublicKey<L>,
}

impl<L: Level> Clone for PublicKey<L> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<L: Level> core::fmt::Debug for PublicKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "PublicKey<{}>(", L::NAME)?;
        for b in self.to_bytes().iter() {
            write!(f, "{b:02x}")?;
        }
        f.write_str(")")
    }
}

impl<L: Level> PublicKey<L> {
    /// Wrap a protocol-level key.
    pub fn from_inner(inner: sqisign::PublicKey<L>) -> Self {
        Self { inner }
    }

    /// The protocol-level key.
    pub fn inner(&self) -> &sqisign::PublicKey<L> {
        &self.inner
    }

    /// Decode `A_pk || hint`; the coefficient must be canonical.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::PkLen::USIZE {
            return Err(Error::InvalidLength);
        }
        sqisign::public_key_from_bytes::<L>(&L::PARAMS, bytes)
            .map(|inner| Self { inner })
            .ok_or(Error::MalformedInput)
    }

    /// Encode as `A_pk || hint`.
    pub fn to_bytes(&self) -> Array<u8, L::PkLen> {
        let mut buf = [0u8; MAX_PUBLICKEY_BYTES];
        let n = sqisign::public_key_to_bytes(&self.inner, &mut buf);
        Array::try_from(&buf[..n]).expect("invariant: the level's encoded length")
    }

    /// Verify `sig` on `msg`.
    pub fn verify(&self, msg: &[u8], sig: &Signature<L>) -> Result<(), Error> {
        if sqisign::verify(&L::PARAMS, &self.inner, &sig.inner, msg) {
            Ok(())
        } else {
            Err(Error::InvalidSignature)
        }
    }

    /// The per-key work done once (the curve normalised, its canonical
    /// basis rebuilt): [`PreparedPublicKey::new`].
    pub fn prepare(&self) -> Result<PreparedPublicKey<L>, Error> {
        PreparedPublicKey::new(&L::PARAMS, &self.inner).ok_or(Error::MalformedInput)
    }
}

impl<L: Level> TryFrom<&[u8]> for PublicKey<L> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Error> {
        Self::from_bytes(bytes)
    }
}

/// A signature in the round-3 wire format: `A_aux`, the basis-change
/// matrix, the challenge and two hints.
pub struct Signature<L: Level> {
    inner: sqisign::Signature<L>,
}

impl<L: Level> Clone for Signature<L> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<L: Level> core::fmt::Debug for Signature<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Signature<{}>(", L::NAME)?;
        for b in self.to_bytes().iter() {
            write!(f, "{b:02x}")?;
        }
        f.write_str(")")
    }
}

impl<L: Level> Signature<L> {
    /// Wrap a protocol-level signature.
    pub fn from_inner(inner: sqisign::Signature<L>) -> Self {
        Self { inner }
    }

    /// The protocol-level signature.
    pub fn inner(&self) -> &sqisign::Signature<L> {
        &self.inner
    }

    /// Decode the wire format; `A_aux` must be canonical. The matrix is
    /// checked for canonical form by the verifier, not here.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::SigLen::USIZE {
            return Err(Error::InvalidLength);
        }
        sqisign::signature_from_bytes::<L>(&L::PARAMS, bytes)
            .map(|inner| Self { inner })
            .ok_or(Error::MalformedInput)
    }

    /// Encode to the wire format.
    pub fn to_bytes(&self) -> Array<u8, L::SigLen> {
        let mut buf = [0u8; MAX_SIGNATURE_BYTES];
        let n = sqisign::signature_to_bytes(&L::PARAMS, &self.inner, &mut buf)
            .expect("invariant: the buffer holds every level's signature");
        Array::try_from(&buf[..n]).expect("invariant: the level's encoded length")
    }
}

impl<L: Level> TryFrom<&[u8]> for Signature<L> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: Level> From<Signature<L>> for Array<u8, L::SigLen> {
    fn from(sig: Signature<L>) -> Self {
        sig.to_bytes()
    }
}

impl<L: Level> signature::SignatureEncoding for Signature<L> {
    type Repr = Array<u8, L::SigLen>;
}

impl<L: Level> signature::Verifier<Signature<L>> for PublicKey<L> {
    fn verify(&self, msg: &[u8], sig: &Signature<L>) -> Result<(), signature::Error> {
        PublicKey::verify(self, msg, sig).map_err(|_| signature::Error::new())
    }
}

/// Verify an encoded signature under a prepared key, both from bytes.
pub fn verify_bytes_prepared<L: Level>(
    key: &PreparedPublicKey<L>,
    msg: &[u8],
    signature: &[u8],
) -> Result<(), Error> {
    let sig = Signature::<L>::from_bytes(signature)?;
    if sqisign::verify_prepared(&L::PARAMS, key, &sig.inner, msg) {
        Ok(())
    } else {
        Err(Error::InvalidSignature)
    }
}

/// A compressed signature (three matrix entries, the fourth recovered from
/// two pairings; [`crate::compressed`]), 176 / 269 / 353 bytes at levels
/// I / III / V.
pub struct CompressedSignature<L: Level> {
    inner: compressed::CompressedSignature<L>,
}

impl<L: Level> Clone for CompressedSignature<L> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<L: Level> core::fmt::Debug for CompressedSignature<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "CompressedSignature<{}>(", L::NAME)?;
        for b in self.to_bytes().iter() {
            write!(f, "{b:02x}")?;
        }
        f.write_str(")")
    }
}

impl<L: Level> CompressedSignature<L> {
    /// The protocol-level value.
    pub fn inner(&self) -> &compressed::CompressedSignature<L> {
        &self.inner
    }

    /// Decode the wire format (strict).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::CompressedLen::USIZE {
            return Err(Error::InvalidLength);
        }
        compressed::compressed_from_bytes::<L>(&L::PARAMS, bytes)
            .map(|inner| Self { inner })
            .ok_or(Error::MalformedInput)
    }

    /// Encode to the wire format.
    pub fn to_bytes(&self) -> Array<u8, L::CompressedLen> {
        let mut buf = [0u8; MAX_COMPRESSED_BYTES];
        let n = compressed::compressed_to_bytes(&L::PARAMS, &self.inner, &mut buf)
            .expect("invariant: the buffer holds every level's compressed signature");
        Array::try_from(&buf[..n]).expect("invariant: the level's encoded length")
    }

    /// The standard signature this one stands for, recovered under `pk`
    /// (equal to the signer's up to the free top bit of each entry).
    pub fn decompress(&self, pk: &PublicKey<L>) -> Result<Signature<L>, Error> {
        compressed::decompress(&L::PARAMS, &pk.inner, &self.inner)
            .map(Signature::from_inner)
            .ok_or(Error::InvalidSignature)
    }
}

impl<L: Level> Signature<L> {
    /// The compressed form: three matrix entries and four bits.
    pub fn compress(&self) -> CompressedSignature<L> {
        CompressedSignature {
            inner: compressed::compress(&L::PARAMS, &self.inner)
                .expect("invariant: a signature's matrix has an odd determinant"),
        }
    }
}

impl<L: Level> PublicKey<L> {
    /// Verify a compressed signature on `msg`.
    pub fn verify_compressed(&self, msg: &[u8], sig: &CompressedSignature<L>) -> Result<(), Error> {
        if compressed::verify_compressed(&L::PARAMS, &self.inner, &sig.inner, msg) {
            Ok(())
        } else {
            Err(Error::InvalidSignature)
        }
    }
}

impl<L: Level> TryFrom<&[u8]> for CompressedSignature<L> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: Level> From<CompressedSignature<L>> for Array<u8, L::CompressedLen> {
    fn from(sig: CompressedSignature<L>) -> Self {
        sig.to_bytes()
    }
}

impl<L: Level> signature::SignatureEncoding for CompressedSignature<L> {
    type Repr = Array<u8, L::CompressedLen>;
}

impl<L: Level> signature::Verifier<CompressedSignature<L>> for PublicKey<L> {
    fn verify(&self, msg: &[u8], sig: &CompressedSignature<L>) -> Result<(), signature::Error> {
        self.verify_compressed(msg, sig)
            .map_err(|_| signature::Error::new())
    }
}
