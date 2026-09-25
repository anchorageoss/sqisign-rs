//! The compact (dimension-4) format's typed API, level I, experimental
//! (feature `compact`; the protocol is in [`crate::hd`], the parameters in
//! `docs/COMPACT_R3.md`). A compact key is not a round-3 key: its basis
//! convention is the SQIsignHD library's, and its signatures verify only
//! with a [`CompactPublicKey`].

use core::marker::PhantomData;

use hybrid_array::typenum::{Unsigned, U142, U84};
use hybrid_array::{Array, ArraySize};

use crate::hd::{hd_verify_bytes, HdLevel};
use crate::params::P324_3;
use crate::Error;

/// A level with compact-format wire sizes.
pub trait CompactLevel: HdLevel {
    /// Encoded public-key length.
    type PkLen: ArraySize;
    /// Encoded signature length.
    type SigLen: ArraySize;
}

impl CompactLevel for P324_3 {
    type PkLen = U84;
    type SigLen = U142;
}

/// A compact public key (84 bytes at level I).
#[derive(Clone, Debug)]
pub struct CompactPublicKey<L: CompactLevel> {
    bytes: Array<u8, L::PkLen>,
    _level: PhantomData<L>,
}

/// A compact signature (142 bytes at level I).
#[derive(Clone, Debug)]
pub struct CompactSignature<L: CompactLevel> {
    bytes: Array<u8, L::SigLen>,
    _level: PhantomData<L>,
}

impl<L: CompactLevel> CompactPublicKey<L> {
    /// Decode (strict: the length and a canonical curve coefficient).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::PkLen::USIZE {
            return Err(Error::InvalidLength);
        }
        crate::hd::parse_public_key::<L>(bytes).map_err(|_| Error::MalformedInput)?;
        Ok(Self {
            bytes: Array::try_from(bytes).expect("length checked"),
            _level: PhantomData,
        })
    }

    /// The encoding.
    pub fn to_bytes(&self) -> Array<u8, L::PkLen> {
        self.bytes.clone()
    }

    /// The encoding as a slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Verify a compact signature on `msg`.
    pub fn verify_compact(&self, msg: &[u8], sig: &CompactSignature<L>) -> Result<(), Error> {
        hd_verify_bytes::<L>(&sig.bytes, &self.bytes, msg).map_err(|_| Error::InvalidSignature)
    }
}

impl<L: CompactLevel> CompactSignature<L> {
    /// Decode (strict: the length, a canonical curve coefficient, the
    /// scalars' and the degree's ranges).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != L::SigLen::USIZE {
            return Err(Error::InvalidLength);
        }
        crate::hd::parse_signature::<L>(bytes).map_err(|_| Error::MalformedInput)?;
        Ok(Self {
            bytes: Array::try_from(bytes).expect("length checked"),
            _level: PhantomData,
        })
    }

    /// The encoding.
    pub fn to_bytes(&self) -> Array<u8, L::SigLen> {
        self.bytes.clone()
    }

    /// The encoding as a slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }
}

impl<L: CompactLevel> TryFrom<&[u8]> for CompactPublicKey<L> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: CompactLevel> TryFrom<&[u8]> for CompactSignature<L> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: CompactLevel> From<CompactSignature<L>> for Array<u8, L::SigLen> {
    fn from(sig: CompactSignature<L>) -> Self {
        sig.bytes
    }
}

impl<L: CompactLevel> signature::SignatureEncoding for CompactSignature<L> {
    type Repr = Array<u8, L::SigLen>;
}

impl<L: CompactLevel> signature::Verifier<CompactSignature<L>> for CompactPublicKey<L> {
    fn verify(&self, msg: &[u8], sig: &CompactSignature<L>) -> Result<(), signature::Error> {
        self.verify_compact(msg, sig)
            .map_err(|_| signature::Error::new())
    }
}
