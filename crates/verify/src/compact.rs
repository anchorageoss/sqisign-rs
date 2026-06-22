//! The compact signature scheme's verification-side types.
//!
//! [`CompactSignature`] and [`CompactPublicKey`] are the public types of the
//! 108-byte-signature ("compact") scheme. They are deliberately **distinct**
//! from the dim-2 [`crate::Signature`] / [`crate::PublicKey`]: the two schemes
//! use different torsion-basis conventions and their keys are not
//! interchangeable.
//!
//! # Cross-type verification rules
//!
//! * a [`CompactPublicKey`] verifies a [`CompactSignature`] (its native format)
//!   and the compact arm of [`crate::AnySignature`];
//! * a [`CompactPublicKey`] does **not** verify dim-2
//!   standard/compressed/expanded signatures;
//! * the dim-2 [`crate::PublicKey`] does **not** verify a [`CompactSignature`].
//!
//! [`crate::AnySignature`] autodetects the wire *format* from its length, but
//! the caller must hold the public key of the matching *scheme*.
//!
//! # Levels
//!
//! The compact scheme is implemented at Level 1 only. These types are generic
//! over the security level for API symmetry with the dim-2 keys and to slot
//! into [`crate::AnySignature`]`<L>`, but the functional impls
//! (`SignatureEncoding`, `Verifier`, (de)serialization) exist for [`Level1`].

use core::marker::PhantomData;

use hybrid_array::typenum::{U108, U64};
use hybrid_array::Array;

use crate::fp::Fp2;
use crate::hd::{
    encode_public_key, encode_signature, hd_verify_l1_parsed, parse_public_key, parse_signature,
    ParsedPublicKey, ParsedSignature,
};
use crate::params::{Level1, SecurityLevel};
use crate::Error;

/// A compact signature - the 108-byte (Level 1) wire format.
///
/// Verify it with a [`CompactPublicKey`] (see the [module docs](self) for the
/// cross-type rules).
#[derive(Clone, Debug)]
pub struct CompactSignature<L: SecurityLevel = Level1> {
    pub(crate) inner: ParsedSignature,
    pub(crate) _marker: PhantomData<L>,
}

impl<L: SecurityLevel> CompactSignature<L> {
    /// Wrap an already-parsed signature (crate-internal; used by
    /// `AnySignature::from_bytes` and the signing-key wrapper).
    pub(crate) fn from_parsed(inner: ParsedSignature) -> Self {
        Self {
            inner,
            _marker: PhantomData,
        }
    }
}

impl CompactSignature<Level1> {
    /// The 108-byte wire encoding (commitment basis hints packed into `A_com`).
    pub fn to_bytes(&self) -> Array<u8, U108> {
        let bytes = encode_signature(
            &self.inner.a_com,
            self.inner.a,
            self.inner.b,
            self.inner.c_or_d,
            &self.inner.q,
            self.inner.hint_com_p,
            self.inner.hint_com_q,
        )
        .expect("invariant: a parsed compact signature always re-encodes");
        let mut out = Array::<u8, U108>::default();
        out.copy_from_slice(&bytes);
        out
    }

    /// Parse a compact signature from its 108 wire bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let inner = parse_signature(bytes).map_err(|_| Error::MalformedInput)?;
        Ok(Self::from_parsed(inner))
    }
}

impl TryFrom<&[u8]> for CompactSignature<Level1> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl From<CompactSignature<Level1>> for Array<u8, U108> {
    fn from(sig: CompactSignature<Level1>) -> Self {
        sig.to_bytes()
    }
}

impl signature::SignatureEncoding for CompactSignature<Level1> {
    type Repr = Array<u8, U108>;
}

impl<L: SecurityLevel> core::fmt::Display for CompactSignature<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Debug::fmt(&self.inner, f)
    }
}

/// A compact public key - verifies [`CompactSignature`]s.
///
/// Distinct from the dim-2 [`crate::PublicKey`]; see the [module docs](self).
#[derive(Clone, Debug)]
pub struct CompactPublicKey<L: SecurityLevel = Level1> {
    /// Commitment/public curve Montgomery coefficient `A_pk` (Level 1 field).
    pub(crate) a_pk: Fp2<Level1>,
    /// Canonical `2^f`-torsion basis hints.
    pub(crate) hint_pk_p: u32,
    pub(crate) hint_pk_q: u32,
    pub(crate) _marker: PhantomData<L>,
}

impl<L: SecurityLevel> CompactPublicKey<L> {
    /// Construct from a curve coefficient and its canonical basis hints
    /// (crate-internal; the signing-key wrapper builds public keys via the
    /// 64-byte encoding instead).
    pub(crate) fn from_parts(a_pk: Fp2<Level1>, hint_pk_p: u32, hint_pk_q: u32) -> Self {
        Self {
            a_pk,
            hint_pk_p,
            hint_pk_q,
            _marker: PhantomData,
        }
    }
}

impl CompactPublicKey<Level1> {
    /// The 64-byte wire encoding (basis hints packed into `A_pk`).
    pub fn to_bytes(&self) -> Array<u8, U64> {
        let bytes = encode_public_key(&self.a_pk, self.hint_pk_p, self.hint_pk_q)
            .expect("invariant: a valid compact public key always re-encodes");
        let mut out = Array::<u8, U64>::default();
        out.copy_from_slice(&bytes);
        out
    }

    /// Parse a compact public key from its 64 wire bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let p = parse_public_key(bytes).map_err(|_| Error::MalformedInput)?;
        Ok(Self::from_parts(p.a_pk, p.hint_pk_p, p.hint_pk_q))
    }
}

impl TryFrom<&[u8]> for CompactPublicKey<Level1> {
    type Error = Error;
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl signature::Verifier<CompactSignature<Level1>> for CompactPublicKey<Level1> {
    fn verify(&self, msg: &[u8], sig: &CompactSignature<Level1>) -> Result<(), signature::Error> {
        let pk = ParsedPublicKey {
            a_pk: self.a_pk.clone(),
            hint_pk_p: self.hint_pk_p,
            hint_pk_q: self.hint_pk_q,
        };
        hd_verify_l1_parsed(&sig.inner, &pk, msg).map_err(|_| signature::Error::new())
    }
}

impl signature::Verifier<crate::formats::AnySignature<Level1>> for CompactPublicKey<Level1> {
    fn verify(
        &self,
        msg: &[u8],
        sig: &crate::formats::AnySignature<Level1>,
    ) -> Result<(), signature::Error> {
        match sig {
            // The compact arm verifies natively.
            crate::formats::AnySignature::Compact(s) => {
                <Self as signature::Verifier<CompactSignature<Level1>>>::verify(self, msg, s)
            }
            // A compact key does not verify dim-2 signatures (wrong scheme).
            crate::formats::AnySignature::Standard(_)
            | crate::formats::AnySignature::Expanded(_)
            | crate::formats::AnySignature::Compressed(_) => Err(signature::Error::new()),
        }
    }
}
