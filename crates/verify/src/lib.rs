//!
//! SQIsign signature verification in pure Rust.
//!
//! This crate is `no_std`-compatible, heap-free, and independent of the
//! quaternion algebra stack. It contains all the arithmetic layers needed
//! for verification: field arithmetic (params, fp), elliptic curves (ec),
//! theta model (theta), precomputed constants (precomp), and the
//! verification protocol itself.
//!
//! # Verify a signature
//!
//! All verification goes through [`pk.verify(msg, &sig)`](Verifier::verify)
//! via the RustCrypto [`Verifier`] trait. It accepts any signature type:
//! [`Signature`], [`ExpandedSignature`], [`CompressedSignature`], or
//! [`AnySignature`](formats::AnySignature) (auto-detected from raw bytes).
//!
//! ```
//! use hex_literal::hex;
//! use sqisign_verify::{PublicKey, Signature, Verifier};
//!
//! # fn main() -> Result<(), sqisign_verify::Error> {
//! let pk_bytes = hex!(
//!     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
//!     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
//!     "0B"
//! );
//! let sig_bytes = hex!(
//!     "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
//!     "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
//!     "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
//!     "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674"
//!     "71D86EFF3447018ADB0A6551EE8322AB30010202"
//! );
//! let msg = hex!(
//!     "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556A"
//!     "C8"
//! );
//!
//! let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
//! let sig: Signature = Signature::from_bytes(&sig_bytes)?;
//! pk.verify(&msg, &sig)?;
//! # Ok(())
//! # }
//! ```
//!
//! For raw bytes where the format is unknown, parse into
//! [`AnySignature`](formats::AnySignature) first:
//!
//! ```
//! use hex_literal::hex;
//! use sqisign_verify::{formats::AnySignature, PublicKey, Verifier};
//!
//! # fn main() -> Result<(), sqisign_verify::Error> {
//! # let pk_bytes = hex!(
//! #     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
//! #     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
//! #     "0B"
//! # );
//! # let sig_bytes = hex!(
//! #     "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
//! #     "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
//! #     "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
//! #     "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674"
//! #     "71D86EFF3447018ADB0A6551EE8322AB30010202"
//! # );
//! # let msg = hex!(
//! #     "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556A"
//! #     "C8"
//! # );
//! let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
//! let sig = AnySignature::from_bytes(&sig_bytes)?;
//! pk.verify(&msg, &sig)?;
//! # Ok(())
//! # }
//! ```

#![no_std]
#![forbid(unsafe_code)]

pub mod ec;
pub mod fp;
pub mod params;
pub mod precomp;
pub mod theta;

pub mod formats;
pub mod hash;
pub mod types;
pub mod verify;

pub use formats::{CompressedSignature, ExpandedSignature};
pub use hash::hash_to_challenge;
pub use types::{PublicKey, Scalar, Signature};

pub use fp::{Fp, Fp2, FpBackend};
pub use params::{Level1, Level3, Level5, SecurityLevel};
pub use precomp::LevelPrecomp;
pub use signature::{self, SignatureEncoding, Verifier};

/// Error type for verification failures.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// The signature is cryptographically invalid.
    InvalidSignature,
    /// The input bytes could not be deserialized.
    MalformedInput,
    /// The input length does not match the expected encoding size.
    InvalidLength,
    /// An internal computation failed (e.g. hash-to-challenge buffer conversion).
    InternalError,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::InvalidSignature => f.write_str("invalid signature"),
            Error::MalformedInput => f.write_str("malformed input"),
            Error::InvalidLength => f.write_str("invalid length"),
            Error::InternalError => f.write_str("internal error"),
        }
    }
}

impl From<signature::Error> for Error {
    #[inline]
    fn from(_: signature::Error) -> Self {
        Error::InvalidSignature
    }
}

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "std")]
impl std::error::Error for Error {}
