//! SQIsign round-3 signature verification in pure Rust.
//!
//! SQIsign as submitted to the third round of NIST's additional-signatures
//! process (`the-sqisign` at tag `nist-v3`): the primes `3 * 2^324 - 1`,
//! `27 * 2^500 - 1` and `17 * 2^664 - 1`, 83 / 129 / 169-byte public keys,
//! 200 / 306 / 406-byte signatures. The 300 known-answer vectors of the
//! reference are reproduced byte for byte (`crates/kat`).
//!
//! This crate is `no_std`, allocates nothing and links no quaternion code:
//! it is the crate a verifier depends on. It holds the parameter sets
//! ([`params`]), the field layer ([`fp`]: on x86-64 inline assembly in the
//! reference's own `mulx`/`adcx`/`adox` schedule, selected at run time from
//! CPUID with a portable fallback of the same layout; elsewhere generated
//! radix backends), the Kummer-line curve layer ([`ec`], with a bounded
//! entangled-basis search that fails closed on crafted curves and a
//! fixed-base two-adic discrete logarithm), the `(2^n, 2^n)`-isogeny chains
//! in the theta model ([`theta`]), the precomputed constants ([`precomp`])
//! and the protocol ([`sqisign`]: verification, the encodings, a prepared
//! key that does the per-key work once, a batch verifier), plus the
//! compressed format of [`compressed`] (three matrix entries, the fourth
//! recovered from two Weil pairings; not in the specification). [`types`] is the
//! typed surface: [`PublicKey`], [`Signature`], the [`Level1`] /
//! [`Level3`] / [`Level5`] markers and the RustCrypto [`Verifier`] trait.
//!
//! # Verify a signature
//!
//! ```
//! use sqisign_verify::{Level1, PublicKey, Signature, Verifier};
//!
//! # fn main() -> Result<(), sqisign_verify::Error> {
//! # let pk_bytes = include_bytes!("../fuzz/fuzz_targets/l1_pk.bin");
//! # let kat = include_str!("../../kat/kat/PQCsignKAT_270_SQIsign_p324_3.rsp");
//! # let field = |k: &str| kat.lines().find_map(|l| l.strip_prefix(k)).unwrap().trim().to_string();
//! # let hex = |s: &str| (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect::<Vec<u8>>();
//! # let msg = hex(&field("msg = "));
//! # let sm = hex(&field("sm = "));
//! # let sig_bytes = &sm[..200];
//! let pk: PublicKey<Level1> = PublicKey::from_bytes(pk_bytes)?;
//! let sig: Signature<Level1> = Signature::from_bytes(sig_bytes)?;
//! pk.verify(&msg, &sig)?;
//! # Ok(())
//! # }
//! ```
//!
//! Round 2 is not supported: round-2 keys and signatures are not accepted
//! and there is no conversion (0.5.0 removed the round-2 implementation).
#![no_std]
// No unsafe code, except on x86-64: the field backends' `asm!` blocks and
// the `cpuid` that decides whether to run them (`fp::dispatch`), the only
// modules allowed to opt out.
#![cfg_attr(not(target_arch = "x86_64"), forbid(unsafe_code))]
#![cfg_attr(target_arch = "x86_64", deny(unsafe_code))]
#![warn(missing_docs)]

#[cfg(feature = "compact")]
extern crate alloc;

#[cfg(feature = "compact")]
pub mod compact;
pub mod compressed;
pub mod ec;
pub mod fp;
#[cfg(feature = "compact")]
pub mod hd;
pub mod params;
pub mod precomp;
pub mod rng;
pub mod sqisign;
pub mod theta;
pub mod types;

#[cfg(feature = "compact")]
pub use compact::{CompactLevel, CompactPublicKey, CompactSignature};
pub use compressed::MAX_COMPRESSED_BYTES;
pub use fp::{Fp, Fp2, FpBackend};
pub use params::{Prime, P324_3, P500_27, P664_17};
pub use precomp::PrimePrecomp;
pub use signature::{self, SignatureEncoding, Verifier};
pub use sqisign::{
    hash_to_challenge, verify_batch, BatchItem, PreparedPublicKey, VerifyParams,
    MAX_PUBLICKEY_BYTES, MAX_SIGNATURE_BYTES,
};
pub use types::{
    verify_bytes_prepared, CompressedSignature, Level, Level1, Level3, Level5, PublicKey, Signature,
};

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
