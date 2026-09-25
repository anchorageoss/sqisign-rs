//! The known-answer tests of the round-3 reference (`kat/*.rsp`, 100
//! entries per level) and the helpers that reproduce them: the NIST
//! AES-256 CTR-DRBG and the kernel-path switch. Not published.

pub mod nist_drbg;

/// Choose the field kernels from the environment: `SQISIGN_FORCE_PORTABLE=1`
/// forces the portable kernels of the x86-64 backends; anything else keeps
/// the run-time detection. Returns the name of the kernels in use.
pub fn arithmetic_from_env() -> &'static str {
    if std::env::var("SQISIGN_FORCE_PORTABLE").is_ok_and(|v| v == "1") {
        sqisign_verify::fp::force_portable_arithmetic();
    }
    sqisign_verify::fp::arithmetic_backend()
}
