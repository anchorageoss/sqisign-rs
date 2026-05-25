//! Shared test helpers: a deterministic PRNG and helpers that draw
//! random Fp and Fp2 elements for the property-based tests.
//!
//! Test inputs need to be deterministic (so failures are reproducible)
//! and cover a wide spread of field values. SHAKE256 with a per-test
//! label gives an independent, reproducible byte stream for each test
//! without dragging in any additional dependencies.

#![allow(dead_code)]

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use sqisign_verify::fp::{Fp, Fp2};
use sqisign_verify::params::{Level1, Level3, Level5, SecurityLevel};

/// A SHAKE256-based deterministic random byte source.
pub struct DetRng {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl DetRng {
    pub fn new(label: &[u8]) -> Self {
        let mut hasher = Shake256::default();
        hasher.update(b"sqisign-fp-test-rng/");
        hasher.update(label);
        Self {
            reader: hasher.finalize_xof(),
        }
    }

    pub fn fill(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }

    pub fn random_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill(&mut buf);
        u32::from_le_bytes(buf)
    }

    /// Draw a uniformly-distributed Level 1 Fp element by filling 32
    /// random bytes and reducing them modulo `p`.
    pub fn random_fp_level1(&mut self) -> Fp<Level1> {
        let mut tmp = [0u8; 32];
        debug_assert_eq!(<Level1 as SecurityLevel>::FpEncodedBytes::USIZE, 32);
        self.fill(&mut tmp);
        Fp::<Level1>::decode_reduce(&tmp)
    }

    /// Draw a Level 1 Fp2 element via two independent `random_fp`
    /// draws (one for the real part, one for the imaginary part).
    pub fn random_fp2_level1(&mut self) -> Fp2<Level1> {
        Fp2 {
            re: self.random_fp_level1(),
            im: self.random_fp_level1(),
        }
    }

    /// Draw a uniformly-distributed Level 3 Fp element by filling 48
    /// random bytes and reducing them modulo `p`.
    pub fn random_fp_level3(&mut self) -> Fp<Level3> {
        let mut tmp = [0u8; 48];
        debug_assert_eq!(<Level3 as SecurityLevel>::FpEncodedBytes::USIZE, 48);
        self.fill(&mut tmp);
        Fp::<Level3>::decode_reduce(&tmp)
    }

    /// Draw a Level 3 Fp2 element.
    pub fn random_fp2_level3(&mut self) -> Fp2<Level3> {
        Fp2 {
            re: self.random_fp_level3(),
            im: self.random_fp_level3(),
        }
    }

    /// Draw a uniformly-distributed Level 5 Fp element by filling 64
    /// random bytes and reducing them modulo `p`.
    pub fn random_fp_level5(&mut self) -> Fp<Level5> {
        let mut tmp = [0u8; 64];
        debug_assert_eq!(<Level5 as SecurityLevel>::FpEncodedBytes::USIZE, 64);
        self.fill(&mut tmp);
        Fp::<Level5>::decode_reduce(&tmp)
    }

    /// Draw a Level 5 Fp2 element.
    pub fn random_fp2_level5(&mut self) -> Fp2<Level5> {
        Fp2 {
            re: self.random_fp_level5(),
            im: self.random_fp_level5(),
        }
    }
}

// Re-export typenum::Unsigned for FpEncodedBytes::USIZE.
use typenum::Unsigned as _;

/// Number of test-loop iterations per property. Small enough that
/// `cargo test` is fast in CI; large enough to catch routine errors
/// in any reasonable random subset of the field.
pub const ITER: usize = 64;
