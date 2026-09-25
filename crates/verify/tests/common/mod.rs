//! Shared test helpers: a deterministic byte source and uniform field
//! element sampling, generic over the prime.
//!
//! Inputs are drawn from SHAKE256 seeded with a per-test label so every
//! failure is reproducible without pulling in a random-number crate.

#![allow(dead_code)]

use hybrid_array::typenum::Unsigned;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use sqisign_verify::fp::{Fp, Fp2, FpBackend};

/// Iterations per property test.
pub const ITER: usize = 64;

/// A SHAKE256-based deterministic byte source.
pub struct DetRng {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl DetRng {
    pub fn new(label: &[u8]) -> Self {
        let mut hasher = Shake256::default();
        hasher.update(b"prism-fp-test-rng/");
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

    /// A field element that is uniform up to negligible bias: 16 bytes
    /// more than the encoding length are drawn and reduced modulo `p`.
    pub fn random_fp<L: FpBackend>(&mut self) -> Fp<L> {
        let mut buf = [0u8; 128];
        let n = L::FpEncodedBytes::USIZE + 16;
        self.fill(&mut buf[..n]);
        Fp::<L>::decode_reduce(&buf[..n])
    }

    pub fn random_fp2<L: FpBackend>(&mut self) -> Fp2<L> {
        Fp2 {
            re: self.random_fp::<L>(),
            im: self.random_fp::<L>(),
        }
    }
}

pub fn eq<L: FpBackend>(a: &Fp<L>, b: &Fp<L>) -> bool {
    bool::from(a.ct_equal(b))
}

pub fn eq2<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}
