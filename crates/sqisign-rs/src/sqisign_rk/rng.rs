//! A deterministic SHAKE256-seeded RNG, so the public `RandSK` / `VerKey`
//! operations need no caller-supplied RNG while the internal ideal-to-isogeny
//! and prime-norm-reduction steps (which sample) still get randomness.
//!
//! Seeding from `(domain ‖ A(E_pk) ‖ rr)` makes `RandSK` deterministic in
//! `(sk, pk, rr)`; the same inputs always yield the same derived key; and
//! keeps the module `no_std`-friendly (no OS RNG, no `std_rng`).

use rand::RngCore;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// An `RngCore` backed by an unbounded SHAKE256 XOF stream.
pub(crate) struct ShakeRng<R: XofReader> {
    reader: R,
}

/// Seed a `ShakeRng` from a domain separator and context bytes.
pub(crate) fn seed_rng(domain: &[u8], context: &[&[u8]]) -> ShakeRng<impl XofReader> {
    let mut hasher = Shake256::default();
    hasher.update(domain);
    for part in context {
        hasher.update(part);
    }
    ShakeRng {
        reader: hasher.finalize_xof(),
    }
}

impl<R: XofReader> RngCore for ShakeRng<R> {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.reader.read(&mut b);
        u32::from_le_bytes(b)
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.reader.read(&mut b);
        u64::from_le_bytes(b)
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.reader.read(dest);
    }

    #[inline]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        self.reader.read(dest);
        Ok(())
    }
}
