//! Randomness: the byte-source trait the randomised routines sample from,
//! a SHAKE256 stream matching the reference's `prng_domain` construction,
//! and the two-stream split the reference's protocol code uses.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// A source of random bytes. `fill` returns `false` on failure, which the
/// callers propagate as `None`.
///
/// The reference draws from two streams at once: the domain a protocol
/// step was given (`"key"`, `"com"`, `"res"`) for the draws that define
/// the step's output, and the process-wide default domain (`"def"`) for
/// the draws inside shared routines (Miller-Rabin bases, square-root
/// searches, Qlapoty, the enumeration inside SmallestEquivalentIdeal).
/// `fill_default` is the second stream; a single-stream generator leaves
/// it at its default, which is the same stream.
pub trait Rng {
    /// Fill `out` with random bytes.
    fn fill(&mut self, out: &mut [u8]) -> bool;

    /// Fill `out` from the default domain (see the trait documentation).
    #[inline]
    fn fill_default(&mut self, out: &mut [u8]) -> bool {
        self.fill(out)
    }
}

/// A view of `R` whose primary stream is `R`'s default domain: passed to
/// the routines the reference draws from `PRNG_default_domain`.
pub struct DefaultDomain<'a, R: ?Sized>(pub &'a mut R);

impl<R: Rng + ?Sized> Rng for DefaultDomain<'_, R> {
    #[inline]
    fn fill(&mut self, out: &mut [u8]) -> bool {
        self.0.fill_default(out)
    }
    #[inline]
    fn fill_default(&mut self, out: &mut [u8]) -> bool {
        self.0.fill_default(out)
    }
}

/// Two streams: `domain` for the caller's draws, `default` for the
/// reference's default domain (see [`Rng`]).
pub struct SplitRng<D, F> {
    /// The stream of the protocol step.
    pub domain: D,
    /// The default domain.
    pub default: F,
}

impl<D: Rng, F: Rng> Rng for SplitRng<D, F> {
    #[inline]
    fn fill(&mut self, out: &mut [u8]) -> bool {
        self.domain.fill(out)
    }
    #[inline]
    fn fill_default(&mut self, out: &mut [u8]) -> bool {
        self.default.fill(out)
    }
}

/// A SHAKE256 output stream seeded as the reference seeds a PRNG domain:
/// `SHAKE256(seed || domain)`. Deterministic, so tests can reproduce the
/// reference's draws.
pub struct ShakeRng {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl ShakeRng {
    /// Derive the domain `domain` from `seed` (the reference uses a 48-byte
    /// root seed and a 3-byte domain separator).
    pub fn new(seed: &[u8], domain: &[u8]) -> Self {
        Self::from_parts(&[seed, domain])
    }

    /// `SHAKE256(parts[0] || parts[1] || ...)` as a stream. The caller
    /// keeps the concatenation unambiguous (fixed-length parts first).
    pub fn from_parts(parts: &[&[u8]]) -> Self {
        let mut h = Shake256::default();
        for p in parts {
            h.update(p);
        }
        Self {
            reader: h.finalize_xof(),
        }
    }
}

impl Rng for ShakeRng {
    #[inline]
    fn fill(&mut self, out: &mut [u8]) -> bool {
        self.reader.read(out);
        true
    }
}

impl<R: Rng + ?Sized> Rng for &mut R {
    #[inline]
    fn fill(&mut self, out: &mut [u8]) -> bool {
        (**self).fill(out)
    }
    #[inline]
    fn fill_default(&mut self, out: &mut [u8]) -> bool {
        (**self).fill_default(out)
    }
}
