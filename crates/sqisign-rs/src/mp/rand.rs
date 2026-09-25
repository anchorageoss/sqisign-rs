//! Uniform sampling of integers in an interval; the byte-source trait and
//! the SHAKE256 streams live in [`sqisign_verify::rng`] and are re-exported
//! from [`crate::mp`].

use super::Ibz;
pub use sqisign_verify::rng::{DefaultDomain, Rng, ShakeRng, SplitRng};

impl<const N: usize> Ibz<N> {
    /// Uniform value in `[a, b]` (`a <= b`), using `b.bitlen + 64` random
    /// bits to make the bias negligible. `b` must be at least 64 bits
    /// below the container. `None` on randomness failure.
    pub fn rand_interval(a: &Self, b: &Self, rng: &mut impl Rng) -> Option<Self> {
        debug_assert!(*b >= *a);
        let mut bma = b.sub(a);
        debug_assert!(bma.is_positive());
        bma = bma.add(&Self::one());
        let numbits = b.bitlen + 64;
        debug_assert!(numbits < Self::MAX_BITS);
        let numbytes = ((numbits + 7) / 8) as usize;
        let numwords = (numbits / 64 + 1) as usize;
        debug_assert!(numwords >= 1 && numwords <= N);
        let padlen = numbits - (numwords as i32 - 1) * 64;
        debug_assert!(padlen >= 0 && padlen < 64);
        let mut t = Self::zero();
        let mut bytes = [0u8; 8 * super::IBZ_MAX_LIMBS];
        if !rng.fill(&mut bytes[..numbytes]) {
            return None;
        }
        for (i, chunk) in bytes[..numwords * 8].chunks(8).enumerate() {
            t.limbs[i] = u64::from_le_bytes(chunk.try_into().expect("8 bytes"));
        }
        t.limbs[numwords - 1] &= (1u64 << padlen) - 1;
        t.bitlen = numbits + 1;
        t = t.add(&bma);
        debug_assert!(t.bitsize() <= numbits);
        t = t.modulo(&bma);
        let mut r = a.add(&t);
        r.set_bound(b.bitlen);
        debug_assert!(r >= *a && r <= *b);
        Some(r)
    }

    /// Uniform value in `[-m, m]` for a small positive `m`.
    pub fn rand_interval_minm_m(m: i32, rng: &mut impl Rng) -> Option<Self> {
        debug_assert!(m > 0);
        let mm = Self::set(m as i64, 32);
        let two_m = mm.add(&mm);
        let r = Self::rand_interval(&Self::zero(), &two_m, rng)?;
        let mut r = r.sub(&mm);
        r.set_bound(32);
        Some(r)
    }
}
