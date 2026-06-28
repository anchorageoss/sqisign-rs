//! Differential/simple power-analysis countermeasures for the secret-dependent
//! elliptic-curve arithmetic in the signing path.
//!
//! Gated behind the `dpa-protect` feature (off by default). When the feature is
//! off, `maybe_randomize_basis` is a no-op and the signing computation is
//! byte-identical to the unprotected path, so the NIST KAT vectors still
//! reproduce exactly. When on, every secret-curve torsion basis that feeds a
//! Montgomery ladder is re-expressed in a fresh random projective representative
//! `(X:Z) -> (rX:rZ)` per point and per signature, so the early ladder steps
//! carry no reproducible power/EM signature across signatures.
//!
//! ## Why this is mathematically transparent
//!
//! The Montgomery x-only differential addition (`xADD`) and doubling (`xDBL`)
//! are homogeneous in each input point's projective coordinates: scaling
//! `(X:Z)` of an input by a nonzero `r` scales both output coordinates by a
//! common factor, leaving the projective (hence affine) point unchanged. The
//! 2-D ladder consumes `P`, `Q`, and the difference `P-Q`; each enters its
//! differential additions homogeneously, so independently scaling all three
//! yields a valid representative of the same `s0*P + s1*Q`. The affine output,
//! and therefore every signature byte, is unchanged — only the intermediate
//! representatives differ from one run to the next.

#[cfg(feature = "dpa-protect")]
use hybrid_array::typenum::Unsigned as _;
use rand::Rng;
use sqisign_verify::ec::EcBasis;
#[cfg(feature = "dpa-protect")]
use sqisign_verify::ec::EcPoint;
use sqisign_verify::fp::FpBackend;
#[cfg(feature = "dpa-protect")]
use sqisign_verify::fp::{Fp, Fp2};

/// A uniformly random nonzero `Fp2` element, used as a projective scaling
/// factor. Each `Fp` component is drawn by reducing fresh random bytes modulo
/// `p` (`decode_reduce`, so there is no rejection at the byte boundary); the
/// negligible-probability zero element is rejected.
#[cfg(feature = "dpa-protect")]
fn fp2_random_nonzero<L: FpBackend>(rng: &mut impl Rng) -> Fp2<L> {
    // One `Fp`-encoded length per component.
    let n = L::FpEncodedBytes::USIZE;
    loop {
        let mut bytes = alloc::vec![0u8; 2 * n];
        rng.fill_bytes(&mut bytes);
        let r = Fp2::<L> {
            re: Fp::<L>::decode_reduce(&bytes[..n]),
            im: Fp::<L>::decode_reduce(&bytes[n..]),
        };
        if !bool::from(r.ct_is_zero()) {
            return r;
        }
    }
}

/// Re-express a point in a fresh random projective representative
/// `(X:Z) -> (rX:rZ)`.
#[cfg(feature = "dpa-protect")]
fn randomize_point_projective<L: FpBackend>(p: &mut EcPoint<L>, rng: &mut impl Rng) {
    let r = fp2_random_nonzero::<L>(rng);
    p.x = p.x.mul(&r);
    p.z = p.z.mul(&r);
}

/// Randomize the projective representatives of a torsion basis `(P, Q, P-Q)`
/// before it enters a Montgomery ladder. Each point gets an independent scaling
/// factor; the ladder's homogeneity makes this transparent to the affine
/// result (see module docs).
///
/// With `dpa-protect` off this is a no-op: the basis is left untouched and no
/// randomness is consumed, so the signing path stays byte-exact against the
/// KAT vectors. The `&mut` borrows of `bas` and `rng` are still taken, so the
/// `mut` bindings at call sites remain justified in both configurations.
#[cfg(feature = "dpa-protect")]
#[inline]
pub(crate) fn maybe_randomize_basis<L: FpBackend>(bas: &mut EcBasis<L>, rng: &mut impl Rng) {
    randomize_point_projective(&mut bas.p, rng);
    randomize_point_projective(&mut bas.q, rng);
    randomize_point_projective(&mut bas.pmq, rng);
}

#[cfg(not(feature = "dpa-protect"))]
#[inline]
pub(crate) fn maybe_randomize_basis<L: FpBackend>(_bas: &mut EcBasis<L>, _rng: &mut impl Rng) {}
