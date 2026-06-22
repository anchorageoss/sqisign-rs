//! Verify-only inversion routing for the dimension-4 theta verifier.
//!
//! SQIsignHD **verification** runs entirely on *public* data: the verifier's
//! inputs are a signature and a public key, with no secret material anywhere
//! in the computation. Constant-time field arithmetic therefore buys no
//! security on this path - there is nothing to leak through timing. This
//! module provides an optional **variable-time** inversion path, selected by
//! the crate's default-off `vartime` Cargo feature, that replaces the shared
//! constant-time inverse with a faster binary extended-GCD inverse.
//!
//! # Why this is sound and cannot affect other paths
//!
//! The vartime inversion is reached **only** from the dimension-4 verification
//! chain (`hd::{isogeny, gluing_chain, structure, self_contained, ...}`). It
//! does not modify the shared constant-time `Fp`/`Fp2` backend in
//! `sqisign-verify`; it only *calls* the public `Fp`/`Fp2` API. So the dim-2
//! verification path and every other constant-time consumer are untouched. The
//! signing path never invokes the dim-4 verifier chain (the compact signer uses
//! only the basis-recovery and encoding helpers in `hd`, which route through the
//! shared constant-time inverse), so enabling `vartime` does not change signing,
//! which is inherently variable-time regardless.
//!
//! # What is (and isn't) variable-time
//!
//! Phase 8d profiling showed that *inversion* dominates the removable cost.
//! The Level-1 `modmul` is a special-form straight-line reduction (it exploits
//! `p = 5·2^248 - 1`) with no constant-time padding to strip away, so a
//! "vartime multiply" would be identical to the constant-time one. The
//! constant-time inverse, by contrast, is a Fermat exponentiation (~250
//! squarings); the binary extended-GCD inverse below is data-dependent but
//! measured ~3.3× faster. Only inversion is routed here; every other field
//! operation stays on the shared constant-time backend.
//!
//! # Correctness / determinism
//!
//! The multiplicative inverse of a field element is unique, so the value
//! produced here is congruent (mod `p`) to the constant-time inverse. Every
//! downstream comparison and serialization in the `hd` verifier fully reduces
//! mod `p` first (`ct_equal`, `encode`, projective cross-multiplication), so the
//! verifier's observable results are **bit-identical** whether the feature is
//! on or off. The accompanying test suite asserts this by running green in
//! both configurations.
//!
//! Only Level 1 (`p = 5·2^248 - 1`, 251 bits ⇒ a 32-byte / 4-limb canonical
//! form) has a binary-GCD path; every other security level transparently
//! falls back to the shared constant-time inverse, so correctness is
//! preserved for all levels even with the feature on.

use crate::{Fp2, FpBackend};

// Public routing surface
//
// Exactly one definition of each function is compiled, chosen by the
// `vartime` feature. With the feature OFF the bodies are literally the shared
// constant-time calls, so a default `sqisign-verify` build is byte-for-byte
// identical to the constant-time path.

/// Montgomery batch inversion (the "Montgomery trick"): invert every element
/// of `x` in place using a single field inversion plus `3·(len-1)`
/// multiplications. `t1`/`t2` are caller-provided scratch slices of the same
/// length (keeps the crate heap-free / `no_std`).
///
/// Drop-in replacement for [`crate::Fp2::batched_inv`]: with the
/// `vartime` feature off it delegates to exactly that; with it on, the single
/// inner inversion uses the variable-time inverse.
#[cfg(not(feature = "vartime"))]
#[inline]
pub fn batched_inv<L: FpBackend>(x: &mut [Fp2<L>], t1: &mut [Fp2<L>], t2: &mut [Fp2<L>]) {
    Fp2::<L>::batched_inv(x, t1, t2);
}

#[cfg(feature = "vartime")]
#[inline]
pub fn batched_inv<L: FpBackend>(x: &mut [Fp2<L>], t1: &mut [Fp2<L>], t2: &mut [Fp2<L>]) {
    let len = x.len();
    debug_assert_eq!(t1.len(), len);
    debug_assert_eq!(t2.len(), len);
    if len == 0 {
        return;
    }
    // t1[i] = x[0] * x[1] * ... * x[i]
    t1[0] = x[0].clone();
    for i in 1..len {
        t1[i] = t1[i - 1].mul(&x[i]);
    }
    // The lone inversion of the whole product - the only place the trick
    // spends an inverse, and hence the only place vartime helps.
    let inverse = vartime_fp2_inv(&t1[len - 1]);
    // t2[i] = 1 / (x[0] * ... * x[len-1-i])
    t2[0] = inverse;
    for i in 1..len {
        t2[i] = t2[i - 1].mul(&x[len - i]);
    }
    // x[0] = 1 / x[0]; x[i] = (x[0]*..*x[i-1]) * (1/(x[0]*..*x[i]))
    x[0] = t2[len - 1].clone();
    for i in 1..len {
        x[i] = t1[i - 1].mul(&t2[len - i - 1]);
    }
}

/// Single `Fp2` inversion. With the `vartime` feature off this is exactly
/// [`crate::Fp2::inv`]; with it on it uses the variable-time inverse
/// (falling back to the constant-time one for any non-Level-1 backend).
#[cfg(not(feature = "vartime"))]
#[inline]
pub fn inv<L: FpBackend>(x: &Fp2<L>) -> Fp2<L> {
    x.inv()
}

#[cfg(feature = "vartime")]
#[inline]
pub fn inv<L: FpBackend>(x: &Fp2<L>) -> Fp2<L> {
    vartime_fp2_inv(x)
}

// Variable-time implementation (compiled only with `--features vartime`)

/// Variable-time `Fp2` inverse: `(a + bi)^-1 = (a - bi) / (a² + b²)`.
///
/// Mirrors the structure of the constant-time `Fp2::inv` exactly, but the one
/// `Fp` inversion of the norm `a² + b²` uses the binary extended-GCD routine
/// below. Only Level 1 (32-byte canonical form) takes that path; every other
/// level falls back to the shared constant-time `Fp2::inv`.
#[cfg(feature = "vartime")]
#[inline]
fn vartime_fp2_inv<L: FpBackend>(x: &Fp2<L>) -> Fp2<L> {
    use typenum::Unsigned;
    if L::FpEncodedBytes::USIZE != 32 {
        return x.inv();
    }
    // norm = re² + im²  (an Fp element: the Fp2 norm form, since i² = -1).
    let norm = x.re.sqr().add(&x.im.sqr());
    let n_inv = vartime_fp_inv_l1(&norm);
    let new_re = x.re.mul(&n_inv);
    let new_im = x.im.mul(&n_inv).neg();
    // `re`/`im` are public fields of the shared `Fp2`; this constructs the
    // result without touching any private state of `sqisign-verify`.
    Fp2 {
        re: new_re,
        im: new_im,
    }
}

/// The Level-1 prime `p = 5·2^248 - 1` as four little-endian `u64` limbs.
/// This is a fixed public parameter of the Level-1 scheme (the `sqisign-verify`
/// field docs state `p = 5 * 2^248 - 1`). The unit test
/// `hardcoded_modulus_matches_backend` asserts this equals the value *derived*
/// from the backend at runtime, so the constant cannot silently drift from the
/// real modulus. Keeping it a constant avoids an `encode` per inversion.
#[cfg(feature = "vartime")]
const P_L1: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0x04FF_FFFF_FFFF_FFFF,
];

/// Variable-time `Fp` inverse for Level 1, via canonical-bytes ↔ 4×u64 limbs
/// and the binary extended-GCD routine. The caller guarantees a 32-byte
/// canonical form (Level 1 only).
///
/// The Montgomery ↔ canonical conversions (`encode`/`decode`) are unavoidable:
/// the shared `sqisign-verify` backend exposes limbs only through them, and
/// this module deliberately does not modify that crate. They bound how much
/// faster than the constant-time Fermat inverse this path can be.
#[cfg(feature = "vartime")]
#[inline]
fn vartime_fp_inv_l1<L: FpBackend>(a: &crate::Fp<L>) -> crate::Fp<L> {
    let limbs = bytes_to_limbs4(a.encode().as_ref());
    let inv_limbs = inv_mod_p(&limbs, &P_L1);
    let out_bytes = limbs4_to_bytes(&inv_limbs);
    // The binary-GCD result is reduced into [0, p), so decode always succeeds.
    crate::Fp::<L>::decode(&out_bytes).expect("vartime Fp inverse is canonical (< p)")
}

/// Variable-time modular inverse modulo an odd prime `p`, on 4×u64
/// little-endian operands, via the binary extended-GCD (Stein) algorithm.
/// Returns `0` for input `0`, matching the shared backend's `inv` convention.
///
/// **Not constant-time:** the iteration count and branch pattern depend on the
/// operand. Used only on public verification data.
#[cfg(feature = "vartime")]
fn inv_mod_p(a: &[u64; 4], p: &[u64; 4]) -> [u64; 4] {
    if is_zero4(a) {
        return [0; 4];
    }
    let mut u = *a;
    let mut v = *p;
    let mut b = [1u64, 0, 0, 0];
    let mut c = [0u64; 4];
    while !is_one4(&u) && !is_one4(&v) {
        while is_even4(&u) {
            u = shr1_4(&u);
            b = if is_even4(&b) {
                shr1_4(&b)
            } else {
                add_p_shr1(&b, p)
            };
        }
        while is_even4(&v) {
            v = shr1_4(&v);
            c = if is_even4(&c) {
                shr1_4(&c)
            } else {
                add_p_shr1(&c, p)
            };
        }
        if ge4(&u, &v) {
            u = sub4(&u, &v);
            b = submod(&b, &c, p);
        } else {
            v = sub4(&v, &u);
            c = submod(&c, &b, p);
        }
    }
    if is_one4(&u) {
        b
    } else {
        c
    }
}

// 4-limb (256-bit) integer helpers for the binary-GCD inverse

#[cfg(feature = "vartime")]
#[inline]
fn is_zero4(a: &[u64; 4]) -> bool {
    (a[0] | a[1] | a[2] | a[3]) == 0
}

#[cfg(feature = "vartime")]
#[inline]
fn is_one4(a: &[u64; 4]) -> bool {
    a[0] == 1 && (a[1] | a[2] | a[3]) == 0
}

#[cfg(feature = "vartime")]
#[inline]
fn is_even4(a: &[u64; 4]) -> bool {
    a[0] & 1 == 0
}

/// `a >> 1` (logical, 256-bit).
#[cfg(feature = "vartime")]
#[inline]
fn shr1_4(a: &[u64; 4]) -> [u64; 4] {
    [
        (a[0] >> 1) | (a[1] << 63),
        (a[1] >> 1) | (a[2] << 63),
        (a[2] >> 1) | (a[3] << 63),
        a[3] >> 1,
    ]
}

/// `a >= b` over 256-bit unsigned integers.
#[cfg(feature = "vartime")]
#[inline]
fn ge4(a: &[u64; 4], b: &[u64; 4]) -> bool {
    for i in (0..4).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

/// `a - b` (256-bit, wrapping); callers only subtract when `a >= b`.
#[cfg(feature = "vartime")]
#[inline]
fn sub4(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut r = [0u64; 4];
    let mut borrow = 0u128;
    for i in 0..4 {
        let t = (a[i] as u128)
            .wrapping_sub(b[i] as u128)
            .wrapping_sub(borrow);
        r[i] = t as u64;
        borrow = (t >> 127) & 1;
    }
    r
}

/// `(a + p) >> 1`, used to halve an odd `Fp` accumulator: `a` is odd and
/// `< p`, so `a + p` is even and `< 2p`, and `(a + p)/2 < p`. The transient
/// carry into bit 256 is preserved via a fifth limb before the shift.
#[cfg(feature = "vartime")]
#[inline]
fn add_p_shr1(a: &[u64; 4], p: &[u64; 4]) -> [u64; 4] {
    let mut r = [0u64; 5];
    let mut c = 0u128;
    for i in 0..4 {
        let t = (a[i] as u128) + (p[i] as u128) + c;
        r[i] = t as u64;
        c = t >> 64;
    }
    r[4] = c as u64;
    [
        (r[0] >> 1) | (r[1] << 63),
        (r[1] >> 1) | (r[2] << 63),
        (r[2] >> 1) | (r[3] << 63),
        (r[3] >> 1) | (r[4] << 63),
    ]
}

/// `(a - b) mod p` with `a, b ∈ [0, p)`.
#[cfg(feature = "vartime")]
#[inline]
fn submod(a: &[u64; 4], b: &[u64; 4], p: &[u64; 4]) -> [u64; 4] {
    if ge4(a, b) {
        sub4(a, b)
    } else {
        let d = sub4(b, a);
        sub4(p, &d)
    }
}

/// Pack 32 little-endian bytes into four `u64` limbs (LSB-first). The caller
/// guarantees `b.len() >= 32`.
#[cfg(feature = "vartime")]
#[inline]
fn bytes_to_limbs4(b: &[u8]) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for (i, limb) in limbs.iter_mut().enumerate() {
        let mut v = 0u64;
        for j in 0..8 {
            v |= (b[i * 8 + j] as u64) << (8 * j);
        }
        *limb = v;
    }
    limbs
}

/// Unpack four `u64` limbs into 32 little-endian bytes (LSB-first).
#[cfg(feature = "vartime")]
#[inline]
fn limbs4_to_bytes(limbs: &[u64; 4]) -> [u8; 32] {
    let mut out = [0u8; 32];
    for (i, &limb) in limbs.iter().enumerate() {
        for j in 0..8 {
            out[i * 8 + j] = (limb >> (8 * j)) as u8;
        }
    }
    out
}

#[cfg(all(test, feature = "vartime"))]
mod tests {
    use super::*;
    use crate::{Fp2, Level1};

    /// Derive `p` from the backend (`-1 ≡ p - 1`, then `+1`) and confirm the
    /// hardcoded [`P_L1`] constant matches it exactly - guards against the
    /// constant drifting from the actual Level-1 modulus.
    #[test]
    fn hardcoded_modulus_matches_backend() {
        let neg_one = crate::Fp::<Level1>::zero().sub(&crate::Fp::<Level1>::one());
        let mut p = bytes_to_limbs4(neg_one.encode().as_ref());
        let mut carry = 1u64;
        for limb in p.iter_mut() {
            let (s, c) = limb.overflowing_add(carry);
            *limb = s;
            carry = c as u64;
        }
        assert_eq!(carry, 0, "p must fit in 4 limbs");
        assert_eq!(
            P_L1, p,
            "hardcoded P_L1 disagrees with backend-derived modulus"
        );
    }

    /// The vartime inverse agrees with the constant-time inverse for a spread
    /// of `Fp2` values (the values are public, so generating them simply).
    #[test]
    fn vartime_matches_constant_time_inv() {
        let mut x = Fp2::<Level1>::from_small(3);
        for _ in 0..200 {
            x = x.mul(&Fp2::from_small(7)).add(&Fp2::one());
            // Mix the imaginary part too.
            let xi = x.mul(&Fp2::<Level1>::from_small(11));
            let probe = x.add(&xi.mul(&Fp2::i_element()));
            if bool::from(probe.ct_is_zero()) {
                continue;
            }
            let ct = probe.inv();
            let vt = vartime_fp2_inv(&probe);
            assert!(
                bool::from(ct.ct_equal(&vt)),
                "vartime inverse disagrees with constant-time inverse"
            );
            // And it is a genuine inverse: probe * inv == 1.
            assert!(bool::from(probe.mul(&vt).ct_equal(&Fp2::one())));
        }
    }

    /// Batch inversion via the router matches per-element vartime inverses.
    #[test]
    fn batched_inv_matches_elementwise() {
        let mut elems: [Fp2<Level1>; 8] = core::array::from_fn(|_| Fp2::one());
        let mut acc = Fp2::<Level1>::from_small(2);
        for e in elems.iter_mut() {
            acc = acc.mul(&Fp2::from_small(5)).add(&Fp2::one());
            *e = acc.clone();
        }
        let expected: [Fp2<Level1>; 8] = core::array::from_fn(|i| vartime_fp2_inv(&elems[i]));
        let mut got = elems.clone();
        let mut t1: [Fp2<Level1>; 8] = core::array::from_fn(|_| Fp2::one());
        let mut t2: [Fp2<Level1>; 8] = core::array::from_fn(|_| Fp2::one());
        batched_inv(&mut got, &mut t1, &mut t2);
        for i in 0..8 {
            assert!(bool::from(got[i].ct_equal(&expected[i])));
        }
    }

    /// inv(inv(x)) == x.
    #[test]
    fn vartime_inv_roundtrip() {
        let mut x = Fp2::<Level1>::from_small(123_456_789);
        for _ in 0..50 {
            x = x.mul(&Fp2::from_small(31)).add(&Fp2::from_small(17));
            let r = vartime_fp2_inv(&vartime_fp2_inv(&x));
            assert!(bool::from(x.ct_equal(&r)));
        }
    }
}
