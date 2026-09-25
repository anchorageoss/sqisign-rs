//! Chains of 4-isogenies (spec Algorithms 4.97 to 4.99). Round 3 computes
//! every `2^e`-isogeny, `e` even, as `e/2 - 1` explicit 4-isogeny steps
//! followed by `NormalizeMontgomery` on the last order-4 kernel point, which
//! yields the codomain in its canonical Montgomery model.

use super::point::xdbl_a24;
use super::{EcCurve, EcIsogEven, EcPoint};
use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;

/// Constants of a 4-isogeny for evaluating it on further points.
#[derive(Clone, Debug)]
pub struct EcKps4<L: Prime> {
    /// `(4 Z_P^2, X_P - Z_P, X_P + Z_P)`.
    pub k: [Fp2<L>; 3],
}

/// `FourIsogenyCodomain`: the 4-isogeny with kernel `<P>`, `[2] P != (0 : 1)`.
/// Returns the constants and the codomain as `(A24 : C24)`.
#[inline]
pub fn iso_xisog_4<L: FpBackend>(p: &EcPoint<L>) -> (EcKps4<L>, EcPoint<L>) {
    let t1 = p.x.sqr();
    let t2 = p.z.sqr();
    let t3 = t2.add(&t1);
    let t4 = t2.sub(&t1);
    let a24 = EcPoint::new(t3.mul(&t4), t2.sqr());
    let k2 = p.x.add(&p.z);
    let k1 = p.x.sub(&p.z);
    let t2 = t2.add(&t2);
    let k0 = t2.add(&t2);
    (EcKps4 { k: [k0, k1, k2] }, a24)
}

/// `FourIsogenyEval` on one point.
#[inline]
pub fn iso_xeval_4<L: FpBackend>(q: &EcPoint<L>, kps: &EcKps4<L>) -> EcPoint<L> {
    let t0 = q.x.add(&q.z);
    let t1 = q.x.sub(&q.z);
    let rx = t0.mul(&kps.k[1]);
    let rz = t1.mul(&kps.k[2]);
    let t0 = t0.mul(&t1).mul(&kps.k[0]);
    let t1 = rx.add(&rz).sqr();
    let rz = rx.sub(&rz).sqr();
    let rx = t0.add(&t1);
    let t0 = t0.sub(&rz);
    EcPoint::new(rx.mul(&t1), rz.mul(&t0))
}

/// Maximum depth of the strategy stack: `ceil(log2(length)) + 1` for any
/// supported length.
const STACK: usize = 16;

/// `TwoIsogenyChain`: the codomain, in canonical Montgomery form, of the
/// `2^length`-isogeny with kernel `<kernel>` on `curve` (`length` even).
/// Returns `None` if the kernel does not have the claimed order, lies above
/// `(0, 0)`, or the final normalisation fails.
///
/// The last 4-isogeny is not computed explicitly: its kernel point, taken as
/// a theta null point, yields the codomain directly (`NormalizeMontgomery`),
/// which is why no points can be pushed through this chain.
#[inline]
pub fn iso_isogeny_2chain_from<L: FpBackend>(
    curve: &EcCurve<L>,
    kernel: &EcPoint<L>,
    length: u32,
) -> Option<EcCurve<L>> {
    if length % 2 != 0 || length == 0 {
        return None;
    }
    let mut domain = curve.clone();
    domain.normalize_a24();
    let mut a24 = domain.a24.clone();

    let mut splits: [EcPoint<L>; STACK] = core::array::from_fn(|_| EcPoint::identity());
    let mut todo = [0u16; STACK];
    splits[0] = kernel.clone();
    todo[0] = length as u16;
    let mut current = 0usize;

    for j in 0..(length / 2 - 1) {
        while todo[current] != 2 {
            debug_assert!(todo[current] >= 3);
            current += 1;
            if current >= STACK {
                return None;
            }
            splits[current] = splits[current - 1].clone();
            let num_dbls = todo[current - 1] / 4 * 2 + todo[current - 1] % 2;
            todo[current] = todo[current - 1] - num_dbls;
            for _ in 0..num_dbls {
                splits[current] = xdbl_a24(&splits[current], &a24, false);
            }
        }
        if j == 0 {
            // input validation: order exactly 4, and not above (0, 0)
            if !bool::from(splits[current].is_four_torsion(&domain)) {
                return None;
            }
            let t = xdbl_a24(&splits[current], &a24, false);
            if bool::from(t.x.ct_is_zero()) {
                return None;
            }
        }
        let (kps, next_a24) = iso_xisog_4(&splits[current]);
        a24 = next_a24;
        for i in 0..current {
            splits[i] = iso_xeval_4(&splits[i], &kps);
            todo[i] -= 2;
        }
        current -= 1;
    }
    debug_assert_eq!(current, 0);
    debug_assert_eq!(todo[0], 2);
    // The last kernel point must have order 4 and not lie above (0, 0); the
    // reference only checks this at the first step, which leaves a chain of
    // length 2 unchecked.
    let last = &splits[0];
    let t = xdbl_a24(last, &a24, false);
    let cur = EcCurve::from_a24(&a24);
    if !bool::from(t.is_two_torsion(&cur)) || bool::from(t.x.ct_is_zero()) {
        return None;
    }
    let (codomain, _) = super::normalize::ec_theta_to_montgomery(last)?;
    Some(codomain)
}

/// [`iso_isogeny_2chain_from`] on an [`EcIsogEven`].
#[inline]
pub fn iso_isogeny_2chain<L: FpBackend>(phi: &EcIsogEven<L>) -> Option<EcCurve<L>> {
    iso_isogeny_2chain_from(&phi.curve, &phi.kernel, phi.length)
}

/// [`iso_isogeny_2chain_from`] that also pushes `points` through the
/// isogeny; the codomain comes in the Montgomery model the 4-isogeny
/// formulas produce, the one the pushed points are on (the compact format's
/// challenge evaluation, feature `compact`).
#[cfg(feature = "compact")]
#[inline]
pub fn iso_isogeny_2chain_eval<L: FpBackend>(
    curve: &EcCurve<L>,
    kernel: &EcPoint<L>,
    length: u32,
    points: &mut [EcPoint<L>],
) -> Option<EcCurve<L>> {
    if length % 2 != 0 || length == 0 {
        return None;
    }
    let mut domain = curve.clone();
    domain.normalize_a24();
    let mut a24 = domain.a24.clone();

    let mut splits: [EcPoint<L>; STACK] = core::array::from_fn(|_| EcPoint::identity());
    let mut todo = [0u16; STACK];
    splits[0] = kernel.clone();
    todo[0] = length as u16;
    let mut current = 0usize;

    for j in 0..(length / 2 - 1) {
        while todo[current] != 2 {
            current += 1;
            if current >= STACK {
                return None;
            }
            splits[current] = splits[current - 1].clone();
            let num_dbls = todo[current - 1] / 4 * 2 + todo[current - 1] % 2;
            todo[current] = todo[current - 1] - num_dbls;
            for _ in 0..num_dbls {
                splits[current] = xdbl_a24(&splits[current], &a24, false);
            }
        }
        if j == 0 {
            if !bool::from(splits[current].is_four_torsion(&domain)) {
                return None;
            }
            let t = xdbl_a24(&splits[current], &a24, false);
            if bool::from(t.x.ct_is_zero()) {
                return None;
            }
        }
        let (kps, next_a24) = iso_xisog_4(&splits[current]);
        a24 = next_a24;
        for i in 0..current {
            splits[i] = iso_xeval_4(&splits[i], &kps);
            todo[i] -= 2;
        }
        for pt in points.iter_mut() {
            *pt = iso_xeval_4(pt, &kps);
        }
        current -= 1;
    }
    let last = &splits[0];
    let t = xdbl_a24(last, &a24, false);
    let mut cur = EcCurve::from_a24(&a24);
    if !bool::from(t.is_two_torsion(&cur)) || bool::from(t.x.ct_is_zero()) {
        return None;
    }
    // The codomain in the model the 4-isogeny formulas produce, which is the
    // model the pushed points live on (the canonical model of
    // `ec_theta_to_montgomery` is for theta coordinates; the compact verifier
    // does not need it).
    cur.normalize();
    cur.is_a24_computed_and_normalized = false;
    Some(cur)
}
