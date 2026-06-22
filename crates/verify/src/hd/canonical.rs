//! Phase 5b.6 (front half) - the dim-1 canonical 4-torsion basis
//! (`basis_change/canonical_basis_dim1.py::make_canonical`).
//!
//! Each dim-1 theta structure in the gluing chain is keyed to a *canonical*
//! 4-torsion basis `(U1, U2)` with `U2[0] = -1` and a fixed Weil pairing. This
//! ports `make_canonical`: from any `E[4]` basis `(P, Q)` it computes the
//! canonical `(U1, U2)` and the `2×2` change matrix `M` over `Z/4`.
//!
//! # Convention note (why this matches the oracle despite the `weil` inverse)
//!
//! The reference uses PARI's Weil pairing; our biextension [`weil`] is its
//! inverse (5b.2). In `make_canonical` the discrete logs `a1, b1` to base `i`
//! therefore come out negated mod 4 - but the basis is built from
//! `d1 = c1·b1 = inv(a1)·b1`, and the two negations cancel (`inv(-a1)·(-b1) =
//! inv(a1)·b1`), and the pairing-fixup compares `e(U1,U2)` to `e(V1,V2)` (a
//! ratio, convention-independent). So `make_canonical` reproduces the oracle's
//! canonical basis and change matrix exactly.

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::weil;
use crate::ec::{EcCurve, JacPoint};
use crate::{Fp2, FpBackend};

/// `last_four_torsion(E)`: the 4-torsion point `(-1, √(A-2))` on
/// `y² = x³ + A x² + x`. Its sign is immaterial to `make_canonical` (it cancels
/// in `d1`), so the canonical `Fp2::sqrt` is used.
fn last_four_torsion<L: FpBackend>(a: &Fp2<L>) -> JacPoint<L> {
    let y = a.sub(&Fp2::<L>::from_small(2)).sqrt();
    JacPoint::new(Fp2::<L>::one().neg(), y, Fp2::<L>::one())
}

/// Discrete log to base `i` modulo 4: the `e ∈ {0,1,2,3}` with `iᵉ = w`.
fn dlog_base_i<L: FpBackend>(w: &Fp2<L>) -> Option<u8> {
    let i = Fp2::<L>::i_element();
    let pows = [Fp2::<L>::one(), i.clone(), Fp2::<L>::one().neg(), i.neg()];
    pows.iter()
        .position(|p| bool::from(w.ct_equal(p)))
        .map(|e| e as u8)
}

/// `[d]·P` for a small `d ∈ {0,1,2,3}` (square-and-add by hand).
fn jac_small_mul<L: FpBackend>(p: &JacPoint<L>, d: u8, curve: &EcCurve<L>) -> JacPoint<L> {
    match d & 3 {
        0 => JacPoint::identity(),
        1 => p.clone(),
        2 => jac_dbl(p, curve),
        _ => jac_add(&jac_dbl(p, curve), p, curve),
    }
}

#[inline]
fn jac_sub<L: FpBackend>(a: &JacPoint<L>, b: &JacPoint<L>, curve: &EcCurve<L>) -> JacPoint<L> {
    jac_add(a, &b.neg(), curve)
}

/// `e_4(U, V)` (biextension `weil`, = PARI's inverse) for full points `U, V`.
fn weil4<L: FpBackend>(u: &JacPoint<L>, v: &JacPoint<L>, curve: &mut EcCurve<L>) -> Fp2<L> {
    let uv = jac_add(u, &v.neg(), curve);
    weil(2, &u.to_xz(), &v.to_xz(), &uv.to_xz(), curve)
}

/// Canonical 4-torsion basis from an `E[4]` basis `(p, q)`
/// (`make_canonical(p, q, 4, preserve_pairing=True)`).
///
/// `curve` must be the Montgomery curve in `(A : 1)` form (`curve.a` the affine
/// coefficient); it is left with `a24` normalised. Returns `(U1, U2, M)` where
/// `(U1, U2)` is the canonical basis (`U2` above `x = -1`) and `M` is the `2×2`
/// change matrix over `Z/4` (`M·(U1,U2)ᵀ = (p,q)ᵀ`). Returns `None` if a Weil
/// pairing is not a 4th root of unity (degenerate input).
#[allow(clippy::type_complexity)] // (U1, U2, change-matrix) - a small, clear tuple.
pub fn make_canonical<L: FpBackend>(
    p: &JacPoint<L>,
    q: &JacPoint<L>,
    curve: &mut EcCurve<L>,
) -> Option<(JacPoint<L>, JacPoint<L>, [[u8; 2]; 2])> {
    curve.normalize_a24();
    let t2 = last_four_torsion(&curve.a);

    // A = 4, so V1 = P, V2 = Q and U1 = P, U2 = Q initially.
    let a1 = dlog_base_i(&weil4(p, &t2, curve))?;
    let b1 = dlog_base_i(&weil4(q, &t2, curve))?;

    let (u1, mut u2, mut m): (JacPoint<L>, JacPoint<L>, [[u8; 2]; 2]);
    if a1 % 2 != 0 {
        let c1 = a1 & 3; // inverse of an odd residue mod 4 is itself
        let d1 = (c1 * b1) & 3;
        // U1 = P, U2 = Q - d1·P ; M = [[1,0],[d1,1]]
        u1 = p.clone();
        u2 = jac_sub(q, &jac_small_mul(p, d1, curve), curve);
        m = [[1, 0], [d1, 1]];
    } else {
        let c1 = b1 & 3;
        let d1 = (c1 * a1) & 3;
        // U1 = Q, U2 = P - d1·Q ; M = [[d1,1],[1,0]]
        u1 = q.clone();
        u2 = jac_sub(p, &jac_small_mul(q, d1, curve), curve);
        m = [[d1, 1], [1, 0]];
    }

    // preserve_pairing: fix the sign so e(U1, U2) = e(V1, V2) = e(P, Q).
    let e4 = weil4(p, q, curve);
    if !bool::from(weil4(&u1, &u2, curve).ct_equal(&e4)) {
        u2 = u2.neg();
        m[0][1] = (4 - m[0][1]) & 3;
        m[1][1] = (4 - m[1][1]) & 3;
    }

    Some((u1, u2, m))
}
