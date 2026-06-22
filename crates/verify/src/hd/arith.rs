//! Coordinate-level theta operations: the Hadamard transform, pointwise
//! squaring, and the 2-torsion translation action. All operate on a bare
//! `[Fp2<L>; 16]` so they can be reused by points, structures, and (later)
//! isogenies without committing to a wrapper type.

use crate::{Fp2, FpBackend};

use crate::hd::point::THETA_DIM4_N as N;

/// The (unnormalised) dimension-4 Hadamard transform of a level-2 theta point.
///
/// `H[χ] = Σⱼ (-1)^(popcount(χ & j) mod 2) · P[j]`, indexed by
/// `k = i₀ + 2·i₁ + 4·i₂ + 8·i₃`.
///
/// Implemented as the in-place radix-2 fast Walsh-Hadamard transform with
/// strides 1, 2, 4, 8 (the butterfly `(x, y) ↦ (x+y, x-y)`), which is exactly
/// the `hadamard16` recursion in the sage reference. The transform is its own
/// inverse up to the scalar `16`; since theta points are projective that
/// factor is irrelevant, so we never divide it out.
#[inline]
pub fn hadamard<L: FpBackend>(p: &[Fp2<L>; N]) -> [Fp2<L>; N] {
    let mut a = p.clone();
    for stride in [1usize, 2, 4, 8] {
        let mut block = 0;
        while block < N {
            for j in block..block + stride {
                let u = a[j].clone();
                let v = a[j + stride].clone();
                a[j] = u.add(&v);
                a[j + stride] = u.sub(&v);
            }
            block += 2 * stride;
        }
    }
    a
}

/// Coordinate-wise squaring: `S[k] = P[k]²`.
#[inline]
pub fn pointwise_square<L: FpBackend>(p: &[Fp2<L>; N]) -> [Fp2<L>; N] {
    core::array::from_fn(|k| p[k].sqr())
}

/// Square every coordinate, then apply the Hadamard transform:
/// `to_squared_theta(P) = H(S(P))`. These are the "dual squared" coordinates
/// used by the doubling formula.
#[inline]
pub fn to_squared_theta<L: FpBackend>(p: &[Fp2<L>; N]) -> [Fp2<L>; N] {
    hadamard(&pointwise_square(p))
}

/// The translation-by-2-torsion action `(I, χ_J)` on a theta point.
///
/// Matches `theta_helpers_dim4.act_point`:
/// `Q[k] = (-1)^(popcount((i ⊕ k) & j) mod 2) · P[i ⊕ k]`, where `i`, `j` are
/// the 4-bit indices of the multi-indices `I`, `J`. Applied to a theta null
/// point this produces the canonical 2-torsion points used to characterise a
/// theta structure.
#[inline]
pub fn act_point<L: FpBackend>(p: &[Fp2<L>; N], i: usize, j: usize) -> [Fp2<L>; N] {
    core::array::from_fn(|k| {
        let ipk = i ^ k;
        let v = p[ipk].clone();
        if ((ipk & j).count_ones() & 1) == 1 {
            v.neg()
        } else {
            v
        }
    })
}
