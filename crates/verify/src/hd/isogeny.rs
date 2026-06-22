//! A single dimension-4 `(2,2,2,2)`-isogeny step between theta structures.
//!
//! This ports the *plain* `IsogenyDim4` of the sage `Theta_dim4` reference
//! (`isogenies/isogeny_dim4.py`): given a kernel described by four 8-torsion
//! theta points on the domain, it computes the codomain theta structure and
//! evaluates the isogeny on domain points. It does **not** cover the
//! gluing/splitting isogenies (those have a different kernel description and the
//! `special_image` evaluation; they belong to a later phase).
//!
//! # Algorithm
//!
//! Let `K₈ = (K₀, K₁, K₂, K₃)` be the four 8-torsion kernel points and
//! `HSKₖ = H(S(Kₖ))` their squared-dual coordinates (`S` = coordinate-wise
//! square, `H` = Hadamard). The codomain **dual** theta null point `O` is the
//! unique-up-to-scale vector with
//!
//! ```text
//!   O[j ⊕ 2ᵏ] / O[j] = HSKₖ[j ⊕ 2ᵏ] / HSKₖ[j]
//! ```
//!
//! for every index `j` and bit `k` with `HSKₖ[j] ≠ 0`. We pin `O[j₀] = 1` at a
//! root `j₀` and propagate these ratios over a spanning tree of the 4-cube
//! (edges = single-bit flips, the flip of bit `k` measured by kernel point
//! `k`). The codomain (standard) theta null point is then `H(O)`, and the
//! isogeny image is `image(P) = H( H(S(P)) ⊙ O⁻¹ )`.
//!
//! Because the kernel comes from a genuine isogeny, the ratios are consistent
//! (path-independent around every 2-face), so the spanning tree chosen is
//! irrelevant: `O` - hence the codomain - is deterministic, matching the
//! oracle. Comparisons are nevertheless projective (see
//! [`ThetaPointDim4::projective_eq`]).

use crate::{Fp2, FpBackend};

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::point::{ThetaPointDim4, THETA_DIM4_N as N};
use crate::hd::structure::ThetaStructureDim4;

/// A computed dimension-4 `(2,2,2,2)`-isogeny: its codomain structure plus the
/// precomputed inverse dual theta-null constants used for image evaluation.
#[derive(Clone, Debug)]
pub struct IsogenyDim4<L: FpBackend> {
    codomain: ThetaStructureDim4<L>,
    /// `1 / O[i]` for each non-zero dual null coordinate; `0` where `O[i] = 0`.
    inv_dual_null: [Fp2<L>; N],
    /// `true` for every index where the codomain dual null coordinate is zero.
    dual_null_zero: [bool; N],
    /// `true` if any dual null coordinate is zero (then [`Self::image`] needs
    /// the un-ported `special_image`; see crate notes).
    has_zero_dual: bool,
}

impl<L: FpBackend> IsogenyDim4<L> {
    /// Compute the isogeny from its kernel: four 8-torsion theta points on the
    /// domain (`4·Kₖ` is a basis of the `(2,2,2,2)` kernel).
    ///
    /// Returns `None` if the codomain dual null point cannot be reconstructed
    /// (no root spans the 4-cube). The sage reference falls back to a random
    /// symplectic base change here; that fallback is **not** ported (it is
    /// documented as exceptional in large characteristic). At Level 1 this has
    /// not been observed.
    pub fn from_kernel(k8: &[ThetaPointDim4<L>; 4]) -> Option<Self> {
        // Squared-dual coordinates of each kernel point.
        let hsk: [[Fp2<L>; N]; 4] =
            core::array::from_fn(|k| hadamard(&pointwise_square(k8[k].coords())));

        // Plain step: four single-bit kernel directions. One batched inversion
        // yields both the dual null `O` and its inverse `1/O` (image precomp).
        let (dual_null, inv_dual_null, dual_null_zero, has_zero_dual) =
            codomain_dual_null_and_inverse(&hsk, &[1, 2, 4, 8])?;

        // Codomain (standard) theta null point is the Hadamard of the dual.
        let codomain_null = ThetaPointDim4::new(hadamard(&dual_null));

        Some(Self {
            codomain: ThetaStructureDim4::new(codomain_null),
            inv_dual_null,
            dual_null_zero,
            has_zero_dual,
        })
    }

    /// The codomain theta structure.
    #[inline]
    pub fn codomain(&self) -> &ThetaStructureDim4<L> {
        &self.codomain
    }

    /// The codomain theta null point.
    #[inline]
    pub fn codomain_null(&self) -> &ThetaPointDim4<L> {
        self.codomain.null_point()
    }

    /// `true` if the codomain has a zero dual theta-null coordinate, in which
    /// case [`Self::image`] returns `None` (the generic image formula does not
    /// apply; `special_image` would be required).
    #[inline]
    pub fn has_zero_dual_null(&self) -> bool {
        self.has_zero_dual
    }

    /// Per-index flags: `true` where the codomain dual null coordinate is zero.
    #[inline]
    pub fn dual_null_zero_mask(&self) -> &[bool; N] {
        &self.dual_null_zero
    }

    /// The image precomputation `1/O` (for reuse without rebuilding the isogeny).
    #[inline]
    pub fn inv_dual_null(&self) -> &[Fp2<L>; N] {
        &self.inv_dual_null
    }

    /// Evaluate the isogeny on a domain theta point:
    /// `image(P) = H( H(S(P)) ⊙ O⁻¹ )`.
    ///
    /// Returns `None` when the codomain has a zero dual null coordinate (see
    /// [`Self::has_zero_dual_null`]); the result is otherwise a theta point on
    /// the codomain.
    #[inline]
    pub fn image(&self, p: &ThetaPointDim4<L>) -> Option<ThetaPointDim4<L>> {
        if self.has_zero_dual {
            return None;
        }
        let mut hs = hadamard(&pointwise_square(p.coords()));
        for (h, inv) in hs.iter_mut().zip(self.inv_dual_null.iter()) {
            *h = h.mul(inv);
        }
        Some(ThetaPointDim4::new(hadamard(&hs)))
    }
}

/// Apply a *plain*-step image map given its precomputed `1/O`:
/// `image(P) = H( H(S(P)) ⊙ O⁻¹ )`. Lets a caller that already holds the image
/// precomputation (e.g. the strategy chain) push a point through a plain step
/// without rebuilding the isogeny. Only valid for zero-free `O` (plain steps).
#[inline]
pub fn apply_plain_image<L: FpBackend>(
    inv_dual_null: &[Fp2<L>; N],
    p: &ThetaPointDim4<L>,
) -> ThetaPointDim4<L> {
    let mut hs = hadamard(&pointwise_square(p.coords()));
    for (h, inv) in hs.iter_mut().zip(inv_dual_null.iter()) {
        *h = h.mul(inv);
    }
    ThetaPointDim4::new(hadamard(&hs))
}

/// Reconstruct a codomain dual theta null point `O` from kernel points'
/// squared-dual coordinates and their index-flip directions.
///
/// `hsk[k] = H(S(K_k))` for each kernel point, and `dirs[k]` is the index XOR
/// that translation by `K_k` induces (`2ᵏ` for a single-bit kernel generator,
/// or e.g. `3 = 2⁰⊕2¹` for a diagonal generator used by the gluing). Pin
/// `O[j₀] = 1` at a root and propagate `O[j⊕dirs[k]] = O[j]·HSKₖ[j⊕dirs[k]]/
/// HSKₖ[j]` over a spanning tree (an edge is usable iff its denominator
/// `HSKₖ[j] ≠ 0`; the numerator may be zero, yielding a zero `O` coordinate).
/// Tries successive roots until one spans all 16 indices; `None` if none does.
///
/// Shared by the plain step (`dirs = [1,2,4,8]`) and the gluing step
/// (`dirs = [1,2,4,8,3]`, whose diagonal edge spans the cube despite the zero
/// dual coordinates).
#[allow(clippy::type_complexity)] // (O, 1/O, zero-mask, has-zero) - a small tuple.
pub(crate) fn codomain_dual_null_and_inverse<L: FpBackend>(
    hsk: &[[Fp2<L>; N]],
    dirs: &[usize],
) -> Option<([Fp2<L>; N], [Fp2<L>; N], [bool; N], bool)> {
    debug_assert_eq!(hsk.len(), dirs.len());
    for j0 in 0..N {
        // The root needs at least one usable outgoing edge.
        if (0..hsk.len()).all(|k| bool::from(hsk[k][j0].ct_is_zero())) {
            continue;
        }

        // Pass 1 - BFS over the 4-cube using only the (value-independent) zero
        // pattern, recording the spanning-tree edges in discovery order with
        // their numerators/denominators. No inversions yet.
        let mut covered = [false; N];
        let mut queue = [0usize; N];
        let (mut head, mut tail) = (0usize, 0usize);
        let mut edge_j = [0usize; N];
        let mut edge_jpk = [0usize; N];
        let mut num: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        let mut den: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        let mut ne = 0usize;

        covered[j0] = true;
        queue[tail] = j0;
        tail += 1;
        while head < tail {
            let j = queue[head];
            head += 1;
            for (k, hsk_k) in hsk.iter().enumerate() {
                let jpk = j ^ dirs[k];
                if !covered[jpk] && !bool::from(hsk_k[j].ct_is_zero()) {
                    edge_j[ne] = j;
                    edge_jpk[ne] = jpk;
                    num[ne] = hsk_k[jpk].clone();
                    den[ne] = hsk_k[j].clone();
                    ne += 1;
                    covered[jpk] = true;
                    queue[tail] = jpk;
                    tail += 1;
                }
            }
        }
        if !covered.iter().all(|&c| c) {
            continue;
        }

        // Pass 2 - propagate the path products (no inversions):
        // `O[j] = p_num[j] / p_den[j]` where `p_num`/`p_den` are the running
        // products of the edge numerators/denominators from the root.
        let mut p_num: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        let mut p_den: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        p_num[j0] = Fp2::one();
        p_den[j0] = Fp2::one();
        for e in 0..ne {
            let (j, jpk) = (edge_j[e], edge_jpk[e]);
            p_num[jpk] = p_num[j].mul(&num[e]);
            p_den[jpk] = p_den[j].mul(&den[e]);
        }

        // A single batched inversion over all `p_den` (always non-zero) plus the
        // non-zero `p_num`. Then `O = p_num·(1/p_den)` and `1/O = p_den·(1/p_num)`.
        let mut batch: [Fp2<L>; 2 * N] = core::array::from_fn(|_| Fp2::zero());
        batch[..N].clone_from_slice(&p_den[..N]);
        let mut idx = [0usize; N];
        let mut cnt = 0;
        let mut zero_mask = [false; N];
        let mut has_zero = false;
        for j in 0..N {
            // `p_num[j] == 0` ⇔ a zero edge numerator on the path ⇔ `O[j] = 0`.
            if bool::from(p_num[j].ct_is_zero()) {
                zero_mask[j] = true;
                has_zero = true;
            } else {
                batch[N + cnt] = p_num[j].clone();
                idx[cnt] = j;
                cnt += 1;
            }
        }
        let nb = N + cnt;
        let mut t1: [Fp2<L>; 2 * N] = core::array::from_fn(|_| Fp2::zero());
        let mut t2: [Fp2<L>; 2 * N] = core::array::from_fn(|_| Fp2::zero());
        crate::hd::field::batched_inv(&mut batch[..nb], &mut t1[..nb], &mut t2[..nb]);

        let mut o: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        let mut inv_o: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
        for j in 0..N {
            // `O[j] = p_num[j] / p_den[j]` (zero where `p_num[j] = 0`).
            o[j] = p_num[j].mul(&batch[j]);
        }
        for e in 0..cnt {
            let j = idx[e];
            inv_o[j] = p_den[j].mul(&batch[N + e]); // `1/O[j] = p_den[j] / p_num[j]`
        }
        return Some((o, inv_o, zero_mask, has_zero));
    }
    None
}
