//! Phase 5b.6 (front half) - the dimension-4 symplectic→theta base change
//! `N_dim4` (`basis_change/base_change_dim4.py::base_change_theta_dim4`) and its
//! application to a 16-coordinate theta point.
//!
//! This is the dim-4 analogue of [`crate::hd::base_change_theta_dim2`]: given a
//! symplectic `M ∈ Sp₈(Z/4)` (the gluing change of basis on `Am²`, derived from
//! the Kani integers) and the 4th root of unity `e4 = e₄(T1,T2)`, it builds the
//! `16×16` matrix `N` mapping product-theta coordinates on `Am²` to the theta
//! coordinates adapted to the dim-4 gluing. The index convention is the
//! `Theta_dim4` one (`k = i₀ + 2i₁ + 4i₂ + 8i₃`).

use crate::{Fp2, FpBackend};

/// `[k]` from a 4-bit multi-index reduced mod 2 (`multindex_to_index∘red_mod_2`).
#[inline]
fn red2_index(x: &[i64; 4]) -> usize {
    (0..4).map(|k| ((x[k].rem_euclid(2)) as usize) << k).sum()
}

/// `A·I` for a `4×4` integer block and a 4-vector (`mat_prod_vect`).
#[inline]
fn mat_prod_vect(a: &[[i64; 4]; 4], i: &[i64; 4]) -> [i64; 4] {
    core::array::from_fn(|r| (0..4).map(|k| a[r][k] * i[k]).sum())
}

#[inline]
fn scal_prod(i: &[i64; 4], j: &[i64; 4]) -> i64 {
    (0..4).map(|k| i[k] * j[k]).sum()
}

#[inline]
fn add4(i: &[i64; 4], j: &[i64; 4]) -> [i64; 4] {
    core::array::from_fn(|k| i[k] + j[k])
}

/// `base_change_theta_dim4(M, zeta)`: the `16×16` theta-coordinate change `N`
/// induced by a symplectic `M ∈ Sp₈(Z/4)` (stored as an `8×8` array of `i64`,
/// entries taken mod 4) and a primitive 4th root of unity `zeta`. `(Nᵢ) = N·(θᵢ)`.
pub fn base_change_theta_dim4<L: FpBackend>(m: &[[i64; 8]; 8], zeta: &Fp2<L>) -> [[Fp2<L>; 16]; 16] {
    // Blocks mod 4: A = M[0..4][0..4], B = M[4..8][0..4], C = M[0..4][4..8],
    // D = M[4..8][4..8].
    let blk = |r0: usize, c0: usize| -> [[i64; 4]; 4] {
        core::array::from_fn(|i| core::array::from_fn(|j| m[r0 + i][c0 + j].rem_euclid(4)))
    };
    let a = blk(0, 0);
    let b = blk(4, 0);
    let c = blk(0, 4);
    let d = blk(4, 4);

    let zpow = |e: i64| -> Fp2<L> {
        match e.rem_euclid(4) {
            0 => Fp2::<L>::one(),
            1 => zeta.clone(),
            2 => zeta.sqr(),
            _ => zeta.sqr().mul(zeta),
        }
    };

    let zero16 = || core::array::from_fn::<Fp2<L>, 16, _>(|_| Fp2::zero());

    // choose_non_vanishing_index: first I0 ∈ {0,1}^4 whose row L0 is non-zero.
    let mut i0_ref = [0i64; 4];
    let mut l0 = zero16();
    'outer: for i0bits in 0..16u32 {
        let i0: [i64; 4] = core::array::from_fn(|k| ((i0bits >> k) & 1) as i64);
        let mut l = zero16();
        for jbits in 0..16u32 {
            let jj: [i64; 4] = core::array::from_fn(|k| ((jbits >> k) & 1) as i64);
            let cj = mat_prod_vect(&c, &jj);
            let dj = mat_prod_vect(&d, &jj);
            let e = -scal_prod(&cj, &dj) - 2 * scal_prod(&i0, &dj);
            let idx = red2_index(&add4(&i0, &cj));
            l[idx] = l[idx].add(&zpow(e));
        }
        if l.iter().any(|x| !bool::from(x.ct_is_zero())) {
            i0_ref = i0;
            l0 = l;
            break 'outer;
        }
    }

    let mut n: [[Fp2<L>; 16]; 16] = core::array::from_fn(|_| zero16());
    n[0] = l0; // row for I = (0,0,0,0)
    for ibits in 1..16u32 {
        let ii: [i64; 4] = core::array::from_fn(|k| ((ibits >> k) & 1) as i64);
        let ai = mat_prod_vect(&a, &ii);
        let bi = mat_prod_vect(&b, &ii);
        for jbits in 0..16u32 {
            let jj: [i64; 4] = core::array::from_fn(|k| ((jbits >> k) & 1) as i64);
            let cj = mat_prod_vect(&c, &jj);
            let dj = mat_prod_vect(&d, &jj);
            let aicj = add4(&ai, &cj);
            let bidj = add4(&bi, &dj);
            let e =
                scal_prod(&ii, &jj) - scal_prod(&aicj, &bidj) - 2 * scal_prod(&i0_ref, &bidj);
            let col = red2_index(&add4(&aicj, &i0_ref));
            n[ibits as usize][col] = n[ibits as usize][col].add(&zpow(e));
        }
    }
    n
}

/// Split a dim-4 **product** theta point (16 coords) into its two dim-2 factors
/// (`theta_helpers_dim4.product_to_theta_points_dim4_dim2`, the inverse of
/// [`crate::hd::product_theta_dim2to4`]). Finds the first nonzero index
/// `k0 = i0 + 4·j0`, then `P1[i] = P[i+4·j0]/P[k0]` (`P1[i0]=1`) and
/// `P2[j] = P[i0+4·j]/P[k0]` (`P2[j0]=1`). Returns `None` if `P = 0`.
#[allow(clippy::type_complexity)] // the two dim-2 factors of a dim-4 product.
pub fn product_to_theta_points_dim4_dim2<L: FpBackend>(
    p: &[Fp2<L>; 16],
) -> Option<([Fp2<L>; 4], [Fp2<L>; 4])> {
    let k0 = (0..16).find(|&k| !bool::from(p[k].ct_is_zero()))?;
    let (i0, j0) = (k0 % 4, k0 / 4);
    let inv = crate::hd::field::inv(&p[k0]);
    let p1 = core::array::from_fn(|i| {
        if i == i0 {
            Fp2::one()
        } else {
            p[i + 4 * j0].mul(&inv)
        }
    });
    let p2 = core::array::from_fn(|j| {
        if j == j0 {
            Fp2::one()
        } else {
            p[i0 + 4 * j].mul(&inv)
        }
    });
    Some((p1, p2))
}

/// Apply a `16×16` theta base-change matrix to a 16-coordinate theta point.
#[inline]
pub fn apply_base_change_theta_dim4<L: FpBackend>(
    n: &[[Fp2<L>; 16]; 16],
    p: &[Fp2<L>; 16],
) -> [Fp2<L>; 16] {
    core::array::from_fn(|i| {
        let mut acc = n[i][0].mul(&p[0]);
        for j in 1..16 {
            acc = acc.add(&n[i][j].mul(&p[j]));
        }
        acc
    })
}
