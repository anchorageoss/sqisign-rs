//! Phase 5b.6 (front half) - the dim-4 gluing chain `F_{m+1}∘…∘F_1`
//! (`isogenies/Kani_gluing_isogeny_chain_dim4.py::KaniGluingIsogenyChainDim4Half`).
//!
//! This is the assembly that produces the *gluing output* the strategy loop
//! (5b.6 back half) consumes: the gluing-chain codomain (a dim-4 theta structure)
//! and - via [`KaniGluingChainHalf::evaluate`] - the post-gluing kernel basis.
//! It glues `m+1` levels: an `m`-step dim-2 `(2,2)`-isogeny chain
//! ([`crate::hd::IsogenyChainDim2`]) on each of the two `E1×E2` factors, then one
//! dim-4 gluing ([`crate::hd::GluingIsogenyDim4`], Phase 4) of the resulting product
//! `Am²` into `B`.
//!
//! # Completion dependence
//!
//! The dim-4 base change `N_dim4 = base_change_theta_dim4(M_gluing_dim4, e4)` is
//! built from the dim-4 symplectic matrix `M_gluing_dim4`, which derives from the
//! *completed* starting matrix `M1`/`M2` (`bloc_decomposition(M1)`). The
//! symplectic completion is non-unique (5b.4); a different completion yields a
//! different `N_dim4` and hence a different (but isomorphic) gluing codomain in a
//! different theta coordinate system. The completion-*independent* inner data
//! (the dim-2 codomain) matches the oracle exactly; the completion-*dependent*
//! gluing output need only yield the invariant middle-codomain match (the F1/F2
//! meet-in-the-middle), which it does because `M1` and `M2` share one completion
//! propagated through the closed-form `matrix_F`/`matrix_F_dual`.

use alloc::vec::Vec;

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::{EcCurve, JacPoint};
use crate::{Fp2, FpBackend};

use crate::hd::arith::{hadamard, pointwise_square};
use crate::hd::dim2::{hadamard2, IsogenyChainDim2, TuplePoint};
use crate::hd::dim4::{
    apply_base_change_theta_dim4, base_change_theta_dim4, product_to_theta_points_dim4_dim2,
};
use crate::hd::gluing::{GluingIsogenyDim4, GLUING_KERNEL_DIRS};
use crate::hd::kani::{gluing_dim2_f1, gluing_dim2_f2};
use crate::hd::point::ThetaPointDim4;
use crate::hd::product_theta::product_theta_dim2to4;

/// `[k]·P` on a Montgomery/Weierstrass curve for an arbitrary `u128` scalar,
/// MSB-first double-and-add. `jac_add`/`jac_dbl` handle identity / `P=±Q`.
pub fn jac_mul_u128<L: FpBackend>(p: &JacPoint<L>, k: u128, curve: &EcCurve<L>) -> JacPoint<L> {
    if k == 0 {
        return JacPoint::identity();
    }
    let top = 127 - k.leading_zeros();
    let mut acc = p.clone();
    for i in (0..top).rev() {
        acc = jac_dbl(&acc, curve);
        if (k >> i) & 1 == 1 {
            acc = jac_add(&acc, p, curve);
        }
    }
    acc
}

/// A point on the product fourfold `E1² × E2²` (components `[E1, E1, E2, E2]`).
#[derive(Clone, Debug)]
pub struct TuplePoint4<L: FpBackend> {
    pub c: [JacPoint<L>; 4],
}

impl<L: FpBackend> TuplePoint4<L> {
    #[inline]
    pub fn new(c0: JacPoint<L>, c1: JacPoint<L>, c2: JacPoint<L>, c3: JacPoint<L>) -> Self {
        Self {
            c: [c0, c1, c2, c3],
        }
    }
    #[inline]
    pub fn add(&self, o: &Self, e1: &EcCurve<L>, e2: &EcCurve<L>) -> Self {
        Self {
            c: [
                jac_add(&self.c[0], &o.c[0], e1),
                jac_add(&self.c[1], &o.c[1], e1),
                jac_add(&self.c[2], &o.c[2], e2),
                jac_add(&self.c[3], &o.c[3], e2),
            ],
        }
    }
    #[inline]
    pub fn double(&self, e1: &EcCurve<L>, e2: &EcCurve<L>) -> Self {
        Self {
            c: [
                jac_dbl(&self.c[0], e1),
                jac_dbl(&self.c[1], e1),
                jac_dbl(&self.c[2], e2),
                jac_dbl(&self.c[3], e2),
            ],
        }
    }
}

/// `point_matrix_product(M, [P1,P2,R1,R2], [4,5,6,7], modulus)`
/// (`basis_change/kani_base_change.py`): the four kernel `TuplePoint4`s read
/// off columns 4..7 of an `8×8` matrix `M` (entries reduced mod `mask+1`), where
/// the rows act on `(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1),(P2,…),…`.
#[allow(clippy::too_many_arguments)]
pub fn point_matrix_product_k<L: FpBackend>(
    m: &[[u128; 8]; 8],
    p1: &JacPoint<L>,
    p2: &JacPoint<L>,
    r1: &JacPoint<L>,
    r2: &JacPoint<L>,
    mask: u128,
    e1: &EcCurve<L>,
    e2: &EcCurve<L>,
) -> [TuplePoint4<L>; 4] {
    core::array::from_fn(|jj| {
        let j = 4 + jj;
        let comp0 = jac_add(
            &jac_mul_u128(p1, m[0][j] & mask, e1),
            &jac_mul_u128(p2, m[4][j] & mask, e1),
            e1,
        );
        let comp1 = jac_add(
            &jac_mul_u128(p1, m[1][j] & mask, e1),
            &jac_mul_u128(p2, m[5][j] & mask, e1),
            e1,
        );
        let comp2 = jac_add(
            &jac_mul_u128(r1, m[2][j] & mask, e2),
            &jac_mul_u128(r2, m[6][j] & mask, e2),
            e2,
        );
        let comp3 = jac_add(
            &jac_mul_u128(r1, m[3][j] & mask, e2),
            &jac_mul_u128(r2, m[7][j] & mask, e2),
            e2,
        );
        TuplePoint4::new(comp0, comp1, comp2, comp3)
    })
}

/// `4×4 · 4×4` over `Z/4`, returned as `i64` (for the dim-2 base change).
fn m4_mul(a: &[[u8; 4]; 4], b: &[[u8; 4]; 4]) -> [[i64; 4]; 4] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut s = 0u32;
            for k in 0..4 {
                s += (a[i][k] as u32) * (b[k][j] as u32);
            }
            (s % 4) as i64
        })
    })
}

/// A computed dim-4 gluing chain (one half-chain), with the data needed to
/// evaluate it on the full kernel basis.
pub struct KaniGluingChainHalf<L: FpBackend> {
    isogenies_dim2: IsogenyChainDim2<L>,
    n_dim4: [[Fp2<L>; 16]; 16],
    gluing_dim4: GluingIsogenyDim4<L>,
    l_trans: [TuplePoint4<L>; 2],
    l_trans_ind: [usize; 2],
    e1: EcCurve<L>,
    e2: EcCurve<L>,
}

impl<L: FpBackend> KaniGluingChainHalf<L> {
    /// Build the gluing chain.
    ///
    /// * `points_m = [P1_m, Q1_m, R2_m, S2_m]` of order `2^{m+3}`.
    /// * `zero12`: the product dim-2 theta null (domain of the dim-2 chain).
    /// * `m0`: `M_product_dim2` (`4×4` over `Z/4`, from the canonical bases).
    /// * `e4`: `e₄(T1,T2)` (PARI convention).
    /// * `m1_full`: the starting symplectic matrix `M1`/`M2` (`8×8` over `Z/2^f`),
    ///   used for the dim-4 kernel basis `B_K_dim4`.
    /// * `m_gluing_dim4`: the dim-4 gluing change of basis (`8×8`, entries mod 4).
    /// * `dual`: `false` for F1, `true` for F2_dual.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        points_m: &[JacPoint<L>; 4],
        zero12: &[Fp2<L>; 4],
        m0: &[[u8; 4]; 4],
        e4: &Fp2<L>,
        a1: u128,
        a2: u128,
        q: u128,
        m: usize,
        m1_full: &[[u128; 8]; 8],
        m_gluing_dim4: &[[i64; 8]; 8],
        dual: bool,
        e1: &EcCurve<L>,
        e2: &EcCurve<L>,
    ) -> Option<Self> {
        let [p1_m, q1_m, r2_m, s2_m] = points_m;

        // dim-2 gluing base change: N_dim2 = base_change(M0·M1_glue, e4)
        let m1_glue = if !dual {
            gluing_dim2_f1(a1, a2, q)
        } else {
            gluing_dim2_f2(a1, a2, q)
        };
        let m10 = m4_mul(m0, &m1_glue);
        let n_dim2 = crate::hd::dim2::base_change_theta_dim2(&m10, e4);

        // dim-4 gluing base change
        let n_dim4 = base_change_theta_dim4(m_gluing_dim4, e4);

        // the dim-2 (2,2)-isogeny chain kernel B_K_dim2
        let two_mp2 = 1u128 << (m + 2);
        let a1r = a1 % two_mp2;
        let a2r = a2 % two_mp2;
        let (s1c, s2c) = (2 * a1r, 2 * a2r);
        let s1p = jac_mul_u128(p1_m, s1c, e1);
        let s2p = jac_mul_u128(p1_m, s2c, e1);
        let s1q = jac_mul_u128(q1_m, s1c, e1);
        let s2q = jac_mul_u128(q1_m, s2c, e1);
        let two_r2 = jac_dbl(r2_m, e2);
        let two_s2 = jac_dbl(s2_m, e2);
        let (tp0, tp1) = if !dual {
            (
                TuplePoint::new(jac_add(&s1p, &s2q.neg(), e1), two_r2.clone()),
                TuplePoint::new(jac_add(&s1q, &s2p, e1), two_s2.clone()),
            )
        } else {
            (
                TuplePoint::new(jac_add(&s1p, &s2q, e1), two_r2.neg()),
                TuplePoint::new(jac_add(&s1q, &s2p.neg(), e1), two_s2.neg()),
            )
        };
        let isogenies_dim2 = IsogenyChainDim2::new(&tp0, &tp1, zero12, &n_dim2, m, e1, e2)?;

        // the dim-4 gluing kernel B_K_dim4 = kernel_basis(M1, m+1, …)
        let two_mp3 = 1u128 << (m + 3);
        let mask_mp3 = two_mp3 - 1;
        let lamb = crate::hd::kani::inverse_mod_pow2(q, mask_mp3); // inverse_mod(q, 2^{m+3})
        let lamb_s2 = jac_mul_u128(s2_m, lamb, e2);
        let b_k_dim4 =
            point_matrix_product_k(m1_full, p1_m, q1_m, r2_m, &lamb_s2, mask_mp3, e1, e2);

        // L_K_dim4 = B_K_dim4 + [B_K_dim4[0] + B_K_dim4[1]] (5 kernel directions).
        let mut l_k: Vec<TuplePoint4<L>> = b_k_dim4.to_vec();
        l_k.push(b_k_dim4[0].add(&b_k_dim4[1], e1, e2));

        // Push each through the two dim-2 chains, product, base change → dim-4.
        let mut k8: Vec<ThetaPointDim4<L>> = Vec::with_capacity(5);
        for t in &l_k {
            k8.push(self_eval_to_dim4(&isogenies_dim2, &n_dim4, t, e1, e2));
        }
        let k8_arr: [ThetaPointDim4<L>; 5] = core::array::from_fn(|i| k8[i].clone());
        let gluing_dim4 = GluingIsogenyDim4::from_kernel(&k8_arr, &GLUING_KERNEL_DIRS)?;

        // Translates for special_image: 2·B_K_dim4[0], 2·B_K_dim4[1], indices 1,2.
        let l_trans = [b_k_dim4[0].double(e1, e2), b_k_dim4[1].double(e1, e2)];

        Some(Self {
            isogenies_dim2,
            n_dim4,
            gluing_dim4,
            l_trans,
            l_trans_ind: [1, 2],
            e1: e1.clone(),
            e2: e2.clone(),
        })
    }

    /// The gluing-chain codomain theta null point (`glue_codomain_null`).
    #[inline]
    pub fn codomain_null(&self) -> &ThetaPointDim4<L> {
        self.gluing_dim4.codomain_null()
    }

    /// The inner dim-2 isogeny chain (used by the dual/splitting evaluation).
    #[inline]
    pub fn dim2_chain(&self) -> &IsogenyChainDim2<L> {
        &self.isogenies_dim2
    }

    /// The splitting (dual gluing): a dim-4 theta point on `C₀.hadamard()` back
    /// to a `TuplePoint4` on `E1²×E2²` (`KaniSplittingIsogenyChainDim4`). Used by
    /// the dual chain `F2 = F2_dual.dual()` for the stage-6 HD-image check.
    ///
    /// `gluing_dim4.dual()` is `DualIsogenyDim4(C₀, domain_base_change,
    /// hadamard=false)`: image `= H(S(P)) ⊙ inv(domain_base_change null)`. Then
    /// `N_dim4⁻¹` maps to the product `Am²`, the product is split into two dim-2
    /// theta points, each Hadamard'd and pushed through the dim-2 dual chain.
    pub fn splitting_eval(&self, y: &ThetaPointDim4<L>) -> Option<TuplePoint4<L>> {
        // domain_base_change null = N_dim4 · (dim2_cod ⊗ dim2_cod).
        let dim2_cod = self.isogenies_dim2.codomain_null();
        let product_null = product_theta_dim2to4(dim2_cod, dim2_cod);
        let dbc_null = apply_base_change_theta_dim4(&self.n_dim4, &product_null);
        if dbc_null.iter().any(|x| bool::from(x.ct_is_zero())) {
            return None;
        }
        let inv_dbc: [Fp2<L>; 16] = core::array::from_fn(|i| crate::hd::field::inv(&dbc_null[i]));
        let hs = hadamard(&pointwise_square(y.coords()));
        let q: [Fp2<L>; 16] = core::array::from_fn(|i| hs[i].mul(&inv_dbc[i])); // hadamard=false

        let n_inv = crate::hd::dim2::mat_inverse(&self.n_dim4)?;
        let q = apply_base_change_theta_dim4(&n_inv, &q); // product Am² coords

        let (q1, q2) = product_to_theta_points_dim4_dim2(&q)?;
        let q1h = hadamard2(&q1);
        let q2h = hadamard2(&q2);
        let tp1 = self.isogenies_dim2.dual_eval(&q1h)?;
        let tp2 = self.isogenies_dim2.dual_eval(&q2h)?;
        Some(TuplePoint4::new(tp1.p1, tp2.p1, tp1.p2, tp2.p2))
    }

    /// The gluing-chain codomain theta structure.
    #[inline]
    pub fn codomain(&self) -> &crate::hd::structure::ThetaStructureDim4<L> {
        self.gluing_dim4.codomain()
    }

    /// Evaluate the gluing chain on a `TuplePoint4` of the domain `E1²×E2²`,
    /// returning the dim-4 theta image on the codomain.
    pub fn evaluate(&self, p: &TuplePoint4<L>) -> ThetaPointDim4<L> {
        let eval_p = self.split_eval(p);
        let l_p_trans: [ThetaPointDim4<L>; 2] = core::array::from_fn(|k| {
            let q = p.add(&self.l_trans[k], &self.e1, &self.e2);
            self.split_eval(&q)
        });
        self.gluing_dim4
            .special_image(&eval_p, &l_p_trans, &self.l_trans_ind)
    }

    /// dim4 → (dim2 × dim2) split, push through the dim-2 chain, product, base change.
    fn split_eval(&self, p: &TuplePoint4<L>) -> ThetaPointDim4<L> {
        self_eval_to_dim4(&self.isogenies_dim2, &self.n_dim4, p, &self.e1, &self.e2)
    }
}

/// Shared splitting/evaluation: `TuplePoint4 → [Φ(P0,P2), Φ(P1,P3)] → product →
/// N_dim4 base change`, where `Φ` is the dim-2 chain.
fn self_eval_to_dim4<L: FpBackend>(
    chain: &IsogenyChainDim2<L>,
    n_dim4: &[[Fp2<L>; 16]; 16],
    p: &TuplePoint4<L>,
    _e1: &EcCurve<L>,
    _e2: &EcCurve<L>,
) -> ThetaPointDim4<L> {
    let a = TuplePoint::new(p.c[0].clone(), p.c[2].clone());
    let b = TuplePoint::new(p.c[1].clone(), p.c[3].clone());
    let ta = chain.eval(&a);
    let tb = chain.eval(&b);
    let prod = product_theta_dim2to4(&ta, &tb);
    ThetaPointDim4::new(apply_base_change_theta_dim4(n_dim4, &prod))
}
