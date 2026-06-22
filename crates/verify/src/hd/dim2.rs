//! Phase 5b.6 (front half) - the dimension-2 `(2,2)`-isogeny chain in the theta
//! model, the gluing chain's inner engine.
//!
//! Ported from the SQIsignHD reference (`isogenies_dim2/` and
//! `basis_change/base_change_dim2.py`), which is the ThetaIsogenies/
//! two-isogenies library - the same object as `sqisign-verify::theta`. But that
//! crate's dim-2 theta is built over **elliptic products** with its own
//! theta-null convention, whereas HD works over a **product theta structure**
//! keyed to a canonical 4-torsion basis with an explicit symplectic→theta base
//! change `N` (5b.5 found these conventions differ). So this is a fresh port in
//! the dim-4 theta idiom; the reuse is `Fp2`, the EC Jacobian arithmetic (for the
//! `TuplePoint`s on `E1×E2`), and [`crate::hd::optimised_strategy`].
//!
//! The chain is short at Level 1 (`m ≤ 6` steps): one gluing
//! ([`GluingThetaIsogenyDim2`], elliptic product → theta) then `m-1` plain
//! steps ([`ThetaIsogenyDim2`]). Everything is projective (4 coords).

use alloc::vec;
use alloc::vec::Vec;

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::{EcCurve, JacPoint};
use crate::{Fp2, FpBackend};

use crate::hd::strategy::optimised_strategy;

/// A pair of dim-2 theta points (two kernel generators living on one codomain).
type ThetaPair<L> = ([Fp2<L>; 4], [Fp2<L>; 4]);

/// Inverse of an `N×N` matrix over `Fp2` by Gauss-Jordan elimination with
/// pivoting. Returns `None` if singular. Used for the dual base changes `N⁻¹`
/// (4×4 dim-2 gluing split, 16×16 dim-4 split).
pub fn mat_inverse<L: FpBackend, const N: usize>(
    m: &[[Fp2<L>; N]; N],
) -> Option<[[Fp2<L>; N]; N]> {
    let mut a: [[Fp2<L>; N]; N] =
        core::array::from_fn(|i| core::array::from_fn(|j| m[i][j].clone()));
    let mut inv: [[Fp2<L>; N]; N] =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { Fp2::one() } else { Fp2::zero() }));
    for col in 0..N {
        // Find a nonzero pivot in this column at or below `col`.
        let piv = (col..N).find(|&r| !bool::from(a[r][col].ct_is_zero()))?;
        a.swap(col, piv);
        inv.swap(col, piv);
        let s = crate::hd::field::inv(&a[col][col]);
        for j in 0..N {
            a[col][j] = a[col][j].mul(&s);
            inv[col][j] = inv[col][j].mul(&s);
        }
        for r in 0..N {
            if r != col && !bool::from(a[r][col].ct_is_zero()) {
                let f = a[r][col].clone();
                for j in 0..N {
                    a[r][j] = a[r][j].sub(&f.mul(&a[col][j]));
                    inv[r][j] = inv[r][j].sub(&f.mul(&inv[col][j]));
                }
            }
        }
    }
    Some(inv)
}

/// `[x+1, …]` per-coordinate inverse of a 4-vector; `None` if any coord is zero.
fn inv4<L: FpBackend>(v: &[Fp2<L>; 4]) -> Option<[Fp2<L>; 4]> {
    if v.iter().any(|x| bool::from(x.ct_is_zero())) {
        return None;
    }
    Some(core::array::from_fn(|i| crate::hd::field::inv(&v[i])))
}

/// Lift a Montgomery `(X : Z)` Kummer coordinate to a Jacobian point on
/// `y² = x³ + A x² + x`. The `y`-sign is arbitrary (the HD-image check compares
/// only `x`, i.e. up to `±`). `Z = 0` is the identity.
fn lift_kummer<L: FpBackend>(curve: &EcCurve<L>, x: &Fp2<L>, z: &Fp2<L>) -> JacPoint<L> {
    if bool::from(z.ct_is_zero()) {
        return JacPoint::identity();
    }
    let xa = x.mul(&crate::hd::field::inv(z));
    let y2 = xa.mul(&xa.sqr().add(&curve.a.mul(&xa)).add(&Fp2::one()));
    JacPoint::new(xa, y2.sqrt(), Fp2::one())
}

// 4-coordinate theta arithmetic.

/// Hadamard transform of 4 coordinates (`ThetaPointDim2.to_hadamard`):
/// `(a+b+c+d, a-b+c-d, a+b-c-d, a-b-c+d)`.
#[inline]
pub fn hadamard2<L: FpBackend>(c: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
    let x00 = c[0].add(&c[1]);
    let x10 = c[0].sub(&c[1]);
    let x01 = c[2].add(&c[3]);
    let x11 = c[2].sub(&c[3]);
    [x00.add(&x01), x10.add(&x11), x00.sub(&x01), x10.sub(&x11)]
}

/// Squared-theta transform: Hadamard of the squared coordinates.
#[inline]
pub fn squared_theta2<L: FpBackend>(c: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
    let sq: [Fp2<L>; 4] = core::array::from_fn(|i| c[i].sqr());
    hadamard2(&sq)
}

/// A dim-2 theta structure: its 4-coordinate null point and the 6 precomputed
/// arithmetic constants `(y0,z0,t0, Y0,Z0,T0)` used by [`Self::double`].
#[derive(Clone, Debug)]
pub struct ThetaStructureDim2<L: FpBackend> {
    null: [Fp2<L>; 4],
    prec: Option<[Fp2<L>; 6]>,
}

impl<L: FpBackend> ThetaStructureDim2<L> {
    #[inline]
    pub fn new(null: [Fp2<L>; 4]) -> Self {
        Self { null, prec: None }
    }

    #[inline]
    pub fn null(&self) -> &[Fp2<L>; 4] {
        &self.null
    }

    /// Precompute `y0=a/b, z0=a/c, t0=a/d, Y0=AA/BB, Z0=AA/CC, T0=AA/DD`
    /// (`ThetaStructureDim2._arithmetic_precomputation`), where `(AA,BB,CC,DD)`
    /// is the squared-theta of the null. Idempotent.
    pub fn precompute(&mut self) {
        if self.prec.is_some() {
            return;
        }
        let [a, b, c, d] = &self.null;
        let st = squared_theta2(&self.null);
        // One batched inversion of (b, c, d, BB, CC, DD) instead of six.
        let mut inv = [
            b.clone(),
            c.clone(),
            d.clone(),
            st[1].clone(),
            st[2].clone(),
            st[3].clone(),
        ];
        let mut t1: [Fp2<L>; 6] = core::array::from_fn(|_| Fp2::zero());
        let mut t2: [Fp2<L>; 6] = core::array::from_fn(|_| Fp2::zero());
        crate::hd::field::batched_inv(&mut inv, &mut t1, &mut t2);
        let y0 = a.mul(&inv[0]);
        let z0 = a.mul(&inv[1]);
        let t0 = a.mul(&inv[2]);
        let big_y0 = st[0].mul(&inv[3]);
        let big_z0 = st[0].mul(&inv[4]);
        let big_t0 = st[0].mul(&inv[5]);
        self.prec = Some([y0, z0, t0, big_y0, big_z0, big_t0]);
    }

    /// `2·P` (`ThetaPointDim2.double`); requires [`Self::precompute`].
    pub fn double(&self, p: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
        let [y0, z0, t0, big_y0, big_z0, big_t0] =
            self.prec.as_ref().expect("ThetaStructureDim2::double needs precompute()");
        let st = squared_theta2(p);
        let xp = st[0].sqr();
        let yp = big_y0.mul(&st[1].sqr());
        let zp = big_z0.mul(&st[2].sqr());
        let tp = big_t0.mul(&st[3].sqr());
        let h = hadamard2(&[xp, yp, zp, tp]);
        [h[0].clone(), y0.mul(&h[1]), z0.mul(&h[2]), t0.mul(&h[3])]
    }

    /// `2ⁿ·P`.
    pub fn double_iter(&self, p: &[Fp2<L>; 4], n: u32) -> [Fp2<L>; 4] {
        let mut acc = p.clone();
        for _ in 0..n {
            acc = self.double(&acc);
        }
        acc
    }
}

// Symplectic → theta base change.

/// `base_change_theta_dim2(M, zeta)`: the `4×4` theta-coordinate change `N`
/// induced by a symplectic `M ∈ Sp₄(Z/4)`, with `zeta` a primitive 4th root of
/// unity (`e₄` of the symplectic basis). `M` is the four `2×2` blocks
/// `[[A,C],[B,D]]` stored as a `4×4` array of `i64` (entries mod 4).
pub fn base_change_theta_dim2<L: FpBackend>(m: &[[i64; 4]; 4], zeta: &Fp2<L>) -> [[Fp2<L>; 4]; 4] {
    // Powers zeta^e for e mod 4.
    let zpow = |e: i64| -> Fp2<L> {
        match e.rem_euclid(4) {
            0 => Fp2::<L>::one(),
            1 => zeta.clone(),
            2 => zeta.sqr(),
            _ => zeta.sqr().mul(zeta),
        }
    };
    // Blocks: A = M[0..2][0..2], B = M[2..4][0..2], C = M[0..2][2..4], D = M[2..4][2..4].
    let a = [[m[0][0], m[0][1]], [m[1][0], m[1][1]]];
    let b = [[m[2][0], m[2][1]], [m[3][0], m[3][1]]];
    let c = [[m[0][2], m[0][3]], [m[1][2], m[1][3]]];
    let d = [[m[2][2], m[2][3]], [m[3][2], m[3][3]]];

    // choose_non_vanishing_index over (ir0, ir1).
    let mut chosen = (0i64, 0i64);
    'outer: for ir0 in 0..2i64 {
        for ir1 in 0..2i64 {
            let mut l = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
            for j0 in 0..2i64 {
                for j1 in 0..2i64 {
                    let k0 = c[0][0] * j0 + c[0][1] * j1;
                    let k1 = c[1][0] * j0 + c[1][1] * j1;
                    let l0 = d[0][0] * j0 + d[0][1] * j1;
                    let l1 = d[1][0] * j0 + d[1][1] * j1;
                    let e = -(k0 + 2 * ir0) * l0 - (k1 + 2 * ir1) * l1;
                    let idx = (((k0 + ir0) % 2 + 2) % 2 + 2 * (((k1 + ir1) % 2 + 2) % 2)) as usize;
                    l[idx] = l[idx].add(&zpow(e));
                }
            }
            if l.iter().any(|x| !bool::from(x.ct_is_zero())) {
                chosen = (ir0, ir1);
                break 'outer;
            }
        }
    }
    let (ir0, ir1) = chosen;

    let mut n: [[Fp2<L>; 4]; 4] = core::array::from_fn(|_| core::array::from_fn(|_| Fp2::zero()));
    for i0 in 0..2i64 {
        for i1 in 0..2i64 {
            for j0 in 0..2i64 {
                for j1 in 0..2i64 {
                    let k0 = a[0][0] * i0 + a[0][1] * i1 + c[0][0] * j0 + c[0][1] * j1;
                    let k1 = a[1][0] * i0 + a[1][1] * i1 + c[1][0] * j0 + c[1][1] * j1;
                    let l0 = b[0][0] * i0 + b[0][1] * i1 + d[0][0] * j0 + d[0][1] * j1;
                    let l1 = b[1][0] * i0 + b[1][1] * i1 + d[1][0] * j0 + d[1][1] * j1;
                    let e = i0 * j0 + i1 * j1 - (k0 + 2 * ir0) * l0 - (k1 + 2 * ir1) * l1;
                    let row = (i0 + 2 * i1) as usize;
                    let col = (((k0 + ir0) % 2 + 2) % 2 + 2 * (((k1 + ir1) % 2 + 2) % 2)) as usize;
                    n[row][col] = n[row][col].add(&zpow(e));
                }
            }
        }
    }
    n
}

/// Apply a `4×4` matrix to a 4-vector.
#[inline]
pub fn apply_mat4<L: FpBackend>(n: &[[Fp2<L>; 4]; 4], p: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
    core::array::from_fn(|i| {
        let mut acc = n[i][0].mul(&p[0]);
        for j in 1..4 {
            acc = acc.add(&n[i][j].mul(&p[j]));
        }
        acc
    })
}

/// `montgomery_to_theta_matrix_dim2(zero12, N)`: the matrix mapping product
/// Montgomery coordinates `(X1X2, X1Z2, Z1X2, Z1Z2)` to theta coordinates under
/// the base change `N`, and the resulting theta null point. `zero12` is the
/// product theta-null of `E1×E2[2]`.
pub fn montgomery_to_theta_matrix_dim2<L: FpBackend>(
    zero12: &[Fp2<L>; 4],
    n: &[[Fp2<L>; 4]; 4],
) -> ([[Fp2<L>; 4]; 4], [Fp2<L>; 4]) {
    // M[i,j] = N[i,j] * zero12[j]
    let mm: [[Fp2<L>; 4]; 4] = core::array::from_fn(|i| core::array::from_fn(|j| n[i][j].mul(&zero12[j])));
    // M2[i] = (m0+m1+m2+m3, -m0-m1+m2+m3, -m0+m1-m2+m3, m0-m1-m2+m3)
    let m2: [[Fp2<L>; 4]; 4] = core::array::from_fn(|i| {
        let r = &mm[i];
        [
            r[0].add(&r[1]).add(&r[2]).add(&r[3]),
            r[2].add(&r[3]).sub(&r[0]).sub(&r[1]),
            r[1].add(&r[3]).sub(&r[0]).sub(&r[2]),
            r[0].add(&r[3]).sub(&r[1]).sub(&r[2]),
        ]
    });
    let null: [Fp2<L>; 4] = core::array::from_fn(|i| {
        mm[i][0].add(&mm[i][1]).add(&mm[i][2]).add(&mm[i][3])
    });
    (m2, null)
}

// TuplePoint on E1 × E2.

/// A point on the product `E1 × E2` (a pair of Jacobian points).
#[derive(Clone, Debug)]
pub struct TuplePoint<L: FpBackend> {
    pub p1: JacPoint<L>,
    pub p2: JacPoint<L>,
}

impl<L: FpBackend> TuplePoint<L> {
    #[inline]
    pub fn new(p1: JacPoint<L>, p2: JacPoint<L>) -> Self {
        Self { p1, p2 }
    }
    #[inline]
    pub fn double(&self, e1: &EcCurve<L>, e2: &EcCurve<L>) -> Self {
        Self::new(jac_dbl(&self.p1, e1), jac_dbl(&self.p2, e2))
    }
    #[inline]
    pub fn double_iter(&self, n: u32, e1: &EcCurve<L>, e2: &EcCurve<L>) -> Self {
        let mut acc = self.clone();
        for _ in 0..n {
            acc = acc.double(e1, e2);
        }
        acc
    }
    #[inline]
    pub fn add(&self, other: &Self, e1: &EcCurve<L>, e2: &EcCurve<L>) -> Self {
        Self::new(jac_add(&self.p1, &other.p1, e1), jac_add(&self.p2, &other.p2, e2))
    }
}

/// `(X1·X2, X1·Z2, Z1·X2, Z1·Z2)` of a `TuplePoint`, with `(0:0)` components
/// normalised to `(1:0)` (`GluingThetaIsogenyDim2.base_change`'s prelude).
fn tuple_product_xz<L: FpBackend>(p: &TuplePoint<L>) -> [Fp2<L>; 4] {
    // (X : Z) of each component, from Jacobian (X : Y : Z) → (X/Z² : 1) - but we
    // only need projective (X:Z), and to_xz gives (X : Z²).
    let q1 = p.p1.to_xz();
    let q2 = p.p2.to_xz();
    let (mut x1, mut z1) = (q1.x, q1.z);
    let (mut x2, mut z2) = (q2.x, q2.z);
    if bool::from(x1.ct_is_zero()) && bool::from(z1.ct_is_zero()) {
        x1 = Fp2::one();
        z1 = Fp2::zero();
    }
    if bool::from(x2.ct_is_zero()) && bool::from(z2.ct_is_zero()) {
        x2 = Fp2::one();
        z2 = Fp2::zero();
    }
    [x1.mul(&x2), x1.mul(&z2), z1.mul(&x2), z1.mul(&z2)]
}

// The (2,2)-isogeny gluing step (elliptic product → theta).

/// A computed dim-2 gluing isogeny (`GluingThetaIsogenyDim2`).
pub struct GluingThetaIsogenyDim2<L: FpBackend> {
    base_change_matrix: [[Fp2<L>; 4]; 4],
    /// The product theta-null after the base change `N` (`domain_bc` null); kept
    /// for the dual (splitting) evaluation.
    null_bc: [Fp2<L>; 4],
    codomain: ThetaStructureDim2<L>,
    precomp: [Fp2<L>; 4],
    zero_idx: usize,
    t_shift: TuplePoint<L>,
    e1: EcCurve<L>,
    e2: EcCurve<L>,
}

impl<L: FpBackend> GluingThetaIsogenyDim2<L> {
    /// Build the gluing from the 8-torsion `(k1_8, k2_8)` above the kernel, the
    /// product structure null `zero12`, and the theta base change `n`.
    pub fn new(
        k1_8: &TuplePoint<L>,
        k2_8: &TuplePoint<L>,
        zero12: &[Fp2<L>; 4],
        n: &[[Fp2<L>; 4]; 4],
        e1: &EcCurve<L>,
        e2: &EcCurve<L>,
    ) -> Self {
        let (base_change_matrix, null_bc) = montgomery_to_theta_matrix_dim2(zero12, n);
        let t_shift = k1_8.double(e1, e2); // 2·K1_8 (4-torsion shift)
        let t1 = apply_mat4(&base_change_matrix, &tuple_product_xz(k1_8));
        let t2 = apply_mat4(&base_change_matrix, &tuple_product_xz(k2_8));
        let (codomain, precomp, zero_idx) = Self::special_compute_codomain(&t1, &t2);
        Self {
            base_change_matrix,
            null_bc,
            codomain,
            precomp,
            zero_idx,
            t_shift,
            e1: e1.clone(),
            e2: e2.clone(),
        }
    }

    pub fn codomain(&self) -> &ThetaStructureDim2<L> {
        &self.codomain
    }

    // The `i ^ zero_idx` index permutation mirrors the reference exactly; the
    // `0 ^ zero_idx` term is kept for that visual symmetry.
    #[allow(clippy::identity_op)]
    fn special_compute_codomain(
        t1: &[Fp2<L>; 4],
        t2: &[Fp2<L>; 4],
    ) -> (ThetaStructureDim2<L>, [Fp2<L>; 4], usize) {
        let xa = squared_theta2(t1);
        let za = squared_theta2(t2);
        let zero_idx = xa
            .iter()
            .position(|x| bool::from(x.ct_is_zero()))
            .expect("gluing: a vanishing index must exist");

        let num1 = za[1 ^ zero_idx].clone();
        let num2 = xa[2 ^ zero_idx].clone();
        let num3 = za[3 ^ zero_idx].clone();
        let num4 = xa[3 ^ zero_idx].clone();
        let den1 = crate::hd::field::inv(&num1);
        let den2 = crate::hd::field::inv(&num2);
        let den3 = crate::hd::field::inv(&num3);
        let den4 = crate::hd::field::inv(&num4);

        let mut abcd = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
        abcd[0 ^ zero_idx] = Fp2::zero();
        abcd[1 ^ zero_idx] = num1.mul(&den3);
        abcd[2 ^ zero_idx] = num2.mul(&den4);
        abcd[3 ^ zero_idx] = Fp2::one();

        let mut precomp = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
        precomp[0 ^ zero_idx] = Fp2::zero();
        precomp[1 ^ zero_idx] = den1.mul(&num3);
        precomp[2 ^ zero_idx] = den2.mul(&num4);
        precomp[3 ^ zero_idx] = Fp2::one();

        let null = hadamard2(&abcd);
        (ThetaStructureDim2::new(null), precomp, zero_idx)
    }

    #[allow(clippy::identity_op)] // `0 ^ z0` kept for symmetry with `i ^ z0`.
    fn special_image(&self, p: &[Fp2<L>; 4], translate: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
        let axby = squared_theta2(p);
        let aybx = squared_theta2(translate);
        let z0 = self.zero_idx;

        let y = axby[1 ^ z0].mul(&self.precomp[1 ^ z0]);
        let z = axby[2 ^ z0].mul(&self.precomp[2 ^ z0]);
        let t = axby[3 ^ z0].clone();

        let lam = if !bool::from(z.ct_is_zero()) {
            let zb = aybx[3 ^ z0].clone();
            z.mul(&crate::hd::field::inv(&zb))
        } else {
            let tb = aybx[2 ^ z0].mul(&self.precomp[2 ^ z0]);
            t.mul(&crate::hd::field::inv(&tb))
        };
        let xb = aybx[1 ^ z0].mul(&self.precomp[1 ^ z0]);
        let x = xb.mul(&lam);

        let mut xyzt = [Fp2::<L>::zero(), Fp2::zero(), Fp2::zero(), Fp2::zero()];
        xyzt[0 ^ z0] = x;
        xyzt[1 ^ z0] = y;
        xyzt[2 ^ z0] = z;
        xyzt[3 ^ z0] = t;
        hadamard2(&xyzt)
    }

    /// Evaluate on a `TuplePoint` → codomain theta point.
    pub fn eval(&self, p: &TuplePoint<L>) -> [Fp2<L>; 4] {
        let p_sum = p.add(&self.t_shift, &self.e1, &self.e2);
        let iso_p = apply_mat4(&self.base_change_matrix, &tuple_product_xz(p));
        let iso_p_sum = apply_mat4(&self.base_change_matrix, &tuple_product_xz(&p_sum));
        self.special_image(&iso_p, &iso_p_sum)
    }

    /// Dual (splitting) evaluation (`DualGluingThetaIsogenyDim2`): a dim-2 theta
    /// point on the codomain back to a `TuplePoint` on `E1 × E2`, via the inverse
    /// base change and a Kummer lift. The `precomputation` is the inverse of the
    /// `domain_bc` null (`inv(codomain_bc.null_dual)` up to the global Hadamard
    /// scalar, which cancels in the Kummer ratios).
    pub fn dual_eval(&self, coords: &[Fp2<L>; 4]) -> Option<TuplePoint<L>> {
        let inv_bc = inv4(&self.null_bc)?;
        let st = squared_theta2(coords);
        let img: [Fp2<L>; 4] = core::array::from_fn(|k| st[k].mul(&inv_bc[k]));
        let n_split = mat_inverse(&self.base_change_matrix)?;
        let bc = apply_mat4(&n_split, &img); // (X1X2, X1Z2, Z1X2, Z1Z2)
        let (x1x2, x1z2, z1x2, z1z2) = (bc[0].clone(), bc[1].clone(), bc[2].clone(), bc[3].clone());

        let (p1, p2) = if !bool::from(z1z2.ct_is_zero()) {
            let z2_inv = crate::hd::field::inv(&z1z2);
            let x2 = z1x2.mul(&z2_inv);
            let x1 = x1z2.mul(&z2_inv);
            (
                lift_kummer(&self.e1, &x1, &Fp2::one()),
                lift_kummer(&self.e2, &x2, &Fp2::one()),
            )
        } else if bool::from(z1x2.ct_is_zero()) && !bool::from(x1z2.ct_is_zero()) {
            let x2 = x1x2.mul(&crate::hd::field::inv(&x1z2));
            (JacPoint::identity(), lift_kummer(&self.e2, &x2, &Fp2::one()))
        } else if !bool::from(z1x2.ct_is_zero()) && bool::from(x1z2.ct_is_zero()) {
            let x1 = x1x2.mul(&crate::hd::field::inv(&z1x2));
            (lift_kummer(&self.e1, &x1, &Fp2::one()), JacPoint::identity())
        } else {
            (JacPoint::identity(), JacPoint::identity())
        };
        Some(TuplePoint::new(p1, p2))
    }
}

// The plain (2,2)-isogeny step (theta → theta), hadamard = (false, true).

/// A computed plain dim-2 `(2,2)`-isogeny (`ThetaIsogenyDim2`, the default
/// `hadamard=(False, True)` mode used throughout the chain).
pub struct ThetaIsogenyDim2<L: FpBackend> {
    codomain: ThetaStructureDim2<L>,
    precomp: [Fp2<L>; 3], // (B_inv, C_inv, D_inv)
}

impl<L: FpBackend> ThetaIsogenyDim2<L> {
    /// Build from the 8-torsion `(t1_8, t2_8)` above the kernel. Uses the full
    /// inversion path (image precomp `1/B,1/C,1/D`), which is independent of any
    /// domain precomputation.
    pub fn new(t1_8: &[Fp2<L>; 4], t2_8: &[Fp2<L>; 4]) -> Self {
        let st1 = squared_theta2(t1_8);
        let st2 = squared_theta2(t2_8);
        let (xa, xb) = (st1[0].clone(), st1[1].clone());
        let (za, tb, zc, td) = (st2[0].clone(), st2[1].clone(), st2[2].clone(), st2[3].clone());

        // One batched inversion of (xA, zA, tB, xB, zC, tD) instead of six.
        let mut inv = [
            xa.clone(),
            za.clone(),
            tb.clone(),
            xb.clone(),
            zc.clone(),
            td.clone(),
        ];
        let mut s1: [Fp2<L>; 6] = core::array::from_fn(|_| Fp2::zero());
        let mut s2: [Fp2<L>; 6] = core::array::from_fn(|_| Fp2::zero());
        crate::hd::field::batched_inv(&mut inv, &mut s1, &mut s2);
        let (xa_inv, za_inv, tb_inv, xb_inv, zc_inv, td_inv) = (
            inv[0].clone(),
            inv[1].clone(),
            inv[2].clone(),
            inv[3].clone(),
            inv[4].clone(),
            inv[5].clone(),
        );

        let big_b = xb.mul(&xa_inv);
        let big_c = zc.mul(&za_inv);
        let big_d = td.mul(&tb_inv).mul(&big_b);
        let b_inv = xb_inv.mul(&xa);
        let c_inv = zc_inv.mul(&za);
        let d_inv = td_inv.mul(&tb).mul(&b_inv);

        // codomain dual null (A,B,C,D) with A=1, then Hadamard to standard.
        let abcd = [Fp2::<L>::one(), big_b, big_c, big_d];
        let null = hadamard2(&abcd);
        Self {
            codomain: ThetaStructureDim2::new(null),
            precomp: [b_inv, c_inv, d_inv],
        }
    }

    pub fn codomain(&self) -> &ThetaStructureDim2<L> {
        &self.codomain
    }

    /// Evaluate on a theta point.
    pub fn eval(&self, p: &[Fp2<L>; 4]) -> [Fp2<L>; 4] {
        let st = squared_theta2(p);
        let img = [
            st[0].clone(),
            st[1].mul(&self.precomp[0]),
            st[2].mul(&self.precomp[1]),
            st[3].mul(&self.precomp[2]),
        ];
        hadamard2(&img)
    }
}

// The chain.

/// A computed dim-2 `(2,2)`-isogeny chain of `m` steps (`IsogenyChainDim2`):
/// gluing then `m-1` plain steps, via the optimal strategy. Evaluates
/// `TuplePoint`s on `E1×E2` to theta points on the codomain.
pub struct IsogenyChainDim2<L: FpBackend> {
    gluing: GluingThetaIsogenyDim2<L>,
    plain: Vec<ThetaIsogenyDim2<L>>,
    codomain_null: [Fp2<L>; 4],
}

impl<L: FpBackend> IsogenyChainDim2<L> {
    /// Build the chain from the kernel `(tp1, tp2)` (order `2^{m+2}`
    /// `TuplePoint`s), the product null `zero12`, the theta base change `n`, and
    /// the length `m ≥ 1`.
    pub fn new(
        tp1: &TuplePoint<L>,
        tp2: &TuplePoint<L>,
        zero12: &[Fp2<L>; 4],
        n: &[[Fp2<L>; 4]; 4],
        m: usize,
        e1: &EcCurve<L>,
        e2: &EcCurve<L>,
    ) -> Self {
        let strategy = optimised_strategy(m, 1.0);
        let mut plain: Vec<ThetaIsogenyDim2<L>> = Vec::with_capacity(m.saturating_sub(1));

        let mut strat_idx = 0usize;
        let mut level: Vec<u32> = vec![0];
        // Before the gluing, kernel elements are TuplePoints; after, they are
        // theta points. We keep one active stack and switch at the gluing.
        let mut tuple_stack: Vec<(TuplePoint<L>, TuplePoint<L>)> = vec![(tp1.clone(), tp2.clone())];
        let mut theta_stack: Vec<ThetaPair<L>> = Vec::new();
        let mut gluing: Option<GluingThetaIsogenyDim2<L>> = None;
        let mut cur: Option<ThetaStructureDim2<L>> = None;
        let mut codomain_null = zero12.clone();

        for k in 0..m {
            let mut prev: u32 = level.iter().sum();
            let target = (m - 1 - k) as u32;
            while prev != target {
                let s = strategy[strat_idx];
                level.push(s);
                prev += s;
                if gluing.is_none() {
                    let last = tuple_stack.last().unwrap();
                    let nx = (last.0.double_iter(s, e1, e2), last.1.double_iter(s, e1, e2));
                    tuple_stack.push(nx);
                } else {
                    let c = cur.as_ref().unwrap();
                    let last = theta_stack.last().unwrap();
                    theta_stack.push((c.double_iter(&last.0, s), c.double_iter(&last.1, s)));
                }
                strat_idx += 1;
            }

            if k == 0 {
                let (a, b) = tuple_stack.last().unwrap().clone();
                let g = GluingThetaIsogenyDim2::new(&a, &b, zero12, n, e1, e2);
                let mut cc = g.codomain().clone();
                cc.precompute();
                codomain_null = cc.null().clone();
                tuple_stack.pop();
                level.pop();
                theta_stack = tuple_stack.iter().map(|(x, y)| (g.eval(x), g.eval(y))).collect();
                cur = Some(cc);
                gluing = Some(g);
            } else {
                let (a, b) = theta_stack.last().unwrap().clone();
                let iso = ThetaIsogenyDim2::new(&a, &b);
                let mut cc = iso.codomain().clone();
                cc.precompute();
                codomain_null = cc.null().clone();
                theta_stack.pop();
                level.pop();
                theta_stack = theta_stack.iter().map(|(x, y)| (iso.eval(x), iso.eval(y))).collect();
                cur = Some(cc);
                plain.push(iso);
            }
        }

        IsogenyChainDim2 {
            gluing: gluing.expect("m ≥ 1, so the gluing was built"),
            plain,
            codomain_null,
        }
    }

    pub fn codomain_null(&self) -> &[Fp2<L>; 4] {
        &self.codomain_null
    }

    /// Evaluate a `TuplePoint` through the whole chain → codomain theta point.
    pub fn eval(&self, p: &TuplePoint<L>) -> [Fp2<L>; 4] {
        let mut cur = self.gluing.eval(p);
        for iso in &self.plain {
            cur = iso.eval(&cur);
        }
        cur
    }

    /// Dual evaluation (`DualIsogenyChainDim2`): a dim-2 theta point on the chain
    /// codomain back to a `TuplePoint` on `E1 × E2`. Applies the plain duals
    /// (reverse order, each `precomp = inv(domain null)`) then the gluing dual.
    pub fn dual_eval(&self, coords: &[Fp2<L>; 4]) -> Option<TuplePoint<L>> {
        let mut cur = coords.clone();
        for j in (0..self.plain.len()).rev() {
            // Forward plain[j] maps D_j -> D_{j+1}; its dual uses precomp = inv(D_j).
            let dj = if j == 0 {
                self.gluing.codomain().null()
            } else {
                self.plain[j - 1].codomain().null()
            };
            let inv = inv4(dj)?;
            let st = squared_theta2(&cur);
            let prod: [Fp2<L>; 4] = core::array::from_fn(|k| st[k].mul(&inv[k]));
            cur = hadamard2(&prod);
        }
        self.gluing.dual_eval(&cur)
    }
}
