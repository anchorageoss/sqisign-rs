//! The gluing step: from a product `E1 x E2` with a kernel given by couple
//! points of order 8, the theta structure adapted to the kernel, the
//! codomain of the first (2, 2)-isogeny, and its evaluation (spec Section
//! 4.5.6). Points on the product are evaluated through barycentric
//! coordinates, since level-2 theta coordinates of a sum and a difference
//! are only available together.

use super::structure::{hadamard, invert_point, pointwise_square, to_squared_theta};
use super::{ThetaCoupleCurve, ThetaCouplePoint, ThetaKernelCouplePoints, ThetaPoint};
use crate::ec::point::{ec_points_to_bary_coordinates, xdbladd};
use crate::ec::{EcBaryCoordinates, EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;
use subtle::Choice;

/// The three entries of a symmetric 2x2 element `((s00, s01), (s10, -s00))`
/// acting on theta coordinates of one factor.
#[derive(Clone)]
struct SymmetricElement<L: FpBackend> {
    s00: Fp2<L>,
    s01: Fp2<L>,
    s10: Fp2<L>,
}

struct SymmetricElementBasis<L: FpBackend> {
    delta: Fp2<L>,
    g1: SymmetricElement<L>,
    g2: SymmetricElement<L>,
    g3: SymmetricElement<L>,
}

/// The 4x4 change of theta coordinates induced by the kernel.
#[derive(Clone, Debug)]
pub struct GluingChangeCoordMatrix<L: FpBackend> {
    /// Row-major entries.
    pub m: [[Fp2<L>; 4]; 4],
}

/// A gluing isogeny `E1 x E2 -> A`.
#[derive(Clone, Debug)]
pub struct ThetaGluing<L: FpBackend> {
    /// Inverse of the image of the first kernel generator (order 8).
    pub inv_image_k1_8: ThetaPoint<L>,
    /// Inverse of the image of the second kernel generator (order 8).
    pub inv_image_k2_8: ThetaPoint<L>,
    /// Change of coordinates from the product theta structure.
    pub matrix: GluingChangeCoordMatrix<L>,
    /// Inverse of the dual theta null point of the codomain (`t = 0`).
    pub inv_dual_theta_null: ThetaPoint<L>,
    /// Theta null point of the codomain.
    pub codomain: ThetaPoint<L>,
}

/// Whether the 2-torsion below the kernel is what a well-formed kernel gives:
/// four non-zero points of order 2 with independent components.
fn verify_two_torsion<L: FpBackend>(
    k1_2: &ThetaCouplePoint<L>,
    k2_2: &ThetaCouplePoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> bool {
    let all_two_torsion = k1_2.p1.is_two_torsion(&e12.e1)
        & k1_2.p2.is_two_torsion(&e12.e2)
        & k2_2.p1.is_two_torsion(&e12.e1)
        & k2_2.p2.is_two_torsion(&e12.e2);
    let independent = !k1_2.p1.ct_equal(&k2_2.p1) & !k1_2.p2.ct_equal(&k2_2.p2);
    bool::from(all_two_torsion & independent)
}

/// `xDBLADD` specialised to a normalised `A24` and a normalised difference
/// `P - Q = (x : 1)`.
#[inline]
fn special_xdbladd<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    pq_x: &Fp2<L>,
    a24: &EcPoint<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    let pq = EcPoint::new(pq_x.clone(), Fp2::one());
    xdbladd(p, q, &pq, a24, true)
}

/// From `P` and `[k] P - P`... precisely: given `(P, R)` with `P - R` equal to
/// `diff` (a normalised point), return `([2^pow2] P, [2^pow2] P - diff)`,
/// which for `R = 0` is `([2^pow2] P, [2^pow2 - 1] P)`. The ladder keeps the
/// difference invariant.
#[inline]
fn gluing_ladder<L: FpBackend>(
    p: &EcPoint<L>,
    r: &EcPoint<L>,
    diff_x: &Fp2<L>,
    pow2: u16,
    curve: &EcCurve<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    debug_assert!(curve.is_a24_computed_and_normalized);
    let mut r1 = p.clone();
    let mut r0 = r.clone();
    for _ in 0..pow2 {
        let (d, s) = special_xdbladd(&r1, &r0, diff_x, &curve.a24);
        r1 = d;
        r0 = s;
    }
    (r1, r0)
}

impl<L: FpBackend> GluingChangeCoordMatrix<L> {
    /// Apply the matrix to `P`, skipping the fourth column when `P.t = 0`.
    #[inline]
    fn apply(&self, p: &ThetaPoint<L>, t_not_zero: bool) -> ThetaPoint<L> {
        let row = |r: usize| {
            let mut acc = p.x.mul(&self.m[r][0]);
            acc = acc.add(&p.y.mul(&self.m[r][1]));
            acc = acc.add(&p.z.mul(&self.m[r][2]));
            if t_not_zero {
                acc = acc.add(&p.t.mul(&self.m[r][3]));
            }
            acc
        };
        ThetaPoint {
            x: row(0),
            y: row(1),
            z: row(2),
            t: row(3),
        }
    }

    /// Theta coordinates of the couple point `T` on the product, in the
    /// adapted structure: `(x1 x2 : x1 z2 : z1 x2 : z1 z2)` transformed.
    #[inline]
    pub fn couple_point_to_theta(&self, t: &ThetaCouplePoint<L>) -> ThetaPoint<L> {
        let prod = ThetaPoint {
            x: t.p1.x.mul(&t.p2.x),
            y: t.p1.x.mul(&t.p2.z),
            z: t.p2.x.mul(&t.p1.z),
            t: t.p1.z.mul(&t.p2.z),
        };
        self.apply(&prod, true)
    }
}

/// All-ones when `a == b`, as a `Choice`.
#[inline]
fn eq_choice(a: u32, b: u32) -> Choice {
    Choice::from((a == b) as u8)
}

/// The symmetric elements above `T1`, `T2` and `T1 + T2`, a basis of `E[4]`,
/// that define the theta structure on one factor (spec Algorithm 4.118).
/// Constant time in which of the six special positions the points occupy.
fn compute_symmetric_element<L: FpBackend>(
    t1_4: &EcPoint<L>,
    t2_4: &EcPoint<L>,
) -> SymmetricElementBasis<L> {
    let mut pos: u32 = 255;
    let mask = |c: Choice| 0u32.wrapping_sub(c.unwrap_u8() as u32);

    let r1 = t1_4.x.add(&t1_4.z);
    let s1 = t1_4.x.sub(&t1_4.z);
    let r2 = t2_4.x.add(&t2_4.z);
    let s2 = t2_4.x.sub(&t2_4.z);

    pos ^= pos & mask(s1.ct_is_zero());
    pos ^= (pos ^ 1) & mask(r1.ct_is_zero());
    pos ^= (pos ^ 2) & mask(s2.ct_is_zero());
    pos ^= (pos ^ 3) & mask(r2.ct_is_zero());

    // is T1 + T2 a special point?
    let s2 = r1.mul(&s2); // (x1 + z1)(x2 - z2)
    let r1 = s2.mul_by_i(Choice::from(1)); // i (x1 + z1)(x2 - z2)
    let r2 = r2.mul(&s1); // (x1 - z1)(x2 + z2)
    let s1 = r1.add(&r2);
    let s2 = r1.sub(&r2);
    pos ^= (pos ^ 4) & mask(s1.ct_is_zero());
    pos ^= (pos ^ 5) & mask(s2.ct_is_zero());

    let pos01 = eq_choice(pos >> 1, 0); // pos in {0, 1}
    let tmp1 = Fp2::select(&t1_4.x, &t2_4.x, pos01);
    let tmp2 = Fp2::select(&t1_4.z, &t2_4.z, pos01);
    let s1 = tmp1.sqr();
    let s2 = tmp2.sqr();
    let delta = s1.sub(&s2);
    let r1 = s1.add(&s2);
    let r2 = tmp1.mul(&tmp2);
    let r2 = r2.add(&r2);

    // g1 = ((-r1, r2), (-r2, r1))
    let mut g1 = SymmetricElement {
        s00: r1.neg(),
        s01: r2.clone(),
        s10: r2.neg(),
    };
    // g2 = ((0, -d), (-d, 0)), or with +d when pos is odd
    let d_signed = Fp2::select(&delta.neg(), &delta, Choice::from((pos & 1) as u8));
    let mut g2 = SymmetricElement {
        s00: Fp2::zero(),
        s01: d_signed.clone(),
        s10: d_signed,
    };
    // swap g1 and g2 when pos in {0, 1}
    g1.s00.cswap(&mut g2.s00, pos01);
    g1.s01.cswap(&mut g2.s01, pos01);
    g1.s10.cswap(&mut g2.s10, pos01);

    // g3 = ((i r2, i r1), (-i r1, -i r2)) or its negative
    let flg1 = eq_choice(pos, 0) | eq_choice(pos, 3) | eq_choice(pos, 4);
    let flg2 = eq_choice(pos, 1) | eq_choice(pos, 2) | eq_choice(pos, 5);
    let mut g3 = SymmetricElement {
        s00: r2.mul_by_i(flg1),
        s01: r1.mul_by_i(flg2),
        s10: r1.mul_by_i(flg1),
    };
    // swap g2 and g3 when pos in {4, 5}
    let pos45 = eq_choice(pos >> 1, 2);
    g2.s00.cswap(&mut g3.s00, pos45);
    g2.s01.cswap(&mut g3.s01, pos45);
    g2.s10.cswap(&mut g3.s10, pos45);

    SymmetricElementBasis { delta, g1, g2, g3 }
}

/// The change of theta coordinates from the kernel's 4-torsion.
pub fn basis_compute<L: FpBackend>(
    k1_4: &ThetaCouplePoint<L>,
    k2_4: &ThetaCouplePoint<L>,
) -> GluingChangeCoordMatrix<L> {
    let b1 = compute_symmetric_element(&k1_4.p1, &k2_4.p1);
    let b2 = compute_symmetric_element(&k1_4.p2, &k2_4.p2);
    let delta12 = b1.delta.mul(&b2.delta);
    let mut m: [[Fp2<L>; 4]; 4] = core::array::from_fn(|_| core::array::from_fn(|_| Fp2::zero()));

    // first row
    m[0][0] = delta12
        .add(&b1.g1.s00.mul(&b2.g1.s00))
        .add(&b1.g2.s00.mul(&b2.g2.s00))
        .sub(&b1.g3.s00.mul(&b2.g3.s00));
    m[0][1] = b1
        .g1
        .s00
        .mul(&b2.g1.s01)
        .add(&b1.g2.s00.mul(&b2.g2.s01))
        .sub(&b1.g3.s00.mul(&b2.g3.s01));
    m[0][2] = b1
        .g1
        .s01
        .mul(&b2.g1.s00)
        .add(&b1.g2.s01.mul(&b2.g2.s00))
        .sub(&b1.g3.s01.mul(&b2.g3.s00));
    m[0][3] = b1
        .g1
        .s01
        .mul(&b2.g1.s01)
        .add(&b1.g2.s01.mul(&b2.g2.s01))
        .sub(&b1.g3.s01.mul(&b2.g3.s01));

    // second row: action of g2 of the second factor on the first row
    m[1][0] = b2.g2.s00.mul(&m[0][0]).add(&b2.g2.s10.mul(&m[0][1]));
    m[1][1] = b2.g2.s01.mul(&m[0][0]).sub(&b2.g2.s00.mul(&m[0][1]));
    m[1][2] = b2.g2.s00.mul(&m[0][2]).add(&b2.g2.s10.mul(&m[0][3]));
    m[1][3] = b2.g2.s01.mul(&m[0][2]).sub(&b2.g2.s00.mul(&m[0][3]));

    // third row: action of g1 of the first factor on the first row
    m[2][0] = b1.g1.s00.mul(&m[0][0]).add(&b1.g1.s10.mul(&m[0][2]));
    m[2][1] = b1.g1.s00.mul(&m[0][1]).add(&b1.g1.s10.mul(&m[0][3]));
    m[2][2] = b1.g1.s01.mul(&m[0][0]).sub(&b1.g1.s00.mul(&m[0][2]));
    m[2][3] = b1.g1.s01.mul(&m[0][1]).sub(&b1.g1.s00.mul(&m[0][3]));

    // last row: action of g1 of the first factor on the second row
    m[3][0] = b1.g1.s00.mul(&m[1][0]).add(&b1.g1.s10.mul(&m[1][2]));
    m[3][1] = b1.g1.s00.mul(&m[1][1]).add(&b1.g1.s10.mul(&m[1][3]));
    m[3][2] = b1.g1.s01.mul(&m[1][0]).sub(&b1.g1.s00.mul(&m[1][2]));
    m[3][3] = b1.g1.s01.mul(&m[1][1]).sub(&b1.g1.s00.mul(&m[1][3]));

    // scalar realignment of the rows
    let (row0, rest) = m.split_first_mut().expect("four rows");
    let (row1, rest) = rest.split_first_mut().expect("three rows");
    let row2 = &mut rest[0];
    for (r0, (r1, r2)) in row0.iter_mut().zip(row1.iter_mut().zip(row2.iter_mut())) {
        *r0 = r0.mul(&delta12);
        *r1 = r1.mul(&b1.delta);
        *r2 = r2.mul(&b2.delta);
    }
    GluingChangeCoordMatrix { m }
}

/// The gluing isogeny `E1 x E2 -> A` with kernel `[4] <K1_8, K2_8>`. With
/// `dual_codomain` the codomain null point is left in dual coordinates.
/// `verify` adds the 2-torsion sanity check. `None` on a malformed kernel.
pub fn gluing_compute<L: FpBackend>(
    e12: &ThetaCoupleCurve<L>,
    k1_8: &ThetaCouplePoint<L>,
    k2_8: &ThetaCouplePoint<L>,
    dual_codomain: bool,
    verify: bool,
) -> Option<ThetaGluing<L>> {
    let k1_4 = k1_8.double(e12);
    let k2_4 = k2_8.double(e12);
    let k1_2 = k1_4.double(e12);
    let k2_2 = k2_4.double(e12);
    if verify && !verify_two_torsion(&k1_2, &k2_2, e12) {
        return None;
    }
    let matrix = basis_compute(&k1_4, &k2_4);

    let tt1 = to_squared_theta(&matrix.couple_point_to_theta(k1_8));
    let tt2 = to_squared_theta(&matrix.couple_point_to_theta(k2_8));

    // a well-formed kernel has zero t coordinates here
    if !bool::from(tt1.t.ct_is_zero() & tt2.t.ct_is_zero()) {
        return None;
    }
    let proj_zero = tt1.x.ct_is_zero()
        | tt2.x.ct_is_zero()
        | tt1.y.ct_is_zero()
        | tt2.z.ct_is_zero()
        | tt1.z.ct_is_zero();
    if bool::from(proj_zero) {
        return None;
    }

    // codomain with projective factor Ax
    let codomain = ThetaPoint {
        x: tt1.x.mul(&tt2.x),
        y: tt1.y.mul(&tt2.x),
        z: tt1.x.mul(&tt2.z),
        t: Fp2::zero(),
    };
    // inverse dual null point with projective factor ABCxz
    let inv_dual_theta_null = ThetaPoint {
        x: tt1.y.mul(&tt2.z),
        y: codomain.z.clone(),
        z: codomain.y.clone(),
        t: Fp2::zero(),
    };
    // phi(K1_8) = (x : x : y : y), inverse (y : y : x : x)
    let ik1x = tt1.z.mul(&inv_dual_theta_null.z);
    let ik1z = tt1.x.mul(&inv_dual_theta_null.x);
    let inv_image_k1_8 = ThetaPoint {
        x: ik1x.clone(),
        y: ik1x,
        z: ik1z.clone(),
        t: ik1z,
    };
    // phi(K2_8) = (z : w : z : w), inverse (w : z : w : z)
    let ik2x = tt2.y.mul(&inv_dual_theta_null.y);
    let ik2y = tt2.x.mul(&inv_dual_theta_null.x);
    let inv_image_k2_8 = ThetaPoint {
        x: ik2x.clone(),
        y: ik2y.clone(),
        z: ik2x,
        t: ik2y,
    };
    let codomain = if dual_codomain {
        codomain
    } else {
        hadamard(&codomain)
    };
    Some(ThetaGluing {
        inv_image_k1_8,
        inv_image_k2_8,
        matrix,
        inv_dual_theta_null,
        codomain,
    })
}

impl<L: FpBackend> ThetaGluing<L> {
    /// Image of a couple point given through the barycentric coordinates of
    /// `(P_i, R_i, P_i - R_i)` on each factor, divided by the image of the
    /// auxiliary point (`aux_inv` is that image's inverse), so that the
    /// result is the image of `P`.
    pub fn eval_point_bary(
        &self,
        c1: &EcBaryCoordinates<L>,
        c2: &EcBaryCoordinates<L>,
        aux_inv: &ThetaPoint<L>,
        dual_codomain: bool,
    ) -> ThetaPoint<L> {
        let mut t1 = ThetaPoint {
            x: c1.u.mul(&c2.u).add(&c1.v.mul(&c2.v)),
            y: c1.u.mul(&c2.w),
            z: c1.w.mul(&c2.u),
            t: c1.w.mul(&c2.w),
        };
        let mut t2 = ThetaPoint {
            x: c1.u.add(&c1.v).mul(&c2.u.add(&c2.v)).sub(&t1.x),
            y: c1.v.mul(&c2.w),
            z: c1.w.mul(&c2.v),
            t: Fp2::zero(),
        };
        // theta(P + Q) = M T1 - M T2, theta(P - Q) = M T1 + M T2
        t1 = pointwise_square(&self.matrix.apply(&t1, true));
        t2 = pointwise_square(&self.matrix.apply(&t2, false));
        let diff = ThetaPoint {
            x: t1.x.sub(&t2.x),
            y: t1.y.sub(&t2.y),
            z: t1.z.sub(&t2.z),
            t: t1.t.sub(&t2.t),
        };
        let h = hadamard(&diff);
        let image = ThetaPoint {
            x: h.x.mul(&aux_inv.x),
            y: h.y.mul(&aux_inv.y),
            z: h.z.mul(&aux_inv.z),
            t: h.t.mul(&aux_inv.t),
        };
        if dual_codomain {
            image
        } else {
            hadamard(&image)
        }
    }

    /// Image of a couple point with one zero component, where the dual null
    /// point's zero coordinate causes no trouble. `None` if the point is not
    /// of that shape.
    pub fn eval_point_special_case(
        &self,
        p: &ThetaCouplePoint<L>,
        dual_codomain: bool,
    ) -> Option<ThetaPoint<L>> {
        let t = to_squared_theta(&self.matrix.couple_point_to_theta(p));
        if !bool::from(t.t.ct_is_zero()) {
            return None;
        }
        let image = ThetaPoint {
            x: t.x.mul(&self.inv_dual_theta_null.x),
            y: t.y.mul(&self.inv_dual_theta_null.y),
            z: t.z.mul(&self.inv_dual_theta_null.z),
            t: Fp2::zero(),
        };
        Some(if dual_codomain {
            image
        } else {
            hadamard(&image)
        })
    }
}

/// Maximum stack depth of the chain strategy.
pub(crate) const STACK: usize = 16;

/// A couple point with both components non-zero, to be pushed through a
/// chain: `x`-only data cannot tell `(P1, P2)` from `(P1, -P2)`, so the
/// point comes with `diff = K1_8 - point` componentwise, where `K1_8 =
/// [2^(n-1)] t1` is the first kernel generator's 8-torsion multiple of an
/// `n`-step chain, computed by the caller with full coordinates from one
/// consistent choice of signs.
#[derive(Clone, Debug)]
pub struct GeneralCouplePoint<L: Prime> {
    /// The point `X = (X1, X2)`.
    pub point: ThetaCouplePoint<L>,
    /// `K1_8 - X`, componentwise, `x`-only.
    pub diff: ThetaCouplePoint<L>,
}

/// Output of the gluing stage: the gluing itself, the images of the kernel
/// stack entries, the images of the extra points, and the stack pointer
/// (the index of the entry consumed by the gluing; the caller drops it).
pub(crate) struct GluingStage<L: FpBackend> {
    pub(crate) first_step: ThetaGluing<L>,
    pub(crate) theta_q1: [ThetaPoint<L>; STACK],
    pub(crate) theta_q2: [ThetaPoint<L>; STACK],
    pub(crate) current: usize,
}

/// The gluing stage of a chain: double the kernel down to the 8-torsion
/// along a strategy, compute the gluing, and push the strategy stack and the
/// extra points through it. `todo[0]` must hold `n`; the entries are updated
/// in place. Extra points must have one zero component.
#[allow(clippy::too_many_arguments)]
pub(crate) fn start_chain<L: FpBackend>(
    todo: &mut [u16; STACK],
    e12: &ThetaCoupleCurve<L>,
    ker: &ThetaKernelCouplePoints<L>,
    points: &[ThetaCouplePoint<L>],
    general: &[GeneralCouplePoint<L>],
    pts_out: &mut [ThetaPoint<L>],
    verify: bool,
) -> Option<GluingStage<L>> {
    // the four kernel components, normalised, each paired with its running
    // "minus one" companion for the gluing ladder: (P, [2^k - 1] P)
    let mut ker_pts = [
        ker.t1.p1.clone(),
        ker.t1.p2.clone(),
        ker.t2.p1.clone(),
        ker.t2.p2.clone(),
    ];
    {
        let mut zs = [
            ker_pts[0].z.clone(),
            ker_pts[1].z.clone(),
            ker_pts[2].z.clone(),
            ker_pts[3].z.clone(),
        ];
        let mut s1: [Fp2<L>; 4] = core::array::from_fn(|_| Fp2::zero());
        let mut s2: [Fp2<L>; 4] = core::array::from_fn(|_| Fp2::zero());
        Fp2::batched_inv(&mut zs, &mut s1, &mut s2);
        for (p, zi) in ker_pts.iter_mut().zip(zs.iter()) {
            p.x = p.x.mul(zi);
            p.z = Fp2::one();
        }
    }
    // t_11, t_21 on E1; t_12, t_22 on E2; entry k holds ([2^s] P, [2^s - 1] P)
    let mut t11: [(EcPoint<L>, EcPoint<L>); STACK] =
        core::array::from_fn(|_| (EcPoint::identity(), EcPoint::identity()));
    let mut t12 = t11.clone();
    let mut t21 = t11.clone();
    let mut t22 = t11.clone();
    t11[0].0 = ker_pts[0].clone();
    t12[0].0 = ker_pts[1].clone();
    t21[0].0 = ker_pts[2].clone();
    t22[0].0 = ker_pts[3].clone();

    let mut current = 0usize;
    while todo[current] != 1 {
        debug_assert!(todo[current] >= 2);
        current += 1;
        if current >= STACK {
            return None;
        }
        // the gluing is expensive, so near the end recompute doublings
        // rather than pushing intermediate points
        let prev = todo[current - 1];
        let num_dbls = if prev >= 16 { prev / 2 } else { prev - 1 };
        let step = |arr: &mut [(EcPoint<L>, EcPoint<L>); STACK], e: &EcCurve<L>| {
            let (p, r) = (arr[current - 1].0.clone(), arr[current - 1].1.clone());
            let diff_x = arr[0].0.x.clone();
            arr[current] = gluing_ladder(&p, &r, &diff_x, num_dbls, e);
        };
        step(&mut t11, &e12.e1);
        step(&mut t12, &e12.e2);
        step(&mut t21, &e12.e1);
        step(&mut t22, &e12.e2);
        todo[current] = prev - num_dbls;
    }
    debug_assert_eq!(todo[current], 1);

    let k1_8 = ThetaCouplePoint::new(t11[current].0.clone(), t12[current].0.clone());
    let k2_8 = ThetaCouplePoint::new(t21[current].0.clone(), t22[current].0.clone());
    let first_step = gluing_compute(e12, &k1_8, &k2_8, false, verify)?;

    for (out, p) in pts_out.iter_mut().zip(points.iter()) {
        if !bool::from(p.p1.is_zero() | p.p2.is_zero()) {
            return None;
        }
        *out = first_step.eval_point_special_case(p, false)?;
    }
    // general couple points: through the barycentric coordinates of
    // (K1_8, X, K1_8 - X) on each factor, exactly as the kernel generator
    // itself is pushed below; the caller supplies K1_8 - X, which carries the
    // relative sign of the components that x-only data cannot
    for (out, g) in pts_out[points.len()..].iter_mut().zip(general.iter()) {
        if bool::from(g.point.p1.is_zero() | g.point.p2.is_zero()) {
            return None;
        }
        let b1 = ec_points_to_bary_coordinates(&k1_8.p1, &g.point.p1, &g.diff.p1);
        let b2 = ec_points_to_bary_coordinates(&k1_8.p2, &g.point.p2, &g.diff.p2);
        let img = first_step.eval_point_bary(&b1, &b2, &first_step.inv_image_k1_8, true);
        *out = hadamard(&img);
    }

    let mut theta_q1: [ThetaPoint<L>; STACK] = core::array::from_fn(|_| ThetaPoint::zero());
    let mut theta_q2: [ThetaPoint<L>; STACK] = core::array::from_fn(|_| ThetaPoint::zero());

    // the kernel generators themselves (index 0): via [2^(k)] P, P, [2^k - 1] P
    let bary = |arr: &[(EcPoint<L>, EcPoint<L>); STACK], k: usize| {
        ec_points_to_bary_coordinates(&arr[k].0, &arr[0].0, &arr[k].1)
    };
    theta_q1[0] = first_step.eval_point_bary(
        &bary(&t11, current),
        &bary(&t12, current),
        &first_step.inv_image_k1_8,
        true,
    );
    theta_q2[0] = first_step.eval_point_bary(
        &bary(&t21, current),
        &bary(&t22, current),
        &first_step.inv_image_k2_8,
        true,
    );
    todo[0] -= 1;
    let inv_phi_k1 = invert_point(&theta_q1[0]);
    let inv_phi_k2 = invert_point(&theta_q2[0]);
    theta_q1[0] = hadamard(&theta_q1[0]);
    theta_q2[0] = hadamard(&theta_q2[0]);

    for j in 1..current {
        theta_q1[j] =
            first_step.eval_point_bary(&bary(&t11, j), &bary(&t12, j), &inv_phi_k1, false);
        theta_q2[j] =
            first_step.eval_point_bary(&bary(&t21, j), &bary(&t22, j), &inv_phi_k2, false);
        todo[j] -= 1;
    }
    Some(GluingStage {
        first_step,
        theta_q1,
        theta_q2,
        current,
    })
}
