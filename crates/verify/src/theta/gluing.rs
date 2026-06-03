use crate::ec::jacobian::jac_to_xz_add_components;
use crate::fp::{Fp2, FpBackend};

use super::basis_change::{apply_isomorphism, apply_isomorphism_general};
use super::couple::{couple_jac_to_xz, double_couple_jac_point, double_couple_point};
use super::theta_structure::{hadamard, pointwise_square, to_squared_theta};
use super::{
    BasisChangeMatrix, ThetaCoupleCurve, ThetaCoupleJacPoint, ThetaCouplePoint, ThetaGluing,
    ThetaPoint, TranslationMatrix,
};

/// Compute the theta point on the product from a couple point via the
/// Segre embedding `(P1.x*P2.x : P1.x*P2.z : P1.z*P2.x : P1.z*P2.z)`
/// followed by the gluing basis change.
#[inline]
fn base_change<L: FpBackend>(phi: &ThetaGluing<L>, t: &ThetaCouplePoint<L>) -> ThetaPoint<L> {
    let null_point = ThetaPoint {
        x: t.p1.x.mul(&t.p2.x),
        y: t.p1.x.mul(&t.p2.z),
        z: t.p1.z.mul(&t.p2.x),
        t: t.p1.z.mul(&t.p2.z),
    };
    apply_isomorphism(&phi.basis_change, &null_point)
}

/// Collect the Z-coordinate and the determinant for batched inversion
/// in the action-by-translation computation.
#[inline]
fn action_by_translation_z_and_det<L: FpBackend>(
    p4: &crate::ec::EcPoint<L>,
    p2: &crate::ec::EcPoint<L>,
) -> (Fp2<L>, Fp2<L>) {
    let z_inv = p4.z.clone();
    let det = p4.x.mul(&p2.z).sub(&p4.z.mul(&p2.x));
    (z_inv, det)
}

/// Build the 2×2 translation matrix from the 4-torsion and 2-torsion
/// points after their Z-coordinate and determinant have been inverted.
#[inline]
fn action_by_translation_compute_matrix<L: FpBackend>(
    p4: &crate::ec::EcPoint<L>,
    p2: &crate::ec::EcPoint<L>,
    z_inv: &Fp2<L>,
    det_inv: &Fp2<L>,
) -> TranslationMatrix<L> {
    // g10 = P4.x * P2.x / det - P4.x / P4.z
    let tmp = p4.x.mul(z_inv);
    let g10 = p4.x.mul(&p2.x).mul(det_inv).sub(&tmp);

    // g11 = P2.x * P4.z * det_inv
    let g11 = p2.x.mul(det_inv).mul(&p4.z);

    // g00 = -g11
    let g00 = g11.neg();

    // g01 = -(P2.z * P4.z * det_inv)
    let g01 = p2.z.mul(det_inv).mul(&p4.z).neg();

    TranslationMatrix { g00, g01, g10, g11 }
}

/// Verify that the 2-torsion points form a valid basis.
///
/// Returns `true` if K1_2 and K2_2 are independent order-2 points
/// on E1 x E2 (none zero, none equal within each component, and
/// doubling gives the identity).
#[inline]
pub fn verify_two_torsion<L: FpBackend>(
    k1_2: &ThetaCouplePoint<L>,
    k2_2: &ThetaCouplePoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> bool {
    if bool::from(k1_2.p1.is_zero() | k1_2.p2.is_zero() | k2_2.p1.is_zero() | k2_2.p2.is_zero()) {
        return false;
    }
    if bool::from(k1_2.p1.ct_equal(&k2_2.p1) | k1_2.p2.ct_equal(&k2_2.p2)) {
        return false;
    }
    let o1 = double_couple_point(k1_2, e12);
    let o2 = double_couple_point(k2_2, e12);
    if !bool::from(o1.p1.is_zero() & o1.p2.is_zero() & o2.p1.is_zero() & o2.p2.is_zero()) {
        return false;
    }
    true
}

/// Compute the four action-by-translation matrices from the 4-torsion
/// kernel generators. Returns `None` if the input does not have the
/// expected order.
#[inline]
fn action_by_translation<L: FpBackend>(
    k1_4: &ThetaCouplePoint<L>,
    k2_4: &ThetaCouplePoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> Option<[TranslationMatrix<L>; 4]> {
    let k1_2 = double_couple_point(k1_4, e12);
    let k2_2 = double_couple_point(k2_4, e12);

    if !verify_two_torsion(&k1_2, &k2_2, e12) {
        return None;
    }

    let (z0, d0) = action_by_translation_z_and_det(&k1_4.p1, &k1_2.p1);
    let (z1, d1) = action_by_translation_z_and_det(&k1_4.p2, &k1_2.p2);
    let (z2, d2) = action_by_translation_z_and_det(&k2_4.p1, &k2_2.p1);
    let (z3, d3) = action_by_translation_z_and_det(&k2_4.p2, &k2_2.p2);

    let mut inverses = [z0, z1, z2, z3, d0, d1, d2, d3];
    let mut t1 = [
        Fp2::<L>::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
    ];
    let mut t2 = [
        Fp2::<L>::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
        Fp2::zero(),
    ];
    Fp2::batched_inv(&mut inverses, &mut t1, &mut t2);

    if bool::from(inverses[0].ct_is_zero()) {
        return None;
    }

    let g0 = action_by_translation_compute_matrix(&k1_4.p1, &k1_2.p1, &inverses[0], &inverses[4]);
    let g1 = action_by_translation_compute_matrix(&k1_4.p2, &k1_2.p2, &inverses[1], &inverses[5]);
    let g2 = action_by_translation_compute_matrix(&k2_4.p1, &k2_2.p1, &inverses[2], &inverses[6]);
    let g3 = action_by_translation_compute_matrix(&k2_4.p2, &k2_2.p2, &inverses[3], &inverses[7]);

    Some([g0, g1, g2, g3])
}

/// Compute the 4×4 basis change matrix for the gluing isogeny from
/// 4-torsion kernel generators on E1 x E2.
///
/// Returns `None` if the kernel does not have the expected order.
#[inline]
pub fn gluing_change_of_basis<L: FpBackend>(
    k1_4: &ThetaCouplePoint<L>,
    k2_4: &ThetaCouplePoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> Option<BasisChangeMatrix<L>> {
    let gi = action_by_translation(k1_4, k2_4, e12)?;

    let mut m = BasisChangeMatrix::default();

    // Products of first columns: M11*M21 and M12*M22
    let t001 = gi[0].g00.mul(&gi[2].g00).add(&gi[0].g01.mul(&gi[2].g10));
    let t101 = gi[0].g10.mul(&gi[2].g00).add(&gi[0].g11.mul(&gi[2].g10));
    let t002 = gi[1].g00.mul(&gi[3].g00).add(&gi[1].g01.mul(&gi[3].g10));
    let t102 = gi[1].g10.mul(&gi[3].g00).add(&gi[1].g11.mul(&gi[3].g10));

    // Row 0: trace
    m.m[0][0] = Fp2::one()
        .add(&t001.mul(&t002))
        .add(&gi[2].g00.mul(&gi[3].g00))
        .add(&gi[0].g00.mul(&gi[1].g00));

    m.m[0][1] = t001
        .mul(&t102)
        .add(&gi[2].g00.mul(&gi[3].g10))
        .add(&gi[0].g00.mul(&gi[1].g10));

    m.m[0][2] = t101
        .mul(&t002)
        .add(&gi[2].g10.mul(&gi[3].g00))
        .add(&gi[0].g10.mul(&gi[1].g00));

    m.m[0][3] = t101
        .mul(&t102)
        .add(&gi[2].g10.mul(&gi[3].g10))
        .add(&gi[0].g10.mul(&gi[1].g10));

    // Row 1: action of (0, K2_4.P2)
    m.m[1][0] = gi[3].g00.mul(&m.m[0][0]).add(&gi[3].g01.mul(&m.m[0][1]));
    m.m[1][1] = gi[3].g10.mul(&m.m[0][0]).add(&gi[3].g11.mul(&m.m[0][1]));
    m.m[1][2] = gi[3].g00.mul(&m.m[0][2]).add(&gi[3].g01.mul(&m.m[0][3]));
    m.m[1][3] = gi[3].g10.mul(&m.m[0][2]).add(&gi[3].g11.mul(&m.m[0][3]));

    // Row 2: action of (K1_4.P1, 0)
    m.m[2][0] = gi[0].g00.mul(&m.m[0][0]).add(&gi[0].g01.mul(&m.m[0][2]));
    m.m[2][1] = gi[0].g00.mul(&m.m[0][1]).add(&gi[0].g01.mul(&m.m[0][3]));
    m.m[2][2] = gi[0].g10.mul(&m.m[0][0]).add(&gi[0].g11.mul(&m.m[0][2]));
    m.m[2][3] = gi[0].g10.mul(&m.m[0][1]).add(&gi[0].g11.mul(&m.m[0][3]));

    // Row 3: action of (K1_4.P1, K2_4.P2)
    m.m[3][0] = gi[0].g00.mul(&m.m[1][0]).add(&gi[0].g01.mul(&m.m[1][2]));
    m.m[3][1] = gi[0].g00.mul(&m.m[1][1]).add(&gi[0].g01.mul(&m.m[1][3]));
    m.m[3][2] = gi[0].g10.mul(&m.m[1][0]).add(&gi[0].g11.mul(&m.m[1][2]));
    m.m[3][3] = gi[0].g10.mul(&m.m[1][1]).add(&gi[0].g11.mul(&m.m[1][3]));

    Some(m)
}

/// Compute the gluing isogeny from an elliptic product.
///
/// Given 8-torsion generators `xy_k1_8`, `xy_k2_8` in Jacobian
/// coordinates on E1 x E2, computes the (2,2)-isogeny to a theta
/// structure. The kernel is `[4](K1_8, K2_8)`.
///
/// When `verify` is true, extra checks ensure the 4-torsion is isotropic.
/// Returns `None` if the kernel has incorrect order or the gluing is
/// malformed.
#[inline]
pub fn gluing_compute<L: FpBackend>(
    e12: &ThetaCoupleCurve<L>,
    xy_k1_8: &ThetaCoupleJacPoint<L>,
    xy_k2_8: &ThetaCoupleJacPoint<L>,
    verify: bool,
) -> Option<ThetaGluing<L>> {
    let xy_k1_4 = double_couple_jac_point(xy_k1_8, e12);
    let xy_k2_4 = double_couple_jac_point(xy_k2_8, e12);

    let k1_8 = couple_jac_to_xz(xy_k1_8);
    let k2_8 = couple_jac_to_xz(xy_k2_8);
    let k1_4 = couple_jac_to_xz(&xy_k1_4);
    let k2_4 = couple_jac_to_xz(&xy_k2_4);

    let basis_change = gluing_change_of_basis(&k1_4, &k2_4, e12)?;

    let mut out = ThetaGluing {
        domain: e12.clone(),
        xy_k1_8: xy_k1_8.clone(),
        image_k1_8: Default::default(),
        basis_change,
        precomputation: ThetaPoint::default(),
        codomain: ThetaPoint::default(),
    };

    let tt1 = {
        let t = base_change(&out, &k1_8);
        to_squared_theta(&t)
    };
    let tt2 = {
        let t = base_change(&out, &k2_8);
        to_squared_theta(&t)
    };

    // Kernel is well-formed only if TT1.t and TT2.t are zero
    if !bool::from(tt1.t.ct_is_zero() & tt2.t.ct_is_zero()) {
        return None;
    }
    // Projective factors must be nonzero
    if bool::from(
        tt1.x.ct_is_zero()
            | tt2.x.ct_is_zero()
            | tt1.y.ct_is_zero()
            | tt2.z.ct_is_zero()
            | tt1.z.ct_is_zero(),
    ) {
        return None;
    }

    // Codomain: projective factor Ax
    out.codomain.x = tt1.x.mul(&tt2.x);
    out.codomain.y = tt1.y.mul(&tt2.x);
    out.codomain.z = tt1.x.mul(&tt2.z);
    out.codomain.t = Fp2::zero();

    // Precomputation: projective factor ABCxz
    out.precomputation.x = tt1.y.mul(&tt2.z);
    out.precomputation.y = out.codomain.z.clone();
    out.precomputation.z = out.codomain.y.clone();
    out.precomputation.t = Fp2::zero();

    // Image of K1_8: phi(K1_8) = (x:x:y:y), store compact (x, y)
    out.image_k1_8.x = tt1.x.mul(&out.precomputation.x);
    out.image_k1_8.y = tt1.z.mul(&out.precomputation.z);

    if verify {
        let t1 = tt1.y.mul(&out.precomputation.y);
        if !bool::from(t1.ct_equal(&out.image_k1_8.x)) {
            return None;
        }
        let t1 = tt2.x.mul(&out.precomputation.x);
        let t2 = tt2.z.mul(&out.precomputation.z);
        if !bool::from(t2.ct_equal(&t1)) {
            return None;
        }
    }

    out.codomain = hadamard(&out.codomain);
    Some(out)
}

/// Evaluate the gluing isogeny on a point given in Jacobian coordinates.
///
/// Computes the cross-addition components of P+K1_8 and applies
/// the Segre embedding plus basis change to produce the image theta point.
#[inline]
pub fn gluing_eval_point<L: FpBackend>(
    p: &ThetaCoupleJacPoint<L>,
    phi: &ThetaGluing<L>,
) -> ThetaPoint<L> {
    let add_comp1 = jac_to_xz_add_components(&p.p1, &phi.xy_k1_8.p1, &phi.domain.e1);
    let add_comp2 = jac_to_xz_add_components(&p.p2, &phi.xy_k1_8.p2, &phi.domain.e2);

    // Build T1 and T2 from cross products of (u, v, w) components
    let u1u2 = add_comp1.u.mul(&add_comp2.u);
    let v1v2 = add_comp1.v.mul(&add_comp2.v);
    let t1_x = u1u2.add(&v1v2);
    let t1_y = add_comp1.u.mul(&add_comp2.w);
    let t1_z = add_comp1.w.mul(&add_comp2.u);
    let t1_t = add_comp1.w.mul(&add_comp2.w);

    let u1_plus_v1 = add_comp1.u.add(&add_comp1.v);
    let u2_plus_v2 = add_comp2.u.add(&add_comp2.v);
    let t2_x = u1_plus_v1.mul(&u2_plus_v2).sub(&t1_x);
    let t2_y = add_comp1.v.mul(&add_comp2.w);
    let t2_z = add_comp1.w.mul(&add_comp2.v);

    let t1 = ThetaPoint {
        x: t1_x,
        y: t1_y,
        z: t1_z,
        t: t1_t,
    };
    let t2 = ThetaPoint {
        x: t2_x,
        y: t2_y,
        z: t2_z,
        t: Fp2::zero(),
    };

    // Apply basis change, then square
    let t1 = apply_isomorphism_general(&phi.basis_change, &t1, true);
    let t2 = apply_isomorphism_general(&phi.basis_change, &t2, false);
    let t1 = pointwise_square(&t1);
    let t2 = pointwise_square(&t2);

    // Difference = theta(P+Q)*theta(P-Q)
    let diff = ThetaPoint {
        x: t1.x.sub(&t2.x),
        y: t1.y.sub(&t2.y),
        z: t1.z.sub(&t2.z),
        t: t1.t.sub(&t2.t),
    };
    let diff = hadamard(&diff);

    // Scale by inverse of imageK1_8 = (x:x:y:y), so inverse ~ (y:y:x:x)
    let image = ThetaPoint {
        x: diff.x.mul(&phi.image_k1_8.y),
        y: diff.y.mul(&phi.image_k1_8.y),
        z: diff.z.mul(&phi.image_k1_8.x),
        t: diff.t.mul(&phi.image_k1_8.x),
    };

    hadamard(&image)
}

/// Evaluate the gluing isogeny in the special case where the point
/// is known to produce a zero t-coordinate after `to_squared_theta`.
///
/// Returns `None` if the t-coordinate is unexpectedly nonzero.
#[inline]
pub fn gluing_eval_point_special_case<L: FpBackend>(
    p: &ThetaCouplePoint<L>,
    phi: &ThetaGluing<L>,
) -> Option<ThetaPoint<L>> {
    let t = base_change(phi, p);
    let t = to_squared_theta(&t);

    if !bool::from(t.t.ct_is_zero()) {
        return None;
    }

    let image = ThetaPoint {
        x: t.x.mul(&phi.precomputation.x),
        y: t.y.mul(&phi.precomputation.y),
        z: t.z.mul(&phi.precomputation.z),
        t: Fp2::zero(),
    };

    Some(hadamard(&image))
}

/// Evaluate the gluing isogeny on a pair of Jacobian-coordinate points.
#[inline]
pub fn gluing_eval_basis<L: FpBackend>(
    xy_t1: &ThetaCoupleJacPoint<L>,
    xy_t2: &ThetaCoupleJacPoint<L>,
    phi: &ThetaGluing<L>,
) -> (ThetaPoint<L>, ThetaPoint<L>) {
    let image1 = gluing_eval_point(xy_t1, phi);
    let image2 = gluing_eval_point(xy_t2, phi);
    (image1, image2)
}
