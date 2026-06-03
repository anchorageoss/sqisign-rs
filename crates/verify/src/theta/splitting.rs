//! structure that represents a product of two elliptic curves.

use crate::fp::{Fp2, FpBackend};
use crate::precomp::LevelPrecomp;
use rand_core::RngCore;
use subtle::{Choice, ConstantTimeEq};

use super::basis_change::{
    apply_isomorphism, base_change_matrix_multiplication, choose_index_theta_point,
    select_base_change_matrix, set_base_change_matrix_from_precomp,
};
use super::theta_structure::is_product_theta_point;
use super::{
    BasisChangeMatrix, EcCurve, PrecompBasisChangeMatrix, ThetaCouplePoint, ThetaPoint,
    ThetaSplitting, ThetaStructure,
};

#[inline]
fn fp2_constants<L: FpBackend>() -> [Fp2<L>; 5] {
    [
        Fp2::zero(),
        Fp2::one(),
        Fp2::i_element(),
        Fp2::<L>::one().neg(),
        Fp2::<L>::i_element().neg(),
    ]
}

/// Find a splitting isomorphism from a theta structure.
///
/// Searches the 10 candidate splitting transforms for one that sends the
/// null point to product form (`x*t == y*z`). Returns `None` if
/// not exactly one candidate produces a zero (indicating a valid
/// splitting), or if `zero_index >= 0` and the expected index did
/// not produce a zero.
///
/// When `randomize` is true (signing only), applies a random normalization
/// transform to hide the splitting choice.
#[inline]
pub fn splitting_compute<L: FpBackend + LevelPrecomp>(
    a: &ThetaStructure<L>,
    zero_index: i32,
    randomize: bool,
    rng: Option<&mut dyn RngCore>,
) -> Option<ThetaSplitting<L>> {
    let splitting_transforms = L::splitting_transforms();
    let normalization_transforms = L::normalization_transforms();
    let chi_eval = L::chi_eval();
    let even_index = L::even_index();

    let fc = fp2_constants::<L>();
    let mut out_m = BasisChangeMatrix::<L>::default();
    let mut count: u32 = 0;

    for i in 0..10 {
        let mut u_cst = Fp2::<L>::zero();

        #[allow(clippy::needless_range_loop)]
        for t in 0..4usize {
            let t2_val = choose_index_theta_point(t, &a.null_point);
            let t1_val =
                choose_index_theta_point((t ^ even_index[i][1] as usize) & 3, &a.null_point);
            let prod = t1_val.mul(&t2_val);

            // CHI_EVAL value is +1 or -1; ctl is clear for +1, set for -1
            let ctl = Choice::from(((chi_eval[even_index[i][0] as usize][t] >> 1) & 1) as u8);

            let neg_prod = prod.neg();
            let selected = Fp2::select(&prod, &neg_prod, ctl);
            u_cst = u_cst.add(&selected);
        }

        let ctl = u_cst.ct_is_zero();
        count = count.wrapping_add(ctl.unwrap_u8() as u32);

        let precomp = PrecompBasisChangeMatrix {
            m: splitting_transforms[i],
        };
        out_m = select_base_change_matrix(&out_m, &precomp, &fc, ctl);

        if zero_index != -1 && i as i32 == zero_index && !bool::from(ctl) {
            return None;
        }
    }

    if randomize {
        if let Some(rng) = rng {
            let secret_index = sample_random_index(rng);
            let mut m_random = set_base_change_matrix_from_precomp(
                &PrecompBasisChangeMatrix {
                    m: normalization_transforms[0],
                },
                &fc,
            );
            for i in 1u8..6u8 {
                let ctl = i.ct_eq(&secret_index);
                let precomp_i = PrecompBasisChangeMatrix {
                    m: normalization_transforms[i as usize],
                };
                m_random = select_base_change_matrix(&m_random, &precomp_i, &fc, ctl);
            }
            out_m = base_change_matrix_multiplication(&m_random, &out_m);
        }
    }

    let null_point = apply_isomorphism(&out_m, &a.null_point);
    let out = ThetaSplitting {
        basis_change: out_m,
        b: ThetaStructure {
            null_point,
            precomputation: false,
            ..ThetaStructure::default()
        },
    };

    if count == 1 {
        Some(out)
    } else {
        None
    }
}

#[inline]
fn sample_random_index(rng: &mut dyn RngCore) -> u8 {
    let mut seed_arr = [0u8; 4];
    loop {
        rng.fill_bytes(&mut seed_arr);
        let seed = u32::from_le_bytes(seed_arr);
        if seed < 4294967292u32 {
            let secret_index = seed.wrapping_sub(((seed as u64 * 2863311531u64) >> 34) as u32 * 6);
            return secret_index as u8;
        }
    }
}

/// Convert a product theta structure to a pair of Montgomery curves.
///
/// Given a theta structure whose null point is in product form
/// (`x*t == y*z`, after splitting), recovers two Montgomery curves:
///
/// - `E1: A1 = -2(x^4 + z^4) / (x^4 - z^4)`
/// - `E2: A2 = -2(x^4 + y^4) / (x^4 - y^4)`
///
/// Returns `None` if the null point is not in product form or has zero
/// denominators.
#[inline]
pub fn theta_product_structure_to_elliptic_product<L: FpBackend>(
    a: &ThetaStructure<L>,
) -> Option<(EcCurve<L>, EcCurve<L>)> {
    if !bool::from(is_product_theta_point(&a.null_point)) {
        return None;
    }

    let np = &a.null_point;

    if bool::from(np.x.ct_is_zero() | np.y.ct_is_zero() | np.z.ct_is_zero()) {
        return None;
    }

    // E2: A2 = -2(x^4 + y^4) / (x^4 - y^4), C2 = x^4 - y^4
    let xx = np.x.sqr().sqr(); // x^4
    let yy = np.y.sqr().sqr(); // y^4

    let a2 = xx.add(&yy).add(&xx.add(&yy)).neg(); // -2(x^4+y^4)
    let c2 = xx.sub(&yy);

    // E1: A1 = -2(x^4 + z^4) / (x^4 - z^4), C1 = x^4 - z^4
    let xx = np.x.sqr().sqr();
    let zz = np.z.sqr().sqr(); // z^4

    let a1 = xx.add(&zz).add(&xx.add(&zz)).neg(); // -2(x^4+z^4)
    let c1 = xx.sub(&zz);

    if bool::from(c1.ct_is_zero() | c2.ct_is_zero()) {
        return None;
    }

    let e1 = EcCurve::<L> {
        a: a1,
        c: c1,
        ..EcCurve::default()
    };
    let e2 = EcCurve::<L> {
        a: a2,
        c: c2,
        ..EcCurve::default()
    };

    Some((e1, e2))
}

/// Convert a theta point to a pair of Montgomery points.
///
/// Given a theta point `P` in product form and the corresponding theta
/// structure `A`, recovers `(P1, P2)` where `P1` is on `E1` and `P2` is
/// on `E2`.
///
/// Returns `None` if the point is not in product form or both selected
/// coordinate pairs are zero.
#[inline]
pub fn theta_point_to_montgomery_point<L: FpBackend>(
    p: &ThetaPoint<L>,
    a: &ThetaStructure<L>,
) -> Option<ThetaCouplePoint<L>> {
    if !bool::from(is_product_theta_point(p)) {
        return None;
    }

    let (x, z) = if bool::from(p.x.ct_is_zero() & p.y.ct_is_zero()) {
        (&p.z, &p.t)
    } else {
        (&p.x, &p.y)
    };

    if bool::from(x.ct_is_zero() & z.ct_is_zero()) {
        return None;
    }

    // P2.X = a.null.y * x + a.null.x * z
    // P2.Z = -a.null.y * x + a.null.x * z
    let ay_x = a.null_point.y.mul(x);
    let ax_z = a.null_point.x.mul(z);
    let p2_x = ay_x.add(&ax_z);
    let p2_z = ax_z.sub(&ay_x);

    let (x2, z2) = if bool::from(p.x.ct_is_zero() & p.z.ct_is_zero()) {
        (&p.y, &p.t)
    } else {
        (&p.x, &p.z)
    };

    // P1.X = a.null.z * x2 + a.null.x * z2
    // P1.Z = -a.null.z * x2 + a.null.x * z2
    let az_x = a.null_point.z.mul(x2);
    let ax_z2 = a.null_point.x.mul(z2);
    let p1_x = az_x.add(&ax_z2);
    let p1_z = ax_z2.sub(&az_x);

    Some(ThetaCouplePoint {
        p1: crate::ec::EcPoint { x: p1_x, z: p1_z },
        p2: crate::ec::EcPoint { x: p2_x, z: p2_z },
    })
}
