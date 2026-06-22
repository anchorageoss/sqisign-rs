//!
//! Implements the Fiat-Shamir sigma protocol: commit, challenge, response.

use crate::id2iso::sign_precomp::{HasSigningPrecomp, SigningPrecomp};
use crate::id2iso::sign_side::{
    change_of_basis_matrix_tate, change_of_basis_matrix_tate_invert,
    dim2id2iso_arbitrary_isogeny_evaluation, id2iso_ideal_to_kernel_dlogs_even,
    id2iso_kernel_dlogs_to_ideal_even, matrix_application_even_basis,
};
use crate::quaternion::algebra::{quat_alg_conj, quat_alg_make_primitive, quat_alg_norm};
use crate::quaternion::ideal::{quat_lideal_create, quat_lideal_inter};
use crate::quaternion::intbig::{
    ibz_copy_digits, ibz_div, ibz_invmod, ibz_pow, ibz_to_digits, ibz_two_adic, ibz_zeroize, Ibz,
};
use crate::quaternion::lat_ball::quat_lattice_sample_from_ball;
use crate::quaternion::lattice::{quat_lattice_conjugate_without_hnf, quat_lattice_intersect};
use crate::quaternion::lll::quat_lideal_prime_norm_reduced_equivalent;
use crate::quaternion::normeq::quat_sampling_random_ideal_o0_given_norm;
use crate::quaternion::types::{QuatAlgElem, QuatLeftIdeal, QuatRepresentIntegerParams};
use num_traits::{One, Zero};
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::isogeny::{ec_eval_even, ec_eval_small_chain, ec_iso_eval, ec_isomorphism};
use sqisign_verify::ec::point::{
    ec_biscalar_mul, ec_dbl_iter, ec_dbl_iter_basis, ec_ladder3pt, ec_mul,
};
use sqisign_verify::ec::{EcBasis, EcCurve, EcIsogEven, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::hash::hash_to_challenge;
use sqisign_verify::theta::chain::theta_chain_compute_and_eval_randomized;
use sqisign_verify::theta::couple::{copy_bases_to_kernel, double_couple_point_iter};
use sqisign_verify::theta::{ThetaCoupleCurve, ThetaCouplePoint};
use sqisign_verify::types::{PublicKey, Signature};
use zeroize::Zeroize;

use crate::keygen::SecretKey;

const MAX_NWORDS: usize = 8;

/// Convert a `Scalar<L>` (u64 digit array) to a BigInt.
fn scalar_to_ibz(digits: &[u64]) -> Ibz {
    ibz_copy_digits(digits)
}

/// Convert a BigInt to a fixed-size u64 digit array for EC scalar ops.
fn ibz_to_scalar(x: &Ibz, nwords: usize) -> [u64; MAX_NWORDS] {
    let mut out = [0u64; MAX_NWORDS];
    let abs = if x < &Ibz::zero() { -x } else { x.clone() };
    ibz_to_digits(&abs, &mut out[..nwords]);
    out
}

pub(crate) fn commit<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    e_com: &mut EcCurve<L>,
    basis_even_com: &mut EcBasis<L>,
    lideal_com: &mut QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Option<()> {
    let params = QuatRepresentIntegerParams {
        algebra: &precomp.algebra,
        order: &precomp.extremal_orders[0],
        primality_test_iterations: precomp.quat_primality_num_iter,
    };

    let maybe_ideal =
        quat_sampling_random_ideal_o0_given_norm(&precomp.com_degree, true, &params, None, rng);
    *lideal_com = maybe_ideal?;

    if !quat_lideal_prime_norm_reduced_equivalent(
        lideal_com,
        &precomp.algebra,
        precomp.quat_primality_num_iter,
        precomp.quat_equiv_bound_coeff,
        rng,
    ) {
        return None;
    }
    dim2id2iso_arbitrary_isogeny_evaluation(basis_even_com, e_com, lideal_com, precomp, rng)
}

fn compute_challenge_ideal_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &Signature<L>,
    sk: &SecretKey<L>,
    precomp: &SigningPrecomp<L>,
) -> QuatLeftIdeal {
    use crate::quaternion::dim2::ibz_mat_2x2_eval;
    use crate::quaternion::types::IbzVec2;

    let mut vec = IbzVec2::default();
    vec[0] = Ibz::one();
    vec[1] = scalar_to_ibz(sig.chall_coeff().digits());

    vec = ibz_mat_2x2_eval(&sk.mat_ba_can_to_ba0_two, &vec);

    let lideal_chall_two = id2iso_kernel_dlogs_to_ideal_even(&vec, L::TORSION_EVEN_POWER, precomp);
    debug_assert_eq!(lideal_chall_two.norm, precomp.torsion_plus_2power);

    lideal_chall_two
}

fn sample_response(
    lattice: &crate::quaternion::types::QuatLattice,
    lattice_content: &Ibz,
    alg: &crate::quaternion::types::QuatAlg,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
    response_length: u32,
) -> Option<QuatAlgElem> {
    let bound_base = ibz_pow(&Ibz::from(2), response_length);
    let bound = (&bound_base - &Ibz::one()) * lattice_content;

    quat_lattice_sample_from_ball(lattice, alg, &bound, rng)
}

fn compute_response_quat_element<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sk: &SecretKey<L>,
    lideal_chall_two: &QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Option<(QuatAlgElem, Ibz)> {
    let lideal_chall_secret = quat_lideal_inter(lideal_chall_two, &sk.secret_ideal);

    let lat_commit = quat_lattice_conjugate_without_hnf(&lideal_commit.lattice);

    let lattice_hom_chall_to_com =
        quat_lattice_intersect(&lideal_chall_secret.lattice, &lat_commit);

    let lattice_content = &lideal_chall_secret.norm * &lideal_commit.norm;

    let resp_quat = sample_response(
        &lattice_hom_chall_to_com,
        &lattice_content,
        &precomp.algebra,
        rng,
        L::SQISIGN_RESPONSE_LENGTH,
    )?;

    Some((resp_quat, lattice_content))
}

fn compute_backtracking_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &mut Signature<L>,
    resp_quat: &mut QuatAlgElem,
    lattice_content: &mut Ibz,
    remain: &mut Ibz,
    precomp: &SigningPrecomp<L>,
) {
    let maxord_o0 = &precomp.extremal_orders[0].order;

    let (_dummy_coord, content) = quat_alg_make_primitive(resp_quat, maxord_o0);
    resp_quat.denom = &resp_quat.denom * &content;

    let backtracking = ibz_two_adic(&content);
    sig.set_backtracking(backtracking as u8);

    let two_pow = ibz_pow(&Ibz::from(2), backtracking);
    let (q, r) = ibz_div(lattice_content, &two_pow);
    *lattice_content = q;
    *remain = r;
}

#[allow(clippy::too_many_arguments)]
fn compute_random_aux_norm_and_helpers<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &mut Signature<L>,
    random_aux_norm: &mut Ibz,
    degree_resp_inv: &mut Ibz,
    remain: &mut Ibz,
    lattice_content: &Ibz,
    resp_quat: &mut QuatAlgElem,
    lideal_com_resp: &mut QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
) -> Option<u8> {
    let (degree_full_resp, norm_d) = quat_alg_norm(resp_quat, &precomp.algebra);
    debug_assert!(norm_d.is_one());

    let (degree_full_resp, r) = ibz_div(&degree_full_resp, lattice_content);
    debug_assert!(r.is_zero());

    let exp_diadic_val_full_resp = ibz_two_adic(&degree_full_resp);
    sig.set_two_resp_length(exp_diadic_val_full_resp as u8);

    let tmp = ibz_pow(&Ibz::from(2), exp_diadic_val_full_resp);
    let (degree_odd_resp, r) = ibz_div(&degree_full_resp, &tmp);
    debug_assert!(r.is_zero());

    // Conjugate resp_quat for the ideal creation
    *resp_quat = quat_alg_conj(resp_quat);

    // Create lideal_com_resp
    let norm_tmp = &lideal_commit.norm * &degree_odd_resp;
    let maxord_o0 = &precomp.extremal_orders[0].order;
    *lideal_com_resp = quat_lideal_create(resp_quat, &norm_tmp, maxord_o0, &precomp.algebra)
        .expect("invariant: commitment-response ideal has perfect-square index");

    // Compute the random aux ideal norm
    let pow_dim2_deg_resp =
        L::SQISIGN_RESPONSE_LENGTH - exp_diadic_val_full_resp - sig.backtracking() as u32;
    let two_pow_resp = ibz_pow(&Ibz::from(2), pow_dim2_deg_resp);
    *random_aux_norm = &two_pow_resp - &degree_odd_resp;

    // remain = 2^(pow_dim2_deg_resp + HD_extra_torsion)
    let hd_extra = sqisign_verify::theta::HD_EXTRA_TORSION;
    *remain = ibz_pow(&Ibz::from(2), pow_dim2_deg_resp + hd_extra);

    *degree_resp_inv = ibz_invmod(&degree_odd_resp, remain)?;

    Some(pow_dim2_deg_resp as u8)
}

fn evaluate_random_aux_isogeny_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    e_aux: &mut EcCurve<L>,
    b_aux: &mut EcBasis<L>,
    norm: &Ibz,
    lideal_com_resp: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Option<()> {
    let params = QuatRepresentIntegerParams {
        algebra: &precomp.algebra,
        order: &precomp.extremal_orders[0],
        primality_test_iterations: precomp.quat_primality_num_iter,
    };

    let lideal_aux = quat_sampling_random_ideal_o0_given_norm(
        norm,
        false,
        &params,
        Some(&precomp.quat_prime_cofactor),
        rng,
    )?;

    let lideal_aux_resp_com = quat_lideal_inter(lideal_com_resp, &lideal_aux);

    dim2id2iso_arbitrary_isogeny_evaluation(b_aux, e_aux, &lideal_aux_resp_com, precomp, rng)
}

struct ThetaCoupleCurveWithBasis<L: FpBackend + sqisign_verify::precomp::LevelPrecomp> {
    e1: EcCurve<L>,
    e2: EcCurve<L>,
    b1: EcBasis<L>,
    b2: EcBasis<L>,
}

#[allow(clippy::too_many_arguments)]
fn compute_dim2_isogeny_challenge<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    codomain: &mut ThetaCoupleCurveWithBasis<L>,
    domain: &ThetaCoupleCurveWithBasis<L>,
    degree_resp_inv: &Ibz,
    pow_dim2_deg_resp: u8,
    exp_diadic_val_full_resp: u8,
    reduced_order: u32,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Option<()> {
    let nwords = L::NWORDS_ORDER;

    let ecom_x_eaux = ThetaCoupleCurve {
        e1: domain.e1.clone(),
        e2: domain.e2.clone(),
    };

    let mut dim_two_ker = copy_bases_to_kernel(&domain.b1, &domain.b2);

    // Divide aux basis points by the degree of the response
    let scalar = ibz_to_scalar(degree_resp_inv, nwords);
    dim_two_ker.t1.p2 = ec_mul(
        &dim_two_ker.t1.p2,
        &scalar[..nwords],
        reduced_order as usize,
        &mut domain.e2.clone(),
    );
    dim_two_ker.t2.p2 = ec_mul(
        &dim_two_ker.t2.p2,
        &scalar[..nwords],
        reduced_order as usize,
        &mut domain.e2.clone(),
    );
    dim_two_ker.t1m2.p2 = ec_mul(
        &dim_two_ker.t1m2.p2,
        &scalar[..nwords],
        reduced_order as usize,
        &mut domain.e2.clone(),
    );

    // Double by exp_diadic_val_full_resp
    dim_two_ker.t1 = double_couple_point_iter(
        &dim_two_ker.t1,
        exp_diadic_val_full_resp as u32,
        &ecom_x_eaux,
    );
    dim_two_ker.t2 = double_couple_point_iter(
        &dim_two_ker.t2,
        exp_diadic_val_full_resp as u32,
        &ecom_x_eaux,
    );
    dim_two_ker.t1m2 = double_couple_point_iter(
        &dim_two_ker.t1m2,
        exp_diadic_val_full_resp as u32,
        &ecom_x_eaux,
    );

    // Points to push through: commitment basis on E1, zero on E2
    let zero_pt = EcPoint::new(Fp2::one(), Fp2::zero());
    let mut pushed_points = [
        ThetaCouplePoint {
            p1: domain.b1.p.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: domain.b1.q.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: domain.b1.pmq.clone(),
            p2: zero_pt,
        },
    ];

    let mut ecom_x_eaux_mut = ecom_x_eaux;
    let codomain_product = theta_chain_compute_and_eval_randomized(
        pow_dim2_deg_resp as u32,
        &mut ecom_x_eaux_mut,
        &dim_two_ker,
        true,
        &mut pushed_points,
        rng,
    )?;

    // CRITICAL: codomain curves are SWAPPED
    // E_aux_2 = codomain_product.E2, E_chall_2 = codomain_product.E1
    codomain.e1 = codomain_product.e2.clone();
    codomain.e2 = codomain_product.e1.clone();

    // B_aux_2 (on E_aux_2 = codomain.E1): comes from P2 component
    codomain.b1.p = pushed_points[0].p2.clone();
    codomain.b1.q = pushed_points[1].p2.clone();
    codomain.b1.pmq = pushed_points[2].p2.clone();

    // B_chall_2 (on E_chall_2 = codomain.E2): comes from P1 component
    codomain.b2.p = pushed_points[0].p1.clone();
    codomain.b2.q = pushed_points[1].p1.clone();
    codomain.b2.pmq = pushed_points[2].p1.clone();

    Some(())
}

fn compute_small_chain_isogeny_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    e_chall_2: &mut EcCurve<L>,
    b_chall_2: &mut EcBasis<L>,
    resp_quat: &QuatAlgElem,
    pow_dim2_deg_resp: u8,
    length: u8,
    precomp: &SigningPrecomp<L>,
) -> Option<()> {
    let nwords = L::NWORDS_ORDER;
    let hd_extra = sqisign_verify::theta::HD_EXTRA_TORSION;

    let two_pow = ibz_pow(&Ibz::from(2), length as u32);
    let maxord_o0 = &precomp.extremal_orders[0].order;
    let lideal_resp_two = quat_lideal_create(resp_quat, &two_pow, maxord_o0, &precomp.algebra)
        .expect("invariant: response ideal has perfect-square index");

    let vec_resp_two = id2iso_ideal_to_kernel_dlogs_even(&lideal_resp_two, precomp);

    let mut points = [
        b_chall_2.p.clone(),
        b_chall_2.q.clone(),
        b_chall_2.pmq.clone(),
    ];

    // Get down to right order
    let dbl_amount = pow_dim2_deg_resp as u32 + hd_extra;
    *b_chall_2 = ec_dbl_iter_basis(b_chall_2, dbl_amount as usize, e_chall_2);

    // Apply the vector to find the kernel
    let s0 = ibz_to_scalar(&vec_resp_two[0], nwords);
    let s1 = ibz_to_scalar(&vec_resp_two[1], nwords);
    let ker = ec_biscalar_mul(
        &s0[..nwords],
        &s1[..nwords],
        length as usize,
        b_chall_2,
        e_chall_2,
    )?;

    // Compute the isogeny and push the points
    ec_eval_small_chain(e_chall_2, &ker, length as i32, &mut points, true)?;

    b_chall_2.p = points[0].clone();
    b_chall_2.q = points[1].clone();
    b_chall_2.pmq = points[2].clone();

    Some(())
}

fn compute_challenge_codomain_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &Signature<L>,
    sk: &SecretKey<L>,
    e_chall: &mut EcCurve<L>,
    e_chall_2: &EcCurve<L>,
    b_chall_2: &mut EcBasis<L>,
) -> Option<()> {
    let bas_sk = sk.canonical_basis.clone();

    let mut phi_chall = EcIsogEven {
        curve: sk.curve.clone(),
        kernel: EcPoint::identity(),
        length: L::TORSION_EVEN_POWER - sig.backtracking() as u32,
    };

    // Compute the kernel: P + [chall_coeff]*Q
    let chall_digits = sig.chall_coeff().digits();
    let kernel = ec_ladder3pt(chall_digits, &bas_sk.p, &bas_sk.q, &bas_sk.pmq, &sk.curve)?;

    // Double kernel by backtracking amount
    phi_chall.kernel = ec_dbl_iter(&kernel, sig.backtracking() as usize, &mut sk.curve.clone());

    // Compute the codomain
    ec_eval_even(e_chall, &phi_chall, &mut [])?;

    // Apply isomorphism from E_chall_2 to E_chall
    let isom = ec_isomorphism(e_chall_2, e_chall)?;
    ec_iso_eval(&mut b_chall_2.p, &isom);
    ec_iso_eval(&mut b_chall_2.q, &isom);
    ec_iso_eval(&mut b_chall_2.pmq, &isom);

    Some(())
}

fn set_aux_curve_signature<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &mut Signature<L>,
    e_aux: &mut EcCurve<L>,
) {
    e_aux.normalize();
    sig.set_e_aux_a(e_aux.a.clone());
}

#[allow(clippy::too_many_arguments)]
fn compute_and_set_basis_change_matrix<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    sig: &mut Signature<L>,
    b_aux_2: &EcBasis<L>,
    b_chall_2: &mut EcBasis<L>,
    e_aux_2: &mut EcCurve<L>,
    e_chall: &mut EcCurve<L>,
    f: u32,
    precomp: &SigningPrecomp<L>,
) -> Option<()> {
    let nwords = L::NWORDS_ORDER;

    // Compute canonical bases
    let (b_can_chall, hint_chall) = ec_curve_to_basis_2f_to_hint(
        e_chall,
        L::TORSION_EVEN_POWER,
        precomp.basis_e0_px_bytes,
        precomp.basis_e0_qx_bytes,
        precomp.p_cofactor_for_2f,
        precomp.p_cofactor_for_2f_bitlength,
        L::TORSION_EVEN_POWER,
    )
    .expect("invariant: signing-generated curve produces a valid basis");
    sig.set_hint_chall(hint_chall);

    let (b_aux_2_can, hint_aux) = ec_curve_to_basis_2f_to_hint(
        e_aux_2,
        L::TORSION_EVEN_POWER,
        precomp.basis_e0_px_bytes,
        precomp.basis_e0_qx_bytes,
        precomp.p_cofactor_for_2f,
        precomp.p_cofactor_for_2f_bitlength,
        L::TORSION_EVEN_POWER,
    )
    .expect("invariant: signing-generated curve produces a valid basis");
    sig.set_hint_aux(hint_aux);

    // Compute change of basis: B_aux_2 -> B_aux_2_can
    let mat_baux2_to_baux2_can =
        change_of_basis_matrix_tate_invert(&b_aux_2_can, b_aux_2, e_aux_2, f, precomp)?;

    // Apply change of basis to B_chall_2
    let mut mat_clone = mat_baux2_to_baux2_can;
    matrix_application_even_basis(b_chall_2, e_chall, &mut mat_clone, f);

    // Compute change of basis: B_chall_can -> B_chall_2
    let mat_bchall_can_to_bchall =
        change_of_basis_matrix_tate(b_chall_2, &b_can_chall, e_chall, f, precomp)?;

    // Set the basis change matrix to signature
    let mat = sig.mat_mut();
    for (i, row) in mat.iter_mut().enumerate().take(2) {
        for (j, entry) in row.iter_mut().enumerate().take(2) {
            let mut digits = [0u64; MAX_NWORDS];
            ibz_to_digits(&mat_bchall_can_to_bchall.0[i][j], &mut digits[..nwords]);
            Signature::<L>::scalar_digits_mut(entry)[..nwords].copy_from_slice(&digits[..nwords]);
        }
    }
    Some(())
}

/// Execute the SQIsign signing protocol.
///
/// Produces a signature for `message` under secret key `sk` / public key `pk`.
/// The caller provides a cryptographic RNG; use the NIST DRBG for KAT testing
/// or `OsRng` for production.
pub fn protocols_sign<L: FpBackend + HasSigningPrecomp + sqisign_verify::precomp::LevelPrecomp>(
    pk: &PublicKey<L>,
    sk: &SecretKey<L>,
    message: &[u8],
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Result<Signature<L>, sqisign_verify::Error> {
    let precomp = L::signing_precomp();
    let mut sig = Signature::default();
    let hd_extra = sqisign_verify::theta::HD_EXTRA_TORSION;

    let mut lideal_commit = QuatLeftIdeal::default();
    let mut lideal_com_resp = QuatLeftIdeal::default();
    let mut resp_quat;
    let mut remain = Ibz::zero();
    let mut lattice_content;
    let mut random_aux_norm = Ibz::zero();
    let mut degree_resp_inv = Ibz::zero();

    let mut ecom_eaux = ThetaCoupleCurveWithBasis {
        e1: EcCurve::default(),
        e2: EcCurve::default(),
        b1: EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        ),
        b2: EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        ),
    };

    let mut eaux2_echall2 = ThetaCoupleCurveWithBasis {
        e1: EcCurve::default(),
        e2: EcCurve::default(),
        b1: EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        ),
        b2: EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        ),
    };

    let mut e_chall = sk.curve.clone();
    let mut pow_dim2_deg_resp: u8;
    let mut reduced_order: u32;

    loop {
        // 1. Commitment
        if commit(
            &mut ecom_eaux.e1,
            &mut ecom_eaux.b1,
            &mut lideal_commit,
            &precomp,
            rng,
        )
        .is_none()
        {
            continue;
        }

        // 2. Challenge hash
        *sig.chall_coeff_mut() = hash_to_challenge(pk, &ecom_eaux.e1, message)?;

        // 3. Challenge ideal + response quaternion
        {
            let lideal_chall_two = compute_challenge_ideal_signature(&sig, sk, &precomp);
            match compute_response_quat_element(
                sk,
                &lideal_chall_two,
                &lideal_commit,
                &precomp,
                rng,
            ) {
                Some((rq, lc)) => {
                    resp_quat = rq;
                    lattice_content = lc;
                }
                None => continue,
            }
        }

        // 4. Backtracking
        compute_backtracking_signature(
            &mut sig,
            &mut resp_quat,
            &mut lattice_content,
            &mut remain,
            &precomp,
        );

        // 5. Random aux norm and helpers
        pow_dim2_deg_resp = match compute_random_aux_norm_and_helpers(
            &mut sig,
            &mut random_aux_norm,
            &mut degree_resp_inv,
            &mut remain,
            &lattice_content,
            &mut resp_quat,
            &mut lideal_com_resp,
            &lideal_commit,
            &precomp,
        ) {
            Some(val) => val,
            None => continue,
        };

        // SECURITY: the verifier rejects pow_dim2_deg_resp == 0 because the
        // auxiliary curve is unbound in that case (breaks SUF-CMA). Retry.
        if pow_dim2_deg_resp == 0 {
            continue;
        }

        // 6. Evaluate random aux isogeny
        if evaluate_random_aux_isogeny_signature(
            &mut ecom_eaux.e2,
            &mut ecom_eaux.b2,
            &random_aux_norm,
            &lideal_com_resp,
            &precomp,
            rng,
        )
        .is_none()
        {
            continue;
        }

        // Reduce bases to the relevant order
        reduced_order = pow_dim2_deg_resp as u32 + hd_extra + sig.two_resp_length() as u32;
        ecom_eaux.b1 = ec_dbl_iter_basis(
            &ecom_eaux.b1,
            (L::TORSION_EVEN_POWER - reduced_order) as usize,
            &mut ecom_eaux.e1,
        );
        ecom_eaux.b2 = ec_dbl_iter_basis(
            &ecom_eaux.b2,
            (L::TORSION_EVEN_POWER - reduced_order) as usize,
            &mut ecom_eaux.e2,
        );

        // 7. Dim2 isogeny
        if compute_dim2_isogeny_challenge(
            &mut eaux2_echall2,
            &ecom_eaux,
            &degree_resp_inv,
            pow_dim2_deg_resp,
            sig.two_resp_length(),
            reduced_order,
            rng,
        )
        .is_none()
        {
            continue;
        }

        // 8. Small chain isogeny for remaining 2-power
        if sig.two_resp_length() > 0 {
            compute_small_chain_isogeny_signature(
                &mut eaux2_echall2.e2,
                &mut eaux2_echall2.b2,
                &resp_quat,
                pow_dim2_deg_resp,
                sig.two_resp_length(),
                &precomp,
            )
            .expect("invariant: small chain isogeny must succeed");
        }

        // 9. Challenge codomain
        compute_challenge_codomain_signature(
            &sig,
            sk,
            &mut e_chall,
            &eaux2_echall2.e2,
            &mut eaux2_echall2.b2,
        )
        .expect("invariant: challenge codomain must succeed");

        break;
    }

    // 10. Set aux curve to signature
    set_aux_curve_signature(&mut sig, &mut eaux2_echall2.e1);

    // 11. Set basis change matrix
    compute_and_set_basis_change_matrix(
        &mut sig,
        &eaux2_echall2.b1,
        &mut eaux2_echall2.b2,
        &mut eaux2_echall2.e1,
        &mut e_chall,
        reduced_order,
        &precomp,
    )
    .ok_or(sqisign_verify::Error::InternalError)?;

    // Scrub secret intermediates. BigInt heap allocations are NOT zeroed
    // here (num-bigint limitation); use ZeroizingAllocator for that.
    lideal_commit.zeroize();
    lideal_com_resp.zeroize();
    resp_quat.zeroize();
    ibz_zeroize(&mut remain);
    ibz_zeroize(&mut lattice_content);
    ibz_zeroize(&mut random_aux_norm);
    ibz_zeroize(&mut degree_resp_inv);
    ecom_eaux.b1.zeroize();
    ecom_eaux.b2.zeroize();
    eaux2_echall2.b1.zeroize();
    eaux2_echall2.b2.zeroize();

    Ok(sig)
}
