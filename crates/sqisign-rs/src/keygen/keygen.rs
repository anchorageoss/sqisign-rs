//!
//! Implements the keygen protocol: sample a random ideal of prescribed
//! norm from O₀, reduce it, translate to an isogeny, then compute the
//! canonical torsion basis and basis-change matrix on the codomain.

use crate::id2iso::sign_precomp::SigningPrecomp;
use crate::id2iso::sign_side::{
    change_of_basis_matrix_tate, dim2id2iso_arbitrary_isogeny_evaluation,
};
use crate::quaternion::lll::quat_lideal_prime_norm_reduced_equivalent;
use crate::quaternion::normeq::quat_sampling_random_ideal_o0_given_norm;
use crate::quaternion::types::QuatRepresentIntegerParams;
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::point::test_basis_order_twof;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::PublicKey;
use zeroize::Zeroize;

use crate::SecretKey;

/// Generate a fresh SQIsign keypair.
///
/// The caller provides a cryptographic RNG. In production use `OsRng`
/// or `thread_rng()`; for KAT testing, pass the NIST AES-256-CTR-DRBG.
pub fn protocols_keygen<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
    precomp: &SigningPrecomp<L>,
) -> (PublicKey<L>, SecretKey<L>) {
    let torsion_even_power = L::TORSION_EVEN_POWER as usize;

    let mut secret_ideal;
    let mut curve = EcCurve::<L>::default();
    let mut b0_two = EcBasis::new(
        EcPoint::identity(),
        EcPoint::identity(),
        EcPoint::identity(),
    );

    let params = QuatRepresentIntegerParams {
        algebra: &precomp.algebra,
        order: &precomp.extremal_orders[0],
        primality_test_iterations: precomp.quat_primality_num_iter,
    };

    loop {
        // 1. Sample a random ideal of O₀ with norm SEC_DEGREE
        secret_ideal = match quat_sampling_random_ideal_o0_given_norm(
            &precomp.sec_degree,
            true,
            &params,
            None,
            rng,
        ) {
            Some(ideal) => ideal,
            None => continue,
        };

        // 2. Replace by a shorter equivalent ideal with prime norm
        if !quat_lideal_prime_norm_reduced_equivalent(
            &mut secret_ideal,
            &precomp.algebra,
            precomp.quat_primality_num_iter,
            precomp.quat_equiv_bound_coeff,
            rng,
        ) {
            continue;
        }

        // 3. Ideal-to-isogeny via clapotis
        if dim2id2iso_arbitrary_isogeny_evaluation(
            &mut b0_two,
            &mut curve,
            &secret_ideal,
            precomp,
            rng,
        )
        .is_none()
        {
            continue;
        }

        break;
    }

    debug_assert!(bool::from(test_basis_order_twof(
        &b0_two,
        &curve,
        torsion_even_power
    )));

    // 4. Compute a deterministic basis with a hint for fast re-derivation
    let (canonical_basis, hint_pk) = ec_curve_to_basis_2f_to_hint(
        &mut curve,
        L::TORSION_EVEN_POWER,
        L::basis_e0_px_bytes(),
        L::basis_e0_qx_bytes(),
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
        L::torsion_even_power(),
    )
    .expect("invariant: keygen-generated curve produces a valid basis");

    debug_assert!(bool::from(test_basis_order_twof(
        &canonical_basis,
        &curve,
        torsion_even_power
    )));

    // 5. Compute the basis-change matrix (Tate pairing DLP)
    let mat_ba_can_to_ba0_two = change_of_basis_matrix_tate(
        &canonical_basis,
        &b0_two,
        &mut curve,
        L::TORSION_EVEN_POWER,
        precomp,
    )
    .expect("invariant: Tate pairing DLP must succeed on valid keygen basis");

    // 6. Build public key
    let mut pk_curve = curve.clone();
    pk_curve.is_a24_computed_and_normalized = false;
    debug_assert!(bool::from(Fp2::ct_is_one(&pk_curve.c)));

    let pk = PublicKey::new(pk_curve, hint_pk);

    // Scrub the secret E0 basis image (no longer needed).
    b0_two.zeroize();

    let sk = SecretKey {
        curve,
        secret_ideal,
        mat_ba_can_to_ba0_two,
        canonical_basis,
    };

    (pk, sk)
}
