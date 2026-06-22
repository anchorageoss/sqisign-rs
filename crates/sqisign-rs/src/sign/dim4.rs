//! Dimension-4 SQIsignHD signer (Level 1, round-trip validated).
//!
//! This is the signing counterpart to the self-contained dim-4 verifier in
//! `sqisign_verify::hd`. It produces the 108-byte wire signature
//! ([`sqisign_verify::hd::SIG_WIRE_BYTES`]) that `hd_verify_*` consumes.
//!
//! # How it relates to the dim-2 signer
//!
//! The dim4 signer shares the entire *front half*
//! with [`crate::sign::sign::protocols_sign`] - commitment, challenge,
//! challenge-ideal, response-lattice - but its *back half* is strictly simpler:
//! it samples a response quaternion and encodes that quaternion's `2×2` action
//! matrix as `(a, b, c_or_d, q)`. There is **no response isogeny**: the dim-4
//! isogeny is the verifier's Kani embedding, never built here.
//!
//! This is a faithful port of the SQIsignHD reference C signer
//! (`Signature/src/sqisignhd/ref/sqisignhdx/{keygen,sign}.c`); because the
//! reference C signer and the sage verifier our `hd` module ports are a matched
//! pair, a faithful port verifies by construction (the primary validation is the
//! round-trip `sign → verify`).
//!
//! # The one place keygen differs from the dim-2 keygen
//!
//! The dim-4 verifier recovers torsion bases with the **HD basis-from-hint
//! convention** (`hd_torsion_basis_l1` / `canonical_hints_l1`, NQR tables),
//! which differs from the dim-2 `ec_curve_to_basis_2f_to_hint`. So the public
//! key's canonical basis - and therefore `mat_BAcan_to_BA0_two` - must be built
//! with the **HD** convention. A dim-4 key is thus *not* interchangeable with a
//! dim-2 key; [`dim4_keygen`] is the dim-2 keygen with the basis step swapped.
//!
//! Level 1 only (the verifier is Level-1 only); the structure generalizes once
//! the verifier does.

extern crate alloc;

use num_bigint::{BigInt, Sign};
use num_traits::{One, Zero};
use zeroize::Zeroize;

use crypto_bigint::U256;

use sqisign_verify::ec::jacobian::jac_add;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::Fp2;
use sqisign_verify::hd::{
    canonical_hints_l1, encode_signature, hd_challenge, hd_torsion_basis_l1, SIG_WIRE_BYTES,
};
use sqisign_verify::{Error, Level1, SecurityLevel};

use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::id2iso::sign_side::{
    change_of_basis_matrix_tate, dim2id2iso_arbitrary_isogeny_evaluation,
    id2iso_kernel_dlogs_to_ideal_even,
};
use crate::quaternion::algebra::{quat_alg_conj, quat_alg_make_primitive};
use crate::quaternion::dim2::{ibz_2x2_mul_mod, ibz_mat_2x2_eval, ibz_mat_2x2_inv_mod};
use crate::quaternion::dim4::ibz_mat_4x4_eval;
use crate::quaternion::ideal::quat_lideal_inter;
use crate::quaternion::intbig::{
    ibz_div, ibz_invmod, ibz_mod, ibz_pow, ibz_probab_prime, ibz_rand_interval, ibz_two_adic,
    ibz_zeroize, Ibz,
};
use crate::quaternion::lattice::{quat_lattice_conjugate_without_hnf, quat_lattice_intersect};
use crate::quaternion::lll::{quat_lattice_lll, quat_lideal_prime_norm_reduced_equivalent};
use crate::quaternion::normeq::quat_sampling_random_ideal_o0_given_norm;
use crate::quaternion::types::{
    IbzMat2x2, IbzVec2, IbzVec4, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal,
    QuatRepresentIntegerParams,
};

use crate::sign::sign::commit;

// Level-1 SQIsignHD parameters (Signature/src/precomp/ref/lvl1)

/// Challenge isogeny degree exponent `λ` (`EXPONENT_CHAL_HD`, `DEGREE_CHAL_HD = 2^λ`).
const LAMBDA: u32 = 128;
/// Response-scalar / half-torsion rescale exponent `r` (`EXPONENT_SIGN_PT_ORDER_HD`,
/// `SIGN_PT_ORDER_HD = 2^r`).
const R: u32 = 70;
/// Embedding-dimension exponent `e` (`DEGREE_RESP_HD = 2^e`); the response degree
/// `q` must satisfy `q < 2^e` and `2^e - q` prime.
const E_EMBED: u32 = 136;
/// Box half-width for the response short-vector enumeration (reference `m = 10`).
const SAMPLE_BOX_M: i64 = 10;
/// Outer retry cap (a correct keypair succeeds in a handful of iterations).
const MAX_SIGN_ITERS: usize = 256;

// key material

/// A dim-4 public key: the Montgomery coefficient `A_pk` plus the HD canonical
/// `2^f`-torsion basis hints (a deterministic function of `A_pk`).
///
/// Internal: the public compact API ([`crate::sign::compact`]) wraps this in
/// [`sqisign_verify::CompactPublicKey`].
#[derive(Clone, Debug)]
pub(crate) struct Dim4PublicKey {
    pub(crate) a_pk: Fp2<Level1>,
    pub(crate) hint_pk_p: u32,
    pub(crate) hint_pk_q: u32,
}

/// A dim-4 secret key. `mat_ba_can_to_ba0_two` is the change-of-basis matrix from
/// the **HD** canonical basis of `E_pk[2^f]` to the image of `E0`'s basis through
/// the secret isogeny (so it is consistent with the verifier's HD basis).
///
/// Internal: wrapped by [`crate::sign::compact::CompactSigningKey`].
#[derive(Clone)]
pub(crate) struct Dim4SecretKey {
    /// The public curve `E_pk`. Retained so the key is self-contained (the
    /// public key is derivable from it); the signer reads only `secret_ideal`
    /// and `mat_ba_can_to_ba0_two`.
    #[allow(dead_code)]
    pub(crate) curve: EcCurve<Level1>,
    pub(crate) secret_ideal: QuatLeftIdeal,
    pub(crate) mat_ba_can_to_ba0_two: IbzMat2x2,
}

// key generation

/// Generate a dim-4 keypair (Level 1).
///
/// Identical to the dim-2 keygen except the public-key canonical basis and hints
/// use the HD convention (see the module docs). The caller supplies a
/// cryptographic RNG.
pub(crate) fn dim4_keygen(
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> (Dim4PublicKey, Dim4SecretKey) {
    let precomp = Level1::signing_precomp();
    let f = <Level1 as SecurityLevel>::TORSION_EVEN_POWER;
    let params = QuatRepresentIntegerParams {
        algebra: &precomp.algebra,
        order: &precomp.extremal_orders[0],
        primality_test_iterations: precomp.quat_primality_num_iter,
    };

    loop {
        // 1. Sample a random secret ideal of O0 with norm SEC_DEGREE.
        let mut secret_ideal = match quat_sampling_random_ideal_o0_given_norm(
            &precomp.sec_degree,
            true,
            &params,
            None,
            rng,
        ) {
            Some(i) => i,
            None => continue,
        };

        // 2. Replace by a shorter prime-norm equivalent.
        if !quat_lideal_prime_norm_reduced_equivalent(
            &mut secret_ideal,
            &precomp.algebra,
            precomp.quat_primality_num_iter,
            precomp.quat_equiv_bound_coeff,
            rng,
        ) {
            continue;
        }

        // 3. Ideal-to-isogeny via clapotis → public curve + image of E0's basis.
        let mut curve = EcCurve::<Level1>::default();
        let mut b0_two = EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        );
        if dim2id2iso_arbitrary_isogeny_evaluation(
            &mut b0_two,
            &mut curve,
            &secret_ideal,
            &precomp,
            rng,
        )
        .is_none()
        {
            continue;
        }

        curve.normalize(); // C = 1 → curve.a is the canonical A_pk
        let a_pk = curve.a.clone();
        curve.normalize_a24();

        // 4. HD canonical basis + hints (the dim-4-specific step).
        let (hp, hq) = match canonical_hints_l1(&a_pk) {
            Some(h) => h,
            None => continue,
        };
        let (p_can, q_can) = match hd_torsion_basis_l1(&a_pk, hp, hq) {
            Some(b) => b,
            None => continue,
        };
        let pmq = jac_add(&p_can, &q_can.neg(), &curve);
        let canonical_basis = EcBasis::new(p_can.to_xz(), q_can.to_xz(), pmq.to_xz());

        // 5. Change-of-basis matrix (HD canonical basis → E0-image basis).
        let mat =
            match change_of_basis_matrix_tate(&canonical_basis, &b0_two, &mut curve, f, &precomp) {
                Some(m) => m,
                None => continue,
            };

        // Scrub the secret E0 basis image (no longer needed).
        let mut b0z = b0_two;
        b0z.zeroize();

        let pk = Dim4PublicKey {
            a_pk,
            hint_pk_p: hp,
            hint_pk_q: hq,
        };
        let sk = Dim4SecretKey {
            curve,
            secret_ideal,
            mat_ba_can_to_ba0_two: mat,
        };
        return (pk, sk);
    }
}

// response sampling (the new core)

/// `is_good_norm` (reference `sign.c`): the response degree `q` must be
/// `≡ 3 (mod 4)`, `q < 2^e`, and `2^e - q` prime (hence `≡ 1 (mod 4)`, a sum of
/// two squares - the data the verifier's Kani embedding needs).
fn is_good_norm(q: &Ibz, degree_resp: &Ibz) -> bool {
    if ibz_mod(q, &BigInt::from(4)) != BigInt::from(3) {
        return false;
    }
    if degree_resp <= q {
        return false;
    }
    let sum = degree_resp - q;
    ibz_probab_prime(&sum, 40) != 0
}

/// Sample a response quaternion `γ` in the hom lattice whose degree
/// `q = n(γ) / lattice_content` passes [`is_good_norm`].
///
/// LLL-reduce the lattice, then enumerate small integer combinations of the
/// reduced basis (box `[-m, m]^4`), computing each candidate's degree directly
/// from its reduced norm form `c0² + c1² + p·(c2² + c3²)` divided by
/// `denom² · lattice_content`. Returns the quaternion, its degree, and the
/// number of candidates tried (for diagnostics), or `None` if no good norm was
/// found within the reference's `2·(2m+1)^4` budget.
fn dim4_sample_response(
    lattice: &QuatLattice,
    lattice_content: &Ibz,
    alg: &QuatAlg,
    degree_resp: &Ibz,
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Option<(QuatAlgElem, Ibz, usize)> {
    let m = SAMPLE_BOX_M;
    let max_tries: u64 = 2 * (2 * m as u64 + 1).pow(4);

    // Reduced basis: small integer combinations are short quaternions.
    let lll = quat_lattice_lll(lattice, alg);

    let p = &alg.p;
    let denom2 = &lattice.denom * &lattice.denom;
    let divisor = &denom2 * lattice_content;
    let two_m = BigInt::from(2 * m);
    let m_big = BigInt::from(m);

    for attempt in 1..=max_tries {
        // vec ∈ [-m, m]^4
        let mut vec = IbzVec4::default();
        for i in 0..4 {
            let raw = ibz_rand_interval(rng, &Ibz::zero(), &two_m);
            vec[i] = &raw - &m_big;
        }

        // γ = lll · vec (numerator coords); degree = n(γ)/lattice_content.
        let coord = ibz_mat_4x4_eval(&lll, &vec);
        let qf = &(&(&coord[0] * &coord[0]) + &(&coord[1] * &coord[1]))
            + &(p * &(&(&coord[2] * &coord[2]) + &(&coord[3] * &coord[3])));
        if qf.is_zero() {
            continue;
        }
        let (degree, rem) = ibz_div(&qf, &divisor);
        if !rem.is_zero() {
            continue;
        }
        if is_good_norm(&degree, degree_resp) {
            let resp = QuatAlgElem {
                coord,
                denom: lattice.denom.clone(),
            };
            return Some((resp, degree, attempt as usize));
        }
    }
    None
}

// small conversions

/// A small non-negative residue (`< 2^r ≤ 2^70`) as `i128`.
fn ibz_to_i128(x: &Ibz) -> i128 {
    let (sign, digits) = x.to_u64_digits();
    let lo = *digits.first().unwrap_or(&0);
    let hi = *digits.get(1).unwrap_or(&0);
    let mag = (lo as u128) | ((hi as u128) << 64);
    let v = mag as i128;
    if sign == Sign::Minus {
        -v
    } else {
        v
    }
}

/// A non-negative value `< 2^136` as `U256`.
fn ibz_to_u256(x: &Ibz) -> U256 {
    let (_sign, digits) = x.to_u64_digits();
    let mut w = [0u64; 4];
    for (i, d) in digits.iter().take(4).enumerate() {
        w[i] = *d;
    }
    U256::from_words(w)
}

// the signer

/// Sign `message` under `(pk, sk)`, producing a 108-byte dim-4 signature.
///
/// Runs the shared front half (commitment, challenge, challenge ideal, response
/// lattice), samples the response quaternion ([`dim4_sample_response`]), then
/// encodes that quaternion's `2×2` action matrix as `(a, b, c_or_d)` plus the
/// degree `q` and the commitment curve `A_com` (with the commitment basis hints
/// packed into its spare bits). The challenge is **not** transmitted - the
/// verifier recomputes it.
pub(crate) fn dim4_sign(
    pk: &Dim4PublicKey,
    sk: &Dim4SecretKey,
    message: &[u8],
    rng: &mut (impl rand_core::RngCore + rand_core::CryptoRng),
) -> Result<[u8; SIG_WIRE_BYTES], Error> {
    let precomp = Level1::signing_precomp();
    let f = <Level1 as SecurityLevel>::TORSION_EVEN_POWER; // 248
    let two = BigInt::from(2);
    let two_f = ibz_pow(&two, f);
    let two_lambda = ibz_pow(&two, LAMBDA);
    let two_r = ibz_pow(&two, R);
    let degree_resp = ibz_pow(&two, E_EMBED); // DEGREE_RESP_HD = 2^e

    let e_pk = EcCurve::<Level1>::from_a(&pk.a_pk).ok_or(Error::InternalError)?;
    let j_pk = e_pk.j_inv();

    for _ in 0..MAX_SIGN_ITERS {
        // (0) Commitment: random secret isogeny E0 → E_com, push E0's basis.
        let mut e_com = EcCurve::<Level1>::default();
        let mut b_com0 = EcBasis::new(
            EcPoint::identity(),
            EcPoint::identity(),
            EcPoint::identity(),
        );
        let mut lideal_commit = QuatLeftIdeal::default();
        if commit(&mut e_com, &mut b_com0, &mut lideal_commit, &precomp, rng).is_none() {
            continue;
        }
        e_com.normalize(); // C = 1 → e_com.a is the canonical A_com
        let a_com = e_com.a.clone();
        e_com.normalize_a24();
        let j_com = e_com.j_inv();

        // (1) Challenge: chal = SHAKE256(j(E_com) ‖ j(E_pk) ‖ message), as int.
        let mut chal_bytes = [0u8; 32];
        hd_challenge::<Level1>(&j_com, &j_pk, message, &mut chal_bytes);
        let k = BigInt::from_bytes_le(Sign::Plus, &chal_bytes);

        // (2) Challenge ideal (norm 2^λ): pull the kernel back through the
        // secret key, reduce mod 2^λ, translate to an ideal of degree 2^λ.
        let chal_vec =
            ibz_mat_2x2_eval(&sk.mat_ba_can_to_ba0_two, &IbzVec2([Ibz::one(), k.clone()]));
        let vec2 = IbzVec2([
            ibz_mod(&chal_vec[0], &two_lambda),
            ibz_mod(&chal_vec[1], &two_lambda),
        ]);
        let lideal_chall = id2iso_kernel_dlogs_to_ideal_even(&vec2, LAMBDA, &precomp);

        // (3) Response lattice: (challenge ∩ secret) ∩ conj(commitment).
        let lideal_chall_secret = quat_lideal_inter(&lideal_chall, &sk.secret_ideal);
        let lat_commit = quat_lattice_conjugate_without_hnf(&lideal_commit.lattice);
        let lattice_hom = quat_lattice_intersect(&lideal_chall_secret.lattice, &lat_commit);
        let lattice_content = &lideal_chall_secret.norm * &lideal_commit.norm;

        // (4) Sample the response quaternion (the new core).
        let (resp_quat, q_big, _attempts) = match dim4_sample_response(
            &lattice_hom,
            &lattice_content,
            &precomp.algebra,
            &degree_resp,
            rng,
        ) {
            Some(x) => x,
            None => continue,
        };
        debug_assert_eq!(ibz_two_adic(&q_big), 0, "response degree q must be odd");

        // (5/6/7) HD canonical commitment basis + change of basis B_com0 → B_com_can.
        let (hp_com, hq_com) = match canonical_hints_l1(&a_com) {
            Some(h) => h,
            None => continue,
        };
        let (p_can, q_can) = match hd_torsion_basis_l1(&a_com, hp_com, hq_com) {
            Some(b) => b,
            None => continue,
        };
        let pmq = jac_add(&p_can, &q_can.neg(), &e_com);
        let b_com_can = EcBasis::new(p_can.to_xz(), q_can.to_xz(), pmq.to_xz());
        let mat_bcom0_to_bcom_can =
            match change_of_basis_matrix_tate(&b_com_can, &b_com0, &mut e_com, f, &precomp) {
                Some(m) => m,
                None => continue,
            };

        // (8) Response action matrix: conjugate γ, express it in O0, then act
        // on the torsion basis via the precomputed endomorphism-action matrices.
        let resp_conj = quat_alg_conj(&resp_quat);
        let (coeffs, content) =
            quat_alg_make_primitive(&resp_conj, &precomp.extremal_orders[0].order);
        debug_assert!(ibz_mod(&content, &two).is_one(), "content must be odd");
        let action = &precomp.action_matrices[0];
        let mut mat = IbzMat2x2::default();
        for i in 0..2 {
            mat.0[i][i] = &mat.0[i][i] + &coeffs[0];
            for j in 0..2 {
                let g = &(&(&action.gen2.0[i][j] * &coeffs[1])
                    + &(&action.gen3.0[i][j] * &coeffs[2]))
                    + &(&action.gen4.0[i][j] * &coeffs[3]);
                mat.0[i][j] = &mat.0[i][j] + &g;
                mat.0[i][j] = &mat.0[i][j] * &content;
                mat.0[i][j] = ibz_mod(&mat.0[i][j], &two_f);
            }
        }

        // (9) M = inv(mat_BAcan_to_BA0) · mat · mat_Bcom0_to_Bcom_can · N_sk⁻¹ (mod 2^f).
        let (inv_sk, ok) = ibz_mat_2x2_inv_mod(&sk.mat_ba_can_to_ba0_two, &two_f);
        if !ok {
            continue;
        }
        let mut m = ibz_2x2_mul_mod(&inv_sk, &mat, &two_f);
        m = ibz_2x2_mul_mod(&m, &mat_bcom0_to_bcom_can, &two_f);
        let nsk_inv = match ibz_invmod(&sk.secret_ideal.norm, &two_f) {
            Some(x) => x,
            None => continue,
        };
        for i in 0..2 {
            for j in 0..2 {
                m.0[i][j] = ibz_mod(&(&m.0[i][j] * &nsk_inv), &two_f);
            }
        }

        // (10) Encode (a, b, c_or_d): the half-torsion rescaling. `a` is the low
        // `r` bits of M[0][0]; `b` strips the challenge component (÷2^λ) from
        // M[1][0] then keeps the low `r` bits; `c_or_d` is M[0][1] (a odd) or the
        // M[1][1] analogue of `b` (a even). All land in [0, 2^r).
        let a_big = ibz_mod(&m.0[0][0], &two_r);
        let v_b = &m.0[1][0] - &(&m.0[0][0] * &k);
        let (b_q, _) = ibz_div(&v_b, &two_lambda);
        let b_big = ibz_mod(&b_q, &two_r);
        let a_odd = ibz_mod(&a_big, &two).is_one();
        let cod_big = if a_odd {
            ibz_mod(&m.0[0][1], &two_r)
        } else {
            let v_d = &m.0[1][1] - &(&m.0[0][1] * &k);
            let (d_q, _) = ibz_div(&v_d, &two_lambda);
            ibz_mod(&d_q, &two_r)
        };

        let a_i = ibz_to_i128(&a_big);
        let b_i = ibz_to_i128(&b_big);
        let cod_i = ibz_to_i128(&cod_big);
        let q_u = ibz_to_u256(&q_big);

        // Scrub secret intermediates (BigInt heap is not scrubbed; see sign.rs).
        let mut resp_z = resp_quat;
        resp_z.zeroize();
        lideal_commit.zeroize();
        let mut lc = lattice_content;
        ibz_zeroize(&mut lc);

        return encode_signature(&a_com, a_i, b_i, cod_i, &q_u, hp_com, hq_com)
            .ok_or(Error::InternalError);
    }

    Err(Error::InternalError)
}
