//!
//! Implements the verification protocol from the v2.0 spec, Section 6.

use crate::ec::basis::{ec_curve_to_basis_2f_from_hint, ec_curve_to_basis_2f_to_hint};
use crate::ec::isogeny::{ec_eval_even, ec_eval_small_chain};
use crate::ec::point::{ec_biscalar_mul, ec_dbl_iter, ec_dbl_iter_basis, ec_ladder3pt};
use crate::ec::{EcBasis, EcCurve, EcIsogEven};
use crate::fp::FpBackend;
use crate::params::SecurityLevel;
use crate::precomp::LevelPrecomp;
use crate::theta::chain::theta_chain_compute_and_eval_verify;
use crate::theta::couple::copy_bases_to_kernel;
use crate::theta::{ThetaCoupleCurve, HD_EXTRA_TORSION};
use hybrid_array::typenum::Unsigned;

use crate::hash::hash_to_challenge;
use crate::types::{PublicKey, Scalar, Signature};

// Multiprecision helpers operating on Scalar digit arrays.
// These deliberately avoid any dependency on num-bigint / GMP.

pub(crate) fn mp_compare<L: SecurityLevel>(a: &Scalar<L>, b: &Scalar<L>) -> i32 {
    let n = L::MpLimbs::USIZE;
    for i in (0..n).rev() {
        if a.digits[i] > b.digits[i] {
            return 1;
        } else if a.digits[i] < b.digits[i] {
            return -1;
        }
    }
    0
}

pub(crate) fn mp_is_even<L: SecurityLevel>(s: &Scalar<L>) -> bool {
    L::MpLimbs::USIZE != 0 && (s.digits[0] & 1) == 0
}

pub(crate) fn mp_sub_digits(c: &mut [u64], a: &[u64], b: &[u64]) {
    let mut borrow: u64 = 0;
    for i in 0..a.len() {
        let (diff, b1) = a[i].overflowing_sub(b[i]);
        let (diff2, b2) = diff.overflowing_sub(borrow);
        c[i] = diff2;
        borrow = (b1 as u64) + (b2 as u64);
    }
}

pub(crate) fn mp_mod_2exp_digits(a: &mut [u64], e: usize) {
    let q = e / 64;
    let r = e % 64;
    if q < a.len() {
        if r != 0 {
            a[q] &= (1u64 << r) - 1;
        } else {
            a[q] = 0;
        }
        for limb in a[q + 1..].iter_mut() {
            *limb = 0;
        }
    }
}

fn mp_shiftl(x: &mut [u64], shift: usize) {
    if shift == 0 || x.is_empty() {
        return;
    }
    let n = x.len();
    for i in (1..n).rev() {
        x[i] = (x[i] << shift) | (x[i - 1] >> (64 - shift));
    }
    x[0] <<= shift;
}

fn multiple_mp_shiftl(x: &mut [u64], total_shift: usize) {
    let mut remaining = total_shift;
    while remaining > 63 {
        mp_shiftl(x, 63);
        remaining -= 63;
    }
    if remaining > 0 {
        mp_shiftl(x, remaining);
    }
}

pub fn basis_from_hint<L: FpBackend + LevelPrecomp>(
    curve: &mut EcCurve<L>,
    f: u32,
    hint: u8,
) -> Option<EcBasis<L>> {
    let (basis, ok) = ec_curve_to_basis_2f_from_hint(
        curve,
        f,
        hint,
        L::basis_e0_px_bytes(),
        L::basis_e0_qx_bytes(),
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
        L::torsion_even_power(),
    )
    .ok()?;
    if ok != 1 {
        None
    } else {
        Some(basis)
    }
}

pub(crate) fn check_canonical_basis_change_matrix<L: FpBackend>(sig: &Signature<L>) -> Option<()> {
    let mut aux = Scalar::<L>::default();
    aux.digits[0] = 1;

    let shift = L::E_RSP as usize + HD_EXTRA_TORSION as usize - sig.backtracking as usize;
    multiple_mp_shiftl(aux.digits.as_mut_slice(), shift);

    for i in 0..2 {
        for j in 0..2 {
            if mp_compare::<L>(&aux, &sig.mat[i][j]) <= 0 {
                return None;
            }
        }
    }
    Some(())
}

/// Compute the challenge curve `E_chall` from the challenge coefficient
/// and backtracking count, evaluated on the public key curve.
pub fn compute_challenge_curve<L: FpBackend + LevelPrecomp>(
    chall_coeff: &Scalar<L>,
    backtracking: u8,
    e_pk: &EcCurve<L>,
    hint_pk: u8,
) -> Option<EcCurve<L>> {
    let mut phi_curve = e_pk.clone();
    let phi_length = L::F_CHR - backtracking as u32;

    let bas_ea = basis_from_hint(&mut phi_curve, L::F_CHR, hint_pk)?;

    let kernel = ec_ladder3pt(
        chall_coeff.digits.as_slice(),
        &bas_ea.p,
        &bas_ea.q,
        &bas_ea.pmq,
        &phi_curve,
    )?;

    let kernel = ec_dbl_iter(&kernel, backtracking as usize, &mut phi_curve);

    let phi = EcIsogEven {
        curve: phi_curve.clone(),
        kernel,
        length: phi_length,
    };

    let mut e_chall = phi_curve;
    ec_eval_even(&mut e_chall, &phi, &mut [])?;
    Some(e_chall)
}

fn compute_challenge_verify<L: FpBackend + LevelPrecomp>(
    sig: &Signature<L>,
    e_pk: &EcCurve<L>,
    hint_pk: u8,
) -> Option<EcCurve<L>> {
    compute_challenge_curve(&sig.chall_coeff, sig.backtracking, e_pk, hint_pk)
}

pub(crate) fn matrix_scalar_application_even_basis<L: FpBackend>(
    bas: &mut EcBasis<L>,
    e: &EcCurve<L>,
    mat: &[[Scalar<L>; 2]; 2],
    f: usize,
) -> Option<()> {
    let tmp_bas = bas.clone();

    // R = [a]P + [b]Q
    let r = ec_biscalar_mul(
        mat[0][0].digits.as_slice(),
        mat[1][0].digits.as_slice(),
        f,
        &tmp_bas,
        e,
    )?;

    // S = [c]P + [d]Q
    let s = ec_biscalar_mul(
        mat[0][1].digits.as_slice(),
        mat[1][1].digits.as_slice(),
        f,
        &tmp_bas,
        e,
    )?;

    // R - S = [a-c]P + [b-d]Q
    let nwords = L::MpLimbs::USIZE;
    let mut scalar0 = [0u64; 8];
    let mut scalar1 = [0u64; 8];
    mp_sub_digits(
        &mut scalar0[..nwords],
        mat[0][0].digits.as_slice(),
        mat[0][1].digits.as_slice(),
    );
    mp_mod_2exp_digits(&mut scalar0[..nwords], f);
    mp_sub_digits(
        &mut scalar1[..nwords],
        mat[1][0].digits.as_slice(),
        mat[1][1].digits.as_slice(),
    );
    mp_mod_2exp_digits(&mut scalar1[..nwords], f);

    let pmq = ec_biscalar_mul(&scalar0[..nwords], &scalar1[..nwords], f, &tmp_bas, e)?;
    bas.p = r;
    bas.q = s;
    bas.pmq = pmq;
    Some(())
}

pub(crate) fn verify_canonical_hint<L: FpBackend + LevelPrecomp>(
    curve: &mut EcCurve<L>,
    hint: u8,
) -> bool {
    let Ok((_, expected)) = ec_curve_to_basis_2f_to_hint(
        curve,
        L::F_CHR,
        L::basis_e0_px_bytes(),
        L::basis_e0_qx_bytes(),
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
        L::torsion_even_power(),
    ) else {
        return false;
    };
    hint == expected
}

fn challenge_and_aux_basis_verify<L: FpBackend + LevelPrecomp>(
    e_chall: &mut EcCurve<L>,
    e_aux: &mut EcCurve<L>,
    sig: &Signature<L>,
    pow_dim2_deg_resp: i32,
) -> Option<(EcBasis<L>, EcBasis<L>)> {
    let mut b_chall_can = basis_from_hint(e_chall, L::F_CHR, sig.hint_chall)?;

    let dbl_chall = L::F_CHR as usize
        - pow_dim2_deg_resp as usize
        - HD_EXTRA_TORSION as usize
        - sig.two_resp_length as usize;
    b_chall_can = ec_dbl_iter_basis(&b_chall_can, dbl_chall, e_chall);

    let b_aux_can = basis_from_hint(e_aux, L::F_CHR, sig.hint_aux)?;

    let dbl_aux = L::F_CHR as usize - pow_dim2_deg_resp as usize - HD_EXTRA_TORSION as usize;
    let b_aux_can = ec_dbl_iter_basis(&b_aux_can, dbl_aux, e_aux);

    let f = pow_dim2_deg_resp as usize + HD_EXTRA_TORSION as usize + sig.two_resp_length as usize;
    matrix_scalar_application_even_basis(&mut b_chall_can, e_chall, &sig.mat, f)?;

    Some((b_chall_can, b_aux_can))
}

/// Evaluate the small `2^r`-isogeny chain for the `two_resp_length`
/// portion of the response.
///
/// `kernel_is_q`: if true, the kernel generator is `b_chall_can.q`;
/// otherwise it is `b_chall_can.p`.
pub(crate) fn two_response_isogeny_verify_inner<L: FpBackend>(
    e_chall: &mut EcCurve<L>,
    b_chall_can: &mut EcBasis<L>,
    kernel_is_q: bool,
    two_resp_length: u8,
    pow_dim2_deg_resp: i32,
) -> Option<()> {
    let ker = if kernel_is_q {
        b_chall_can.q.clone()
    } else {
        b_chall_can.p.clone()
    };

    let mut points = [
        b_chall_can.p.clone(),
        b_chall_can.q.clone(),
        b_chall_can.pmq.clone(),
    ];

    let ker = ec_dbl_iter(
        &ker,
        pow_dim2_deg_resp as usize + HD_EXTRA_TORSION as usize,
        e_chall,
    );

    ec_eval_small_chain(e_chall, &ker, two_resp_length as i32, &mut points, false)?;

    b_chall_can.p = points[0].clone();
    b_chall_can.q = points[1].clone();
    b_chall_can.pmq = points[2].clone();
    Some(())
}

fn two_response_isogeny_verify<L: FpBackend>(
    e_chall: &mut EcCurve<L>,
    b_chall_can: &mut EcBasis<L>,
    sig: &Signature<L>,
    pow_dim2_deg_resp: i32,
) -> Option<()> {
    let kernel_is_q = mp_is_even::<L>(&sig.mat[0][0]) && mp_is_even::<L>(&sig.mat[1][0]);
    two_response_isogeny_verify_inner(
        e_chall,
        b_chall_can,
        kernel_is_q,
        sig.two_resp_length,
        pow_dim2_deg_resp,
    )
}

pub(crate) fn compute_commitment_curve_verify<L: FpBackend + LevelPrecomp>(
    b_chall_can: &EcBasis<L>,
    b_aux_can: &EcBasis<L>,
    e_chall: &EcCurve<L>,
    e_aux: &EcCurve<L>,
    pow_dim2_deg_resp: i32,
) -> Option<EcCurve<L>> {
    let mut echall_x_eaux = ThetaCoupleCurve {
        e1: e_chall.clone(),
        e2: e_aux.clone(),
    };

    let dim_two_ker = copy_bases_to_kernel(b_chall_can, b_aux_can);

    // Unreachable: callers reject pow_dim2_deg_resp <= 1. Kept for defense-in-depth.
    if pow_dim2_deg_resp == 0 {
        let mut e_chall_mut = e_chall.clone();
        e_chall_mut.normalize_curve_and_a24();
        if !bool::from(b_chall_can.is_four_torsion(&e_chall_mut)) {
            return None;
        }
        Some(echall_x_eaux.e1)
    } else {
        let codomain = theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp as u32,
            &mut echall_x_eaux,
            &dim_two_ker,
            true,
            &mut [],
        )?;
        Some(codomain.e1)
    }
}

/// Verify a standard SQIsign signature against a public key and message.
pub fn protocols_verify<L: FpBackend + LevelPrecomp>(
    pk: &PublicKey<L>,
    message: &[u8],
    sig: &Signature<L>,
) -> Result<(), crate::Error> {
    let err = || crate::Error::InvalidSignature;

    let pow_dim2_deg_resp = L::E_RSP as i32 - sig.two_resp_length as i32 - sig.backtracking as i32;

    // SECURITY: reject pow_dim2_deg_resp <= 0 because the auxiliary curve
    // e_aux is not consumed during verification when the commitment chain
    // length is zero, creating a malleability vector (breaks SUF-CMA).
    // Also reject == 1 (theta chain requires length >= 2).
    if pow_dim2_deg_resp <= 1 {
        return Err(err());
    }

    check_canonical_basis_change_matrix(sig).ok_or_else(err)?;

    if !EcCurve::<L>::verify_a(&pk.curve.a) {
        return Err(err());
    }

    let mut e_aux = EcCurve::<L>::from_a(&sig.e_aux_a).ok_or_else(err)?;

    if !verify_canonical_hint::<L>(&mut e_aux, sig.hint_aux) {
        return Err(err());
    }

    let mut e_chall = compute_challenge_verify(sig, &pk.curve, pk.hint_pk).ok_or_else(err)?;

    if !verify_canonical_hint::<L>(&mut e_chall, sig.hint_chall) {
        return Err(err());
    }

    let (mut b_chall_can, b_aux_can) =
        challenge_and_aux_basis_verify(&mut e_chall, &mut e_aux, sig, pow_dim2_deg_resp)
            .ok_or_else(err)?;

    if sig.two_resp_length > 0 {
        two_response_isogeny_verify(&mut e_chall, &mut b_chall_can, sig, pow_dim2_deg_resp)
            .ok_or_else(err)?;
    }

    let e_com = compute_commitment_curve_verify(
        &b_chall_can,
        &b_aux_can,
        &e_chall,
        &e_aux,
        pow_dim2_deg_resp,
    )
    .ok_or_else(err)?;

    let chk_chall = hash_to_challenge(pk, &e_com, message)?;

    if mp_compare::<L>(&sig.chall_coeff, &chk_chall) != 0 {
        return Err(err());
    }

    Ok(())
}
