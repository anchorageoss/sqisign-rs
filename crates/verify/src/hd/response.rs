//! Stage 3 of the compact verification: the response images (the library's
//! `image_response`). With `R_com = 2^(f − r) P_com`, `S_com = 2^(f − r)
//! Q_com` of order `2^r`, `w_com = e(R_com, S_com)`, `k = log_{w_chall}
//! w_com`, the dropped scalar from `a d − b c ≡ k q (mod 2^r)`, and the
//! images `σ(R_com) = a P_resc + b Q_resc`, `σ(S_com) = c P_resc + d Q_resc`.

use crate::ec::jacobian::jac_add;
use crate::ec::pairing::{fp2_dlog_2e, weil};
use crate::ec::{EcCurve, JacPoint, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};

use super::basis::hd_torsion_basis;
use super::challenge::{jac_dbl_iter, jac_scalar_mul, ChallengeRecovery};
use super::hd_verify::recover_response_cd;
use super::params::HdLevel;

/// The response scalars of a signature, `q` reduced modulo `2^128` (only
/// `q mod 2^r` enters here).
#[derive(Clone, Copy, Debug)]
pub struct ResponseScalars {
    pub a: i128,
    pub b: i128,
    pub c_or_d: i128,
    pub q: u128,
}

/// Stage-3 outputs: `c`, `d` modulo `2^r`, the commitment basis of order
/// `2^r` and its images on `E_chall`.
pub struct ResponseRecovery<L: FpBackend> {
    pub k: u128,
    pub c: u128,
    pub d: u128,
    pub w_com: Fp2<L>,
    pub r_com: JacPoint<L>,
    pub s_com: JacPoint<L>,
    pub phi_rsp_r_com: JacPoint<L>,
    pub phi_rsp_s_com: JacPoint<L>,
}

fn jac_signed_mul<L: FpBackend>(p: &JacPoint<L>, s: i128, curve: &EcCurve<L>) -> JacPoint<L> {
    let mag = s.unsigned_abs();
    let limbs = [mag as u64, (mag >> 64) as u64];
    let base = if s < 0 { p.neg() } else { p.clone() };
    jac_scalar_mul(&base, &limbs, curve)
}

/// Recover the response images from the challenge state and the signature.
pub fn recover_response<L: HdLevel>(
    chal: &ChallengeRecovery<L>,
    a_com: &Fp2<L>,
    hcp: u32,
    hcq: u32,
    s: ResponseScalars,
) -> Option<ResponseRecovery<L>> {
    let ResponseScalars { a, b, c_or_d, q } = s;
    let r = L::R;
    let rescale = (L::TWO_ADIC_EXPONENT - r) as usize;
    let (p_com, q_com) = hd_torsion_basis::<L>(a_com, hcp, hcq)?;
    let mut e_com = EcCurve::from_a(a_com)?;
    e_com.normalize_a24();
    let r_com = jac_dbl_iter(&p_com, rescale, &e_com);
    let s_com = jac_dbl_iter(&q_com, rescale, &e_com);

    let r_xz = r_com.to_xz();
    let s_xz = s_com.to_xz();
    let rms = jac_add(&r_com, &s_com.neg(), &e_com).to_xz();
    let w_com = weil(r, &r_xz, &s_xz, &rms, &mut e_com);

    // k with w_chall^k = w_com: the dlog of w_com to the base w_chall,
    // given 1 / w_chall
    let mut scal = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e(&mut scal[..L::ORDER_WORDS], &w_com, &chal.w_chal.inv(), r)?;
    let k = (scal[0] as u128) | ((scal[1] as u128) << 64);

    let (c, d) = recover_response_cd(a, b, c_or_d, q, k, r);

    let (c_img, d_img): (i128, i128) = if a & 1 != 0 {
        (c_or_d, d as i128)
    } else {
        (c as i128, c_or_d)
    };

    let e_chal = &chal.e_chal;
    let phi_rsp_r_com = jac_add(
        &jac_signed_mul(&chal.p_chal_resc, a, e_chal),
        &jac_signed_mul(&chal.q_chal_resc, b, e_chal),
        e_chal,
    );
    let phi_rsp_s_com = jac_add(
        &jac_signed_mul(&chal.p_chal_resc, c_img, e_chal),
        &jac_signed_mul(&chal.q_chal_resc, d_img, e_chal),
        e_chal,
    );

    Some(ResponseRecovery {
        k,
        c,
        d,
        w_com,
        r_com,
        s_com,
        phi_rsp_r_com,
        phi_rsp_s_com,
    })
}
