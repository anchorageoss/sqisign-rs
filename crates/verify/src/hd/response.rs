//! Phase 5b.3 - response image recovery (stage 3).
//!
//! Given the self-derived challenge data from 5b.2 ([`ChallengeRecovery`]) and
//! the commitment curve, recover the discrete log `k`, the response scalars
//! `(c, d)`, and the response isogeny images `φ_rsp(R_com)`, `φ_rsp(S_com)`.
//! Mirrors `Verify.py::image_response`:
//!
//! ```text
//! R_com = 2^(e-r)·P_com ;  S_com = 2^(e-r)·Q_com         # order 2^r on E_com
//! w_com = weil(R_com, S_com, 2^r)
//! k     = dlog(w_com, w_chal, 2^r)                       # w_chal^k = w_com
//! (a odd)  c = c_or_d, d = a⁻¹(k·q + b·c) mod 2^r
//! (a even) d = c_or_d, c = b⁻¹(a·d - k·q) mod 2^r        # a·d - b·c ≡ k·q
//! φ_rsp(R_com) = a·P_chal_resc + b·Q_chal_resc           # on E_chal
//! φ_rsp(S_com) = c·P_chal_resc + d·Q_chal_resc
//! ```
//!
//! # The `w_chal` convention cancels here (the key 5b.3 fact)
//!
//! 5b.2 found that the dim-2 biextension [`weil`] equals the **inverse** of the
//! oracle's PARI Weil pairing, so its `w_chal` is `w_chal_oracle⁻¹`. This stage
//! computes `w_com` with the **same** `weil`, so `w_com = w_com_oracle⁻¹` too.
//! The discrete log then satisfies
//! `w_com⁻¹·1 = (w_chal⁻¹)^k ⇒ w_com_oracle = w_chal_oracle^k`, i.e. `k` is
//! **identical** to the oracle's - the global convention cancels because it
//! divides out of both pairing arguments. This is verified against the
//! oracle's recorded `k` for all 5 vectors.
//!
//! # Reused vs new
//!
//! Reused: the dim-2 [`weil`] and [`fp2_dlog_2e_pub`] (the field-element
//! `2^e`-dlog, exactly `discrete_log_pari`), [`jac_add`], the 5b.1 basis, the
//! 5b.2 challenge data, and the Phase-5 [`recover_response_cd`] (the `(c,d)`
//! determinant solve mod `2^r`). New here: the connective tissue (rescaled
//! commitment basis, `w_com`, the dlog call, the signed full-point linear
//! combinations for the images).

use crate::ec::jacobian::jac_add;
use crate::ec::pairing::{fp2_dlog_2e_pub, weil};
use crate::ec::{EcCurve, JacPoint};
use crate::{Fp2, Level1};

use crate::hd::challenge::{jac_dbl_iter, jac_scalar_mul, ChallengeRecovery};
use crate::hd::hd_torsion_basis_l1;
use crate::hd::hd_verify::recover_response_cd;

/// Level-1 parameters: torsion exponent `e = 248`, response modulus `r = 70`,
/// so the commitment basis is rescaled by `2^(e-r) = 2^178` to land on the
/// `2^r`-torsion.
const E_L1: u32 = 248;
const R_L1: u32 = 70;
const RESCALE_RESP_BITS: usize = (E_L1 - R_L1) as usize;

/// The signed response scalars from the signature: `a`, `b`, the stored
/// `c_or_d`, and the response degree `q` (only `q mod 2^r` is used, so it may
/// be supplied reduced mod any multiple of `2^r`).
#[derive(Clone, Copy, Debug)]
pub struct ResponseScalars {
    pub a: i128,
    pub b: i128,
    pub c_or_d: i128,
    pub q: u128,
}

/// Stage-3 outputs.
///
/// `c`/`d` are reduced mod `2^r`. `r_com`/`s_com` are on `E_com`;
/// `phi_rsp_r_com`/`phi_rsp_s_com` are on `E_chal`. Use [`crate::hd::jac_to_affine`]
/// to compare against affine references.
pub struct ResponseRecovery {
    /// Discrete log `k` with `w_chal^k = w_com` (matches the oracle exactly).
    pub k: u128,
    /// Response scalar `c` mod `2^r`.
    pub c: u128,
    /// Response scalar `d` mod `2^r`.
    pub d: u128,
    /// Commitment-basis Weil pairing `e_{2^r}(R_com, S_com)` (native `weil`
    /// convention - the inverse of the oracle's PARI value).
    pub w_com: Fp2<Level1>,
    /// Rescaled commitment basis on `E_com`.
    pub r_com: JacPoint<Level1>,
    pub s_com: JacPoint<Level1>,
    /// Response isogeny images on `E_chal`.
    pub phi_rsp_r_com: JacPoint<Level1>,
    pub phi_rsp_s_com: JacPoint<Level1>,
}

/// Full-point scalar multiplication by a signed `i128`: `[s]·P` via
/// `[|s|]·(±P)`. The magnitude fits in two limbs at Level 1 (the response
/// scalars are < 2^71).
fn jac_signed_mul(p: &JacPoint<Level1>, s: i128, curve: &EcCurve<Level1>) -> JacPoint<Level1> {
    let mag = s.unsigned_abs();
    let limbs = [mag as u64, (mag >> 64) as u64];
    let base = if s < 0 { p.neg() } else { p.clone() };
    jac_scalar_mul(&base, &limbs, curve)
}

/// Recover the response images (`Verify.py::image_response`) for Level 1,
/// consuming the self-derived challenge data from 5b.2.
///
/// `chal` is the 5b.2 output; `a_com`/`(hcp, hcq)` recover the commitment
/// basis; `s` carries the signed signature scalars. Returns `None` if a
/// curve/dlog is degenerate.
pub fn recover_response_l1(
    chal: &ChallengeRecovery,
    a_com: &Fp2<Level1>,
    hcp: u32,
    hcq: u32,
    s: ResponseScalars,
) -> Option<ResponseRecovery> {
    let ResponseScalars { a, b, c_or_d, q } = s;
    // Commitment-curve 2^248-torsion basis (5b.1), rescaled to the 2^r torsion.
    let (p_com, q_com) = hd_torsion_basis_l1(a_com, hcp, hcq)?;
    let mut e_com = EcCurve::from_a(a_com)?;
    e_com.normalize_a24();
    let r_com = jac_dbl_iter(&p_com, RESCALE_RESP_BITS, &e_com);
    let s_com = jac_dbl_iter(&q_com, RESCALE_RESP_BITS, &e_com);

    // w_com = e_{2^r}(R_com, S_com) with the native biextension weil.
    let r_xz = r_com.to_xz();
    let s_xz = s_com.to_xz();
    let rms = jac_add(&r_com, &s_com.neg(), &e_com).to_xz();
    let w_com = weil(R_L1, &r_xz, &s_xz, &rms, &mut e_com);

    // k with w_chal^k = w_com (fp2_dlog_2e_pub: f = g^scal, g_inverse = g⁻¹).
    // Both pairings share the inverse convention, so it cancels and k matches
    // the oracle.
    let mut scal = [0u64; 4];
    fp2_dlog_2e_pub(&mut scal, &w_com, &chal.w_chal.inv(), R_L1)?;
    let k = (scal[0] as u128) | ((scal[1] as u128) << 64);

    // Response scalars (c, d) mod 2^r via the Phase-5 determinant solve.
    let (c, d) = recover_response_cd(a, b, c_or_d, q, k, R_L1);

    // Image scalars: sage uses the *raw signed* c_or_d in its slot (and the
    // reduced value in the other). Only the c_or_d·P_chal_resc term is
    // sensitive to this - P_chal_resc has order ≫ 2^r - so it must stay raw.
    let (c_img, d_img): (i128, i128) = if a & 1 != 0 {
        (c_or_d, d as i128)
    } else {
        (c as i128, c_or_d)
    };

    // φ_rsp images on E_chal (chal.e_chal is normalised: a = affine A, C = 1).
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
