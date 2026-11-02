//! Stage 3 of the compact verification: the response images (the library's
//! `image_response`). With `R_com = 2^(f − r) P_com`, `S_com = 2^(f − r)
//! Q_com` of order `2^r`, `w_com = e(R_com, S_com)`, `k = log_{w_chall}
//! w_com`, the dropped scalar from `a d − b c ≡ k q (mod 2^r)`, and the
//! images `σ(R_com) = a P_resc + b Q_resc`, `σ(S_com) = c P_resc + d Q_resc`.

use alloc::vec::Vec;

use crate::ec::jacobian::{jac_add, jac_dbl};
use crate::ec::pairing::{fp2_dlog_2e, weil};
use crate::ec::{EcCurve, JacPoint, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};

use super::basis::hd_torsion_basis;
use super::challenge::{jac_dbl_iter, ChallengeRecovery};
use super::hd_verify::recover_response_cd;
use super::params::HdLevel;
use super::profile::{NoProfiler, Profiler};

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
    /// `P_com`, the first point of the hinted basis of `E_com[2^f]` (the
    /// stage-6 input).
    pub p_com: JacPoint<L>,
    pub r_com: JacPoint<L>,
    pub s_com: JacPoint<L>,
    pub phi_rsp_r_com: JacPoint<L>,
    pub phi_rsp_s_com: JacPoint<L>,
}

/// `[a] P + [b] Q` for signed scalars by one joint double-and-add over the
/// bits of `|a|` and `|b|` with `P + Q` precomputed (Straus), variable time
/// (public data).
fn jac_double_mul<L: FpBackend>(
    p: &JacPoint<L>,
    a: i128,
    q: &JacPoint<L>,
    b: i128,
    curve: &EcCurve<L>,
) -> JacPoint<L> {
    let (ma, mb) = (a.unsigned_abs(), b.unsigned_abs());
    let pa = if a < 0 { p.neg() } else { p.clone() };
    let qb = if b < 0 { q.neg() } else { q.clone() };
    if ma | mb == 0 {
        return JacPoint::identity();
    }
    // the joint sparse form of (|a|, |b|): digit pairs in {-1, 0, 1}^2, at
    // most half of the positions nonzero, against three quarters for the
    // joint binary form; the table holds P, Q, P + Q and P - Q and their
    // negatives
    let digits = joint_sparse_form(ma, mb);
    let pq = jac_add(&pa, &qb, curve);
    let pmq = jac_add(&pa, &qb.neg(), curve);
    let table = [
        pa.clone(), // (1, 0)
        qb.clone(), // (0, 1)
        pq,         // (1, 1)
        pmq,        // (1, -1)
        pa.neg(),   // (-1, 0)
        qb.neg(),   // (0, -1)
    ];
    let pick = |u: i8, v: i8| -> Option<JacPoint<L>> {
        match (u, v) {
            (0, 0) => None,
            (1, 0) => Some(table[0].clone()),
            (0, 1) => Some(table[1].clone()),
            (1, 1) => Some(table[2].clone()),
            (1, -1) => Some(table[3].clone()),
            (-1, 0) => Some(table[4].clone()),
            (0, -1) => Some(table[5].clone()),
            (-1, -1) => Some(table[2].neg()),
            _ => Some(table[3].neg()), // (-1, 1)
        }
    };
    let top = digits.len() - 1;
    let mut acc = pick(digits[top].0, digits[top].1).expect("the top digit pair is nonzero");
    for &(u, v) in digits[..top].iter().rev() {
        acc = jac_dbl(&acc, curve);
        if let Some(t) = pick(u, v) {
            acc = jac_add(&acc, &t, curve);
        }
    }
    acc
}

/// The joint sparse form of `(a, b)` (Solinas), least significant pair
/// first, the last pair nonzero.
fn joint_sparse_form(mut a: u128, mut b: u128) -> Vec<(i8, i8)> {
    let mut out = Vec::with_capacity(130);
    let (mut da, mut db) = (0u128, 0u128);
    while a + da > 0 || b + db > 0 {
        let la = (a & 7) + da;
        let lb = (b & 7) + db;
        let digit = |l: u128, other: u128| -> i8 {
            if l & 1 == 0 {
                return 0;
            }
            let mut u: i8 = if l & 3 == 1 { 1 } else { -1 };
            if (l & 7 == 3 || l & 7 == 5) && other & 3 == 2 {
                u = -u;
            }
            u
        };
        let u = digit(la, lb);
        let v = digit(lb, la);
        if 1 + u as i32 == 2 * da as i32 {
            da = 1 - da;
        }
        if 1 + v as i32 == 2 * db as i32 {
            db = 1 - db;
        }
        a >>= 1;
        b >>= 1;
        out.push((u, v));
    }
    out
}

/// Recover the response images from the challenge state and the signature.
pub fn recover_response<L: HdLevel>(
    chal: &ChallengeRecovery<L>,
    a_com: &Fp2<L>,
    hcp: u32,
    hcq: u32,
    s: ResponseScalars,
) -> Option<ResponseRecovery<L>> {
    recover_response_prof::<L, NoProfiler>(chal, a_com, hcp, hcq, s, &mut NoProfiler)
}

/// [`recover_response`] with a [`Profiler`]: sub-slots 6 to 9.
pub fn recover_response_prof<L: HdLevel, P: Profiler>(
    chal: &ChallengeRecovery<L>,
    a_com: &Fp2<L>,
    hcp: u32,
    hcq: u32,
    s: ResponseScalars,
    prof: &mut P,
) -> Option<ResponseRecovery<L>> {
    let ResponseScalars { a, b, c_or_d, q } = s;
    // The rule for the dropped scalar inverts a when a is odd and b
    // otherwise; with a and b both even it is not defined, and the
    // signature is rejected here, before any inversion.
    if a & 1 == 0 && b & 1 == 0 {
        return None;
    }
    let r = L::R;
    let rescale = (L::TWO_ADIC_EXPONENT - r) as usize;
    let t = prof.now();
    let (p_com, q_com) = hd_torsion_basis::<L>(a_com, hcp, hcq)?;
    prof.sub_add(6, prof.now() - t);
    let t = prof.now();
    let mut e_com = EcCurve::from_a(a_com)?;
    e_com.normalize_a24();
    let r_com = jac_dbl_iter(&p_com, rescale, &e_com);
    let s_com = jac_dbl_iter(&q_com, rescale, &e_com);
    prof.sub_add(7, prof.now() - t);
    let t = prof.now();

    let r_xz = r_com.to_xz();
    let s_xz = s_com.to_xz();
    let rms = jac_add(&r_com, &s_com.neg(), &e_com).to_xz();
    let w_com = weil(r, &r_xz, &s_xz, &rms, &mut e_com);

    // k with w_chall^k = w_com: the dlog of w_com to the base w_chall,
    // given 1 / w_chall
    let mut scal = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e(&mut scal[..L::ORDER_WORDS], &w_com, &chal.w_chal.inv(), r)?;
    let k = (scal[0] as u128) | ((scal[1] as u128) << 64);
    prof.sub_add(8, prof.now() - t);
    let t = prof.now();

    let (c, d) = recover_response_cd(a, b, c_or_d, q, k, r);

    let (c_img, d_img): (i128, i128) = if a & 1 != 0 {
        (c_or_d, d as i128)
    } else {
        (c as i128, c_or_d)
    };

    let e_chal = &chal.e_chal;
    let phi_rsp_r_com = jac_double_mul(&chal.p_chal_resc, a, &chal.q_chal_resc, b, e_chal);
    let phi_rsp_s_com = jac_double_mul(&chal.p_chal_resc, c_img, &chal.q_chal_resc, d_img, e_chal);

    prof.sub_add(9, prof.now() - t);
    Some(ResponseRecovery {
        k,
        c,
        d,
        w_com,
        p_com,
        r_com,
        s_com,
        phi_rsp_r_com,
        phi_rsp_s_com,
    })
}
