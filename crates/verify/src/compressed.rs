//! The compressed signature: three of the four entries of the basis-change
//! matrix, the fourth recovered from two Weil pairings the verifier can
//! compute anyway.
//!
//! **This format is not in the SQIsign specification.** It is the P21
//! method of prism-rs carried to SQIsign round 3.
//!
//! # What the pairings give
//!
//! With `t = RESPONSE_BITS + 2`, the verifier computes the canonical bases
//! `(P, Q)` of `E_chall[2^t]` and `(P_aux, Q_aux)` of `E_aux[2^t]` and forms
//! `P' = [m00] P + [m10] Q`, `Q' = [m01] P + [m11] Q`. The
//! `(2^RESPONSE_BITS, 2^RESPONSE_BITS)`-isogeny has kernel `<[4](P', P_aux),
//! [4](Q', Q_aux)>`, maximal isotropic, so `e((P', P_aux), (Q', Q_aux))^4 =
//! 1` at `2^t`, that is `g^(4 det M) = h^-4` for `g = e(P, Q)`, `h = e(P_aux,
//! Q_aux)`: the two pairings and one discrete logarithm in the
//! `2^t`-subgroup of `F_p^2^*` give `det M` modulo `2^RESPONSE_BITS`, and
//! nothing above it. The isogeny depends on the entries modulo
//! `2^RESPONSE_BITS` only, but the chain that computes it does not: of the
//! sixteen values of the entries' bits `RESPONSE_BITS`, it accepts two
//! (the signer's, and the first column scaled by `1 + 2^RESPONSE_BITS`),
//! while bit `RESPONSE_BITS + 1` of every entry is free (measured on
//! generated signatures against this verifier and the reference; see
//! `COMPRESSION.md`). So the format carries the four bits `RESPONSE_BITS`
//! and drops the bits above.
//!
//! # The format
//!
//! ```text
//! A_aux | m00 | m01 | (m10 if m00 is odd, else m11) | challenge | bits | hint_aux | hint_chall
//! ```
//!
//! The three entries are the low `RESPONSE_BITS` bits, `ceil(RESPONSE_BITS
//! / 8)` bytes each; `bits` holds bit `RESPONSE_BITS` of `m00, m01, m10,
//! m11` in its bits 0 to 3. The verifier recovers the dropped entry modulo
//! `2^RESPONSE_BITS` from the determinant and the other three (the pivot is
//! odd, the determinant being odd), sets the four bits, and runs the
//! standard verification on the result. Sizes `2 fp + 3 ceil(rb / 8) + cb +
//! 3`: 176 / 269 / 353 bytes against 200 / 306 / 406. The recovery costs
//! two Weil pairings and one variable-base discrete logarithm at `2^t`
//! (`BENCH.md`).

use crate::ec::pairing::{fp2_dlog_2e, weil};
use crate::ec::{EcBasis, EcCurve, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;
use crate::precomp::PrimePrecomp;
use crate::sqisign::{
    challenge_state, check_canonical_basis_change_matrix, digits_from_bytes, digits_to_bytes,
    finish, level5, ChallengeState, PreparedPublicKey, PublicKey, Scalar, Signature, VerifyParams,
    HD_EXTRA_TORSION,
};

/// Bytes of one transmitted entry: `ceil(RESPONSE_BITS / 8)`.
pub const fn entry_bytes(params: &VerifyParams) -> usize {
    (params.response_bits as usize).div_ceil(8)
}

/// Encoded length of a compressed signature at these parameters.
pub const fn compressed_bytes(params: &VerifyParams) -> usize {
    2 * params.fp_encoded_bytes + 3 * entry_bytes(params) + params.challenge_bytes + 3
}

/// The largest compressed signature of the three levels (`p664_17`).
pub const MAX_COMPRESSED_BYTES: usize = compressed_bytes(&level5());

/// A compressed signature.
#[derive(Clone, Debug)]
pub struct CompressedSignature<L: Prime> {
    /// Montgomery coefficient of the auxiliary curve.
    pub e_aux_a: Fp2<L>,
    /// `M_00 mod 2^RESPONSE_BITS`.
    pub m00: Scalar,
    /// `M_01 mod 2^RESPONSE_BITS`.
    pub m01: Scalar,
    /// `M_10` when `M_00` is odd, `M_11` otherwise, modulo `2^RESPONSE_BITS`.
    pub m_var: Scalar,
    /// The challenge scalar.
    pub chall_coeff: Scalar,
    /// Bit `RESPONSE_BITS` of `M_00, M_01, M_10, M_11` in bits 0 to 3.
    pub bits: u8,
    /// Basis hint of the auxiliary curve.
    pub hint_aux: u8,
    /// Basis hint of the challenge curve.
    pub hint_chall: u8,
}

// ---- arithmetic modulo 2^t on scalars -------------------------------------

fn words(t: u32) -> usize {
    (t as usize).div_ceil(64)
}

/// `x <- x mod 2^t`.
pub(crate) fn mask_2exp(x: &mut Scalar, t: u32) {
    for (i, w) in x.iter_mut().enumerate() {
        let lo = 64 * i as u32;
        if lo >= t {
            *w = 0;
        } else if lo + 64 > t {
            *w &= (1u64 << (t - lo)) - 1;
        }
    }
}

fn reduced(x: &Scalar, t: u32) -> Scalar {
    let mut y = *x;
    mask_2exp(&mut y, t);
    y
}

fn mul_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let n = words(t);
    let mut out = [0u64; MAX_ORDER_WORDS];
    for i in 0..n {
        let mut carry = 0u128;
        for j in 0..n - i {
            let acc = out[i + j] as u128 + (a[i] as u128) * (b[j] as u128) + carry;
            out[i + j] = acc as u64;
            carry = acc >> 64;
        }
    }
    mask_2exp(&mut out, t);
    out
}

fn add_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let mut out = [0u64; MAX_ORDER_WORDS];
    let mut carry = 0u128;
    for i in 0..words(t) {
        let s = a[i] as u128 + b[i] as u128 + carry;
        out[i] = s as u64;
        carry = s >> 64;
    }
    mask_2exp(&mut out, t);
    out
}

fn sub_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let mut out = [0u64; MAX_ORDER_WORDS];
    let mut borrow = 0u64;
    for i in 0..words(t) {
        let (d, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = d.overflowing_sub(borrow);
        out[i] = d;
        borrow = (b1 | b2) as u64;
    }
    mask_2exp(&mut out, t);
    out
}

/// `a^-1 mod 2^t` for odd `a` (Newton from the three bits `a` gets right
/// itself); `None` for even `a`.
fn inv_mod_2exp(a: &Scalar, t: u32) -> Option<Scalar> {
    if a[0] & 1 == 0 {
        return None;
    }
    let mut two = [0u64; MAX_ORDER_WORDS];
    two[0] = 2;
    let mut x = reduced(a, t);
    let mut bits = 3u32;
    while bits < t {
        let ax = mul_mod_2exp(a, &x, t);
        let f = sub_mod_2exp(&two, &ax, t);
        x = mul_mod_2exp(&x, &f, t);
        bits *= 2;
    }
    Some(x)
}

fn bit_at(x: &Scalar, k: u32) -> u8 {
    ((x[(k / 64) as usize] >> (k % 64)) & 1) as u8
}

fn set_bit(x: &mut Scalar, k: u32, v: u8) {
    x[(k / 64) as usize] |= ((v & 1) as u64) << (k % 64);
}

// ---- compression -----------------------------------------------------------

/// Compress a signature: keep the low `RESPONSE_BITS` bits of three
/// entries and bit `RESPONSE_BITS` of all four. `None` if the matrix has an
/// even determinant (not a signature).
pub fn compress<L: FpBackend>(
    params: &VerifyParams,
    sig: &Signature<L>,
) -> Option<CompressedSignature<L>> {
    let r = params.response_bits;
    let m = &sig.mat_bchall_can_to_bchall;
    let det_odd = (m[0][0][0] & m[1][1][0] ^ m[0][1][0] & m[1][0][0]) & 1 == 1;
    if !det_odd {
        return None;
    }
    let m_var = if m[0][0][0] & 1 == 1 {
        &m[1][0]
    } else {
        &m[1][1]
    };
    let bits = bit_at(&m[0][0], r)
        | bit_at(&m[0][1], r) << 1
        | bit_at(&m[1][0], r) << 2
        | bit_at(&m[1][1], r) << 3;
    Some(CompressedSignature {
        e_aux_a: sig.e_aux_a.clone(),
        m00: reduced(&m[0][0], r),
        m01: reduced(&m[0][1], r),
        m_var: reduced(m_var, r),
        chall_coeff: sig.chall_coeff,
        bits,
        hint_aux: sig.hint_aux,
        hint_chall: sig.hint_chall,
    })
}

/// The matrix from the three entries, the determinant modulo
/// `2^RESPONSE_BITS` and the four bits. `None` if the pivot is even.
fn reconstruct(
    params: &VerifyParams,
    c: &CompressedSignature<impl Prime>,
    det: &Scalar,
) -> Option<[[Scalar; 2]; 2]> {
    let r = params.response_bits;
    let (m10, m11) = if c.m00[0] & 1 == 1 {
        // m11 = (det + m01 m10) / m00
        let num = add_mod_2exp(det, &mul_mod_2exp(&c.m01, &c.m_var, r), r);
        (c.m_var, mul_mod_2exp(&num, &inv_mod_2exp(&c.m00, r)?, r))
    } else {
        // m10 = (m00 m11 - det) / m01
        let num = sub_mod_2exp(&mul_mod_2exp(&c.m00, &c.m_var, r), det, r);
        (mul_mod_2exp(&num, &inv_mod_2exp(&c.m01, r)?, r), c.m_var)
    };
    let mut m = [[c.m00, c.m01], [m10, m11]];
    set_bit(&mut m[0][0], r, c.bits);
    set_bit(&mut m[0][1], r, c.bits >> 1);
    set_bit(&mut m[1][0], r, c.bits >> 2);
    set_bit(&mut m[1][1], r, c.bits >> 3);
    Some(m)
}

// ---- encoding ------------------------------------------------------------

/// Encode; returns the length written, `None` if `out` is too short.
pub fn compressed_to_bytes<L: FpBackend>(
    params: &VerifyParams,
    c: &CompressedSignature<L>,
    out: &mut [u8],
) -> Option<usize> {
    let n = compressed_bytes(params);
    if out.len() < n {
        return None;
    }
    let enc = c.e_aux_a.encode();
    let mut pos = enc.len();
    out[..pos].copy_from_slice(&enc[..]);
    let eb = entry_bytes(params);
    for x in [&c.m00, &c.m01, &c.m_var] {
        digits_to_bytes(x, &mut out[pos..pos + eb]);
        pos += eb;
    }
    digits_to_bytes(&c.chall_coeff, &mut out[pos..pos + params.challenge_bytes]);
    pos += params.challenge_bytes;
    out[pos] = c.bits & 0x0f;
    out[pos + 1] = c.hint_aux;
    out[pos + 2] = c.hint_chall;
    debug_assert_eq!(pos + 3, n);
    Some(n)
}

/// Decode: the length, a canonical `A_aux`, entries below
/// `2^RESPONSE_BITS`, the bits byte within its four bits, an odd pivot.
pub fn compressed_from_bytes<L: FpBackend>(
    params: &VerifyParams,
    bytes: &[u8],
) -> Option<CompressedSignature<L>> {
    if bytes.len() != compressed_bytes(params) {
        return None;
    }
    let r = params.response_bits;
    let mut pos = 2 * params.fp_encoded_bytes;
    let e_aux_a = Fp2::<L>::decode(&bytes[..pos])?;
    let eb = entry_bytes(params);
    let mut next = |n: usize| {
        let x = digits_from_bytes(&bytes[pos..pos + n]);
        pos += n;
        x
    };
    let m00 = next(eb);
    let m01 = next(eb);
    let m_var = next(eb);
    let chall_coeff = next(params.challenge_bytes);
    for x in [&m00, &m01, &m_var] {
        if reduced(x, r) != *x {
            return None;
        }
    }
    // the pivot is odd: m00, or m01 when m00 is even (else det is even)
    if m00[0] & 1 == 0 && m01[0] & 1 == 0 {
        return None;
    }
    let bits = bytes[pos];
    if bits & !0x0f != 0 {
        return None;
    }
    Some(CompressedSignature {
        e_aux_a,
        m00,
        m01,
        m_var,
        chall_coeff,
        bits,
        hint_aux: bytes[pos + 1],
        hint_chall: bytes[pos + 2],
    })
}

// ---- verification -----------------------------------------------------------

/// `log_g(h^-1)` in the `2^t`-subgroup for `g = e(P, Q)` on `E_chall` and
/// `h = e(P_aux, Q_aux)` on `E_aux`, `t = RESPONSE_BITS + 2`; `det M` is
/// its residue modulo `2^RESPONSE_BITS`. `None` if `h^-1` is not a power
/// of `g` (a basis of the wrong order).
pub fn pairing_dlog<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    e_chall: &EcCurve<L>,
    b_chall: &EcBasis<L>,
    e_aux: &EcCurve<L>,
    b_aux: &EcBasis<L>,
) -> Option<Scalar> {
    let t = params.response_bits + HD_EXTRA_TORSION;
    let mut e_chall = e_chall.clone();
    let mut e_aux = e_aux.clone();
    let g = weil(t, &b_chall.p, &b_chall.q, &b_chall.pmq, &mut e_chall);
    let h = weil(t, &b_aux.p, &b_aux.q, &b_aux.pmq, &mut e_aux);
    let mut y = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e(&mut y[..L::ORDER_WORDS], &h.inv(), &g.inv(), t)?;
    Some(y)
}

/// The recovery: the challenge state, the pairings, the matrix. `None` on
/// any failure.
fn recover<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    c: &CompressedSignature<L>,
) -> Option<(ChallengeState<L>, Signature<L>)> {
    let st = challenge_state(
        params,
        key,
        &c.e_aux_a,
        &c.chall_coeff,
        c.hint_chall,
        c.hint_aux,
    )?;
    let y = pairing_dlog(
        params,
        &st.e_chall,
        &st.b_chall_can,
        &st.e_aux,
        &st.b_aux_can,
    )?;
    let det = reduced(&y, params.response_bits);
    let mat = reconstruct(params, c, &det)?;
    let sig = Signature {
        e_aux_a: c.e_aux_a.clone(),
        mat_bchall_can_to_bchall: mat,
        chall_coeff: c.chall_coeff,
        hint_aux: c.hint_aux,
        hint_chall: c.hint_chall,
    };
    if !check_canonical_basis_change_matrix(params, &sig) {
        return None;
    }
    Some((st, sig))
}

/// `protocols_verify` on a compressed signature under a prepared key.
pub fn verify_compressed_prepared<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    key: &PreparedPublicKey<L>,
    c: &CompressedSignature<L>,
    msg: &[u8],
) -> bool {
    let Some((st, sig)) = recover(params, key, c) else {
        return false;
    };
    finish(
        params,
        key,
        st,
        &sig.mat_bchall_can_to_bchall,
        &sig.chall_coeff,
        msg,
    )
    .is_some()
}

/// `protocols_verify` on a compressed signature.
pub fn verify_compressed<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    c: &CompressedSignature<L>,
    msg: &[u8],
) -> bool {
    let Some(key) = PreparedPublicKey::new(params, pk) else {
        return false;
    };
    verify_compressed_prepared(params, &key, c, msg)
}

/// The standard signature a compressed one stands for, recovered under
/// `pk` (bit `RESPONSE_BITS + 1` of every entry zero, which the standard
/// verifier accepts). `None` if the recovery fails, as the verification
/// would.
pub fn decompress<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    c: &CompressedSignature<L>,
) -> Option<Signature<L>> {
    let key = PreparedPublicKey::new(params, pk)?;
    recover(params, &key, c).map(|(_, sig)| sig)
}
