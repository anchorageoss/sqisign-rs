//! The compressed signature: three of the four entries of the basis-change
//! matrix, the fourth recovered from two Weil pairings the verifier can
//! compute anyway, every field packed at its bit width. Compressing then
//! decompressing gives back the signer's signature exactly.
//!
//! **This format is not in the SQIsign specification.** It is the P21
//! method of eclipse-rs carried to SQIsign round 3.
//!
//! # What the pairings give
//!
//! With `t = RESPONSE_BITS + 2`, the verifier computes the canonical bases
//! `(P, Q)` of `E_chall[2^t]` and `(P_aux, Q_aux)` of `E_aux[2^t]` and forms
//! `P' = [m00] P + [m10] Q`, `Q' = [m01] P + [m11] Q`. The
//! `(2^RESPONSE_BITS, 2^RESPONSE_BITS)`-isogeny has kernel `<[4](P', P_aux),
//! [4](Q', Q_aux)>`. For a signature the signer produced, the Kani square
//! behind that kernel gives `g^(det M) h = zeta` with `g = e(P, Q)`,
//! `h = e(P_aux, Q_aux)` and `zeta` a primitive fourth root of unity: with
//! `y = log_g(h^-1)` taken in the `2^t`-subgroup of `F_p^2^*`,
//! `det M = y + 2^RESPONSE_BITS u (mod 2^t)` with `u` in `{1, 3}`. So two
//! pairings and one discrete logarithm give `det M` modulo
//! `2^(RESPONSE_BITS + 1)` outright, as `y + 2^RESPONSE_BITS`, and only bit
//! `RESPONSE_BITS + 1` of it is free. The dropped entry is therefore
//! recovered modulo `2^(RESPONSE_BITS + 1)` from the other three, and its
//! top bit, bit `RESPONSE_BITS + 1`, is transmitted.
//!
//! The signer normalises the matrix (the verifier's
//! `check_canonical_basis_change_matrix`): the first entry in wire order
//! that is neither `0` nor `2^(RESPONSE_BITS + 1)`, that is, whose low
//! `RESPONSE_BITS + 1` bits are not all zero, has its top bit zero. That
//! entry is `m00` or `m01` (the pivot is odd), never the dropped one, and
//! the decoder identifies it from the low bits as the encoder does. So the
//! format transmits the low `RESPONSE_BITS + 1` bits of the three entries
//! and three top bits: that of the dropped entry and those of the two
//! transmitted entries other than the normalising one. `compress` returns
//! `None` for a matrix the normal-form check rejects, and on every other
//! signature `decompress` returns the signer's signature.
//!
//! # The format
//!
//! One bit string, bit `i` of the string being bit `i mod 8` of byte
//! `i / 8`, each field the low bits of its little-endian integer:
//!
//! ```text
//! field       bits               content
//! A_aux       2 fp_bits          real then imaginary part, each in [0, p)
//! m00, m01    RESPONSE_BITS + 1  the low bits of the first row
//! m_var       RESPONSE_BITS + 1  the low bits of m10 if m00 is odd, else of m11
//! top         3                  bit RESPONSE_BITS + 1 of the dropped entry, then of
//!                                the two transmitted entries other than the
//!                                normalising one, in wire order
//! chl         CHALLENGE_BITS     the challenge scalar
//! hint_aux    8                  as in the standard format
//! hint_chall  8                  as in the standard format
//! padding     0 to 7             zero, to a byte boundary
//! ```
//!
//! with `fp_bits = ceil(log2 p)`: `2 fp_bits + 3 (RESPONSE_BITS + 1) + 3 +
//! CHALLENGE_BITS + 16` bits, 174 / 267 / 352 bytes against 200 / 306 /
//! 406. The decoder rejects a wrong length, a coordinate at or above `p`,
//! a non-zero padding bit and two even pivots. The verifier recovers the
//! dropped entry, sets its top bit from the string, and runs the standard
//! verification on the result. The recovery costs two Weil pairings and
//! one variable-base discrete logarithm at `2^t`.

use crate::ec::pairing::{fp2_dlog_2e, weil};
use crate::ec::{EcBasis, EcCurve, MAX_ORDER_WORDS};
use crate::fp::{Fp2, FpBackend};
use crate::params::Prime;
use crate::precomp::PrimePrecomp;
use crate::sqisign::{
    challenge_state, check_canonical_basis_change_matrix, digits_from_bytes, digits_to_bytes,
    finish, level5, ChallengeState, PreparedPublicKey, PublicKey, Scalar, Signature, VerifyParams,
    HD_EXTRA_TORSION, MAX_PUBLICKEY_BYTES,
};

/// Bits of one transmitted entry: `RESPONSE_BITS + 1`.
pub const fn entry_bits(params: &VerifyParams) -> usize {
    params.response_bits as usize + 1
}

/// Bits of a compressed signature before the padding.
pub const fn compressed_bits(params: &VerifyParams) -> usize {
    2 * params.fp_bits as usize + 3 * entry_bits(params) + 3 + params.challenge_bits as usize + 16
}

/// Encoded length of a compressed signature at these parameters.
pub const fn compressed_bytes(params: &VerifyParams) -> usize {
    compressed_bits(params).div_ceil(8)
}

/// The largest compressed signature of the three levels (`p664_17`).
pub const MAX_COMPRESSED_BYTES: usize = compressed_bytes(&level5());

/// A compressed signature.
#[derive(Clone, Debug)]
pub struct CompressedSignature<L: Prime> {
    /// Montgomery coefficient of the auxiliary curve.
    pub e_aux_a: Fp2<L>,
    /// `M_00`, below `2^(RESPONSE_BITS + 2)`.
    pub m00: Scalar,
    /// `M_01`, below `2^(RESPONSE_BITS + 2)`.
    pub m01: Scalar,
    /// `M_10` when `M_00` is odd, `M_11` otherwise, below
    /// `2^(RESPONSE_BITS + 2)`.
    pub m_var: Scalar,
    /// Bit `RESPONSE_BITS + 1` of the dropped entry.
    pub top: u8,
    /// The challenge scalar.
    pub chall_coeff: Scalar,
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

fn set_bit(x: &mut Scalar, k: u32, v: u8) {
    x[(k / 64) as usize] |= ((v & 1) as u64) << (k % 64);
}

fn bit_at(x: &Scalar, k: u32) -> u8 {
    ((x[(k / 64) as usize] >> (k % 64)) & 1) as u8
}

/// Whether `m00`, the first entry in wire order, is the normalising one
/// (its low `t = RESPONSE_BITS + 1` bits are not all zero); otherwise it
/// is `m01`.
fn m00_normalises(m00: &Scalar, t: u32) -> bool {
    reduced(m00, t).iter().any(|w| *w != 0)
}

// ---- bit strings -------------------------------------------------------------

#[inline]
fn bit(bytes: &[u8], i: usize) -> u8 {
    (bytes[i / 8] >> (i % 8)) & 1
}

/// Write the low `width` bits of the little-endian `src` at bit `pos` of
/// `out`, whose bits there are zero.
fn put_bits(out: &mut [u8], pos: usize, src: &[u8], width: usize) {
    for i in 0..width {
        out[(pos + i) / 8] |= bit(src, i) << ((pos + i) % 8);
    }
}

/// The `width` bits at `pos` of `inp` as a little-endian integer in `dst`
/// (cleared first; at least `ceil(width / 8)` bytes).
fn get_bits(inp: &[u8], pos: usize, width: usize, dst: &mut [u8]) {
    dst.iter_mut().for_each(|b| *b = 0);
    for i in 0..width {
        dst[i / 8] |= bit(inp, pos + i) << (i % 8);
    }
}

/// A cursor over the bit string of a compressed signature.
struct Bits<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Bits<'_> {
    fn scalar(&mut self, width: usize) -> Scalar {
        let mut buf = [0u8; 8 * MAX_ORDER_WORDS];
        let k = width.div_ceil(8);
        get_bits(self.bytes, self.pos, width, &mut buf[..k]);
        self.pos += width;
        digits_from_bytes(&buf[..k])
    }

    fn byte(&mut self) -> u8 {
        let mut buf = [0u8; 1];
        get_bits(self.bytes, self.pos, 8, &mut buf);
        self.pos += 8;
        buf[0]
    }

    /// Whether every bit from the cursor to the end is zero.
    fn rest_is_zero(&self) -> bool {
        (self.pos..8 * self.bytes.len()).all(|i| bit(self.bytes, i) == 0)
    }
}

// ---- compression -----------------------------------------------------------

/// Compress a signature: three entries and the top bit of the fourth.
/// `None` if the matrix has an even determinant or fails the verifier's
/// normal-form check (not a signature).
pub fn compress<L: FpBackend>(
    params: &VerifyParams,
    sig: &Signature<L>,
) -> Option<CompressedSignature<L>> {
    let r = params.response_bits;
    let m = &sig.mat_bchall_can_to_bchall;
    let det_odd = (m[0][0][0] & m[1][1][0] ^ m[0][1][0] & m[1][0][0]) & 1 == 1;
    if !det_odd || !check_canonical_basis_change_matrix(params, sig) {
        return None;
    }
    let (m_var, dropped) = if m[0][0][0] & 1 == 1 {
        (&m[1][0], &m[1][1])
    } else {
        (&m[1][1], &m[1][0])
    };
    Some(CompressedSignature {
        e_aux_a: sig.e_aux_a.clone(),
        m00: m[0][0],
        m01: m[0][1],
        m_var: *m_var,
        top: bit_at(dropped, r + 1),
        chall_coeff: sig.chall_coeff,
        hint_aux: sig.hint_aux,
        hint_chall: sig.hint_chall,
    })
}

/// The matrix from the three transmitted entries, the determinant modulo
/// `2^(RESPONSE_BITS + 1)` and the dropped entry's top bit. `None` if the
/// pivot is even.
fn reconstruct(
    params: &VerifyParams,
    c: &CompressedSignature<impl Prime>,
    det: &Scalar,
) -> Option<[[Scalar; 2]; 2]> {
    let r = params.response_bits;
    let t = r + 1;
    let (m00, m01, m_var) = (c.m00, c.m01, c.m_var);
    let (l00, l01, l_var) = (reduced(&m00, t), reduced(&m01, t), reduced(&m_var, t));
    let mut dropped = if l00[0] & 1 == 1 {
        // m11 = (det + m01 m10) / m00
        let num = add_mod_2exp(det, &mul_mod_2exp(&l01, &l_var, t), t);
        mul_mod_2exp(&num, &inv_mod_2exp(&l00, t)?, t)
    } else {
        // m10 = (m00 m11 - det) / m01
        let num = sub_mod_2exp(&mul_mod_2exp(&l00, &l_var, t), det, t);
        mul_mod_2exp(&num, &inv_mod_2exp(&l01, t)?, t)
    };
    set_bit(&mut dropped, r + 1, c.top);
    Some(if l00[0] & 1 == 1 {
        [[m00, m01], [m_var, dropped]]
    } else {
        [[m00, m01], [dropped, m_var]]
    })
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
    let out = &mut out[..n];
    out.iter_mut().for_each(|b| *b = 0);
    let fp = params.fp_encoded_bytes;
    let fb = params.fp_bits as usize;
    let enc = c.e_aux_a.encode();
    let mut pos = 0;
    put_bits(out, pos, &enc[..fp], fb);
    pos += fb;
    put_bits(out, pos, &enc[fp..2 * fp], fb);
    pos += fb;
    let eb = entry_bits(params);
    let t = params.response_bits + 1;
    let mut buf = [0u8; 8 * MAX_ORDER_WORDS];
    for x in [&c.m00, &c.m01, &c.m_var] {
        digits_to_bytes(&reduced(x, t), &mut buf[..eb.div_ceil(8)]);
        put_bits(out, pos, &buf, eb);
        pos += eb;
    }
    // the top bits: the dropped entry's, then those of the two transmitted
    // entries other than the normalising one, in wire order
    let (o1, o2) = if m00_normalises(&c.m00, t) {
        (&c.m01, &c.m_var)
    } else {
        (&c.m00, &c.m_var)
    };
    let top = (c.top & 1) | bit_at(o1, t) << 1 | bit_at(o2, t) << 2;
    put_bits(out, pos, &[top], 3);
    pos += 3;
    let cb = params.challenge_bits as usize;
    digits_to_bytes(&c.chall_coeff, &mut buf[..params.challenge_bytes]);
    put_bits(out, pos, &buf, cb);
    pos += cb;
    put_bits(out, pos, &[c.hint_aux], 8);
    pos += 8;
    put_bits(out, pos, &[c.hint_chall], 8);
    pos += 8;
    debug_assert_eq!(pos, compressed_bits(params));
    Some(n)
}

/// Decode: the length, a canonical `A_aux`, zero padding, an odd pivot;
/// the three top bits go to the dropped entry and to the two transmitted
/// entries other than the normalising one.
pub fn compressed_from_bytes<L: FpBackend>(
    params: &VerifyParams,
    bytes: &[u8],
) -> Option<CompressedSignature<L>> {
    if bytes.len() != compressed_bytes(params) {
        return None;
    }
    let fp = params.fp_encoded_bytes;
    let fb = params.fp_bits as usize;
    let mut r = Bits { bytes, pos: 0 };
    // the two coordinates, each fp_bits wide, into the canonical encoding
    let mut enc = [0u8; MAX_PUBLICKEY_BYTES];
    get_bits(bytes, r.pos, fb, &mut enc[..fp]);
    r.pos += fb;
    get_bits(bytes, r.pos, fb, &mut enc[fp..2 * fp]);
    r.pos += fb;
    let e_aux_a = Fp2::<L>::decode(&enc[..2 * fp])?;
    let eb = entry_bits(params);
    let mut m00 = r.scalar(eb);
    let mut m01 = r.scalar(eb);
    let mut m_var = r.scalar(eb);
    let top = r.scalar(3)[0] as u8;
    let chall_coeff = r.scalar(params.challenge_bits as usize);
    let hint_aux = r.byte();
    let hint_chall = r.byte();
    debug_assert_eq!(r.pos, compressed_bits(params));
    if !r.rest_is_zero() {
        return None;
    }
    // the pivot is odd: m00, or m01 when m00 is even (else det is even)
    if m00[0] & 1 == 0 && m01[0] & 1 == 0 {
        return None;
    }
    let t = params.response_bits + 1;
    let (o1, o2) = if m00_normalises(&m00, t) {
        (&mut m01, &mut m_var)
    } else {
        (&mut m00, &mut m_var)
    };
    set_bit(o1, t, top >> 1);
    set_bit(o2, t, top >> 2);
    Some(CompressedSignature {
        e_aux_a,
        m00,
        m01,
        m_var,
        top: top & 1,
        chall_coeff,
        hint_aux,
        hint_chall,
    })
}

// ---- verification -----------------------------------------------------------

/// `y = log_g(h^-1)` in the `2^t`-subgroup for `g = e(P, Q)` on `E_chall`
/// and `h = e(P_aux, Q_aux)` on `E_aux`, `t = RESPONSE_BITS + 2`; for a
/// signature the signer produced, `det M = y + 2^RESPONSE_BITS u (mod 2^t)`
/// with `u` odd. `None` if `h^-1` is not a power of `g` (a basis of the
/// wrong order).
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
    // det M = y + 2^RESPONSE_BITS u (mod 2^(RESPONSE_BITS + 2)) with u odd,
    // so det M = y + 2^RESPONSE_BITS (mod 2^(RESPONSE_BITS + 1)): the
    // dropped entry modulo 2^(RESPONSE_BITS + 1), its top bit transmitted.
    let r = params.response_bits;
    let mut two_r = [0u64; MAX_ORDER_WORDS];
    set_bit(&mut two_r, r, 1);
    let det = add_mod_2exp(&y, &two_r, r + 1);
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
/// `pk`: the signer's signature. `None` if the recovery fails, as the
/// verification would.
pub fn decompress<L: FpBackend + PrimePrecomp>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    c: &CompressedSignature<L>,
) -> Option<Signature<L>> {
    let key = PreparedPublicKey::new(params, pk)?;
    recover(params, &key, c).map(|(_, sig)| sig)
}
