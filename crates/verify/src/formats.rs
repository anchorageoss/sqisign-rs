//!
//! Three wire formats with different size/speed tradeoffs:
//!
//! | Format | Level 1 | Description |
//! |--------|---------|-------------|
//! | Standard | 148 B | Default: 2×2 matrix + hints |
//! | Expanded | 212 B | Pre-evaluated kernel points, faster verify |
//! | Compressed | 129 B | 3 of 4 matrix entries, 4th via determinant |
//!
//! The standard format is the NIST v2.0 wire format and is the only one
//! validated against KAT vectors. The other two are defined by the
//! SQIsign specification as optional compression levels.

use crate::ec::basis::{
    difference_point, difference_point_with_hint, ec_curve_to_basis_2f_to_hint,
};
use crate::ec::pairing::{fp2_dlog_2e_pub, weil};
use crate::ec::point::{ec_dbl_iter_basis, xadd};
use crate::ec::{EcBasis, EcPoint};
use crate::fp::{Fp2, FpBackend};
use crate::params::{Level1, SecurityLevel};
use crate::precomp::LevelPrecomp;
use crate::theta::HD_EXTRA_TORSION;
use hybrid_array::typenum::Unsigned;

use crate::hash::hash_to_challenge;
use crate::types::{decode_digits, encode_digits, PublicKey, Scalar, Signature};
use crate::verify::{
    basis_from_hint, check_canonical_basis_change_matrix, compute_challenge_curve,
    compute_commitment_curve_verify, matrix_scalar_application_even_basis, mp_compare, mp_is_even,
    protocols_verify, two_response_isogeny_verify_inner,
};
use crate::Error;

/// Identifies which wire format a signature uses.
///
/// Format detection is purely length-based: each format has a unique
/// wire size at every security level, so no prefix byte is needed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SignatureFormat {
    Expanded,
    Standard,
    Compressed,
}

/// Expanded signature with pre-evaluated kernel points.
///
/// Instead of storing the 2×2 basis change matrix, stores the resulting
/// x-coordinates of P_chl and Q_chl directly. The difference point
/// P_chl - Q_chl is recomputed during verification via `difference_point`.
/// This may yield faster verification at the cost of a larger wire format;
/// the actual speedup varies depending on the security level and platform.
#[derive(Clone, Debug)]
pub struct ExpandedSignature<L: SecurityLevel = Level1> {
    /// Montgomery A-coefficient of the auxiliary curve E_aux.
    pub(crate) e_aux_a: Fp2<L>,
    /// Number of backtracking steps in the response isogeny.
    /// Wire encoding packs flags into the high bits of this byte:
    /// bit 7 = kernel_is_q, bit 6 = pmq_sign_hint, bits 0-5 = backtracking.
    pub(crate) backtracking: u8,
    /// Length of the initial 2-isogeny chain in the response.
    pub(crate) two_resp_length: u8,
    /// Challenge coefficient (LAMBDA bits).
    pub(crate) chall_coeff: Scalar<L>,
    /// x-coordinate of the first basis image under the matrix action.
    pub(crate) p_chl_x: Fp2<L>,
    /// x-coordinate of the second basis image under the matrix action.
    pub(crate) q_chl_x: Fp2<L>,
    /// If true, the kernel generator for the small 2-isogeny chain is
    /// `Q_chl`; otherwise it is `P_chl`.
    pub(crate) kernel_is_q: bool,
    /// Sign hint for reconstructing P_chl - Q_chl via `difference_point`.
    /// When true, the verifier negates the discriminant to get the correct root.
    pub(crate) pmq_sign_hint: bool,
    /// Torsion basis hint for the auxiliary curve.
    pub(crate) hint_aux: u8,
    /// Torsion basis hint for the challenge curve.
    pub(crate) hint_chall: u8,
}

impl<L: FpBackend> ExpandedSignature<L> {
    /// Wire size in bytes.
    ///
    /// Layout: `Fp2 (e_aux) | backtracking_and_flags | two_resp_length |
    ///          LAMBDA/8 (challenge) | Fp2 (P_chl_x) | Fp2 (Q_chl_x) |
    ///          hint_aux | hint_chall`
    ///
    /// The `kernel_is_q` flag is packed into bit 7 of the backtracking byte.
    pub const WIRE_BYTES: usize = <L as SecurityLevel>::Fp2EncodedBytes::USIZE * 3
        + <L as SecurityLevel>::LAMBDA as usize / 8
        + 4;

    fn chall_coeff_bytes() -> usize {
        L::LAMBDA as usize / 8
    }

    /// Decode from bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != Self::WIRE_BYTES {
            return Err(Error::InvalidLength);
        }

        let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
        let mut pos = 0;

        let e_aux_a = Fp2::<L>::decode(&bytes[pos..pos + fp2_len]).ok_or(Error::MalformedInput)?;
        pos += fp2_len;

        let bt_byte = bytes[pos];
        pos += 1;
        let kernel_is_q = (bt_byte & 0x80) != 0;
        let pmq_sign_hint = (bt_byte & 0x40) != 0;
        let backtracking = bt_byte & 0x3F;

        let two_resp_length = bytes[pos];
        pos += 1;

        let chall_bytes = Self::chall_coeff_bytes();
        let mut chall_coeff = Scalar::<L>::default();
        decode_digits(
            chall_coeff.digits.as_mut_slice(),
            &bytes[pos..],
            chall_bytes,
        );
        pos += chall_bytes;

        let p_chl_x = Fp2::<L>::decode(&bytes[pos..pos + fp2_len]).ok_or(Error::MalformedInput)?;
        pos += fp2_len;

        let q_chl_x = Fp2::<L>::decode(&bytes[pos..pos + fp2_len]).ok_or(Error::MalformedInput)?;
        pos += fp2_len;

        let hint_aux = bytes[pos];
        pos += 1;
        let hint_chall = bytes[pos];
        pos += 1;
        debug_assert_eq!(pos, Self::WIRE_BYTES);

        Ok(Self {
            e_aux_a,
            backtracking,
            two_resp_length,
            chall_coeff,
            p_chl_x,
            q_chl_x,
            kernel_is_q,
            pmq_sign_hint,
            hint_aux,
            hint_chall,
        })
    }

    /// Encode to bytes.
    ///
    /// Returns a fixed-size buffer; the first `WIRE_BYTES` bytes are
    /// the meaningful payload.
    pub fn to_bytes(&self) -> [u8; 420] {
        let mut buf = [0u8; 420];
        let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
        let mut pos = 0;

        let enc = self.e_aux_a.encode();
        buf[pos..pos + fp2_len].copy_from_slice(&enc);
        pos += fp2_len;

        buf[pos] = self.backtracking
            | if self.kernel_is_q { 0x80 } else { 0 }
            | if self.pmq_sign_hint { 0x40 } else { 0 };
        pos += 1;
        buf[pos] = self.two_resp_length;
        pos += 1;

        let chall_bytes = Self::chall_coeff_bytes();
        encode_digits(
            &mut buf[pos..],
            self.chall_coeff.digits.as_slice(),
            chall_bytes,
        );
        pos += chall_bytes;

        let enc = self.p_chl_x.encode();
        buf[pos..pos + fp2_len].copy_from_slice(&enc);
        pos += fp2_len;

        let enc = self.q_chl_x.encode();
        buf[pos..pos + fp2_len].copy_from_slice(&enc);
        pos += fp2_len;

        buf[pos] = self.hint_aux;
        pos += 1;
        buf[pos] = self.hint_chall;
        pos += 1;
        debug_assert_eq!(pos, Self::WIRE_BYTES);

        buf
    }
}

impl<L: FpBackend + LevelPrecomp> ExpandedSignature<L> {
    /// Verify this expanded signature against a public key and message.
    pub fn verify(&self, pk: &PublicKey<L>, msg: &[u8]) -> Result<(), Error> {
        verify_expanded(pk, msg, self)
    }
}

/// Compute a canonical `2^f`-torsion basis and its hint byte from a curve.
fn basis_to_hint<L: FpBackend + LevelPrecomp>(
    curve: &mut crate::ec::EcCurve<L>,
    f: u32,
) -> Result<(EcBasis<L>, u8), Error> {
    ec_curve_to_basis_2f_to_hint(
        curve,
        f,
        L::basis_e0_px_bytes(),
        L::basis_e0_qx_bytes(),
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
        L::torsion_even_power(),
    )
    .map_err(|()| Error::MalformedInput)
}

/// Multiply `a * b`, truncated to `out.len()` words, then mask to `mod_bits`.
fn mp_mul_mod(out: &mut [u64], a: &[u64], b: &[u64], mod_bits: usize) {
    let n = out.len();
    out.fill(0);
    for i in 0..a.len().min(n) {
        let mut carry: u64 = 0;
        for j in 0..b.len().min(n - i) {
            let prod = (a[i] as u128) * (b[j] as u128) + (out[i + j] as u128) + (carry as u128);
            out[i + j] = prod as u64;
            carry = (prod >> 64) as u64;
        }
    }
    mp_mod_2exp(out, mod_bits);
}

/// Right-shift a digit array by `shift` bits.
fn mp_shiftr(a: &mut [u64], shift: usize) {
    let n = a.len();
    let word_shift = shift / 64;
    let bit_shift = shift % 64;

    if word_shift >= n {
        a.fill(0);
        return;
    }

    if bit_shift == 0 {
        for i in 0..n - word_shift {
            a[i] = a[i + word_shift];
        }
    } else {
        for i in 0..n - word_shift - 1 {
            a[i] = (a[i + word_shift] >> bit_shift) | (a[i + word_shift + 1] << (64 - bit_shift));
        }
        a[n - word_shift - 1] = a[n - 1] >> bit_shift;
    }
    a[n - word_shift..n].fill(0);
}

/// Mask a digit array to `e` bits.
fn mp_mod_2exp(a: &mut [u64], e: usize) {
    let q = e / 64;
    let r = e % 64;
    if q < a.len() {
        if r != 0 {
            a[q] &= (1u64 << r) - 1;
        } else {
            a[q] = 0;
        }
        a[q + 1..].fill(0);
    }
}

/// Compute `a^{-1} mod 2^e` via Newton/Hensel lifting.
/// Requires `a` to be odd.
fn hensel_inv_mod_2e(out: &mut [u64], a: &[u64], e: usize) {
    let n = out.len();
    out.fill(0);
    out[0] = 1;

    let mut ax = [0u64; 8];
    let mut factor = [0u64; 8];
    let mut x_copy = [0u64; 8];

    let mut prec = 1usize;
    while prec < e {
        let next = (prec * 2).min(e);

        ax[..n].fill(0);
        for i in 0..n {
            let mut carry: u64 = 0;
            for j in 0..n.min(n - i) {
                if i + j >= n {
                    break;
                }
                let prod =
                    (a[i] as u128) * (out[j] as u128) + (ax[i + j] as u128) + (carry as u128);
                ax[i + j] = prod as u64;
                carry = (prod >> 64) as u64;
            }
        }
        mp_mod_2exp(&mut ax[..n], next);

        // factor = 2 - ax = !ax + 3 (two's complement: -ax + 2 = !ax + 1 + 2)
        let mut carry = 3u128;
        for i in 0..n {
            let val = (!ax[i]) as u128 + carry;
            factor[i] = val as u64;
            carry = val >> 64;
        }
        mp_mod_2exp(&mut factor[..n], next);

        x_copy[..n].copy_from_slice(&out[..n]);
        out.fill(0);
        for i in 0..n {
            let mut carry_mul: u64 = 0;
            for j in 0..n.min(n - i) {
                if i + j >= n {
                    break;
                }
                let prod = (x_copy[i] as u128) * (factor[j] as u128)
                    + (out[i + j] as u128)
                    + (carry_mul as u128);
                out[i + j] = prod as u64;
                carry_mul = (prod >> 64) as u64;
            }
        }
        mp_mod_2exp(out, next);

        prec = next;
    }
}

/// Set bits `[start_bit..start_bit+n_bits)` from a hint byte.
fn set_high_bits(digits: &mut [u64], hint: u8, start_bit: usize, n_bits: usize) {
    for b in 0..n_bits {
        if hint & (1u8 << b) != 0 {
            let pos = start_bit + b;
            let limb = pos / 64;
            let bit = pos % 64;
            if limb < digits.len() {
                digits[limb] |= 1u64 << bit;
            }
        }
    }
}

/// Compressed signature: 3 of 4 matrix coefficients stored.
///
/// One entry from the second row is dropped and recovered during
/// decompression via the Weil pairing determinant relation. Which entry
/// is dropped depends on `M[0][0]` parity:
///
/// - **M\[0\]\[0\] odd**: drop `M[1][1]`, store `M[1][0]` as `mat_var`.
///   Recover `M[1][1] = (det + M[0][1]·M[1][0]) · M[0][0]⁻¹`.
/// - **M\[0\]\[0\] even**: drop `M[1][0]`, store `M[1][1]` as `mat_var`.
///   Recover `M[1][0] = (M[0][0]·M[1][1] − det) · M[0][1]⁻¹`.
///
/// The Weil pairing determinant gives `E_RSP - bt` bits of precision,
/// leaving exactly 2 unknown bits (HD_EXTRA_TORSION). These are packed
/// into bits 2-3 of the backtracking byte on the wire (129 bytes at
/// Level 1).
///
/// The canonical basis hints for E_chall and E_aux are not stored;
/// they are recomputed from the curves during decompression.
#[derive(Clone, Debug)]
pub struct CompressedSignature<L: SecurityLevel = Level1> {
    /// Montgomery A-coefficient of the auxiliary curve E_aux.
    pub(crate) e_aux_a: Fp2<L>,
    /// Number of backtracking steps (0-3).
    pub(crate) backtracking: u8,
    /// Length of the initial 2-isogeny chain in the response.
    pub(crate) two_resp_length: u8,
    /// Matrix entry M[0][0].
    pub(crate) mat_00: Scalar<L>,
    /// Matrix entry M[0][1].
    pub(crate) mat_01: Scalar<L>,
    /// Third matrix entry: M[1][0] when M[0][0] is odd, M[1][1] when even.
    pub(crate) mat_var: Scalar<L>,
    /// Challenge coefficient (LAMBDA bits).
    pub(crate) chall_coeff: Scalar<L>,
    /// 2-bit hint: bits of the dropped entry above what the Weil pairing
    /// determinant recovers.
    pub(crate) det_hint: u8,
}

impl<L: SecurityLevel> CompressedSignature<L> {
    /// Pack backtracking (2 bits), det_hint (2 bits), and
    /// two_resp_length (4 bits) into one metadata byte.
    ///
    /// Layout: `[trl:4 | det_hint:2 | bt:2]` (LSB first).
    fn pack_meta(&self) -> u8 {
        (self.backtracking & 0x03) | ((self.det_hint & 0x03) << 2) | (self.two_resp_length << 4)
    }

    /// Unpack the metadata byte into (backtracking, det_hint, two_resp_length).
    fn unpack_meta(packed: u8) -> (u8, u8, u8) {
        let backtracking = packed & 0x03;
        let det_hint = (packed >> 2) & 0x03;
        let two_resp_length = (packed >> 4) & 0x0F;
        (backtracking, det_hint, two_resp_length)
    }
}

impl<L: FpBackend> CompressedSignature<L> {
    /// Wire size in bytes.
    ///
    /// Layout: `Fp2 (e_aux) | packed_meta |
    ///          3 × matrix_entry_bytes | LAMBDA/8 (challenge)`
    ///
    /// The packed metadata byte holds backtracking (bits 0-1),
    /// det_hint (bits 2-3), and two_resp_length (bits 4-7).
    /// Canonical basis hints are not stored.
    pub const WIRE_BYTES: usize = <L as SecurityLevel>::Fp2EncodedBytes::USIZE
        + 3 * ((<L as SecurityLevel>::E_RSP as usize + 9) / 8)
        + <L as SecurityLevel>::LAMBDA as usize / 8
        + 1;

    fn matrix_entry_bytes() -> usize {
        (L::E_RSP as usize + 9) / 8
    }

    fn chall_coeff_bytes() -> usize {
        L::LAMBDA as usize / 8
    }

    /// Decode from bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != Self::WIRE_BYTES {
            return Err(Error::InvalidLength);
        }

        let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
        let mat_bytes = Self::matrix_entry_bytes();
        let chall_bytes = Self::chall_coeff_bytes();
        let mut pos = 0;

        let e_aux_a = Fp2::<L>::decode(&bytes[pos..pos + fp2_len]).ok_or(Error::MalformedInput)?;
        pos += fp2_len;

        let (backtracking, det_hint, two_resp_length) = Self::unpack_meta(bytes[pos]);
        pos += 1;

        let mut mat_00 = Scalar::<L>::default();
        decode_digits(mat_00.digits.as_mut_slice(), &bytes[pos..], mat_bytes);
        pos += mat_bytes;

        let mut mat_01 = Scalar::<L>::default();
        decode_digits(mat_01.digits.as_mut_slice(), &bytes[pos..], mat_bytes);
        pos += mat_bytes;

        let mut mat_var = Scalar::<L>::default();
        decode_digits(mat_var.digits.as_mut_slice(), &bytes[pos..], mat_bytes);
        pos += mat_bytes;

        let mut chall_coeff = Scalar::<L>::default();
        decode_digits(
            chall_coeff.digits.as_mut_slice(),
            &bytes[pos..],
            chall_bytes,
        );
        pos += chall_bytes;

        debug_assert_eq!(pos, Self::WIRE_BYTES);

        Ok(Self {
            e_aux_a,
            backtracking,
            two_resp_length,
            mat_00,
            mat_01,
            mat_var,
            chall_coeff,
            det_hint,
        })
    }

    /// Encode to bytes.
    ///
    /// Returns a fixed-size buffer; the first `WIRE_BYTES` bytes are
    /// the meaningful payload.
    pub fn to_bytes(&self) -> [u8; 300] {
        let mut buf = [0u8; 300];
        let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
        let mat_bytes = Self::matrix_entry_bytes();
        let chall_bytes = Self::chall_coeff_bytes();
        let mut pos = 0;

        let enc = self.e_aux_a.encode();
        buf[pos..pos + fp2_len].copy_from_slice(&enc);
        pos += fp2_len;

        buf[pos] = self.pack_meta();
        pos += 1;

        encode_digits(&mut buf[pos..], self.mat_00.digits.as_slice(), mat_bytes);
        pos += mat_bytes;
        encode_digits(&mut buf[pos..], self.mat_01.digits.as_slice(), mat_bytes);
        pos += mat_bytes;
        encode_digits(&mut buf[pos..], self.mat_var.digits.as_slice(), mat_bytes);
        pos += mat_bytes;

        encode_digits(
            &mut buf[pos..],
            self.chall_coeff.digits.as_slice(),
            chall_bytes,
        );
        pos += chall_bytes;

        debug_assert_eq!(pos, Self::WIRE_BYTES);

        buf
    }

    /// Decompress into a standard signature by recovering the dropped entry.
    ///
    /// The Weil pairing determinant gives `pow_dim2 + trl` bits of det(M).
    /// The remaining 2 bits (HD_EXTRA_TORSION) come from `det_hint`.
    pub fn decompress(&self, pk: &PublicKey<L>) -> Result<Signature<L>, Error>
    where
        L: LevelPrecomp,
    {
        let pow_dim2 = L::E_RSP as i32 - self.two_resp_length as i32 - self.backtracking as i32;
        // SECURITY: reject pow_dim2 <= 0 (auxiliary curve unbound, breaks SUF-CMA).
        if pow_dim2 <= 1 {
            return Err(Error::InvalidSignature);
        }
        let pow_dim2 = pow_dim2 as usize;
        let trl = self.two_resp_length as usize;
        let det_precision = pow_dim2 + trl;
        let f = det_precision + HD_EXTRA_TORSION as usize;
        let g = pow_dim2 + HD_EXTRA_TORSION as usize;
        let nw = L::MpLimbs::USIZE;

        let m00_odd = self.mat_00.digits[0] & 1 != 0;

        let pivot = if m00_odd {
            &self.mat_00
        } else {
            if self.mat_01.digits[0] & 1 == 0 {
                return Err(Error::InvalidSignature);
            }
            &self.mat_01
        };

        // Compute det(M) mod 2^det_precision via Weil pairing --
        // Basis hints are not stored; recompute them from the curves.

        let mut e_chall =
            compute_challenge_curve(&self.chall_coeff, self.backtracking, &pk.curve, pk.hint_pk)
                .ok_or(Error::InvalidSignature)?;
        let (b_chall, hint_chall) = basis_to_hint(&mut e_chall, L::F_CHR)?;
        let b_chall = ec_dbl_iter_basis(&b_chall, L::F_CHR as usize - f, &mut e_chall);
        let ppq_chall = xadd(&b_chall.p, &b_chall.q, &b_chall.pmq);
        let omega_f = weil::<L>(f as u32, &b_chall.p, &b_chall.q, &ppq_chall, &mut e_chall);

        let mut e_aux =
            crate::ec::EcCurve::<L>::from_a(&self.e_aux_a).ok_or(Error::InvalidSignature)?;
        let (b_aux, hint_aux) = basis_to_hint(&mut e_aux, L::F_CHR)?;
        let b_aux = ec_dbl_iter_basis(&b_aux, L::F_CHR as usize - g, &mut e_aux);
        let ppq_aux = xadd(&b_aux.p, &b_aux.q, &b_aux.pmq);
        let omega_aux = weil::<L>(g as u32, &b_aux.p, &b_aux.q, &ppq_aux, &mut e_aux);

        let omega_f_inv = omega_f.inv();
        let omega_aux_inv = omega_aux.inv();
        let mut det_digits = [0u64; 8];
        fp2_dlog_2e_pub::<L>(
            &mut det_digits[..nw],
            &omega_aux_inv,
            &omega_f_inv,
            f as u32,
        )
        .ok_or(Error::InvalidSignature)?;
        mp_mod_2exp(&mut det_digits[..nw], det_precision);

        // Recover the dropped entry mod 2^det_precision, then set 2 hint bits

        let mut inv_pivot = Scalar::<L>::default();
        hensel_inv_mod_2e(
            inv_pivot.digits.as_mut_slice(),
            pivot.digits.as_slice(),
            det_precision,
        );

        let (mat_10, mat_11) = if m00_odd {
            let mut product = Scalar::<L>::default();
            mp_mul_mod(
                product.digits.as_mut_slice(),
                self.mat_01.digits.as_slice(),
                self.mat_var.digits.as_slice(),
                det_precision,
            );
            let mut numerator = Scalar::<L>::default();
            let mut carry: u64 = 0;
            for (i, &det_limb) in det_digits.iter().enumerate().take(nw) {
                let sum = (det_limb as u128) + (product.digits[i] as u128) + (carry as u128);
                numerator.digits[i] = sum as u64;
                carry = (sum >> 64) as u64;
            }
            mp_mod_2exp(numerator.digits.as_mut_slice(), det_precision);

            let mut recovered = Scalar::<L>::default();
            mp_mul_mod(
                recovered.digits.as_mut_slice(),
                numerator.digits.as_slice(),
                inv_pivot.digits.as_slice(),
                det_precision,
            );
            set_high_bits(
                recovered.digits.as_mut_slice(),
                self.det_hint,
                det_precision,
                HD_EXTRA_TORSION as usize,
            );
            (self.mat_var.clone(), recovered)
        } else {
            let mut product = Scalar::<L>::default();
            mp_mul_mod(
                product.digits.as_mut_slice(),
                self.mat_00.digits.as_slice(),
                self.mat_var.digits.as_slice(),
                det_precision,
            );
            let mut numerator = Scalar::<L>::default();
            let mut borrow: u64 = 0;
            for (i, &det_limb) in det_digits.iter().enumerate().take(nw) {
                let (diff, b1) = product.digits[i].overflowing_sub(det_limb);
                let (diff2, b2) = diff.overflowing_sub(borrow);
                numerator.digits[i] = diff2;
                borrow = (b1 as u64) + (b2 as u64);
            }
            mp_mod_2exp(numerator.digits.as_mut_slice(), det_precision);

            let mut recovered = Scalar::<L>::default();
            mp_mul_mod(
                recovered.digits.as_mut_slice(),
                numerator.digits.as_slice(),
                inv_pivot.digits.as_slice(),
                det_precision,
            );
            set_high_bits(
                recovered.digits.as_mut_slice(),
                self.det_hint,
                det_precision,
                HD_EXTRA_TORSION as usize,
            );
            (recovered, self.mat_var.clone())
        };

        Ok(Signature {
            e_aux_a: self.e_aux_a.clone(),
            backtracking: self.backtracking,
            two_resp_length: self.two_resp_length,
            mat: [[self.mat_00.clone(), self.mat_01.clone()], [mat_10, mat_11]],
            chall_coeff: self.chall_coeff.clone(),
            hint_aux,
            hint_chall,
        })
    }
}

impl<L: FpBackend + LevelPrecomp> CompressedSignature<L> {
    /// Verify this compressed signature against a public key and message.
    pub fn verify(&self, pk: &PublicKey<L>, msg: &[u8]) -> Result<(), Error> {
        verify_compressed(pk, msg, self)
    }
}

/// Any signature format, auto-detected from wire length.
///
/// Each format has a unique wire size at every security level, so no
/// prefix byte is needed. Use [`AnySignature::from_bytes`] to parse a
/// signature of unknown format, then call [`.verify()`](AnySignature::verify).
#[derive(Clone, Debug)]
pub enum AnySignature<L: SecurityLevel = Level1> {
    Expanded(ExpandedSignature<L>),
    Standard(Signature<L>),
    Compressed(CompressedSignature<L>),
}

impl<L: FpBackend> AnySignature<L> {
    /// Parse a signature, detecting the format from its byte length.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let len = bytes.len();
        if len == ExpandedSignature::<L>::WIRE_BYTES {
            Ok(AnySignature::Expanded(ExpandedSignature::from_bytes(
                bytes,
            )?))
        } else if len == L::SigLen::USIZE {
            let sig = Signature::<L>::from_bytes(bytes)?;
            Ok(AnySignature::Standard(sig))
        } else if len == CompressedSignature::<L>::WIRE_BYTES {
            Ok(AnySignature::Compressed(CompressedSignature::from_bytes(
                bytes,
            )?))
        } else {
            Err(Error::MalformedInput)
        }
    }

    /// The format of this signature.
    pub fn format(&self) -> SignatureFormat {
        match self {
            AnySignature::Expanded(_) => SignatureFormat::Expanded,
            AnySignature::Standard(_) => SignatureFormat::Standard,
            AnySignature::Compressed(_) => SignatureFormat::Compressed,
        }
    }
}

impl<L: FpBackend + LevelPrecomp> AnySignature<L> {
    /// Verify this signature against a public key and message.
    pub fn verify(&self, pk: &PublicKey<L>, msg: &[u8]) -> Result<(), Error> {
        match self {
            AnySignature::Standard(s) => s.verify(pk, msg),
            AnySignature::Expanded(s) => s.verify(pk, msg),
            AnySignature::Compressed(s) => s.verify(pk, msg),
        }
    }
}

impl<L: FpBackend + LevelPrecomp> Signature<L> {
    /// Expand by pre-evaluating the matrix action to get kernel points.
    ///
    /// The expanded form may verify faster because the verifier
    /// skips the biscalar multiplication step; the actual speedup
    /// varies depending on the security level and platform.
    ///
    /// Requires the public key because the matrix action is evaluated
    /// on the challenge curve, which is derived from the public key.
    pub fn expand(&self, pk: &PublicKey<L>) -> Result<ExpandedSignature<L>, Error> {
        let pow_dim2_deg_resp =
            L::E_RSP as i32 - self.two_resp_length as i32 - self.backtracking as i32;
        // SECURITY: reject pow_dim2_deg_resp <= 0 (auxiliary curve unbound, breaks SUF-CMA).
        if pow_dim2_deg_resp <= 1 {
            return Err(Error::InvalidSignature);
        }

        check_canonical_basis_change_matrix(self).ok_or(Error::InvalidSignature)?;

        let mut e_chall =
            compute_challenge_curve(&self.chall_coeff, self.backtracking, &pk.curve, pk.hint_pk)
                .ok_or(Error::InvalidSignature)?;

        let mut b_chall_can = basis_from_hint(&mut e_chall, L::F_CHR, self.hint_chall)
            .ok_or(Error::InvalidSignature)?;

        let dbl_chall = L::F_CHR as usize
            - pow_dim2_deg_resp as usize
            - HD_EXTRA_TORSION as usize
            - self.two_resp_length as usize;
        b_chall_can = crate::ec::point::ec_dbl_iter_basis(&b_chall_can, dbl_chall, &mut e_chall);

        let f =
            pow_dim2_deg_resp as usize + HD_EXTRA_TORSION as usize + self.two_resp_length as usize;
        matrix_scalar_application_even_basis(&mut b_chall_can, &e_chall, &self.mat, f)
            .ok_or(Error::InvalidSignature)?;

        let kernel_is_q = self.two_resp_length > 0
            && mp_is_even::<L>(&self.mat[0][0])
            && mp_is_even::<L>(&self.mat[1][0]);

        let p_aff = b_chall_can.p.x.mul(&b_chall_can.p.z.inv());
        let q_aff = b_chall_can.q.x.mul(&b_chall_can.q.z.inv());

        // Determine whether difference_point's default sqrt sign gives
        // the correct P-Q. If not, the verifier must negate the discriminant.
        let p_pt = EcPoint::new(p_aff.clone(), Fp2::one());
        let q_pt = EcPoint::new(q_aff.clone(), Fp2::one());
        let candidate = difference_point(&p_pt, &q_pt, &e_chall);
        let known_pmq = &b_chall_can.pmq;
        let cross1 = candidate.x.mul(&known_pmq.z);
        let cross2 = candidate.z.mul(&known_pmq.x);
        let pmq_sign_hint = !bool::from(cross1.ct_equal(&cross2));

        Ok(ExpandedSignature {
            e_aux_a: self.e_aux_a.clone(),
            backtracking: self.backtracking,
            two_resp_length: self.two_resp_length,
            chall_coeff: self.chall_coeff.clone(),
            p_chl_x: p_aff,
            q_chl_x: q_aff,
            kernel_is_q,
            pmq_sign_hint,
            hint_aux: self.hint_aux,
            hint_chall: self.hint_chall,
        })
    }

    /// Compress by dropping one second-row entry based on `M[0][0]` parity.
    ///
    /// The Weil pairing recovers `det_precision = pow_dim2 + trl` bits of
    /// the dropped entry. The remaining 2 bits are stored as `det_hint`,
    /// packed into the backtracking byte on the wire.
    pub fn compress(&self) -> CompressedSignature<L> {
        let pow_dim2 =
            L::E_RSP as usize - self.two_resp_length as usize - self.backtracking as usize;
        let det_precision = pow_dim2 + self.two_resp_length as usize;

        let m00_odd = self.mat[0][0].digits[0] & 1 != 0;
        let (kept, dropped) = if m00_odd {
            (&self.mat[1][0], &self.mat[1][1])
        } else {
            (&self.mat[1][1], &self.mat[1][0])
        };

        let mut dropped_shifted = dropped.clone();
        mp_shiftr(dropped_shifted.digits.as_mut_slice(), det_precision);
        let det_hint = dropped_shifted.digits[0] as u8 & 0x03;

        CompressedSignature {
            e_aux_a: self.e_aux_a.clone(),
            backtracking: self.backtracking,
            two_resp_length: self.two_resp_length,
            mat_00: self.mat[0][0].clone(),
            mat_01: self.mat[0][1].clone(),
            mat_var: kept.clone(),
            chall_coeff: self.chall_coeff.clone(),
            det_hint,
        }
    }
}

/// Verify an expanded signature (skips matrix action).
///
/// Uses the pre-evaluated kernel point x-coordinates stored in the
/// expanded signature, bypassing the three biscalar multiplications
/// that the standard format requires.
pub(crate) fn verify_expanded<L: FpBackend + LevelPrecomp>(
    pk: &PublicKey<L>,
    msg: &[u8],
    sig: &ExpandedSignature<L>,
) -> Result<(), Error> {
    let pow_dim2_deg_resp = L::E_RSP as i32 - sig.two_resp_length as i32 - sig.backtracking as i32;

    // SECURITY: reject pow_dim2_deg_resp <= 0 (auxiliary curve unbound, breaks SUF-CMA).
    if pow_dim2_deg_resp <= 1 {
        return Err(Error::InvalidSignature);
    }

    if !crate::ec::EcCurve::<L>::verify_a(&pk.curve.a) {
        return Err(Error::InvalidSignature);
    }

    let mut e_aux = crate::ec::EcCurve::<L>::from_a(&sig.e_aux_a).ok_or(Error::InvalidSignature)?;

    if !crate::verify::verify_canonical_hint::<L>(&mut e_aux, sig.hint_aux) {
        return Err(Error::InvalidSignature);
    }

    let mut e_chall =
        compute_challenge_curve(&sig.chall_coeff, sig.backtracking, &pk.curve, pk.hint_pk)
            .ok_or(Error::InvalidSignature)?;

    // Validate that hint_chall matches the canonical hint for this curve.
    // basis_to_hint normalizes e_chall internally, satisfying the precondition
    // for difference_point_with_hint below.
    let (_, expected_hint_chall) = basis_to_hint(&mut e_chall, L::F_CHR)?;
    if sig.hint_chall != expected_hint_chall {
        return Err(Error::InvalidSignature);
    }

    // Reconstruct b_chall_can from stored affine x-coordinates (Z = 1).
    // P - Q is recomputed via the quadratic formula on the challenge curve,
    // using the sign hint to select the correct root.
    let p_pt = EcPoint::new(sig.p_chl_x.clone(), Fp2::one());
    let q_pt = EcPoint::new(sig.q_chl_x.clone(), Fp2::one());
    let pmq_pt = difference_point_with_hint(&p_pt, &q_pt, &e_chall, sig.pmq_sign_hint)
        .ok_or(Error::InvalidSignature)?;
    let mut b_chall_can = EcBasis {
        p: p_pt,
        q: q_pt,
        pmq: pmq_pt,
    };

    // Auxiliary basis (still requires hint-based generation).
    let b_aux_can =
        basis_from_hint(&mut e_aux, L::F_CHR, sig.hint_aux).ok_or(Error::InvalidSignature)?;
    let dbl_aux = L::F_CHR as usize - pow_dim2_deg_resp as usize - HD_EXTRA_TORSION as usize;
    let b_aux_can = crate::ec::point::ec_dbl_iter_basis(&b_aux_can, dbl_aux, &mut e_aux);

    // Canonical encoding: unused hint bits must be zero to prevent malleability.
    if sig.two_resp_length == 0 && sig.kernel_is_q {
        return Err(Error::InvalidSignature);
    }

    if sig.two_resp_length > 0 {
        two_response_isogeny_verify_inner(
            &mut e_chall,
            &mut b_chall_can,
            sig.kernel_is_q,
            sig.two_resp_length,
            pow_dim2_deg_resp,
        )
        .ok_or(Error::InvalidSignature)?;
    }

    // Theta chain and commitment curve.
    let e_com = compute_commitment_curve_verify(
        &b_chall_can,
        &b_aux_can,
        &e_chall,
        &e_aux,
        pow_dim2_deg_resp,
    )
    .ok_or(Error::InvalidSignature)?;

    // Final hash check.
    let chk_chall = hash_to_challenge(pk, &e_com, msg)?;
    if mp_compare::<L>(&sig.chall_coeff, &chk_chall) != 0 {
        return Err(Error::InvalidSignature);
    }

    Ok(())
}

/// Verify a compressed signature (reconstructs the dropped entry then verifies).
pub(crate) fn verify_compressed<L: FpBackend + LevelPrecomp>(
    pk: &PublicKey<L>,
    msg: &[u8],
    sig: &CompressedSignature<L>,
) -> Result<(), Error> {
    let standard = sig.decompress(pk)?;
    protocols_verify(pk, msg, &standard)
}

#[cfg(test)]
mod tests {
    extern crate std;

    use super::*;
    use crate::params::{Level1, Level3, Level5};

    #[test]
    fn standard_sizes() {
        assert_eq!(<Level1 as SecurityLevel>::SigLen::USIZE, 148);
        assert_eq!(<Level3 as SecurityLevel>::SigLen::USIZE, 224);
        assert_eq!(<Level5 as SecurityLevel>::SigLen::USIZE, 292);
    }

    #[test]
    fn expanded_sizes() {
        // 3 × Fp2EncodedBytes + LAMBDA/8 + 4
        assert_eq!(ExpandedSignature::<Level1>::WIRE_BYTES, 212);
        assert_eq!(ExpandedSignature::<Level3>::WIRE_BYTES, 316);
        assert_eq!(ExpandedSignature::<Level5>::WIRE_BYTES, 420);
    }

    #[test]
    fn compressed_sizes() {
        assert_eq!(CompressedSignature::<Level1>::WIRE_BYTES, 129);
        assert_eq!(CompressedSignature::<Level3>::WIRE_BYTES, 196);
        assert_eq!(CompressedSignature::<Level5>::WIRE_BYTES, 257);
    }

    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn size_ordering() {
        assert!(
            CompressedSignature::<Level1>::WIRE_BYTES < <Level1 as SecurityLevel>::SigLen::USIZE
        );
        assert!(<Level1 as SecurityLevel>::SigLen::USIZE < ExpandedSignature::<Level1>::WIRE_BYTES);
    }

    #[test]
    fn standard_signature_backward_compatible() {
        let sig = Signature::<Level1>::default();
        let bytes = sig.to_bytes();
        let decoded = Signature::<Level1>::from_bytes(&bytes);
        assert!(decoded.is_ok());
    }

    #[test]
    fn any_signature_standard_roundtrip() {
        let sig = Signature::<Level1>::default();
        let bytes = sig.to_bytes();
        let decoded = AnySignature::<Level1>::from_bytes(&bytes);
        assert!(decoded.is_ok());
        match decoded.unwrap() {
            AnySignature::Standard(_) => {}
            _ => panic!("expected Standard variant"),
        }
    }

    #[test]
    fn any_signature_rejects_wrong_length() {
        let bad = [0u8; 200];
        assert!(AnySignature::<Level1>::from_bytes(&bad).is_err());
    }

    #[test]
    fn any_signature_rejects_empty() {
        assert!(AnySignature::<Level1>::from_bytes(&[]).is_err());
    }

    #[test]
    fn any_signature_format_accessor() {
        let sig = Signature::<Level1>::default();
        let any = AnySignature::Standard(sig);
        assert_eq!(any.format(), SignatureFormat::Standard);
    }

    #[test]
    fn expanded_from_bytes_too_short() {
        let data = [0u8; 100];
        assert!(matches!(
            ExpandedSignature::<Level1>::from_bytes(&data),
            Err(Error::InvalidLength)
        ));
    }

    #[test]
    fn expanded_serialization_roundtrip() {
        let sig = ExpandedSignature::<Level1> {
            e_aux_a: Fp2::zero(),
            backtracking: 3,
            two_resp_length: 1,
            chall_coeff: Scalar::default(),
            p_chl_x: Fp2::zero(),
            q_chl_x: Fp2::zero(),
            kernel_is_q: true,
            pmq_sign_hint: true,
            hint_aux: 0xAB,
            hint_chall: 0xCD,
        };
        let buf = sig.to_bytes();
        let decoded = ExpandedSignature::<Level1>::from_bytes(
            &buf[..ExpandedSignature::<Level1>::WIRE_BYTES],
        )
        .expect("roundtrip decode failed");
        assert_eq!(decoded.backtracking, 3);
        assert_eq!(decoded.two_resp_length, 1);
        assert!(decoded.kernel_is_q);
        assert!(decoded.pmq_sign_hint);
        assert_eq!(decoded.hint_aux, 0xAB);
        assert_eq!(decoded.hint_chall, 0xCD);
    }

    #[test]
    fn expanded_kernel_flag_packing() {
        let sig = ExpandedSignature::<Level1> {
            e_aux_a: Fp2::zero(),
            backtracking: 5,
            two_resp_length: 0,
            chall_coeff: Scalar::default(),
            p_chl_x: Fp2::zero(),
            q_chl_x: Fp2::zero(),
            kernel_is_q: false,
            pmq_sign_hint: false,
            hint_aux: 0,
            hint_chall: 0,
        };
        let buf = sig.to_bytes();
        let decoded = ExpandedSignature::<Level1>::from_bytes(
            &buf[..ExpandedSignature::<Level1>::WIRE_BYTES],
        )
        .unwrap();
        assert_eq!(decoded.backtracking, 5);
        assert!(!decoded.kernel_is_q);
    }

    #[test]
    fn compressed_from_bytes_too_short() {
        let data = [0u8; 100];
        assert!(matches!(
            CompressedSignature::<Level1>::from_bytes(&data),
            Err(Error::InvalidLength)
        ));
    }

    #[test]
    fn compressed_serialization_roundtrip() {
        let sig = CompressedSignature::<Level1> {
            e_aux_a: Fp2::zero(),
            backtracking: 2,
            two_resp_length: 1,
            mat_00: Scalar::default(),
            mat_01: Scalar::default(),
            mat_var: Scalar::default(),
            chall_coeff: Scalar::default(),
            det_hint: 0x03,
        };
        let buf = sig.to_bytes();
        let decoded = CompressedSignature::<Level1>::from_bytes(
            &buf[..CompressedSignature::<Level1>::WIRE_BYTES],
        )
        .expect("roundtrip decode failed");
        assert_eq!(decoded.backtracking, 2);
        assert_eq!(decoded.two_resp_length, 1);
        assert_eq!(decoded.det_hint, 0x03);
    }

    #[test]
    fn expand_does_not_panic_on_default_inputs() {
        let sig = Signature::<Level1>::default();
        let pk = PublicKey::<Level1>::default();
        let _ = sig.expand(&pk);
    }

    #[test]
    fn compress_default_signature() {
        let sig = Signature::<Level1>::default();
        let compressed = sig.compress();
        assert_eq!(compressed.backtracking, 0);
        assert_eq!(compressed.two_resp_length, 0);
    }

    #[test]
    fn hensel_inverse_basic() {
        let mut out = [0u64; 4];
        let a = [3u64, 0, 0, 0];
        hensel_inv_mod_2e(&mut out, &a, 64);
        // Verify: a * out ≡ 1 mod 2^64
        let product = (3u128) * (out[0] as u128);
        assert_eq!(product as u64, 1);
    }

    #[test]
    fn hensel_inverse_multiword() {
        let mut out = [0u64; 4];
        let a = [7u64, 0, 0, 0];
        hensel_inv_mod_2e(&mut out, &a, 125);
        // Verify: a * out ≡ 1 mod 2^125
        let mut check = [0u64; 4];
        mp_mul_mod(&mut check, &a, &out, 125);
        assert_eq!(check[0], 1);
        assert_eq!(check[1], 0);
    }
}
