//!
//! Wire formats exactly match the v2.0 specification.

use crate::ec::basis::ec_curve_to_basis_2f_to_hint;
use crate::ec::{EcCurve, EcPoint};
use crate::fp::{Fp2, FpBackend};
use crate::params::{Level1, SecurityLevel};
use crate::precomp::LevelPrecomp;
use hybrid_array::typenum::Unsigned;
use hybrid_array::Array;

/// Write a byte slice as lowercase hex into a formatter.
pub(crate) fn fmt_hex(f: &mut core::fmt::Formatter<'_>, bytes: &[u8]) -> core::fmt::Result {
    for &b in bytes {
        write!(f, "{b:02x}")?;
    }
    Ok(())
}

/// Write the little-endian u64 limbs of a [`Scalar`] as a hex byte string.
pub(crate) fn fmt_scalar<L>(f: &mut core::fmt::Formatter<'_>, s: &Scalar<L>) -> core::fmt::Result
where
    L: FpBackend,
{
    for &limb in s.digits.as_slice() {
        for &b in &limb.to_le_bytes() {
            write!(f, "{b:02x}")?;
        }
    }
    Ok(())
}

/// Write an `Fp2` element as a hex byte string of its canonical encoding.
pub(crate) fn fmt_fp2<L: FpBackend>(
    f: &mut core::fmt::Formatter<'_>,
    v: &Fp2<L>,
) -> core::fmt::Result {
    fmt_hex(f, &v.encode())
}

/// A fixed-width multi-precision integer for scalars, matrix entries,
/// and challenge coefficients.
///
/// Represented as `NWORDS_ORDER` little-endian 64-bit limbs.
#[derive(Clone)]
pub struct Scalar<L: SecurityLevel = Level1> {
    pub(crate) digits: Array<u64, L::MpLimbs>,
}

impl<L: FpBackend> core::fmt::Debug for Scalar<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("Scalar(")?;
        fmt_scalar(f, self)?;
        f.write_str(")")
    }
}

impl<L: FpBackend> Scalar<L> {
    #[inline]
    pub fn digits(&self) -> &[u64] {
        self.digits.as_slice()
    }
}

impl<L: FpBackend> Default for Scalar<L> {
    fn default() -> Self {
        Self {
            digits: Array::default(),
        }
    }
}

/// SQIsign public key: a Montgomery curve coefficient plus a torsion hint byte.
///
/// Wire size: 65 bytes (Level 1), 97 bytes (Level 3), 129 bytes (Level 5).
///
/// # Verify a signature
///
/// Use [`pk.verify(msg, &sig)`](signature::Verifier::verify) via the
/// [`Verifier`](signature::Verifier) trait. It accepts any signature type:
/// [`Signature`], [`ExpandedSignature`](crate::ExpandedSignature),
/// [`CompressedSignature`](crate::CompressedSignature), or
/// [`AnySignature`](crate::formats::AnySignature).
///
/// ```
/// use hex_literal::hex;
/// use sqisign_verify::{PublicKey, Signature, Verifier};
///
/// # fn main() -> Result<(), sqisign_verify::Error> {
/// let pk_bytes = hex!(
///     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
///     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
///     "0B"
/// );
/// let sig_bytes = hex!(
///     "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
///     "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
///     "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
///     "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674"
///     "71D86EFF3447018ADB0A6551EE8322AB30010202"
/// );
/// let msg = hex!(
///     "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556A"
///     "C8"
/// );
///
/// let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
/// let sig: Signature = Signature::from_bytes(&sig_bytes)?;
/// pk.verify(&msg, &sig)?;
/// # Ok(())
/// # }
/// ```
///
/// # Decode and re-encode
///
/// ```
/// use hex_literal::hex;
/// use sqisign_verify::PublicKey;
///
/// # fn main() -> Result<(), sqisign_verify::Error> {
/// let pk_bytes = hex!(
///     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
///     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
///     "0B"
/// );
/// let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
/// assert_eq!(pk_bytes, pk.to_bytes().as_slice());
/// # Ok(())
/// # }
/// ```
///
/// # Display
///
/// `Display` prints the wire-format bytes as lowercase hex:
///
/// ```
/// # use hex_literal::hex;
/// # use sqisign_verify::PublicKey;
/// # fn main() -> Result<(), sqisign_verify::Error> {
/// # let pk_bytes = hex!(
/// #     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
/// #     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
/// #     "0B"
/// # );
/// # let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
/// let hex_str = format!("{pk}");
/// assert!(hex_str.starts_with("07ccd214"));
/// # Ok(())
/// # }
/// ```
#[derive(Clone)]
pub struct PublicKey<L: SecurityLevel = Level1> {
    pub(crate) curve: EcCurve<L>,
    pub(crate) hint_pk: u8,
}

impl<L: FpBackend> core::fmt::Debug for PublicKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("PublicKey { curve_a: ")?;
        fmt_fp2(f, &self.curve.a)?;
        write!(f, ", hint: 0x{:02x}", self.hint_pk)?;
        f.write_str(" }")
    }
}

impl<L: FpBackend> core::fmt::Display for PublicKey<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        fmt_hex(f, &self.to_bytes())
    }
}

impl<L: FpBackend> PublicKey<L> {
    #[doc(hidden)]
    #[inline]
    pub fn new(curve: EcCurve<L>, hint_pk: u8) -> Self {
        Self { curve, hint_pk }
    }

    #[inline]
    pub fn curve(&self) -> &EcCurve<L> {
        &self.curve
    }

    #[inline]
    pub fn hint_pk(&self) -> u8 {
        self.hint_pk
    }
}

impl<L: FpBackend> Default for PublicKey<L> {
    fn default() -> Self {
        Self {
            curve: EcCurve::default(),
            hint_pk: 0,
        }
    }
}

/// SQIsign signature (standard wire format).
///
/// Wire size: 148 bytes (Level 1), 224 bytes (Level 3), 292 bytes (Level 5).
///
/// # Verify
///
/// Use [`pk.verify(msg, &sig)`](signature::Verifier::verify) via the
/// [`Verifier`](signature::Verifier) trait:
///
/// ```
/// use hex_literal::hex;
/// use sqisign_verify::{PublicKey, Signature, Verifier};
///
/// # fn main() -> Result<(), sqisign_verify::Error> {
/// let pk_bytes = hex!(
///     "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
///     "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
///     "0B"
/// );
/// let sig_bytes = hex!(
///     "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
///     "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
///     "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
///     "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674"
///     "71D86EFF3447018ADB0A6551EE8322AB30010202"
/// );
/// let msg = hex!(
///     "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556A"
///     "C8"
/// );
///
/// let pk: PublicKey = PublicKey::from_bytes(&pk_bytes)?;
/// let sig: Signature = Signature::from_bytes(&sig_bytes)?;
/// pk.verify(&msg, &sig)?;
/// # Ok(())
/// # }
/// ```
///
/// # Serialize
///
/// ```
/// # use hex_literal::hex;
/// # use sqisign_verify::{PublicKey, Signature};
/// # fn main() -> Result<(), sqisign_verify::Error> {
/// # let sig_bytes = hex!(
/// #     "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
/// #     "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
/// #     "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
/// #     "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674"
/// #     "71D86EFF3447018ADB0A6551EE8322AB30010202"
/// # );
/// # let sig: Signature = Signature::from_bytes(&sig_bytes)?;
/// let wire = sig.to_bytes();
/// assert_eq!(wire.len(), 148); // Level 1 standard signature
/// assert_eq!(wire.as_slice(), &sig_bytes);
/// # Ok(())
/// # }
/// ```
#[derive(Clone)]
pub struct Signature<L: SecurityLevel = Level1> {
    pub(crate) e_aux_a: Fp2<L>,
    pub(crate) backtracking: u8,
    pub(crate) two_resp_length: u8,
    pub(crate) mat: [[Scalar<L>; 2]; 2],
    pub(crate) chall_coeff: Scalar<L>,
    pub(crate) hint_aux: u8,
    pub(crate) hint_chall: u8,
}

impl<L: FpBackend> Signature<L> {
    #[inline]
    pub fn e_aux_a(&self) -> &Fp2<L> {
        &self.e_aux_a
    }

    #[inline]
    pub fn backtracking(&self) -> u8 {
        self.backtracking
    }

    #[inline]
    pub fn two_resp_length(&self) -> u8 {
        self.two_resp_length
    }

    #[inline]
    pub fn mat(&self) -> &[[Scalar<L>; 2]; 2] {
        &self.mat
    }

    #[inline]
    pub fn chall_coeff(&self) -> &Scalar<L> {
        &self.chall_coeff
    }

    #[inline]
    pub fn hint_aux(&self) -> u8 {
        self.hint_aux
    }

    #[inline]
    pub fn hint_chall(&self) -> u8 {
        self.hint_chall
    }
}

#[doc(hidden)]
impl<L: FpBackend> Signature<L> {
    #[inline]
    pub fn set_e_aux_a(&mut self, v: Fp2<L>) {
        self.e_aux_a = v;
    }

    #[inline]
    pub fn set_backtracking(&mut self, v: u8) {
        self.backtracking = v;
    }

    #[inline]
    pub fn set_two_resp_length(&mut self, v: u8) {
        self.two_resp_length = v;
    }

    #[inline]
    pub fn set_hint_aux(&mut self, v: u8) {
        self.hint_aux = v;
    }

    #[inline]
    pub fn set_hint_chall(&mut self, v: u8) {
        self.hint_chall = v;
    }

    #[inline]
    pub fn mat_mut(&mut self) -> &mut [[Scalar<L>; 2]; 2] {
        &mut self.mat
    }

    #[inline]
    pub fn chall_coeff_mut(&mut self) -> &mut Scalar<L> {
        &mut self.chall_coeff
    }

    #[inline]
    pub fn scalar_digits_mut(s: &mut Scalar<L>) -> &mut [u64] {
        s.digits.as_mut_slice()
    }
}

impl<L: FpBackend> Default for Signature<L> {
    fn default() -> Self {
        Self {
            e_aux_a: Fp2::zero(),
            backtracking: 0,
            two_resp_length: 0,
            mat: [
                [Scalar::default(), Scalar::default()],
                [Scalar::default(), Scalar::default()],
            ],
            chall_coeff: Scalar::default(),
            hint_aux: 0,
            hint_chall: 0,
        }
    }
}

impl<L: FpBackend> core::fmt::Debug for Signature<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("Signature { e_aux_a: ")?;
        fmt_fp2(f, &self.e_aux_a)?;
        write!(
            f,
            ", bt: {}, trl: {}, mat: [[",
            self.backtracking, self.two_resp_length
        )?;
        fmt_scalar(f, &self.mat[0][0])?;
        f.write_str(", ")?;
        fmt_scalar(f, &self.mat[0][1])?;
        f.write_str("], [")?;
        fmt_scalar(f, &self.mat[1][0])?;
        f.write_str(", ")?;
        fmt_scalar(f, &self.mat[1][1])?;
        f.write_str("]], chall: ")?;
        fmt_scalar(f, &self.chall_coeff)?;
        write!(
            f,
            ", hint_aux: 0x{:02x}, hint_chall: 0x{:02x} }}",
            self.hint_aux, self.hint_chall
        )
    }
}

impl<L: FpBackend> core::fmt::Display for Signature<L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        fmt_hex(f, &self.to_bytes())
    }
}

pub(crate) fn encode_digits(dst: &mut [u8], src: &[u64], nbytes: usize) {
    let mut pos = 0;
    for &d in src {
        if pos >= nbytes {
            break;
        }
        let bytes = d.to_le_bytes();
        let take = core::cmp::min(8, nbytes - pos);
        dst[pos..pos + take].copy_from_slice(&bytes[..take]);
        pos += take;
    }
}

pub(crate) fn decode_digits(dst: &mut [u64], src: &[u8], nbytes: usize) {
    dst.fill(0);
    for (i, &byte) in src.iter().enumerate().take(nbytes) {
        let digit_idx = i / 8;
        let byte_idx = i % 8;
        if digit_idx < dst.len() {
            dst[digit_idx] |= (byte as u64) << (byte_idx * 8);
        }
    }
}

pub(crate) fn proj_to_bytes<L: FpBackend>(dst: &mut [u8], x: &Fp2<L>, z: &Fp2<L>) -> usize {
    let z_inv = z.inv();
    let affine = x.mul(&z_inv);
    let enc = affine.encode();
    let len = enc.len();
    dst[..len].copy_from_slice(&enc);
    len
}

pub(crate) fn proj_from_bytes<L: FpBackend>(
    x: &mut Fp2<L>,
    z: &mut Fp2<L>,
    src: &[u8],
) -> Result<usize, crate::Error> {
    let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
    *x = Fp2::<L>::decode(&src[..fp2_len]).ok_or(crate::Error::MalformedInput)?;
    *z = Fp2::one();
    Ok(fp2_len)
}

impl<L: FpBackend> PublicKey<L> {
    /// Encode a public key to bytes (wire format).
    pub fn to_bytes(&self) -> Array<u8, L::PkLen> {
        let mut enc = Array::<u8, L::PkLen>::default();
        let mut pos = proj_to_bytes::<L>(&mut enc, &self.curve.a, &self.curve.c);
        enc[pos] = self.hint_pk;
        pos += 1;
        debug_assert_eq!(pos, L::PkLen::USIZE);
        enc
    }

    /// Decode a public key from bytes (wire format).
    ///
    /// Rejects non-canonical `hint_pk` values to ensure each logical
    /// public key has exactly one valid wire encoding.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, crate::Error>
    where
        L: LevelPrecomp,
    {
        // Dimension-4 SQIsignHD public key (Level 1 only): 64 bytes = the curve
        // coefficient `A_pk` with its two torsion-basis hints packed into the
        // spare top bits (see [`crate::hd`]). Detected purely by length - HD
        // exists only at Level 1, identified here by `Fp2EncodedBytes == 64`,
        // and 64 ≠ any dim-2 `PkLen`. The packed hints are not needed by the
        // verifier (it recomputes them), so we keep only the curve and set the
        // dim-2 hint to its canonical value, yielding a `PublicKey` identical to
        // the 65-byte dim-2 encoding of the same curve.
        if L::Fp2EncodedBytes::USIZE == 64 && bytes.len() == crate::hd::PK_WIRE_BYTES {
            let parsed =
                crate::hd::parse_public_key(bytes).map_err(|_| crate::Error::MalformedInput)?;
            let a_bytes = parsed.a_pk.encode();
            let mut pk = PublicKey::<L>::default();
            proj_from_bytes::<L>(&mut pk.curve.a, &mut pk.curve.c, a_bytes.as_ref())?;
            let mut check_curve = pk.curve.clone();
            let (_, canonical_hint) = ec_curve_to_basis_2f_to_hint::<L>(
                &mut check_curve,
                L::F_CHR,
                L::basis_e0_px_bytes(),
                L::basis_e0_qx_bytes(),
                L::p_cofactor_for_2f(),
                L::p_cofactor_for_2f_bitlength() as usize,
                L::torsion_even_power(),
            )
            .map_err(|()| crate::Error::MalformedInput)?;
            pk.hint_pk = canonical_hint;
            return Ok(pk);
        }

        if bytes.len() != L::PkLen::USIZE {
            return Err(crate::Error::InvalidLength);
        }
        let mut pk = PublicKey::<L>::default();
        let mut pos = proj_from_bytes::<L>(&mut pk.curve.a, &mut pk.curve.c, bytes)?;
        pk.hint_pk = bytes[pos];
        pos += 1;
        debug_assert_eq!(pos, L::PkLen::USIZE);

        let mut check_curve = pk.curve.clone();
        let (_, canonical_hint) = ec_curve_to_basis_2f_to_hint::<L>(
            &mut check_curve,
            L::F_CHR,
            L::basis_e0_px_bytes(),
            L::basis_e0_qx_bytes(),
            L::p_cofactor_for_2f(),
            L::p_cofactor_for_2f_bitlength() as usize,
            L::torsion_even_power(),
        )
        .map_err(|()| crate::Error::MalformedInput)?;
        if pk.hint_pk != canonical_hint {
            return Err(crate::Error::MalformedInput);
        }

        Ok(pk)
    }
}

impl<L: FpBackend> Signature<L> {
    /// Number of bytes per matrix entry in the wire format.
    fn matrix_entry_bytes() -> usize {
        (L::E_RSP as usize + 9) / 8
    }

    /// Number of bytes for the challenge coefficient.
    fn chall_coeff_bytes() -> usize {
        L::LAMBDA as usize / 8
    }

    /// Encode a signature to bytes (wire format).
    pub fn to_bytes(&self) -> Array<u8, L::SigLen> {
        let mut enc = Array::<u8, L::SigLen>::default();
        let mut pos = 0;

        // E_aux_A (Fp2 element)
        let fp2_enc = self.e_aux_a.encode();
        enc[pos..pos + fp2_enc.len()].copy_from_slice(&fp2_enc);
        pos += fp2_enc.len();

        // Metadata bytes
        enc[pos] = self.backtracking;
        pos += 1;
        enc[pos] = self.two_resp_length;
        pos += 1;

        // 2x2 scalar matrix
        let mat_bytes = Self::matrix_entry_bytes();
        for row in &self.mat {
            for entry in row {
                encode_digits(&mut enc[pos..], entry.digits.as_slice(), mat_bytes);
                pos += mat_bytes;
            }
        }

        // Challenge coefficient
        let chall_bytes = Self::chall_coeff_bytes();
        encode_digits(
            &mut enc[pos..],
            self.chall_coeff.digits.as_slice(),
            chall_bytes,
        );
        pos += chall_bytes;

        // Hint bytes
        enc[pos] = self.hint_aux;
        pos += 1;
        enc[pos] = self.hint_chall;
        pos += 1;

        debug_assert_eq!(pos, L::SigLen::USIZE);
        enc
    }

    /// Decode a signature from bytes (wire format).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, crate::Error> {
        if bytes.len() != L::SigLen::USIZE {
            return Err(crate::Error::InvalidLength);
        }

        let mut sig = Signature::<L>::default();
        let mut pos = 0;

        // E_aux_A
        let fp2_len = <L as SecurityLevel>::Fp2EncodedBytes::USIZE;
        sig.e_aux_a =
            Fp2::decode(&bytes[pos..pos + fp2_len]).ok_or(crate::Error::MalformedInput)?;
        pos += fp2_len;

        // Metadata
        sig.backtracking = bytes[pos];
        pos += 1;
        sig.two_resp_length = bytes[pos];
        pos += 1;

        // 2x2 scalar matrix
        let mat_bytes = Self::matrix_entry_bytes();
        for row in sig.mat.iter_mut() {
            for entry in row.iter_mut() {
                decode_digits(entry.digits.as_mut_slice(), &bytes[pos..], mat_bytes);
                pos += mat_bytes;
            }
        }

        // Challenge coefficient
        let chall_bytes = Self::chall_coeff_bytes();
        decode_digits(
            sig.chall_coeff.digits.as_mut_slice(),
            &bytes[pos..],
            chall_bytes,
        );
        pos += chall_bytes;

        // Hints
        sig.hint_aux = bytes[pos];
        pos += 1;
        sig.hint_chall = bytes[pos];
        pos += 1;

        debug_assert_eq!(pos, L::SigLen::USIZE);
        Ok(sig)
    }
}

impl<L: FpBackend + LevelPrecomp> TryFrom<&[u8]> for PublicKey<L> {
    type Error = crate::Error;

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: FpBackend> TryFrom<&[u8]> for Signature<L> {
    type Error = crate::Error;

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl<L: FpBackend> From<Signature<L>> for Array<u8, L::SigLen> {
    fn from(sig: Signature<L>) -> Self {
        sig.to_bytes()
    }
}

impl<L: FpBackend> signature::SignatureEncoding for Signature<L>
where
    Array<u8, L::SigLen>: Send + Sync,
{
    type Repr = Array<u8, L::SigLen>;
}

impl<L: FpBackend + crate::precomp::LevelPrecomp> signature::Verifier<Signature<L>>
    for PublicKey<L>
{
    fn verify(&self, msg: &[u8], sig: &Signature<L>) -> Result<(), signature::Error> {
        crate::verify::protocols_verify(self, msg, sig).map_err(|_| signature::Error::new())
    }
}

/// Encode an elliptic curve `(A:C)` to bytes (affine A-coefficient).
pub fn ec_curve_to_bytes<L: FpBackend>(dst: &mut [u8], curve: &EcCurve<L>) -> usize {
    proj_to_bytes::<L>(dst, &curve.a, &curve.c)
}

/// Decode an elliptic curve from bytes (affine A-coefficient, C=1).
pub fn ec_curve_from_bytes<L: FpBackend>(
    curve: &mut EcCurve<L>,
    src: &[u8],
) -> Result<usize, crate::Error> {
    *curve = EcCurve::default();
    proj_from_bytes::<L>(&mut curve.a, &mut curve.c, src)
}

/// Encode a projective point `(X:Z)` to bytes (affine X-coordinate).
pub fn ec_point_to_bytes<L: FpBackend>(dst: &mut [u8], point: &EcPoint<L>) -> usize {
    proj_to_bytes::<L>(dst, &point.x, &point.z)
}

/// Decode a projective point from bytes (affine X-coordinate, Z=1).
pub fn ec_point_from_bytes<L: FpBackend>(
    point: &mut EcPoint<L>,
    src: &[u8],
) -> Result<usize, crate::Error> {
    proj_from_bytes::<L>(&mut point.x, &mut point.z, src)
}
