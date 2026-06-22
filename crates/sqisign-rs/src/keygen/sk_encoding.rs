//! `SecretKey::to_bytes` / `SecretKey::from_bytes` round-trip.
//!
//! The wire format packs: ideal norm ‖ generator coords ‖ matrix.

use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::quaternion::algebra::quat_alg_elem_is_zero;
use crate::quaternion::dim2::ibz_mat_2x2_inv_mod;
use crate::quaternion::ideal::{quat_lideal_create, quat_lideal_generator};
use crate::quaternion::intbig::Ibz;
use crate::quaternion::types::{IbzMat2x2, QuatAlgElem};
use hybrid_array::typenum::Unsigned;
use hybrid_array::Array;
use num_bigint::{BigInt, Sign};
use num_traits::{One, Zero};
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_from_hint;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::PublicKey;

use crate::SecretKey;

/// Encode a `BigInt` as `nbytes` little-endian bytes.
///
/// If `signed` is true, negative values use two's complement.
fn ibz_to_bytes(dst: &mut [u8], x: &Ibz, nbytes: usize, signed: bool) {
    debug_assert!(dst.len() >= nbytes);

    if x >= &Ibz::zero() {
        let (_, le_bytes) = x.to_bytes_le();
        let copy_len = le_bytes.len().min(nbytes);
        dst[..copy_len].copy_from_slice(&le_bytes[..copy_len]);
        dst[copy_len..nbytes].fill(0);
    } else {
        debug_assert!(signed, "negative value in unsigned encoding");
        // Two's complement: encode -(x) - 1 then bitwise invert
        let pos = -x - BigInt::one();
        let (_, le_bytes) = pos.to_bytes_le();
        let copy_len = le_bytes.len().min(nbytes);
        dst[..copy_len].copy_from_slice(&le_bytes[..copy_len]);
        dst[copy_len..nbytes].fill(0);
        for b in dst[..nbytes].iter_mut() {
            *b = !*b;
        }
    }
}

/// Decode `nbytes` little-endian bytes into a `BigInt`.
///
/// If `signed` is true, the MSB of the last byte determines sign
/// (two's complement).
fn ibz_from_bytes(src: &[u8], nbytes: usize, signed: bool) -> Ibz {
    debug_assert!(src.len() >= nbytes);
    debug_assert!(nbytes > 0);

    let is_negative = signed && (src[nbytes - 1] >> 7) != 0;

    if is_negative {
        let mut buf = src[..nbytes].to_vec();
        for b in buf.iter_mut() {
            *b = !*b;
        }
        let pos = BigInt::from_bytes_le(Sign::Plus, &buf);
        -(pos + BigInt::one())
    } else {
        BigInt::from_bytes_le(Sign::Plus, &src[..nbytes])
    }
}

impl<L: HasSigningPrecomp + LevelPrecomp> SecretKey<L> {
    /// Encode secret data: `norm ‖ gen[0..3] ‖ mat[0..3]`.
    pub fn to_bytes(&self) -> Result<Array<u8, L::SkLen>, sqisign_verify::Error> {
        let precomp = L::signing_precomp();
        let mut enc = Array::<u8, L::SkLen>::default();
        let mut pos = 0;
        let fp_bytes = L::FpEncodedBytes::USIZE;

        // 1. Ideal norm (unsigned)
        ibz_to_bytes(&mut enc[pos..], &self.secret_ideal.norm, fp_bytes, false);
        pos += fp_bytes;

        // 2. Generator coordinates (signed two's complement)
        let gen = quat_lideal_generator(&self.secret_ideal, &precomp.algebra)
            .ok_or(sqisign_verify::Error::InternalError)?;

        debug_assert!({
            use crate::quaternion::intbig::ibz_gcd;
            ibz_gcd(&gen.denom, &self.secret_ideal.norm).is_one()
        });

        for i in 0..4 {
            ibz_to_bytes(&mut enc[pos..], &gen.coord[i], fp_bytes, true);
            pos += fp_bytes;
        }

        // 3. Basis-change matrix entries (unsigned, mod 2^TORSION_EVEN_POWER)
        let t2p_bytes = precomp.torsion_2power_bytes;
        for row in &self.mat_ba_can_to_ba0_two.0 {
            for entry in row {
                ibz_to_bytes(&mut enc[pos..], entry, t2p_bytes, false);
                pos += t2p_bytes;
            }
        }

        debug_assert!(pos <= L::SkLen::USIZE);
        Ok(enc)
    }

    /// Decode secret data: `norm ‖ gen[0..3] ‖ mat[0..3]`.
    ///
    /// The returned key has default `curve` and `canonical_basis` fields.
    /// Call `populate_from_pk` to reconstruct them from the public key
    /// before signing.
    pub fn from_bytes(enc: &[u8]) -> Result<Self, sqisign_verify::Error> {
        let precomp = L::signing_precomp();
        if enc.len() != L::SkLen::USIZE {
            return Err(sqisign_verify::Error::InvalidLength);
        }
        let mut pos = 0;
        let fp_bytes = L::FpEncodedBytes::USIZE;

        // 1. Ideal norm + generator → reconstruct the ideal
        let norm = ibz_from_bytes(&enc[pos..], fp_bytes, false);
        pos += fp_bytes;

        if norm.is_zero() || norm.sign() == Sign::Minus {
            return Err(sqisign_verify::Error::MalformedInput);
        }

        let mut gen = QuatAlgElem::default();
        for i in 0..4 {
            gen.coord[i] = ibz_from_bytes(&enc[pos..], fp_bytes, true);
            pos += fp_bytes;
        }

        if quat_alg_elem_is_zero(&gen) {
            return Err(sqisign_verify::Error::MalformedInput);
        }

        let parent_order = precomp.extremal_orders[0].order.clone();
        let secret_ideal = quat_lideal_create(&gen, &norm, &parent_order, &precomp.algebra)
            .ok_or(sqisign_verify::Error::MalformedInput)?;

        // The recomputed ideal norm must match the decoded norm; a mismatch is
        // a malformed key.
        if secret_ideal.norm != norm {
            return Err(sqisign_verify::Error::MalformedInput);
        }

        // 2. Basis-change matrix
        let t2p_bytes = precomp.torsion_2power_bytes;
        let mut mat = IbzMat2x2::default();
        for row in &mut mat.0 {
            for entry in row {
                *entry = ibz_from_bytes(&enc[pos..], t2p_bytes, false);
                pos += t2p_bytes;
            }
        }

        // The basis-change matrix must be invertible mod 2^TORSION_EVEN_POWER
        // (signing inverts it); reject a non-invertible one here rather than
        // letting it surface later in the signing path.
        let (_, invertible) = ibz_mat_2x2_inv_mod(&mat, &precomp.torsion_plus_2power);
        if !invertible {
            return Err(sqisign_verify::Error::MalformedInput);
        }

        debug_assert!(pos <= L::SkLen::USIZE);

        Ok(SecretKey {
            curve: EcCurve::default(),
            secret_ideal,
            mat_ba_can_to_ba0_two: mat,
            canonical_basis: EcBasis::new(
                EcPoint::identity(),
                EcPoint::identity(),
                EcPoint::identity(),
            ),
        })
    }

    /// Reconstruct the curve and canonical torsion basis from the public key.
    pub fn populate_from_pk(&mut self, pk: &PublicKey<L>) {
        self.curve = pk.curve().clone();
        let (canonical_basis, _ok) = ec_curve_to_basis_2f_from_hint(
            &mut self.curve,
            L::TORSION_EVEN_POWER,
            pk.hint_pk(),
            L::basis_e0_px_bytes(),
            L::basis_e0_qx_bytes(),
            L::p_cofactor_for_2f(),
            L::p_cofactor_for_2f_bitlength() as usize,
            L::torsion_even_power(),
        )
        .expect("invariant: validated PK produces a valid basis");
        self.canonical_basis = canonical_basis;
    }
}
