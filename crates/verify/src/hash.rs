use crate::ec::EcCurve;
use crate::fp::FpBackend;
use hybrid_array::typenum::Unsigned;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use crate::types::{PublicKey, Scalar};

/// Compute the challenge scalar as the iterated SHAKE256 hash of the
/// j-invariants of the public key curve and commitment curve, concatenated
/// with the message.
///
/// Uses an iterated re-hash pattern (absorb-squeeze repeated
/// `HASH_ITERATIONS` times) with final modular truncation.
pub fn hash_to_challenge<L: FpBackend>(
    pk: &PublicKey<L>,
    com_curve: &EcCurve<L>,
    message: &[u8],
) -> Result<Scalar<L>, crate::Error> {
    let security_bits = L::LAMBDA as usize;
    let hash_iterations = L::HASH_ITERATIONS as usize;
    let torsion_even_power = L::F_CHR as usize;
    let response_length = L::E_RSP as usize;
    let nwords = L::MpLimbs::USIZE;

    // Encode j-invariants
    let j1 = pk.curve.j_inv();
    let j2 = com_curve.j_inv();
    let j1_bytes = j1.encode();
    let j2_bytes = j2.encode();

    let scalar_byte_len = nwords * 8;
    let mut scalar_buf = [0u8; 64]; // max NWORDS_ORDER=8 → 64 bytes
    debug_assert!(scalar_byte_len <= 64);

    // First pass: absorb j-invariants + message, squeeze
    let mut hash_bytes = (2 * security_bits).div_ceil(8);
    let mut limbs = hash_bytes.div_ceil(8);
    let mut bits = (2 * security_bits) % 64;
    let mut mask = if bits == 0 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    {
        let mut hasher = Shake256::default();
        hasher.update(&j1_bytes);
        hasher.update(&j2_bytes);
        hasher.update(message);
        let mut reader = hasher.finalize_xof();
        reader.read(&mut scalar_buf[..hash_bytes]);
    }
    mask_top_limb(&mut scalar_buf, limbs, mask);

    // Iterations 2..HASH_ITERATIONS-1: re-hash the current scalar
    for _ in 2..hash_iterations {
        let mut hasher = Shake256::default();
        hasher.update(&scalar_buf[..hash_bytes]);
        let mut reader = hasher.finalize_xof();
        reader.read(&mut scalar_buf[..hash_bytes]);
        mask_top_limb(&mut scalar_buf, limbs, mask);
    }

    // Final iteration: absorb current scalar, squeeze with different parameters
    let mut hasher = Shake256::default();
    hasher.update(&scalar_buf[..hash_bytes]);
    let mut reader = hasher.finalize_xof();

    hash_bytes = (torsion_even_power - response_length).div_ceil(8);
    limbs = hash_bytes.div_ceil(8);
    bits = (torsion_even_power - response_length) % 64;
    mask = if bits == 0 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    scalar_buf[..scalar_byte_len].fill(0);
    reader.read(&mut scalar_buf[..hash_bytes]);
    mask_top_limb(&mut scalar_buf, limbs, mask);

    // Truncate to SECURITY_BITS
    mp_mod_2exp(&mut scalar_buf, security_bits, scalar_byte_len);

    // Convert byte buffer to Scalar<L>
    let mut result = Scalar::<L>::default();
    for (i, digit) in result.digits.as_mut_slice().iter_mut().enumerate() {
        let base = i * 8;
        if base + 8 <= scalar_byte_len {
            let bytes: [u8; 8] = scalar_buf[base..base + 8]
                .try_into()
                .map_err(|_| crate::Error::InternalError)?;
            *digit = u64::from_le_bytes(bytes);
        } else if base < scalar_byte_len {
            let mut bytes = [0u8; 8];
            bytes[..scalar_byte_len - base].copy_from_slice(&scalar_buf[base..scalar_byte_len]);
            *digit = u64::from_le_bytes(bytes);
        }
    }

    Ok(result)
}

/// Mask the top limb of a little-endian byte buffer interpreted as u64 limbs.
fn mask_top_limb(buf: &mut [u8], limb_count: usize, mask: u64) {
    if limb_count == 0 {
        return;
    }
    let base = (limb_count - 1) * 8;
    let end = core::cmp::min(base + 8, buf.len());
    if base >= buf.len() {
        return;
    }

    let mut limb_bytes = [0u8; 8];
    limb_bytes[..end - base].copy_from_slice(&buf[base..end]);
    let limb = u64::from_le_bytes(limb_bytes) & mask;
    let masked = limb.to_le_bytes();
    buf[base..end].copy_from_slice(&masked[..end - base]);
}

/// Multiprecision modulo 2^e: zero all bits at position e and above.
fn mp_mod_2exp(buf: &mut [u8], e: usize, total_bytes: usize) {
    let q = e / 8;
    let r = e % 8;

    if q < total_bytes {
        if r != 0 {
            buf[q] &= (1u8 << r) - 1;
        }
        let start = if r != 0 { q + 1 } else { q };
        for b in buf[start..total_bytes].iter_mut() {
            *b = 0;
        }
    }
}
