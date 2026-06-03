//!
//! Contains id2iso primitives needed by `sqisign-verify` that do NOT
//! depend on the quaternion stack.
//!
//! Functions here operate on EC types (EcPoint, EcBasis, EcCurve)
//! and fixed-size u64 digit arrays. No quaternion ideals or orders.
//! No heap allocation.

use hybrid_array::Array;
use sqisign_verify::ec::pairing::ec_dlog_2_tate;
use sqisign_verify::ec::point::ec_biscalar_mul;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::FpBackend;
use sqisign_verify::params::SecurityLevel;

type BasisChangeMat<L> = [[Array<u64, <L as SecurityLevel>::MpLimbs>; 2]; 2];

/// Scalar multiplication `[x]P + [y]Q` where `x` and `y` are given
/// as u64 limb arrays and `P, Q` form a basis of `E[2ᶠ]`.
pub fn ec_biscalar_mul_ibz<L: FpBackend>(
    scalar0: &[u64],
    scalar1: &[u64],
    f: u32,
    basis: &EcBasis<L>,
    curve: &EcCurve<L>,
) -> Option<EcPoint<L>> {
    ec_biscalar_mul(scalar0, scalar1, f as usize, basis, curve)
}

/// Applies a 2×2 scalar matrix to a torsion basis of `E[2ᶠ]`, in place.
///
/// Matrix layout: `mat[row][col]`, each entry is a `Array<u64, L::MpLimbs>`.
/// For matrix `[[a, c], [b, d]]`, computes:
/// - `R = [a]P + [b]Q`
/// - `S = [c]P + [d]Q`
/// - `R-S = [a-c]P + [b-d]Q`
pub fn matrix_application_even_basis<L: FpBackend>(
    bas: &mut EcBasis<L>,
    curve: &EcCurve<L>,
    mat: &[[Array<u64, L::MpLimbs>; 2]; 2],
    f: u32,
) -> Option<()> {
    let tmp_bas = bas.clone();
    let n = mat[0][0].len();

    // R = [a]P + [b]Q
    bas.p = ec_biscalar_mul(&mat[0][0], &mat[1][0], f as usize, &tmp_bas, curve)?;

    // S = [c]P + [d]Q
    bas.q = ec_biscalar_mul(&mat[0][1], &mat[1][1], f as usize, &tmp_bas, curve)?;

    // R - S = [a-c]P + [b-d]Q
    let mut scalar0 = Array::<u64, L::MpLimbs>::default();
    let mut scalar1 = Array::<u64, L::MpLimbs>::default();
    mp_sub(
        scalar0.as_mut_slice(),
        mat[0][0].as_slice(),
        mat[0][1].as_slice(),
        n,
    );
    mp_mod_2exp(scalar0.as_mut_slice(), f as usize);
    mp_sub(
        scalar1.as_mut_slice(),
        mat[1][0].as_slice(),
        mat[1][1].as_slice(),
        n,
    );
    mp_mod_2exp(scalar1.as_mut_slice(), f as usize);

    bas.pmq = ec_biscalar_mul(&scalar0, &scalar1, f as usize, &tmp_bas, curve)?;
    Some(())
}

/// Change-of-basis matrix via Tate pairing.
///
/// Returns 2×2 matrix `[[r1, s1], [r2, s2]]` as `Array<u64, L::MpLimbs>`.
/// Matrix semantics: `(mat * v) . B2 = v . B1`.
#[allow(clippy::too_many_arguments)]
pub fn change_of_basis_matrix_tate<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    torsion_even_power: u32,
    cofactor: &[u64],
) -> Option<BasisChangeMat<L>> {
    change_of_basis_matrix_tate_impl::<L>(b1, b2, curve, f, torsion_even_power, cofactor, false)
}

/// Change-of-basis matrix via Tate pairing, with inversion.
///
/// Returns 2×2 matrix `[[r1, s1], [r2, s2]]` as `Array<u64, L::MpLimbs>`.
/// Matrix semantics: `(mat * v) . B2 = [2^(e−f)] * v . B1`.
#[allow(clippy::too_many_arguments)]
pub fn change_of_basis_matrix_tate_invert<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    torsion_even_power: u32,
    cofactor: &[u64],
) -> Option<BasisChangeMat<L>> {
    change_of_basis_matrix_tate_impl::<L>(b1, b2, curve, f, torsion_even_power, cofactor, true)
}

#[allow(clippy::too_many_arguments)]
fn change_of_basis_matrix_tate_impl<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    torsion_even_power: u32,
    cofactor: &[u64],
    invert: bool,
) -> Option<BasisChangeMat<L>> {
    let mut x1 = Array::<u64, L::MpLimbs>::default();
    let mut x2 = Array::<u64, L::MpLimbs>::default();
    let mut x3 = Array::<u64, L::MpLimbs>::default();
    let mut x4 = Array::<u64, L::MpLimbs>::default();

    if invert {
        ec_dlog_2_tate(
            x1.as_mut_slice(),
            x2.as_mut_slice(),
            x3.as_mut_slice(),
            x4.as_mut_slice(),
            b1,
            b2,
            curve,
            f,
            torsion_even_power,
            cofactor,
        )?;
        let n = x1.len();
        mp_invert_matrix(
            x1.as_mut_slice(),
            x2.as_mut_slice(),
            x3.as_mut_slice(),
            x4.as_mut_slice(),
            f as usize,
            n,
        );
    } else {
        ec_dlog_2_tate(
            x1.as_mut_slice(),
            x2.as_mut_slice(),
            x3.as_mut_slice(),
            x4.as_mut_slice(),
            b2,
            b1,
            curve,
            f,
            torsion_even_power,
            cofactor,
        )?;
    }

    // mat[row][col]: [[x1, x3], [x2, x4]]
    Some([[x1, x3], [x2, x4]])
}

/// Multi-precision subtraction: `out = a − b` (mod 2^(64*nwords)).
fn mp_sub(out: &mut [u64], a: &[u64], b: &[u64], nwords: usize) {
    let mut borrow: u64 = 0;
    for i in 0..nwords {
        let (diff, b1) = a[i].overflowing_sub(b[i]);
        let (diff, b2) = diff.overflowing_sub(borrow);
        out[i] = diff;
        borrow = (b1 as u64) | (b2 as u64);
    }
}

/// Reduce `x` modulo `2ᵉ` in place.
fn mp_mod_2exp(x: &mut [u64], e: usize) {
    let full_words = e / 64;
    let remaining_bits = e % 64;
    for w in x.iter_mut().skip(full_words + 1) {
        *w = 0;
    }
    if full_words < x.len() {
        if remaining_bits > 0 {
            x[full_words] &= (1u64 << remaining_bits) - 1;
        } else {
            x[full_words] = 0;
        }
    }
}

/// Multi-precision multiplication: `out = a * b` (lower nwords only).
fn mp_mul(out: &mut [u64], a: &[u64], b: &[u64], nwords: usize) {
    out[..nwords].fill(0);
    for i in 0..nwords {
        let mut carry: u64 = 0;
        for j in 0..nwords {
            if i + j >= nwords {
                break;
            }
            let product = (a[i] as u128) * (b[j] as u128) + (out[i + j] as u128) + (carry as u128);
            out[i + j] = product as u64;
            carry = (product >> 64) as u64;
        }
    }
}

/// Negate `x` in place (two's complement, mod 2^(64×nwords)).
fn mp_neg(x: &mut [u64], nwords: usize) {
    let mut carry: u64 = 1;
    for limb in x.iter_mut().take(nwords) {
        let (val, c) = (!*limb).overflowing_add(carry);
        *limb = val;
        carry = c as u64;
    }
}

/// Compute `x⁻¹ mod 2ᵉ` using Hensel lifting.
///
/// Requires `x` to be odd (invertible mod 2).
fn mp_inv_2e(out: &mut [u64], x: &[u64], e: usize, nwords: usize) {
    debug_assert!(x[0] & 1 == 1, "mp_inv_2e: input must be odd");

    // Hensel lift: if inv*x ≡ 1 mod 2^k, then inv*(2 - inv*x) ≡ 1 mod 2^(2k).
    // Base: x is odd ⇒ x ≡ 1 mod 2, so inv=1 gives x*1 ≡ 1 mod 2.

    // Use stack buffers of max size. For Level 1, nwords=4 (256 bits).
    // For Level 5, nwords=8 (512 bits). 8 is a safe max.
    const MAX_NWORDS: usize = 8;
    debug_assert!(nwords <= MAX_NWORDS);

    let mut inv = [0u64; MAX_NWORDS];
    inv[0] = 1;

    let mut tmp = [0u64; MAX_NWORDS];
    let mut two_minus = [0u64; MAX_NWORDS];

    let mut k = 1;
    while k < e {
        // tmp = x * inv
        mp_mul(&mut tmp, x, &inv, nwords);
        // two_minus = 2 - tmp
        let mut two_val = [0u64; MAX_NWORDS];
        two_val[0] = 2;
        mp_sub(&mut two_minus, &two_val, &tmp, nwords);
        // inv = inv * two_minus
        mp_mul(&mut tmp, &inv, &two_minus, nwords);
        inv[..nwords].copy_from_slice(&tmp[..nwords]);
        let mask_bits = core::cmp::min(k * 2, e);
        mp_mod_2exp(&mut inv, mask_bits);
        k *= 2;
    }
    mp_mod_2exp(&mut inv, e);
    out[..nwords].copy_from_slice(&inv[..nwords]);
}

/// Invert a 2×2 matrix of u64 digit arrays mod 2ᵉ.
///
/// Given `[[r1, r2], [s1, s2]]` = `[[a, b], [c, d]]`,
/// computes `(1/det) * [[d, −b], [−c, a]]` mod 2ᵉ in place.
pub fn mp_invert_matrix(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    e: usize,
    nwords: usize,
) {
    const MAX_NWORDS: usize = 8;
    debug_assert!(nwords <= MAX_NWORDS);

    // Round e up to next power of 2
    let mut p = 1;
    while (1 << p) < e {
        p += 1;
    }
    let w = 1 << p;

    let mut det = [0u64; MAX_NWORDS];
    let mut tmp = [0u64; MAX_NWORDS];
    let mut det_inv = [0u64; MAX_NWORDS];

    // det = a*d - b*c
    let mut bc = [0u64; MAX_NWORDS];
    mp_mul(&mut tmp, r1, s2, nwords); // tmp = a*d
    mp_mul(&mut bc, r2, s1, nwords); // bc = b*c
    mp_sub(&mut det, &tmp, &bc, nwords);

    mp_inv_2e(&mut det_inv, &det, e, nwords);

    let mut resa = [0u64; MAX_NWORDS];
    let mut resb = [0u64; MAX_NWORDS];
    let mut resc = [0u64; MAX_NWORDS];
    let mut resd = [0u64; MAX_NWORDS];

    mp_mul(&mut resa, &det_inv, s2, nwords);
    mp_mul(&mut resb, &det_inv, r2, nwords);
    mp_mul(&mut resc, &det_inv, s1, nwords);
    mp_mul(&mut resd, &det_inv, r1, nwords);

    mp_neg(&mut resb, nwords);
    mp_neg(&mut resc, nwords);

    mp_mod_2exp(&mut resa, w);
    mp_mod_2exp(&mut resb, w);
    mp_mod_2exp(&mut resc, w);
    mp_mod_2exp(&mut resd, w);

    r1[..nwords].copy_from_slice(&resa[..nwords]);
    r2[..nwords].copy_from_slice(&resb[..nwords]);
    s1[..nwords].copy_from_slice(&resc[..nwords]);
    s2[..nwords].copy_from_slice(&resd[..nwords]);
}
