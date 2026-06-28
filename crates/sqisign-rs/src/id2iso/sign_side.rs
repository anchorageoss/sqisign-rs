//!
//! Contains the id2iso primitives that require quaternion algebra:
//! ideal-to-kernel-dlog translation, endomorphism application,
//! clapotis (find_uv), and fixed-degree isogeny computation.
//!
//! Gated behind the `sign` feature flag.
//!
//! Functions in this module receive level-specific precomputed data via
//! [`SigningPrecomp`] rather than importing it directly.

extern crate alloc;

use crate::quaternion::algebra::{
    quat_alg_conj, quat_alg_make_primitive, quat_alg_mul, quat_alg_normalize,
};
use crate::quaternion::dim2::{ibz_mat_2x2_eval, ibz_mat_2x2_inv_mod};
use crate::quaternion::dim4::{ibz_mat_4x4_eval, quat_qf_eval};
use crate::quaternion::ideal::{
    quat_lideal_conjugate_without_hnf, quat_lideal_copy, quat_lideal_create, quat_lideal_generator,
};
use crate::quaternion::intbig::{
    ibz_bitsize, ibz_div, ibz_gcd, ibz_invmod, ibz_mod, ibz_mod_ui, ibz_pow, ibz_to_digits,
    ibz_two_adic, Ibz,
};
use crate::quaternion::integers::ibz_cornacchia_prime;
use crate::quaternion::lattice::quat_lattice_alg_elem_mul;
use crate::quaternion::lll::{quat_lideal_lideal_mul_reduced, quat_lideal_reduce_basis};
use crate::quaternion::normeq::{quat_change_to_o0_basis, quat_represent_integer};
use crate::quaternion::types::{
    IbzMat2x2, IbzMat4x4, IbzVec2, IbzVec4, QuatAlg, QuatAlgElem, QuatLeftIdeal,
    QuatRepresentIntegerParams,
};
use alloc::vec;
use alloc::vec::Vec;
use num_integer::Integer;
use num_traits::{One, Zero};
use sqisign_verify::ec::pairing::{ec_dlog_2_tate, weil};
use sqisign_verify::ec::point::{ec_biscalar_mul, ec_dbl_iter_basis};
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
#[cfg(not(feature = "dpa-protect"))]
use sqisign_verify::theta::chain::theta_chain_compute_and_eval;
use sqisign_verify::theta::chain::theta_chain_compute_and_eval_randomized;
use sqisign_verify::theta::couple::{copy_bases_to_kernel, double_couple_point_iter};
use sqisign_verify::theta::{ThetaCoupleCurve, ThetaCouplePoint};

use rand::Rng;

use super::sign_precomp::SigningPrecomp;

/// Upper bound on `NWORDS_ORDER` across all security levels (Level 5 = 8).
/// Stack-allocated digit buffers use this size and are sliced to `[..nwords]`.
const MAX_NWORDS: usize = 8;

/// Convert an Ibz to a fixed-size u64 digit array for EC scalar ops.
fn ibz_to_scalar(x: &Ibz, nwords: usize) -> [u64; MAX_NWORDS] {
    let mut out = [0u64; MAX_NWORDS];
    let abs = if x < &Ibz::zero() { -x } else { x.clone() };
    ibz_to_digits(&abs, &mut out[..nwords]);
    out
}

/// Translates an ideal of norm `2ᶠ` into kernel dlog coefficients
/// `[s0, s1]` over the canonical basis of `E0[2ᶠ]`.
pub fn id2iso_ideal_to_kernel_dlogs_even<L: FpBackend>(
    lideal: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
) -> IbzVec2 {
    // SAFETY: the ideal is constructed by the signing algorithm and is always
    // well-formed, so generator extraction cannot fail.
    let alpha = quat_lideal_generator(lideal, &precomp.algebra)
        .expect("invariant: ideal generator extraction must succeed");
    let alpha = quat_alg_conj(&alpha);

    let coeffs = quat_change_to_o0_basis(&alpha);

    let action_gen2 = &precomp.action_matrices[0].gen2;
    let action_gen3 = &precomp.action_matrices[0].gen3;
    let action_gen4 = &precomp.action_matrices[0].gen4;

    let mut mat = IbzMat2x2::default();
    for i in 0..2 {
        mat.0[i][i] = &mat.0[i][i] + &coeffs[0];
        for j in 0..2 {
            let tmp = &action_gen2.0[i][j] * &coeffs[1];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
            let tmp = &action_gen3.0[i][j] * &coeffs[2];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
            let tmp = &action_gen4.0[i][j] * &coeffs[3];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
        }
    }

    let norm = &lideal.norm;
    let mut vec = IbzVec2::default();
    vec[0] = ibz_mod(&mat.0[0][0], norm);
    vec[1] = ibz_mod(&mat.0[1][0], norm);
    let g = ibz_gcd(&vec[0], &vec[1]);
    if g.is_even() {
        vec[0] = ibz_mod(&mat.0[0][1], norm);
        vec[1] = ibz_mod(&mat.0[1][1], norm);
    }

    vec
}

/// Applies a 2×2 BigInt matrix to a torsion basis of `E[2ᶠ]`.
pub fn matrix_application_even_basis<L: FpBackend>(
    bas: &mut EcBasis<L>,
    curve: &EcCurve<L>,
    mat: &mut IbzMat2x2,
    f: u32,
    rng: &mut impl Rng,
) -> Option<()> {
    let nwords = L::NWORDS_ORDER;
    let pow_two = ibz_pow(&Ibz::from(2), f);

    let mut tmp_bas = bas.clone();
    // Side-channel hardening (no-op unless `dpa-protect` is enabled): place the
    // secret-curve basis in a fresh random projective representative before it
    // enters the 2-D Montgomery ladder. Transparent to the affine result.
    crate::id2iso::dpa::maybe_randomize_basis(&mut tmp_bas, rng);

    for i in 0..2 {
        for j in 0..2 {
            mat.0[i][j] = ibz_mod(&mat.0[i][j], &pow_two);
        }
    }

    let s0 = ibz_to_scalar(&mat.0[0][0], nwords);
    let s1 = ibz_to_scalar(&mat.0[1][0], nwords);
    bas.p = ec_biscalar_mul(&s0[..nwords], &s1[..nwords], f as usize, &tmp_bas, curve)?;

    let s0 = ibz_to_scalar(&mat.0[0][1], nwords);
    let s1 = ibz_to_scalar(&mat.0[1][1], nwords);
    bas.q = ec_biscalar_mul(&s0[..nwords], &s1[..nwords], f as usize, &tmp_bas, curve)?;

    let diff0 = ibz_mod(&(&mat.0[0][0] - &mat.0[0][1]), &pow_two);
    let diff1 = ibz_mod(&(&mat.0[1][0] - &mat.0[1][1]), &pow_two);
    let s0 = ibz_to_scalar(&diff0, nwords);
    let s1 = ibz_to_scalar(&diff1, nwords);
    bas.pmq = ec_biscalar_mul(&s0[..nwords], &s1[..nwords], f as usize, &tmp_bas, curve)?;
    Some(())
}

/// Applies an endomorphism of an alternate curve to the even torsion basis.
pub fn endomorphism_application_even_basis<L: FpBackend>(
    bas: &mut EcBasis<L>,
    index_alternate_curve: usize,
    curve: &EcCurve<L>,
    theta: &QuatAlgElem,
    f: u32,
    precomp: &SigningPrecomp<L>,
    rng: &mut impl Rng,
) {
    let order = &precomp.extremal_orders[index_alternate_curve];
    // theta must be contained in the order (precondition).
    let (coeffs, content) = quat_alg_make_primitive(theta, &order.order);
    debug_assert!(content.is_odd());

    let action_gen2 = &precomp.action_matrices[index_alternate_curve].gen2;
    let action_gen3 = &precomp.action_matrices[index_alternate_curve].gen3;
    let action_gen4 = &precomp.action_matrices[index_alternate_curve].gen4;

    let mut mat = IbzMat2x2::default();
    for i in 0..2 {
        mat.0[i][i] = &mat.0[i][i] + &coeffs[0];
        for j in 0..2 {
            let tmp = &action_gen2.0[i][j] * &coeffs[1];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
            let tmp = &action_gen3.0[i][j] * &coeffs[2];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
            let tmp = &action_gen4.0[i][j] * &coeffs[3];
            mat.0[i][j] = &mat.0[i][j] + &tmp;
            mat.0[i][j] = &mat.0[i][j] * &content;
        }
    }

    matrix_application_even_basis(bas, curve, &mut mat, f, rng)
        .expect("invariant: endomorphism matrix application must succeed");
}

/// Computes the ideal whose kernel is generated by
/// `vec2[0]*B0[0] + vec2[1]*B0[1]` where `B0` is the canonical basis of `E0`.
pub fn id2iso_kernel_dlogs_to_ideal_even<L: FpBackend>(
    vec2: &IbzVec2,
    f: u32,
    precomp: &SigningPrecomp<L>,
) -> QuatLeftIdeal {
    let two_pow = if f == L::TORSION_EVEN_POWER {
        precomp.torsion_plus_2power.clone()
    } else {
        ibz_pow(&Ibz::from(2), f)
    };

    let action_i = &precomp.action_matrices[0].i;
    let action_j = &precomp.action_matrices[0].j;
    let action_gen4 = &precomp.action_matrices[0].gen4;

    let mut mat = IbzMat2x2::default();
    mat.0[0][0] = vec2[0].clone();
    mat.0[1][0] = vec2[1].clone();

    let jv = ibz_mat_2x2_eval(action_j, vec2);
    let g4v = ibz_mat_2x2_eval(action_gen4, vec2);
    mat.0[0][1] = ibz_mod(&(&jv[0] + &g4v[0]), &two_pow);
    mat.0[1][1] = ibz_mod(&(&jv[1] + &g4v[1]), &two_pow);

    let (inv, ok) = ibz_mat_2x2_inv_mod(&mat, &two_pow);
    debug_assert!(ok, "matrix inversion failed");

    let iv = ibz_mat_2x2_eval(action_i, vec2);
    let vec = ibz_mat_2x2_eval(&inv, &iv);

    let mut gen = QuatAlgElem {
        denom: Ibz::from(2),
        ..QuatAlgElem::default()
    };
    gen.coord[0] = &(&vec[0] + &vec[0]) + &vec[1]; // 2*a + b
    gen.coord[1] = Ibz::from(-2);
    gen.coord[2] = &vec[1] + &vec[1]; // 2*b
    gen.coord[3] = vec[1].clone();

    let maxord_o0 = precomp.extremal_orders[0].order.clone();

    quat_lideal_create(&gen, &two_pow, &maxord_o0, &precomp.algebra)
        .expect("invariant: kernel dlog ideal has perfect-square index")
}

/// Change-of-basis matrix via Tate pairing (BigInt version).
fn change_of_basis_matrix_tate_impl<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    invert: bool,
    precomp: &SigningPrecomp<L>,
) -> Option<IbzMat2x2> {
    let nwords = L::NWORDS_ORDER;
    let mut x1 = [0u64; MAX_NWORDS];
    let mut x2 = [0u64; MAX_NWORDS];
    let mut x3 = [0u64; MAX_NWORDS];
    let mut x4 = [0u64; MAX_NWORDS];

    if invert {
        ec_dlog_2_tate(
            &mut x1[..nwords],
            &mut x2[..nwords],
            &mut x3[..nwords],
            &mut x4[..nwords],
            b1,
            b2,
            curve,
            f,
            L::TORSION_EVEN_POWER,
            precomp.p_cofactor_for_2f,
        )?;
        mp_invert_matrix(
            &mut x1[..nwords],
            &mut x2[..nwords],
            &mut x3[..nwords],
            &mut x4[..nwords],
            f as usize,
        );
    } else {
        ec_dlog_2_tate(
            &mut x1[..nwords],
            &mut x2[..nwords],
            &mut x3[..nwords],
            &mut x4[..nwords],
            b2,
            b1,
            curve,
            f,
            L::TORSION_EVEN_POWER,
            precomp.p_cofactor_for_2f,
        )?;
    }

    let mut mat = IbzMat2x2::default();
    mat.0[0][0] = crate::quaternion::intbig::ibz_copy_digits(&x1[..nwords]);
    mat.0[1][0] = crate::quaternion::intbig::ibz_copy_digits(&x2[..nwords]);
    mat.0[0][1] = crate::quaternion::intbig::ibz_copy_digits(&x3[..nwords]);
    mat.0[1][1] = crate::quaternion::intbig::ibz_copy_digits(&x4[..nwords]);
    Some(mat)
}

pub fn change_of_basis_matrix_tate<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    precomp: &SigningPrecomp<L>,
) -> Option<IbzMat2x2> {
    change_of_basis_matrix_tate_impl(b1, b2, curve, f, false, precomp)
}

pub fn change_of_basis_matrix_tate_invert<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
    curve: &mut EcCurve<L>,
    f: u32,
    precomp: &SigningPrecomp<L>,
) -> Option<IbzMat2x2> {
    change_of_basis_matrix_tate_impl(b1, b2, curve, f, true, precomp)
}

/// Invert a 2×2 matrix of u64 digit arrays mod 2ᵉ.
fn mp_invert_matrix(r1: &mut [u64], r2: &mut [u64], s1: &mut [u64], s2: &mut [u64], e: usize) {
    let a = crate::quaternion::intbig::ibz_copy_digits(r1);
    let b = crate::quaternion::intbig::ibz_copy_digits(r2);
    let c = crate::quaternion::intbig::ibz_copy_digits(s1);
    let d = crate::quaternion::intbig::ibz_copy_digits(s2);

    let modulus = ibz_pow(&Ibz::from(2), e as u32);

    let det = ibz_mod(&(&(&a * &d) - &(&b * &c)), &modulus);
    // SAFETY: called on change-of-basis matrices from torsion point pairings,
    // whose determinant is always invertible mod 2^e.
    let det_inv = ibz_invmod(&det, &modulus)
        .expect("invariant: torsion basis change matrix must be invertible");

    let ra = ibz_mod(&(&det_inv * &d), &modulus);
    let rb = ibz_mod(&(&(-&det_inv) * &b), &modulus);
    let rc = ibz_mod(&(&(-&det_inv) * &c), &modulus);
    let rd = ibz_mod(&(&det_inv * &a), &modulus);

    ibz_to_digits(&ra, r1);
    ibz_to_digits(&rb, r2);
    ibz_to_digits(&rc, s1);
    ibz_to_digits(&rd, s2);
}

/// Additional bits of 2-power torsion consumed by the (2,2)-isogeny
/// chain beyond the ideal's own degree. The theta-model chain needs
/// `2^(f + HD_EXTRA_TORSION)` torsion to compute a degree-`2ᶠ` isogeny.
const HD_EXTRA_TORSION: u32 = sqisign_verify::theta::HD_EXTRA_TORSION;

/// Swap two elements in a 4×4 matrix using temporary storage.
fn mat_swap(m: &mut IbzMat4x4, r1: usize, c1: usize, r2: usize, c2: usize) {
    let tmp = m.0[r1][c1].clone();
    m.0[r1][c1] = m.0[r2][c2].clone();
    m.0[r2][c2] = tmp;
}

/// Reorders and adjusts signs in an LLL-reduced basis for the special
/// order (j=1728). Ensures the basis has the form `[gamma, i*gamma, beta, i*beta]`.
fn post_lll_basis_treatment(
    gram: &mut IbzMat4x4,
    reduced: &mut IbzMat4x4,
    _norm: &Ibz,
    is_special_order: bool,
) {
    if !is_special_order {
        return;
    }

    if gram.0[0][0] == gram.0[2][2] {
        for i in 0..4 {
            mat_swap(reduced, i, 1, i, 2);
        }
        mat_swap(gram, 0, 2, 0, 1);
        mat_swap(gram, 2, 0, 1, 0);
        mat_swap(gram, 3, 2, 3, 1);
        mat_swap(gram, 2, 3, 1, 3);
        mat_swap(gram, 2, 2, 1, 1);
    } else if gram.0[0][0] == gram.0[3][3] {
        for i in 0..4 {
            mat_swap(reduced, i, 1, i, 3);
        }
        mat_swap(gram, 0, 3, 0, 1);
        mat_swap(gram, 3, 0, 1, 0);
        mat_swap(gram, 2, 3, 2, 1);
        mat_swap(gram, 3, 2, 1, 2);
        mat_swap(gram, 3, 3, 1, 1);
    } else if gram.0[1][1] == gram.0[3][3] {
        for i in 0..4 {
            mat_swap(reduced, i, 1, i, 2);
        }
        mat_swap(gram, 0, 2, 0, 1);
        mat_swap(gram, 2, 0, 1, 0);
        mat_swap(gram, 3, 2, 3, 1);
        mat_swap(gram, 2, 3, 1, 3);
        mat_swap(gram, 2, 2, 1, 1);
    }

    if reduced.0[0][0] != reduced.0[1][1] {
        for i in 0..4 {
            reduced.0[i][1] = -&reduced.0[i][1];
            gram.0[i][1] = -&gram.0[i][1];
            gram.0[1][i] = -&gram.0[1][i];
        }
    }

    if reduced.0[0][2] != reduced.0[1][3] {
        for i in 0..4 {
            reduced.0[i][3] = -&reduced.0[i][3];
            gram.0[i][3] = -&gram.0[i][3];
            gram.0[3][i] = -&gram.0[3][i];
        }
    }
}

struct VecAndNorm {
    vec: IbzVec4,
    norm: Ibz,
}

/// Enumerates lattice vectors in a 4D hypercube of side `2m+1`, filtered
/// to odd-norm vectors. Returns the list of (vector, norm) pairs.
fn enumerate_hypercube(m: i32, gram: &IbzMat4x4, adjusted_norm: &Ibz) -> Vec<VecAndNorm> {
    debug_assert!(m > 0);

    let dim = 2 * m + 1;
    let dim2 = dim * dim;
    let dim3 = dim2 * dim;

    let need_remove_symmetry = gram.0[0][0] == gram.0[1][1] && gram.0[3][3] == gram.0[2][2];

    let mut results = Vec::new();

    let mut x = -m;
    while x <= 0 {
        let mut y = -m;
        while y < m + 1 {
            if x == 0 && y > 0 {
                break;
            }
            let mut z = -m;
            while z < m + 1 {
                if x == 0 && y == 0 && z > 0 {
                    break;
                }
                let mut w = -m;
                while w < m + 1 {
                    if x == 0 && y == 0 && z == 0 && w >= 0 {
                        break;
                    }

                    if (x | y | z | w) & 1 == 0 {
                        w += 1;
                        continue;
                    }
                    if x % 3 == 0 && y % 3 == 0 && z % 3 == 0 && w % 3 == 0 {
                        w += 1;
                        continue;
                    }

                    let check1 = (m + w) + dim * (m + z) + dim2 * (m + y) + dim3 * (m + x);
                    let check2 = (m - z) + dim * (m + w) + dim2 * (m - x) + dim3 * (m + y);
                    let check3 = (m + z) + dim * (m - w) + dim2 * (m + x) + dim3 * (m - y);

                    if !need_remove_symmetry || (check1 <= check2 && check1 <= check3) {
                        let point =
                            IbzVec4([Ibz::from(x), Ibz::from(y), Ibz::from(z), Ibz::from(w)]);

                        let norm_full = quat_qf_eval(gram, &point);
                        let (norm, remain) = ibz_div(&norm_full, adjusted_norm);
                        debug_assert!(remain.is_zero());

                        if ibz_mod_ui(&norm, 2) == 1 {
                            results.push(VecAndNorm { vec: point, norm });
                        }
                    }

                    w += 1;
                }
                z += 1;
            }
            y += 1;
        }
        x += 1;
    }

    results
}

struct FindUvResult {
    u: Ibz,
    v: Ibz,
    index_sol1: usize,
    index_sol2: usize,
}

/// Searches two sorted norm lists to find `d1, d2` and integers `u, v`
/// such that `u*d1 + v*d2 = target`.
#[allow(clippy::too_many_arguments)]
fn find_uv_from_lists(
    target: &Ibz,
    small_norms1: &[Ibz],
    small_norms2: &[Ibz],
    quotients: &[Ibz],
    index1: usize,
    index2: usize,
    is_diagonal: bool,
    number_sum_square: i32,
) -> Option<FindUvResult> {
    let n = target.clone();

    for (i1, norm1) in small_norms1[..index1].iter().enumerate() {
        let adjusted_norm = ibz_mod(&n, norm1);
        let starting_index2 = if is_diagonal { i1 } else { 0 };

        for (i2_offset, norm2) in small_norms2[starting_index2..index2].iter().enumerate() {
            let i2 = starting_index2 + i2_offset;
            let inv = match ibz_invmod(norm2, norm1) {
                Some(inv) => inv,
                None => continue,
            };
            let mut v = ibz_mod(&(&inv * &adjusted_norm), norm1);

            let mut found = false;
            while !found && v < quotients[i2] {
                if number_sum_square > 0 {
                    found = ibz_cornacchia_prime(&Ibz::one(), &v).is_some();
                } else if number_sum_square == 0 {
                    found = true;
                }
                if found {
                    let v_d2 = &v * norm2;
                    let u_times_d1 = &n - &v_d2;
                    debug_assert!(u_times_d1 > Ibz::zero());
                    let (u, remain) = ibz_div(&u_times_d1, norm1);
                    debug_assert!(remain.is_zero());
                    let u_nonzero = u != Ibz::zero();
                    let v_nonzero = v != Ibz::zero();
                    found = found && u_nonzero && v_nonzero;
                    if found && number_sum_square == 2 {
                        found = ibz_cornacchia_prime(&Ibz::one(), &u).is_some();
                    }
                    if found {
                        return Some(FindUvResult {
                            u,
                            v,
                            index_sol1: i1,
                            index_sol2: i2,
                        });
                    }
                }
                if !found {
                    v = &v + norm1;
                }
            }
        }
    }

    None
}

/// Finds `u, v, d1, d2, beta1, beta2` such that `u*d1 + v*d2 = target`,
/// where `beta_i` are elements of the input ideal with norms related to
/// `d_i` (scaled by the ideal norm and any connecting ideal norm when
/// alternate orders are used).
#[allow(clippy::too_many_arguments)]
pub fn find_uv<L: FpBackend>(
    u: &mut Ibz,
    v: &mut Ibz,
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    d1: &mut Ibz,
    d2: &mut Ibz,
    index_alternate_order_1: &mut usize,
    index_alternate_order_2: &mut usize,
    target: &Ibz,
    lideal: &QuatLeftIdeal,
    bpoo: &QuatAlg,
    num_alternate_order: usize,
    precomp: &SigningPrecomp<L>,
) -> Option<()> {
    let num_ideals = num_alternate_order + 1;
    let n = target.clone();

    let mut adjusted_norm = vec![Ibz::zero(); num_ideals];
    let mut gram = vec![IbzMat4x4::default(); num_ideals];
    let mut reduced = vec![IbzMat4x4::default(); num_ideals];
    let mut ideal = vec![QuatLeftIdeal::default(); num_ideals];

    ideal[0] = quat_lideal_copy(lideal);
    let (red0, gram0) = quat_lideal_reduce_basis(&ideal[0], bpoo);
    reduced[0] = red0;
    gram[0] = gram0;

    ideal[0].lattice.basis = reduced[0].clone();
    adjusted_norm[0] = &ideal[0].lattice.denom * &ideal[0].lattice.denom;
    post_lll_basis_treatment(&mut gram[0], &mut reduced[0], &ideal[0].norm, true);

    let mut reduced_id = quat_lideal_copy(&ideal[0]);
    let unit_vec = IbzVec4([Ibz::one(), Ibz::zero(), Ibz::zero(), Ibz::zero()]);
    let mut delta = QuatAlgElem {
        denom: reduced_id.lattice.denom.clone(),
        coord: ibz_mat_4x4_eval(&reduced[0], &unit_vec),
    };

    quat_alg_conj_in_place(&mut delta);
    delta.denom = &delta.denom * &ideal[0].norm;
    reduced_id.lattice = quat_lattice_alg_elem_mul(&reduced_id.lattice, &delta, bpoo);
    let (new_norm, remain) = ibz_div(&gram[0].0[0][0], &adjusted_norm[0]);
    debug_assert!(remain.is_zero());
    reduced_id.norm = new_norm;

    let (conj_ideal, _right_order) = quat_lideal_conjugate_without_hnf(&reduced_id, bpoo);

    for i in 1..num_ideals {
        let alt_ideal = &precomp.connecting_ideals[i];
        let (prod, gram_i) = quat_lideal_lideal_mul_reduced(&conj_ideal, alt_ideal, bpoo);
        ideal[i] = prod;
        reduced[i] = ideal[i].lattice.basis.clone();
        adjusted_norm[i] = &ideal[i].lattice.denom * &ideal[i].lattice.denom;
        gram[i] = gram_i;
        post_lll_basis_treatment(&mut gram[i], &mut reduced[i], &ideal[i].norm, false);
    }

    let m = precomp.finduv_box_size as i32;

    let mut all_vecs: Vec<Vec<VecAndNorm>> = Vec::with_capacity(num_ideals);
    let mut all_quotients: Vec<Vec<Ibz>> = Vec::with_capacity(num_ideals);

    for j in 0..num_ideals {
        let mut vecs = enumerate_hypercube(m, &gram[j], &adjusted_norm[j]);
        vecs.pop();
        vecs.sort_by(|a, b| a.norm.cmp(&b.norm));

        let quotients: Vec<Ibz> = vecs
            .iter()
            .map(|vn| {
                let (q, _) = ibz_div(&n, &vn.norm);
                q
            })
            .collect();

        all_quotients.push(quotients);
        all_vecs.push(vecs);
    }

    for j1 in 0..num_ideals {
        for j2 in j1..num_ideals {
            let is_diago = j1 == j2;
            let norms1: Vec<Ibz> = all_vecs[j1].iter().map(|vn| vn.norm.clone()).collect();
            let norms2: Vec<Ibz> = all_vecs[j2].iter().map(|vn| vn.norm.clone()).collect();

            if let Some(result) = find_uv_from_lists(
                target,
                &norms1,
                &norms2,
                &all_quotients[j2],
                all_vecs[j1].len(),
                all_vecs[j2].len(),
                is_diago,
                0,
            ) {
                *u = result.u;
                *v = result.v;
                *d1 = all_vecs[j1][result.index_sol1].norm.clone();
                *d2 = all_vecs[j2][result.index_sol2].norm.clone();

                beta1.denom = ideal[j1].lattice.denom.clone();
                beta2.denom = ideal[j2].lattice.denom.clone();
                beta1.coord = ibz_mat_4x4_eval(&reduced[j1], &all_vecs[j1][result.index_sol1].vec);
                beta2.coord = ibz_mat_4x4_eval(&reduced[j2], &all_vecs[j2][result.index_sol2].vec);

                if j1 != 0 || j2 != 0 {
                    let (new_denom, remain) = ibz_div(&delta.denom, &lideal.norm);
                    debug_assert!(remain.is_zero());
                    delta.denom = &new_denom * &conj_ideal.norm;
                }

                if j1 != 0 {
                    *beta1 = quat_alg_mul(&delta, beta1, bpoo);
                    quat_alg_normalize(beta1);
                }
                if j2 != 0 {
                    *beta2 = quat_alg_mul(&delta, beta2, bpoo);
                    quat_alg_normalize(beta2);
                }

                if j1 != 0 {
                    quat_alg_conj_in_place(beta1);
                }
                if j2 != 0 {
                    quat_alg_conj_in_place(beta2);
                }

                *index_alternate_order_1 = j1;
                *index_alternate_order_2 = j2;
                return Some(());
            }
        }
    }

    None
}

/// In-place conjugation of a quaternion algebra element.
fn quat_alg_conj_in_place(x: &mut QuatAlgElem) {
    let tmp = quat_alg_conj(x);
    *x = tmp;
}

/// Computes a fixed-degree isogeny from a left ideal of an alternate
/// extremal order, and evaluates a list of points through it.
///
/// Returns the length of the 2-power chain on success, or 0 on failure.
#[allow(clippy::too_many_arguments)]
fn fixed_degree_isogeny_impl<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    lideal: &mut QuatLeftIdeal,
    u_val: &Ibz,
    small: bool,
    codomain: &mut ThetaCoupleCurve<L>,
    points: &mut [ThetaCouplePoint<L>],
    index_alternate_order: usize,
    precomp: &SigningPrecomp<L>,
    rng: &mut impl Rng,
) -> u32 {
    let torsion_even_power = L::TORSION_EVEN_POWER;

    let mut e0 = precomp.endomorphism_curves[index_alternate_order].clone();
    e0.normalize_a24();

    let u_bitsize = ibz_bitsize(u_val);

    let length: u32 = if !small {
        torsion_even_power - HD_EXTRA_TORSION
    } else {
        let p_bitsize = ibz_bitsize(&precomp.algebra.p);
        let bound = precomp.quat_repres_bound_input;
        let l = p_bitsize + bound - u_bitsize;
        debug_assert!(u_bitsize < l);
        debug_assert!(l < torsion_even_power - HD_EXTRA_TORSION);
        l
    };
    debug_assert!(length > 0);

    let two_pow = ibz_pow(&Ibz::from(2), length);
    debug_assert!(&two_pow > u_val);
    debug_assert!(u_val.is_odd());

    let target = &(&two_pow - u_val) * u_val;
    debug_assert!(target.is_odd());
    let order = &precomp.extremal_orders[index_alternate_order];
    let params = QuatRepresentIntegerParams {
        algebra: &precomp.algebra,
        order,
        primality_test_iterations: precomp.quat_primality_num_iter,
    };

    let mut theta = QuatAlgElem::default();
    if !quat_represent_integer(&mut theta, &target, true, &params, rng) {
        return 0;
    }

    *lideal = quat_lideal_create(&theta, u_val, &order.order, &precomp.algebra)
        .expect("invariant: fixed-degree ideal has perfect-square index");

    let b0_two = precomp.endomorphism_bases[index_alternate_order].clone();
    let b0_two = ec_dbl_iter_basis(
        &b0_two,
        (torsion_even_power - length - HD_EXTRA_TORSION) as usize,
        &mut e0,
    );

    let two_pow_extended = ibz_pow(&Ibz::from(2), length + 2);
    let u_inv = ibz_invmod(u_val, &two_pow_extended)
        .expect("invariant: u must be invertible mod 2^(length+2)");

    theta.coord[0] = &theta.coord[0] * &u_inv;
    theta.coord[1] = &theta.coord[1] * &u_inv;
    theta.coord[2] = &theta.coord[2] * &u_inv;
    theta.coord[3] = &theta.coord[3] * &u_inv;

    let mut b0_two_theta = b0_two.clone();
    endomorphism_application_even_basis(
        &mut b0_two_theta,
        index_alternate_order,
        &e0,
        &theta,
        length + HD_EXTRA_TORSION,
        precomp,
        rng,
    );

    let mut e00 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0,
    };

    let ker = copy_bases_to_kernel(&b0_two, &b0_two_theta);

    // The secondary u/v theta chain. Under `dpa-protect`, use the randomized
    // variant so its theta coordinates are in a fresh representative per
    // signature (matching the response chain, which is always randomized).
    #[cfg(not(feature = "dpa-protect"))]
    let chain_result = theta_chain_compute_and_eval(length, &mut e00, &ker, true, points);
    #[cfg(feature = "dpa-protect")]
    let chain_result =
        theta_chain_compute_and_eval_randomized(length, &mut e00, &ker, true, points, rng);
    match chain_result {
        Some(cod) => {
            *codomain = cod;
            length
        }
        None => 0,
    }
}

/// Computes an isogeny of fixed degree `u` starting from an alternate
/// extremal order curve, and evaluates it on a list of points.
///
/// Returns the chain length on success, or 0 on failure.
#[allow(clippy::too_many_arguments)]
pub fn fixed_degree_isogeny_and_eval<L: FpBackend + sqisign_verify::precomp::LevelPrecomp>(
    lideal: &mut QuatLeftIdeal,
    u_val: &Ibz,
    small: bool,
    codomain: &mut ThetaCoupleCurve<L>,
    points: &mut [ThetaCouplePoint<L>],
    index_alternate_order: usize,
    precomp: &SigningPrecomp<L>,
    rng: &mut impl Rng,
) -> u32 {
    fixed_degree_isogeny_impl(
        lideal,
        u_val,
        small,
        codomain,
        points,
        index_alternate_order,
        precomp,
        rng,
    )
}

/// Main ideal-to-isogeny translation via clapotis.
///
/// Finds `u, v, d1, d2` such that `u*d1 + v*d2 = 2^exp` (where `exp`
/// is `TORSION_EVEN_POWER` minus the 2-adic valuation of `gcd(u, v)`),
/// computes two fixed-degree isogenies, constructs a (2,2)-isogeny kernel,
/// and evaluates points through the chain.
#[allow(clippy::too_many_arguments)]
pub fn dim2id2iso_ideal_to_isogeny_clapotis<
    L: FpBackend + sqisign_verify::precomp::LevelPrecomp,
>(
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    u: &mut Ibz,
    v: &mut Ibz,
    d1: &mut Ibz,
    d2: &mut Ibz,
    codomain: &mut EcCurve<L>,
    basis: &mut EcBasis<L>,
    lideal: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
    rng: &mut impl Rng,
) -> Option<()> {
    let torsion_even_power = L::TORSION_EVEN_POWER;

    let mut index_order1: usize = 0;
    let mut index_order2: usize = 0;

    find_uv(
        u,
        v,
        beta1,
        beta2,
        d1,
        d2,
        &mut index_order1,
        &mut index_order2,
        &precomp.torsion_plus_2power,
        lideal,
        &precomp.algebra,
        precomp.num_alternate_extremal_orders,
        precomp,
    )?;

    debug_assert!(d1.is_odd() && d2.is_odd());
    let gcd = ibz_gcd(u, v);
    debug_assert!(gcd != Ibz::zero());
    let exp_gcd = ibz_two_adic(&gcd);
    let exp = torsion_even_power - exp_gcd;
    let (u_div, _) = ibz_div(u, &gcd);
    *u = u_div;
    let (v_div, _) = ibz_div(v, &gcd);
    *v = v_div;
    let mut e1 = precomp.endomorphism_curves[index_order1].clone();
    let bas1 = precomp.endomorphism_bases[index_order1].clone();
    let bas2_orig = precomp.endomorphism_bases[index_order2].clone();

    let mut theta = quat_alg_conj(beta1);
    theta = quat_alg_mul(beta2, &theta, &precomp.algebra);
    theta.denom = &theta.denom * &lideal.norm;

    let mut idealu = QuatLeftIdeal::default();
    let mut fu_codomain = ThetaCoupleCurve::default();

    let zero_pt = EcPoint::new(Fp2::one(), Fp2::zero());
    let mut pushed_points = [
        ThetaCouplePoint {
            p1: bas1.p.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: bas1.q.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: bas1.pmq.clone(),
            p2: zero_pt.clone(),
        },
    ];

    let fu_length = fixed_degree_isogeny_and_eval(
        &mut idealu,
        u,
        true,
        &mut fu_codomain,
        &mut pushed_points,
        index_order1,
        precomp,
        rng,
    );
    if fu_length == 0 {
        return None;
    }

    let bas_u = EcBasis {
        p: pushed_points[0].p1.clone(),
        q: pushed_points[1].p1.clone(),
        pmq: pushed_points[2].p1.clone(),
    };

    let mut e01 = ThetaCoupleCurve {
        e1: fu_codomain.e1.clone(),
        e2: EcCurve::default(),
    };

    let mut idealv = QuatLeftIdeal::default();
    let mut fv_codomain = ThetaCoupleCurve::default();

    let mut pushed_points_v = [
        ThetaCouplePoint {
            p1: bas2_orig.p.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: bas2_orig.q.clone(),
            p2: zero_pt.clone(),
        },
        ThetaCouplePoint {
            p1: bas2_orig.pmq.clone(),
            p2: zero_pt.clone(),
        },
    ];

    let fv_length = fixed_degree_isogeny_and_eval(
        &mut idealv,
        v,
        true,
        &mut fv_codomain,
        &mut pushed_points_v,
        index_order2,
        precomp,
        rng,
    );
    if fv_length == 0 {
        return None;
    }

    let mut bas2 = EcBasis {
        p: pushed_points_v[0].p1.clone(),
        q: pushed_points_v[1].p1.clone(),
        pmq: pushed_points_v[2].p1.clone(),
    };

    let two_pow_full = ibz_pow(&Ibz::from(2), torsion_even_power);
    let mut tmp = d1.clone();
    if index_order2 > 0 {
        tmp = &tmp * &precomp.connecting_ideals[index_order2].norm;
    }
    let tmp_inv =
        ibz_invmod(&tmp, &two_pow_full).expect("invariant: d1 * norm must be invertible mod 2^e");
    theta.coord[0] = &theta.coord[0] * &tmp_inv;
    theta.coord[1] = &theta.coord[1] * &tmp_inv;
    theta.coord[2] = &theta.coord[2] * &tmp_inv;
    theta.coord[3] = &theta.coord[3] * &tmp_inv;

    endomorphism_application_even_basis(
        &mut bas2,
        0,
        &fv_codomain.e1,
        &theta,
        torsion_even_power,
        precomp,
        rng,
    );

    e01.e2 = fv_codomain.e1.clone();

    let ker_t1 = ThetaCouplePoint {
        p1: bas_u.p.clone(),
        p2: bas2.p.clone(),
    };
    let ker_t2 = ThetaCouplePoint {
        p1: bas_u.q.clone(),
        p2: bas2.q.clone(),
    };
    let ker_t1m2 = ThetaCouplePoint {
        p1: bas_u.pmq.clone(),
        p2: bas2.pmq.clone(),
    };

    let dbl_count = torsion_even_power - exp;
    let ker_t1 = double_couple_point_iter(&ker_t1, dbl_count, &e01);
    let ker_t2 = double_couple_point_iter(&ker_t2, dbl_count, &e01);
    let ker_t1m2 = double_couple_point_iter(&ker_t1m2, dbl_count, &e01);

    let zero_pt2 = EcPoint::new(Fp2::one(), Fp2::zero());
    let mut eval_points = [
        ThetaCouplePoint {
            p1: bas_u.p.clone(),
            p2: zero_pt2.clone(),
        },
        ThetaCouplePoint {
            p1: bas_u.q.clone(),
            p2: zero_pt2.clone(),
        },
        ThetaCouplePoint {
            p1: bas_u.pmq.clone(),
            p2: zero_pt2,
        },
    ];

    let ker = sqisign_verify::theta::ThetaKernelCouplePoints {
        t1: ker_t1,
        t2: ker_t2,
        t1m2: ker_t1m2,
    };

    let theta_codomain =
        theta_chain_compute_and_eval_randomized(exp, &mut e01, &ker, false, &mut eval_points, rng)?;

    let t1 = &eval_points[0];
    let t2 = &eval_points[1];
    let t1m2 = &eval_points[2];

    basis.p = t1.p1.clone();
    basis.q = t2.p1.clone();
    basis.pmq = t1m2.p1.clone();
    *codomain = theta_codomain.e1.clone();

    let nwords = L::NWORDS_ORDER;
    let w0 = weil(torsion_even_power, &bas1.p, &bas1.q, &bas1.pmq, &mut e1);
    let w1 = weil(torsion_even_power, &basis.p, &basis.q, &basis.pmq, codomain);

    let mut exp_scalar = &*d1 * &*u;
    exp_scalar = &exp_scalar * &*u;
    exp_scalar = ibz_mod(&exp_scalar, &precomp.torsion_plus_2power);
    let mut digit_d = [0u64; MAX_NWORDS];
    ibz_to_digits(&exp_scalar, &mut digit_d[..nwords]);
    let test_pow = w0.pow_vartime(&digit_d[..nwords]);

    let weil_match = bool::from(w1.ct_equal(&test_pow));

    if !weil_match {
        basis.p = t1.p2.clone();
        basis.q = t2.p2.clone();
        basis.pmq = t1m2.p2.clone();
        *codomain = theta_codomain.e2.clone();
    }

    let mut scalar = &*u * &*d1;
    if index_order1 != 0 {
        scalar = &scalar * &precomp.connecting_ideals[index_order1].norm;
    }
    let scalar_inv = ibz_invmod(&scalar, &precomp.torsion_plus_2power)
        .expect("invariant: scalar must be invertible mod torsion");
    beta1.coord[0] = &beta1.coord[0] * &scalar_inv;
    beta1.coord[1] = &beta1.coord[1] * &scalar_inv;
    beta1.coord[2] = &beta1.coord[2] * &scalar_inv;
    beta1.coord[3] = &beta1.coord[3] * &scalar_inv;

    endomorphism_application_even_basis(
        basis,
        0,
        codomain,
        beta1,
        torsion_even_power,
        precomp,
        rng,
    );

    Some(())
}

/// Wrapper around clapotis for arbitrary isogeny evaluation.
pub fn dim2id2iso_arbitrary_isogeny_evaluation<
    L: FpBackend + sqisign_verify::precomp::LevelPrecomp,
>(
    basis: &mut EcBasis<L>,
    codomain: &mut EcCurve<L>,
    lideal: &QuatLeftIdeal,
    precomp: &SigningPrecomp<L>,
    rng: &mut impl Rng,
) -> Option<()> {
    let mut beta1 = QuatAlgElem::default();
    let mut beta2 = QuatAlgElem::default();
    let mut u = Ibz::zero();
    let mut v = Ibz::zero();
    let mut d1 = Ibz::zero();
    let mut d2 = Ibz::zero();

    dim2id2iso_ideal_to_isogeny_clapotis(
        &mut beta1, &mut beta2, &mut u, &mut v, &mut d1, &mut d2, codomain, basis, lideal, precomp,
        rng,
    )
}
