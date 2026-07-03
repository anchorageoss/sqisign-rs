//! The `Expand` algorithm (Algorithm 1 of ePrint 2026/1169): a deterministic
//! walk of `RRAND` isogenies of degree `2^ucmp` starting from the public-key
//! curve, driven by a hash of `(A(E_pk) ‖ rr)`.
//!
//! For SQIsign-RK, `ucmp = f` (the full even-torsion power) and `RRAND = 2`,
//! so the whole randomization is two `2^f`-isogeny hops (total degree `2^{2f}`,
//! which exceeds the `log p + 2λ` mixing bound). These are hardcoded here and
//! can be parameterized later.
//!
//! This is the public-derivation foundation: it touches only public data
//! (`verify`-crate EC primitives + the public per-level basis/cofactor
//! constants), never the secret key.

use crate::id2iso::sign_precomp::SigningPrecomp;
use alloc::vec::Vec;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::isogeny::ec_eval_even;
use sqisign_verify::ec::point::ec_ladder3pt;
use sqisign_verify::ec::{EcBasis, EcCurve, EcIsogEven};
use sqisign_verify::fp::FpBackend;

/// Number of `2^ucmp`-isogeny hops in the randomization walk. SQIsign-RK uses
/// two hops at every security level.
pub const RRAND: usize = 2;

/// Upper bound on the number of 64-bit limbs in a `2^f` scalar (Level 5 has
/// `f = 500`, i.e. 8 limbs).
pub(crate) const MAX_SCALAR_LIMBS: usize = 8;

/// One kernel scalar `rr_i ∈ [0, 2^f)` as little-endian `u64` limbs, together
/// with the number of significant limbs.
pub(crate) struct Scalar {
    pub(crate) limbs: [u64; MAX_SCALAR_LIMBS],
    pub(crate) nlimbs: usize,
}

/// `H_path`: expand `(A(E_pk) ‖ rr)` into `RRAND` scalars, each in `[0, 2^f)`.
///
/// Deterministic: the same starting curve and `rr` always yield the same
/// scalars, which is what lets `RandPK` be recomputed by anyone from public
/// data alone. Replaces the Sage reference's `randint` draws with a SHAKE256
/// XOF so the derivation is reproducible.
pub(crate) fn h_path_scalars(curve_a_bytes: &[u8], rr: &[u8], f: u32) -> [Scalar; RRAND] {
    let nlimbs = (f as usize).div_ceil(64);
    let bits_in_top = (f % 64) as u64;
    let top_mask = if bits_in_top == 0 {
        u64::MAX
    } else {
        (1u64 << bits_in_top) - 1
    };
    let bytes_per_scalar = nlimbs * 8;

    let mut xof = {
        let mut hasher = Shake256::default();
        hasher.update(b"SQIsign-RK/H_path");
        hasher.update(curve_a_bytes);
        hasher.update(rr);
        hasher.finalize_xof()
    };

    // Fill each scalar from fresh XOF output, then reduce mod 2^f by masking
    // the top limb. Reading contiguous bytes keeps the two scalars independent.
    core::array::from_fn(|_| {
        let mut raw = [0u8; MAX_SCALAR_LIMBS * 8];
        xof.read(&mut raw[..bytes_per_scalar]);
        let mut limbs = [0u64; MAX_SCALAR_LIMBS];
        for (i, limb) in limbs.iter_mut().enumerate().take(nlimbs) {
            let mut word = [0u8; 8];
            word.copy_from_slice(&raw[i * 8..i * 8 + 8]);
            *limb = u64::from_le_bytes(word);
        }
        limbs[nlimbs - 1] &= top_mask;
        Scalar { limbs, nlimbs }
    })
}

/// Result of the Expand walk: the codomain curves `E_1, …, E_RRAND` and, for
/// each, a deterministic `2^f`-torsion basis with its recomputation hint.
pub struct ExpandOutput<L: FpBackend> {
    /// Codomain of each hop; `curves[RRAND-1]` is the randomized curve `E_pk'`.
    pub curves: Vec<EcCurve<L>>,
    /// Deterministic `2^f` basis on each codomain, with its hint byte.
    pub bases: Vec<(EcBasis<L>, u8)>,
    /// The kernel scalars `rr_i` (as big integers), so `RandSK` can rebuild the
    /// same kernels `P_i + [rr_i] Q_i` it must translate to ideals.
    pub scalars: Vec<crate::quaternion::intbig::Ibz>,
}

/// Run the Expand walk from `pk_curve` under randomness `rr`.
///
/// At each hop `i`: form the kernel `K = P_i + [rr_i] Q_i` from the current
/// deterministic basis, push the `2^f`-isogeny with that kernel, and recompute
/// a deterministic basis on the codomain. Returns `None` if any curve fails to
/// yield a canonical basis or the isogeny evaluation fails.
///
/// `precomp` supplies only public per-level constants (the `E0` basis x-coords
/// and `2^f` cofactor used by the deterministic-basis routine).
pub fn expand<L: FpBackend>(
    pk_curve: &EcCurve<L>,
    rr: &[u8],
    precomp: &SigningPrecomp<L>,
) -> Option<ExpandOutput<L>> {
    let f = L::TORSION_EVEN_POWER;

    // Deterministic basis on the starting curve (also normalizes it, giving a
    // canonical `a` to hash and satisfying `ec_ladder3pt`'s a24 precondition).
    let mut curve = pk_curve.clone();
    let (mut basis, _) = ec_curve_to_basis_2f_to_hint(
        &mut curve,
        f,
        precomp.basis_e0_px_bytes,
        precomp.basis_e0_qx_bytes,
        precomp.p_cofactor_for_2f,
        precomp.p_cofactor_for_2f_bitlength,
        f,
    )
    .ok()?;

    let a_bytes = curve.a.encode();
    let scalars = h_path_scalars(a_bytes.as_ref(), rr, f);

    let mut curves = Vec::with_capacity(RRAND);
    let mut bases = Vec::with_capacity(RRAND);
    let mut scalar_ibz = Vec::with_capacity(RRAND);

    for scalar in scalars.iter() {
        scalar_ibz.push(crate::quaternion::intbig::ibz_copy_digits(
            &scalar.limbs[..scalar.nlimbs],
        ));
        // Kernel K = P + [rr_i] Q  (x-only three-point ladder).
        let kernel = ec_ladder3pt(
            &scalar.limbs[..scalar.nlimbs],
            &basis.p,
            &basis.q,
            &basis.pmq,
            &curve,
        )?;

        // The degree-2^f isogeny with kernel <K>.
        let phi = EcIsogEven {
            curve: curve.clone(),
            kernel,
            length: f,
        };
        let mut codomain = curve.clone();
        ec_eval_even(&mut codomain, &phi, &mut [])?;

        // Deterministic basis on the codomain, ready for the next hop.
        let (next_basis, hint) = ec_curve_to_basis_2f_to_hint(
            &mut codomain,
            f,
            precomp.basis_e0_px_bytes,
            precomp.basis_e0_qx_bytes,
            precomp.p_cofactor_for_2f,
            precomp.p_cofactor_for_2f_bitlength,
            f,
        )
        .ok()?;

        curves.push(codomain.clone());
        bases.push((next_basis.clone(), hint));
        curve = codomain;
        basis = next_basis;
    }

    Some(ExpandOutput {
        curves,
        bases,
        scalars: scalar_ibz,
    })
}
