//! `RandSK`; secret key randomization (Algorithm 3 of ePrint 2026/1169,
//! specialized). Transfers the endomorphism ring through the Expand walk.
//!
//! Because Expand walks in `2^f` steps, each hop's ideal translates in one
//! shot, so the transfer is interleaved with the walk (no post-hoc piecewise
//! `IsoToEnd`). This mirrors the SageMath reference `path_to_new_sk`, mapped
//! onto our SQIsign2D-West Deuring layer:
//!
//! at hop `i`, with running secret ideal `I` (norm `N`) and change-of-basis
//! matrix `M`:
//!   1. pull the kernel `K = P_i + [rr_i] Q_i` back through `hat(phi_I)` to
//!      `E0`-basis coordinates: `vec = [N] · M · (1, rr_i) (mod 2^f)`;
//!   2. translate to an `O0`-ideal `IK` via `id2iso_kernel_dlogs_to_ideal_even`;
//!   3. compose with the running secret: `I·IK = I ∩ IK` (as signing composes
//!      the challenge and secret ideals);
//!   4. reduce to a prime-norm equivalent `J`;
//!   5. push the `E0` basis through `phi_J` (`dim2id2iso`), align its codomain
//!      to the Expand curve `E_i` via an isomorphism, and recompute `M`.
//!
//! The final `(E_pk', J, M, basis)` is a `SecretKey` in exactly the layout the
//! base `Sign` consumes.

use super::expand::{expand, RRAND};
use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::id2iso::sign_side::{
    change_of_basis_matrix_tate, dim2id2iso_arbitrary_isogeny_evaluation,
    id2iso_kernel_dlogs_to_ideal_even,
};
use crate::keygen::SecretKey;
use crate::quaternion::dim2::ibz_mat_2x2_eval;
use crate::quaternion::ideal::quat_lideal_inter;
use crate::quaternion::intbig::{ibz_mod, Ibz};
use crate::quaternion::lll::quat_lideal_prime_norm_reduced_equivalent;
use crate::quaternion::types::IbzVec2;
use num_traits::One;
use sqisign_verify::ec::isogeny::{ec_iso_eval, ec_isomorphism};
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::types::PublicKey;

/// Derive a new secret key corresponding to the public key produced by
/// [`rand_pk`](super::rand_pk)`(pk, rr)`.
///
/// Requires the secret key for `pk`. The output is a valid [`SecretKey`] for
/// `rand_pk(pk, rr)`; the base `Sign` consumes it unchanged.
///
/// Deterministic in `(sk, pk, rr)`: the ideal-to-isogeny and prime-norm steps
/// are driven by a SHAKE256 stream seeded from `pk` and `rr`, so no
/// caller-supplied RNG is needed and repeated calls agree.
///
/// ```no_run
/// # #[cfg(feature = "sqisign-rk")] {
/// use sqisign_rs::keygen::keypair;
/// use sqisign_rs::params::Level1;
/// use sqisign_rs::sqisign_rk::{rand_pk, rand_sk};
///
/// let mut rng = rand::thread_rng();
/// let (pk, sk) = keypair::<Level1>(&mut rng);
/// let pk_child = rand_pk(&pk, b"child-0");
/// let sk_child = rand_sk(&sk, &pk, b"child-0");
/// # }
/// ```
pub fn rand_sk<L: HasSigningPrecomp + LevelPrecomp>(
    sk: &SecretKey<L>,
    pk: &PublicKey<L>,
    rr: &[u8],
) -> SecretKey<L> {
    let precomp = L::signing_precomp();

    // Expand is deterministic and RNG-free; compute it once.
    let walk = expand(pk.curve(), rr, &precomp)
        .expect("invariant: Expand walk must succeed on a valid public key");

    // Deterministic internal randomness for the ideal-to-isogeny / prime-norm
    // sampling. The mathematical result is invariant to these draws, so any
    // successful attempt is a valid secret key for the same E_pk'.
    let a_bytes = pk.curve().a.encode();
    let ctx: [&[u8]; 2] = [a_bytes.as_ref(), rr];
    let mut rng = super::rng::seed_rng(b"SQIsign-RK/RandSK", &ctx);

    for _ in 0..64 {
        if let Some(sk_new) = rand_sk_from_walk(sk, &walk, &precomp, &mut rng) {
            return sk_new;
        }
    }
    panic!("invariant: RandSK failed to converge after 64 attempts");
}

/// The shared per-hop ideal transfer, common to dim-2 and dim-4 RandSK.
///
/// Given the running secret ideal `I` and change-of-basis matrix `M`, the
/// kernel scalar `rr_i` for this hop, and the (dimension-specific) target curve
/// `e_i` with its canonical basis `canonical_basis_i`, returns the composed
/// prime-norm ideal `J` and the new change-of-basis matrix `M'` on `e_i`.
///
/// This is dimension-agnostic: it operates only on `(QuatLeftIdeal, IbzMat2x2)`
/// plus the caller-supplied `(curve, canonical basis)`. The only thing that
/// differs between dim-2 and dim-4 is how that curve/basis (and its hint) is
/// derived and wrapped; see the two `expand` variants.
///
/// Returns `None` if the prime-norm reduction, ideal-to-isogeny, or codomain
/// alignment fails; the caller retries with fresh randomness.
pub(crate) fn ideal_transfer_step<L: HasSigningPrecomp + LevelPrecomp>(
    secret_ideal: &crate::quaternion::types::QuatLeftIdeal,
    mat: &crate::quaternion::types::IbzMat2x2,
    kernel_scalar: &Ibz,
    e_i: &EcCurve<L>,
    canonical_basis_i: &EcBasis<L>,
    precomp: &crate::id2iso::sign_precomp::SigningPrecomp<L>,
    rng: &mut impl rand::Rng,
) -> Option<(
    crate::quaternion::types::QuatLeftIdeal,
    crate::quaternion::types::IbzMat2x2,
)> {
    let f = L::TORSION_EVEN_POWER;
    let two_pow_f = precomp.torsion_plus_2power.clone();
    let n = secret_ideal.norm.clone();

    // 1. Pull the kernel (1, rr_i) back through hat(phi_I): [N] · M · (1, rr_i) mod 2^f.
    let pulled = ibz_mat_2x2_eval(mat, &IbzVec2([Ibz::one(), kernel_scalar.clone()]));
    let c1 = ibz_mod(&(&n * &pulled[0]), &two_pow_f);
    let c2 = ibz_mod(&(&n * &pulled[1]), &two_pow_f);

    // 2. Kernel coordinates → O0-ideal of norm 2^f.
    let pullback_ik = id2iso_kernel_dlogs_to_ideal_even(&IbzVec2([c1, c2]), f, precomp);

    // 3. Compose with the running secret ideal (intersection, as signing
    //    composes challenge ∩ secret).
    let composed = quat_lideal_inter(&pullback_ik, secret_ideal);

    // 4. Reduce to a prime-norm equivalent ideal.
    let mut j = composed;
    if !quat_lideal_prime_norm_reduced_equivalent(
        &mut j,
        &precomp.algebra,
        precomp.quat_primality_num_iter,
        precomp.quat_equiv_bound_coeff,
        rng,
    ) {
        return None;
    }

    // 5. Ideal → isogeny: push the E0 basis through phi_J, then align the
    //    (possibly different) Montgomery model onto the target curve e_i.
    let mut basis_ej = EcBasis::new(
        EcPoint::identity(),
        EcPoint::identity(),
        EcPoint::identity(),
    );
    let mut codomain_ej = EcCurve::<L>::default();
    dim2id2iso_arbitrary_isogeny_evaluation(&mut basis_ej, &mut codomain_ej, &j, precomp, rng)?;

    let eta = ec_isomorphism(&codomain_ej, e_i)?;
    ec_iso_eval(&mut basis_ej.p, &eta);
    ec_iso_eval(&mut basis_ej.q, &eta);
    ec_iso_eval(&mut basis_ej.pmq, &eta);

    // Recompute the change-of-basis matrix on e_i (canonical basis → pushed-E0
    // basis), matching the keygen convention.
    let mut curve = e_i.clone();
    let m_new = change_of_basis_matrix_tate(canonical_basis_i, &basis_ej, &mut curve, f, precomp)?;

    Some((j, m_new))
}

/// One full pass of the per-hop ideal transfer over a precomputed (dim-2)
/// Expand walk. Returns `None` if any hop fails (caller retries).
fn rand_sk_from_walk<L: HasSigningPrecomp + LevelPrecomp>(
    sk: &SecretKey<L>,
    walk: &super::expand::ExpandOutput<L>,
    precomp: &crate::id2iso::sign_precomp::SigningPrecomp<L>,
    rng: &mut impl rand::Rng,
) -> Option<SecretKey<L>> {
    let mut i_ideal = sk.secret_ideal.clone();
    let mut m = sk.mat_ba_can_to_ba0_two.clone();

    // Overwritten on the first hop (RRAND >= 1); initialized for the type.
    let mut final_curve = sk.curve.clone();
    let mut final_basis = sk.canonical_basis.clone();
    let mut final_ideal = i_ideal.clone();

    for hop in 0..RRAND {
        let e_i = &walk.curves[hop];
        let basis_i = &walk.bases[hop].0;
        let (j, m_new) =
            ideal_transfer_step(&i_ideal, &m, &walk.scalars[hop], e_i, basis_i, precomp, rng)?;

        m = m_new;
        i_ideal = j.clone();
        final_curve = e_i.clone();
        final_basis = basis_i.clone();
        final_ideal = j;
    }

    Some(SecretKey {
        curve: final_curve,
        secret_ideal: final_ideal,
        mat_ba_can_to_ba0_two: m,
        canonical_basis: final_basis,
    })
}
