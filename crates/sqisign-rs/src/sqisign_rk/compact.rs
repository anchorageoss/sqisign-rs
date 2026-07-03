//! SQIsign-RK for the **compact** (dimension-4, 108-byte) scheme; Level 1 only
//! (the compact types are hardcoded to Level 1).
//!
//! Reuses the shared per-hop ideal transfer
//! ([`super::iso_to_end::ideal_transfer_step`]); the only difference from the
//! dim-2 path is the Expand walk's per-hop basis, which uses the HD convention
//! (`canonical_hints_l1` + `hd_torsion_basis_l1`, two `u32` hints) instead of
//! `ec_curve_to_basis_2f_to_hint` (one `u8`).

use super::expand::{h_path_scalars, RRAND};
use super::iso_to_end::ideal_transfer_step;
use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::id2iso::sign_side::dim2id2iso_arbitrary_isogeny_evaluation;
use crate::quaternion::intbig::{ibz_copy_digits, Ibz};
use crate::sign::compact::CompactSigningKey;
use crate::sign::dim4::{Dim4PublicKey, Dim4SecretKey};
use alloc::vec::Vec;
use sqisign_verify::ec::isogeny::ec_eval_even;
use sqisign_verify::ec::jacobian::jac_add;
use sqisign_verify::ec::point::ec_ladder3pt;
use sqisign_verify::ec::{EcBasis, EcCurve, EcIsogEven, EcPoint};
use sqisign_verify::fp::Fp2;
use sqisign_verify::hd::{canonical_hints_l1, encode_public_key, hd_torsion_basis_l1};
use sqisign_verify::{CompactPublicKey, Level1, SecurityLevel};

/// HD-convention Expand walk output (Level 1).
struct CompactExpand {
    /// Codomain of each hop; the last is the randomized curve.
    curves: Vec<EcCurve<Level1>>,
    /// HD canonical `2^f` basis on each codomain.
    bases: Vec<EcBasis<Level1>>,
    /// HD basis hints `(hint_pk_p, hint_pk_q)` for each codomain.
    hints: Vec<(u32, u32)>,
    /// Kernel scalars `rr_i`.
    scalars: Vec<Ibz>,
}

/// HD canonical `2^f`-torsion basis on a curve with (normalized) coefficient
/// `a`, plus its `(hint_pk_p, hint_pk_q)`.
fn hd_basis(curve: &EcCurve<Level1>, a: &Fp2<Level1>) -> Option<(EcBasis<Level1>, (u32, u32))> {
    let (hp, hq) = canonical_hints_l1(a)?;
    let (pj, qj) = hd_torsion_basis_l1(a, hp, hq)?;
    let pmq = jac_add(&pj, &qj.neg(), curve);
    let basis = EcBasis::new(pj.to_xz(), qj.to_xz(), pmq.to_xz());
    Some((basis, (hp, hq)))
}

/// The compact (HD) Expand walk from a public curve coefficient `a_pk`.
fn expand_compact(a_pk: &Fp2<Level1>, rr: &[u8]) -> Option<CompactExpand> {
    let f = <Level1 as SecurityLevel>::TORSION_EVEN_POWER;

    let mut curve = EcCurve::<Level1>::from_a(a_pk)?;
    curve.normalize_curve_and_a24();
    let a0 = curve.a.clone();
    let (mut basis, _) = hd_basis(&curve, &a0)?;

    let scalars = h_path_scalars(a0.encode().as_ref(), rr, f);

    let mut curves = Vec::with_capacity(RRAND);
    let mut bases = Vec::with_capacity(RRAND);
    let mut hints = Vec::with_capacity(RRAND);
    let mut scalar_ibz = Vec::with_capacity(RRAND);

    for scalar in scalars.iter() {
        scalar_ibz.push(ibz_copy_digits(&scalar.limbs[..scalar.nlimbs]));

        // Kernel K = P + [rr_i] Q, then the degree-2^f isogeny.
        let kernel = ec_ladder3pt(
            &scalar.limbs[..scalar.nlimbs],
            &basis.p,
            &basis.q,
            &basis.pmq,
            &curve,
        )?;
        let phi = EcIsogEven {
            curve: curve.clone(),
            kernel,
            length: f,
        };
        let mut codomain = curve.clone();
        ec_eval_even(&mut codomain, &phi, &mut [])?;
        codomain.normalize_curve_and_a24();

        let a_cod = codomain.a.clone();
        let (next_basis, hint) = hd_basis(&codomain, &a_cod)?;

        curves.push(codomain.clone());
        bases.push(next_basis.clone());
        hints.push(hint);
        curve = codomain;
        basis = next_basis;
    }

    Some(CompactExpand {
        curves,
        bases,
        hints,
        scalars: scalar_ibz,
    })
}

/// Build a `CompactPublicKey` from a curve coefficient and HD hints (via the
/// public 64-byte encoding, as `generate_compact` does).
fn compact_pk_from_parts(a: &Fp2<Level1>, hp: u32, hq: u32) -> Option<CompactPublicKey<Level1>> {
    let bytes = encode_public_key(a, hp, hq)?;
    CompactPublicKey::from_bytes(&bytes).ok()
}

/// Derive a new compact public key from an existing one and randomness `rr`.
///
/// Deterministic, public, no secret key. Compact / Level-1 analogue of
/// [`super::rand_pk`].
pub fn rand_pk_compact(pk: &CompactPublicKey<Level1>, rr: &[u8]) -> CompactPublicKey<Level1> {
    let a_pk = pk.a_pk();
    let walk = expand_compact(&a_pk, rr).expect("invariant: compact Expand walk must succeed");
    let a_final = walk.curves[RRAND - 1].a.clone();
    let (hp, hq) = walk.hints[RRAND - 1];
    compact_pk_from_parts(&a_final, hp, hq).expect("invariant: derived compact public key encodes")
}

/// Derive a new compact signing key corresponding to `rand_pk_compact(pk, rr)`.
///
/// Deterministic in `(sk, pk, rr)`. Compact / Level-1 analogue of
/// [`super::rand_sk`].
pub fn rand_sk_compact(
    sk: &CompactSigningKey<Level1>,
    _pk: &CompactPublicKey<Level1>,
    rr: &[u8],
) -> CompactSigningKey<Level1> {
    let precomp = Level1::signing_precomp();
    let a_pk = sk.dim4_public().a_pk.clone();
    let walk = expand_compact(&a_pk, rr).expect("invariant: compact Expand walk must succeed");

    let a_bytes = a_pk.encode();
    let ctx: [&[u8]; 2] = [a_bytes.as_ref(), rr];
    let mut rng = super::rng::seed_rng(b"SQIsign-RK/RandSK-compact", &ctx);

    for _ in 0..64 {
        if let Some(sk_new) = rand_sk_compact_once(sk, &walk, &precomp, &mut rng) {
            return sk_new;
        }
    }
    panic!("invariant: compact RandSK failed to converge after 64 attempts");
}

fn rand_sk_compact_once(
    sk: &CompactSigningKey<Level1>,
    walk: &CompactExpand,
    precomp: &crate::id2iso::sign_precomp::SigningPrecomp<Level1>,
    rng: &mut impl rand::Rng,
) -> Option<CompactSigningKey<Level1>> {
    let mut i_ideal = sk.dim4_secret().secret_ideal.clone();
    let mut m = sk.dim4_secret().mat_ba_can_to_ba0_two.clone();
    let mut final_ideal = i_ideal.clone();

    for hop in 0..RRAND {
        let (j, m_new) = ideal_transfer_step(
            &i_ideal,
            &m,
            &walk.scalars[hop],
            &walk.curves[hop],
            &walk.bases[hop],
            precomp,
            rng,
        )?;
        m = m_new;
        i_ideal = j.clone();
        final_ideal = j;
    }

    let a_final = walk.curves[RRAND - 1].a.clone();
    let (hp, hq) = walk.hints[RRAND - 1];
    let pk_dim4 = Dim4PublicKey {
        a_pk: a_final.clone(),
        hint_pk_p: hp,
        hint_pk_q: hq,
    };
    let sk_dim4 = Dim4SecretKey {
        secret_ideal: final_ideal,
        mat_ba_can_to_ba0_two: m,
    };
    let cpk = compact_pk_from_parts(&a_final, hp, hq)?;
    Some(CompactSigningKey::from_parts(sk_dim4, pk_dim4, cpk))
}

/// Verify that a compact secret key is valid for a compact public key.
///
/// Compact / Level-1 analogue of [`super::ver_key`].
pub fn ver_key_compact(pk: &CompactPublicKey<Level1>, sk: &CompactSigningKey<Level1>) -> bool {
    let precomp = Level1::signing_precomp();
    let a_bytes = sk.dim4_public().a_pk.encode();
    let ctx: [&[u8]; 1] = [a_bytes.as_ref()];
    let mut rng = super::rng::seed_rng(b"SQIsign-RK/VerKey-compact", &ctx);

    let mut basis = EcBasis::new(
        EcPoint::identity(),
        EcPoint::identity(),
        EcPoint::identity(),
    );
    let mut codomain = EcCurve::<Level1>::default();
    if dim2id2iso_arbitrary_isogeny_evaluation(
        &mut basis,
        &mut codomain,
        &sk.dim4_secret().secret_ideal,
        &precomp,
        &mut rng,
    )
    .is_none()
    {
        return false;
    }
    let Some(e_pk) = EcCurve::<Level1>::from_a(&pk.a_pk()) else {
        return false;
    };
    codomain.j_inv().encode() == e_pk.j_inv().encode()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sign::generate_compact;
    use crate::Verifier;

    #[test]
    fn rand_compact_roundtrip_level1() {
        let mut rng = rand::thread_rng();
        let (pk, sk) = generate_compact(&mut rng);
        let rr = b"test-randomness";

        let pk_prime = rand_pk_compact(&pk, rr);
        let sk_prime = rand_sk_compact(&sk, &pk, rr);

        // Determinism of public derivation.
        let pk_prime2 = rand_pk_compact(&pk, rr);
        assert_eq!(
            pk_prime.to_bytes(),
            pk_prime2.to_bytes(),
            "rand_pk_compact must be deterministic"
        );

        // The derived pk must differ from the original.
        assert_ne!(
            pk_prime.to_bytes(),
            pk.to_bytes(),
            "derived compact pk must differ from the original"
        );

        // The derived key is valid.
        assert!(
            ver_key_compact(&pk_prime, &sk_prime),
            "ver_key_compact(pk', sk') must hold"
        );

        // The derived signing key's own public key matches the derived pk.
        assert_eq!(
            sk_prime.public_key().to_bytes(),
            pk_prime.to_bytes(),
            "derived sk's public key must equal the derived pk"
        );

        // Sign with the derived sk, verify with the derived pk.
        let msg = b"SQIsign-RK compact round-trip";
        let sig = sk_prime
            .sign(msg, &mut rng)
            .expect("compact sign with derived sk");
        assert!(
            pk_prime.verify(msg, &sig).is_ok(),
            "derived compact signature must verify under derived pk"
        );

        // Cross-key rejection.
        assert!(
            pk.verify(msg, &sig).is_err(),
            "derived compact signature must NOT verify under the original pk"
        );
    }
}
