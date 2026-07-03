//! SQIsign-RK: signatures with randomizable keys (ePrint 2026/1169).
//!
//! Wraps the base SQIsign implementation with key-randomization operations.
//! The base `KeyGen`/`Sign`/`Verify` are unchanged.
//!
//! The public API is exactly three functions, over the same `PublicKey` /
//! `SecretKey` types that keygen/sign/verify already use:
//!
//! - [`rand_pk`]; public key randomization. Deterministic, needs no secret;
//!   walks a chain of `2^f`-isogenies from `E_pk` to a new curve `E_pk'`.
//! - [`rand_sk`]; secret key randomization, transferring the endomorphism ring
//!   through that same walk (the Deuring side).
//! - [`ver_key`]; check a secret key is valid for a public key.
//!
//! Gated behind the `sqisign-rk` feature so it does not affect the default
//! build. See `SQISIGN_RK.md` for the design and the mapping to the SageMath
//! reference.

pub mod compact;
pub mod expand;
pub mod iso_to_end;
mod rng;

pub use compact::{rand_pk_compact, rand_sk_compact, ver_key_compact};
pub use expand::{expand, ExpandOutput, RRAND};
pub use iso_to_end::rand_sk;

use crate::id2iso::sign_precomp::HasSigningPrecomp;
use crate::id2iso::sign_side::dim2id2iso_arbitrary_isogeny_evaluation;
use crate::keygen::SecretKey;
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::types::PublicKey;

/// Derive a new public key from an existing public key and arbitrary
/// randomness.
///
/// This is a deterministic, public operation. No secret key is required. The
/// same `(pk, rr)` always produces the same output, and the derived key is
/// unlinkable to the input key.
///
/// ```no_run
/// # #[cfg(feature = "sqisign-rk")] {
/// use sqisign_rs::keygen::keypair;
/// use sqisign_rs::params::Level1;
/// use sqisign_rs::sqisign_rk::rand_pk;
///
/// let mut rng = rand::thread_rng();
/// let (pk, _sk) = keypair::<Level1>(&mut rng);
/// let pk_child = rand_pk(&pk, b"child-0");
/// # }
/// ```
pub fn rand_pk<L: HasSigningPrecomp + LevelPrecomp>(pk: &PublicKey<L>, rr: &[u8]) -> PublicKey<L> {
    let precomp = L::signing_precomp();
    let out = expand(pk.curve(), rr, &precomp)
        .expect("invariant: Expand walk must succeed on a valid public key");
    // The randomized public key is the final codomain plus its canonical-basis
    // recomputation hint.
    let curve = out.curves[RRAND - 1].clone();
    let hint = out.bases[RRAND - 1].1;
    PublicKey::new(curve, hint)
}

/// Verify that a secret key is valid for a given public key.
///
/// Recomputes the isogeny `E0 → E'` from the secret key's ideal and checks that
/// its codomain matches `pk`'s curve (by `j`-invariant, so it is independent of
/// the Montgomery model).
///
/// ```no_run
/// # #[cfg(feature = "sqisign-rk")] {
/// use sqisign_rs::keygen::keypair;
/// use sqisign_rs::params::Level1;
/// use sqisign_rs::sqisign_rk::ver_key;
///
/// let mut rng = rand::thread_rng();
/// let (pk, sk) = keypair::<Level1>(&mut rng);
/// assert!(ver_key(&pk, &sk));
/// # }
/// ```
pub fn ver_key<L: HasSigningPrecomp + LevelPrecomp>(pk: &PublicKey<L>, sk: &SecretKey<L>) -> bool {
    let precomp = L::signing_precomp();
    // Deterministic internal randomness (the codomain is invariant to it).
    let a_bytes = sk.curve.a.encode();
    let ctx: [&[u8]; 1] = [a_bytes.as_ref()];
    let mut rng = rng::seed_rng(b"SQIsign-RK/VerKey", &ctx);

    let mut basis = EcBasis::new(
        EcPoint::identity(),
        EcPoint::identity(),
        EcPoint::identity(),
    );
    let mut codomain = EcCurve::<L>::default();
    if dim2id2iso_arbitrary_isogeny_evaluation(
        &mut basis,
        &mut codomain,
        &sk.secret_ideal,
        &precomp,
        &mut rng,
    )
    .is_none()
    {
        return false;
    }
    codomain.j_inv().encode() == pk.curve().j_inv().encode()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Verifier;
    use sqisign_verify::params::{Level1, Level3, Level5};

    /// A fresh keypair's public key for the given level.
    fn keygen_pk<L: HasSigningPrecomp + LevelPrecomp>() -> PublicKey<L> {
        let mut rng = rand::rngs::StdRng::from_seed([7u8; 32]);
        let (pk, _sk) = crate::generate::<L>(&mut rng);
        pk
    }

    fn run_level<L: HasSigningPrecomp + LevelPrecomp>() {
        let pk = keygen_pk::<L>();
        let rr_a = [0x11u8; 32];
        let rr_b = [0x22u8; 32];

        // Determinism: same (pk, rr) -> same pk'.
        let pk1 = rand_pk(&pk, &rr_a);
        let pk2 = rand_pk(&pk, &rr_a);
        assert_eq!(
            pk1.curve().a.encode(),
            pk2.curve().a.encode(),
            "RandPK must be deterministic in (pk, rr)"
        );

        // Diversity: different rr -> different pk'.
        let pk3 = rand_pk(&pk, &rr_b);
        assert_ne!(
            pk1.curve().a.encode(),
            pk3.curve().a.encode(),
            "different rr must give a different pk'"
        );

        // The randomized curve must differ from the original (a non-trivial walk).
        assert_ne!(
            pk1.curve().a.encode(),
            pk.curve().a.encode(),
            "randomization must move off the original curve"
        );

        // Validity: the output must be a supersingular curve. We check the
        // strongest cheap invariant available here; that a canonical 2^f basis
        // exists on it (the deterministic-basis routine only succeeds on a
        // valid Montgomery curve of the right order), by round-tripping through
        // rand_pk once more from pk1.
        let pk1_again = rand_pk(&pk1, &rr_a);
        assert_ne!(
            pk1_again.curve().a.encode(),
            pk1.curve().a.encode(),
            "sequential randomization pk' -> pk'' must move again"
        );
    }

    /// Full RandSK round-trip: derive both keys, check the derived key is valid,
    /// sign with the derived sk and verify with the derived pk, and confirm the
    /// derived signature does not verify under the original pk.
    fn run_roundtrip<L: HasSigningPrecomp + LevelPrecomp>() {
        let mut rng = rand::thread_rng();
        let (pk, sk) = crate::keygen::keypair::<L>(&mut rng);
        let rr = [0x5Au8; 32];

        let pk_prime = rand_pk(&pk, &rr);
        let sk_prime = rand_sk(&sk, &pk, &rr);

        // pk' and sk' must agree: sk'.curve is the Expand curve = pk'.curve.
        assert_eq!(
            sk_prime.curve.a.encode(),
            pk_prime.curve().a.encode(),
            "derived sk and pk must share the same curve model"
        );

        // The derived key is valid.
        assert!(ver_key(&pk_prime, &sk_prime), "ver_key(pk', sk') must hold");

        // Sign with the derived sk, verify with the derived pk.
        let msg = b"SQIsign-RK round-trip";
        let sig = crate::sign::sign::<L>(&sk_prime, &pk_prime, msg, &mut rng)
            .expect("sign with derived sk");
        assert!(
            pk_prime.verify(msg, &sig).is_ok(),
            "derived signature must verify under derived pk"
        );

        // Cross-key rejection: derived signature must not verify under the
        // original pk.
        assert!(
            pk.verify(msg, &sig).is_err(),
            "derived signature must NOT verify under the original pk"
        );
    }

    #[test]
    fn rand_sk_roundtrip_level1() {
        run_roundtrip::<Level1>();
    }

    #[test]
    fn rand_sk_roundtrip_level3() {
        run_roundtrip::<Level3>();
    }

    #[test]
    fn rand_sk_roundtrip_level5() {
        run_roundtrip::<Level5>();
    }

    /// Timing comparison (run with `--ignored --nocapture`): RandPK and RandSK
    /// relative to KeyGen. The paper's Table 4 reports ≈0.4× (RandPK) and
    /// ≈2.5× (RandSK) vs KeyGen.
    fn time_level<L: HasSigningPrecomp + LevelPrecomp>(name: &str) {
        use std::time::Instant;
        let mut rng = rand::thread_rng();
        let rr = [0x5Au8; 32];

        let t = Instant::now();
        let (pk, sk) = crate::keygen::keypair::<L>(&mut rng);
        let keygen = t.elapsed().as_secs_f64();

        let t = Instant::now();
        let _ = rand_pk(&pk, &rr);
        let randpk = t.elapsed().as_secs_f64();

        let t = Instant::now();
        let _ = rand_sk(&sk, &pk, &rr);
        let randsk = t.elapsed().as_secs_f64();

        eprintln!(
            "[SQIsign-RK {name}] keygen={keygen:.3}s  rand_pk={randpk:.3}s ({:.2}x)  rand_sk={randsk:.3}s ({:.2}x)",
            randpk / keygen,
            randsk / keygen
        );
    }

    #[test]
    #[ignore]
    fn timings() {
        time_level::<Level1>("L1");
        time_level::<Level3>("L3");
        time_level::<Level5>("L5");
    }

    #[test]
    fn rand_pk_level1() {
        run_level::<Level1>();
    }

    #[test]
    fn rand_pk_level3() {
        run_level::<Level3>();
    }

    #[test]
    fn rand_pk_level5() {
        run_level::<Level5>();
    }

    use rand::SeedableRng;
}
