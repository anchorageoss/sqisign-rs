//! The two-point verifier against the four-point verifier at level I
//! (SECENG-2422): the same verdicts on the shipped vectors and on their
//! tampered variants, the same codomain null points at every step of both
//! half-chains, the gluing images and the fall-back flags, and the
//! fall-back branch on perturbed gluing images. Feature `compact`.
//!
//! The vectors are those of `compact_vectors.rs` (five keys, one signature
//! each, the seeds `"compact level I vector <i>"`), rebuilt here from the
//! signer; the tampered variants are those of `compact_roundtrip.rs` (the
//! wrong message, one bit of `A_com`, of `q` and of `a`, a truncated
//! signature, the wrong key) and 45 single-bit flips over the scalar and
//! degree bytes.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::level1;
use sqisign_verify::hd::{
    challenge_limbs, gluing_images_paired, hd_codomains, hd_verify_bytes, hd_verify_bytes_with,
    parse_public_key, parse_signature, run_strategy_chain, run_strategy_chain_two_point_checked,
    HdReject, HdSignature, HdVerifyOptions, KernelPoints, ThetaStructureDim4,
};
use sqisign_verify::P324_3;

const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
/// The number of shipped level-I vectors (`compact_vectors.rs`).
const COUNT: usize = 5;
const PK_BYTES: usize = 84;
const SIG_BYTES: usize = 142;

const FOUR: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::Four,
};
const TWO_CHECK: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::TwoWithImageCheck,
};
const TWO_NO_CHECK: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::TwoWithoutImageCheck,
};

/// A shipped vector: the public key, the signature and the message.
struct Vector {
    pk: [u8; PK_BYTES],
    sig: [u8; SIG_BYTES],
    msg: Vec<u8>,
    /// A second key of the same seed stream, for the wrong-key negative.
    other_pk: [u8; PK_BYTES],
}

/// The five vectors of `compact_vectors.rs`, rebuilt with the same seeds.
fn vectors() -> Vec<Vector> {
    let params = level1();
    let mut out = Vec::with_capacity(COUNT);
    for i in 0..COUNT {
        let label = format!("compact level I vector {i}");
        let mut rng = ShakeRng::new(label.as_bytes(), b"kat");
        let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen");
        let mut pkb = [0u8; PK_BYTES];
        compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
        let msg = format!("message {i}").into_bytes();
        let mut sig = [0u8; SIG_BYTES];
        compact_sign::<P324_3, N>(&params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
        let (pk2, _) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen 2");
        let mut other_pk = [0u8; PK_BYTES];
        compact_public_key_to_bytes(&pk2, &mut other_pk).unwrap();
        out.push(Vector {
            pk: pkb,
            sig,
            msg,
            other_pk,
        });
    }
    out
}

/// The three verdicts of one input: Four, then the two two-point modes.
fn verdicts(sig: &[u8], pk: &[u8], msg: &[u8]) -> [Result<(), HdReject>; 3] {
    let four = hd_verify_bytes::<P324_3>(sig, pk, msg);
    let (two_check, _) = hd_verify_bytes_with::<P324_3>(sig, pk, msg, &TWO_CHECK);
    let (two_no_check, _) = hd_verify_bytes_with::<P324_3>(sig, pk, msg, &TWO_NO_CHECK);
    [four, two_check, two_no_check]
}

/// A tampered variant: its label, the signature bytes, the public key and
/// the message.
type Tampered = (String, Vec<u8>, [u8; PK_BYTES], Vec<u8>);

/// The tampered variants of `compact_roundtrip.rs` and the 45 single-bit
/// flips over the five vectors, each with a label.
fn tampered(v: &Vector) -> Vec<Tampered> {
    let mut out = Vec::new();
    out.push((
        "wrong message".to_string(),
        v.sig.to_vec(),
        v.pk,
        b"another message".to_vec(),
    ));
    for (name, pos) in [("A_com bit", 0usize), ("q bit", 82), ("a bit", 82 + 22)] {
        let mut t = v.sig;
        t[pos] ^= 1;
        out.push((name.to_string(), t.to_vec(), v.pk, v.msg.clone()));
    }
    out.push((
        "truncated".to_string(),
        v.sig[..SIG_BYTES - 1].to_vec(),
        v.pk,
        v.msg.clone(),
    ));
    out.push((
        "wrong key".to_string(),
        v.sig.to_vec(),
        v.other_pk,
        v.msg.clone(),
    ));
    for pos in (SIG_BYTES - 60..SIG_BYTES).step_by(7) {
        let mut t = v.sig;
        t[pos] ^= 1;
        out.push((format!("flip byte {pos}"), t.to_vec(), v.pk, v.msg.clone()));
    }
    out
}

/// (1) The same verdict, and the same reason on a rejection, in the three
/// modes, on every shipped vector and every tampered variant.
#[test]
fn verdicts_agree_on_vectors_and_tampered_variants() {
    let vs = vectors();
    let mut accepted = 0usize;
    let mut rejected = 0usize;
    for (i, v) in vs.iter().enumerate() {
        let r = verdicts(&v.sig, &v.pk, &v.msg);
        assert_eq!(r, [Ok(()), Ok(()), Ok(())], "vector {i}");
        accepted += 1;
        for (label, sig, pk, msg) in tampered(v) {
            let r = verdicts(&sig, &pk, &msg);
            assert!(r[0].is_err(), "vector {i}, {label}: the four-point verdict");
            assert_eq!(r[1], r[0], "vector {i}, {label}: two-point with check");
            assert_eq!(r[2], r[0], "vector {i}, {label}: two-point without check");
            rejected += 1;
        }
    }
    eprintln!("verdicts: {accepted} vectors accepted in the three modes, {rejected} tampered variants rejected with the same reason in the three modes");
    assert_eq!(accepted, COUNT);
    assert_eq!(rejected, COUNT * (6 + 9));
}

/// The parsed signature of a vector as the `hd` layer's input.
fn with_hd_signature<T>(v: &Vector, f: impl FnOnce(&HdSignature<P324_3>) -> T) -> T {
    let s = parse_signature::<P324_3>(&v.sig).unwrap();
    let p = parse_public_key::<P324_3>(&v.pk).unwrap();
    let chal = challenge_limbs::<P324_3>(&s, &p, &v.msg).unwrap();
    let hdsig = HdSignature {
        a_pk: p.a_pk.clone(),
        a_com: s.a_com.clone(),
        hint_pk_p: p.hint_pk_p,
        hint_pk_q: p.hint_pk_q,
        hint_com_p: s.hint_com_p,
        hint_com_q: s.hint_com_q,
        chal_limbs: &chal,
        resp_a: s.a,
        resp_b: s.b,
        resp_c_or_d: s.c_or_d,
        q: s.q,
    };
    f(&hdsig)
}

/// (2) At every step of both half-chains the two-point modes give the
/// codomain null point of the four-point mode, projectively; the gluing
/// images are the same four points; no half falls back.
#[test]
fn codomains_agree_at_every_step() {
    let mut steps = 0usize;
    for (i, v) in vectors().iter().enumerate() {
        with_hd_signature(v, |hdsig| {
            let c4 = hd_codomains::<P324_3>(hdsig, &FOUR).expect("four-point chains");
            assert_eq!(c4.fell_back, [false, false]);
            for (name, opts) in [("with check", &TWO_CHECK), ("without check", &TWO_NO_CHECK)] {
                let c2 = hd_codomains::<P324_3>(hdsig, opts).expect("two-point chains");
                assert_eq!(c2.fell_back, [false, false], "vector {i} {name}");
                for h in 0..2 {
                    assert!(
                        c4.gluing[h].projective_eq(&c2.gluing[h]),
                        "vector {i} {name} half {h}: gluing codomain"
                    );
                    for (j, (a, b)) in c4.gluing_images[h]
                        .iter()
                        .zip(c2.gluing_images[h].iter())
                        .enumerate()
                    {
                        assert!(a.projective_eq(b), "vector {i} {name} half {h}: image {j}");
                    }
                    assert!(
                        gluing_images_paired(&c2.gluing_images[h]),
                        "vector {i} {name} half {h}: the images pair"
                    );
                    assert_eq!(
                        c4.plain[h].len(),
                        c2.plain[h].len(),
                        "vector {i} {name} half {h}: step count"
                    );
                    for (k, (a, b)) in c4.plain[h].iter().zip(c2.plain[h].iter()).enumerate() {
                        assert!(a.projective_eq(b), "vector {i} {name} half {h}: step {k}");
                    }
                    steps += c2.plain[h].len();
                }
            }
            eprintln!(
                "vector {i}: {} and {} plain steps per half compared in both modes",
                c4.plain[0].len(),
                c4.plain[1].len()
            );
        });
    }
    eprintln!(
        "codomains: {steps} plain-step codomains compared, two modes, two halves, {COUNT} vectors"
    );
}

/// The fall-back branch: with the gluing images of `g2` and `g3` replaced
/// so that they do not pair with those of `g0` and `g1`, the checked
/// two-point run sets its flag and gives the codomains of the four-point
/// chain on the same images; on the honest images it does not fall back.
#[test]
fn fall_back_on_unpaired_gluing_images() {
    let v = &vectors()[0];
    with_hd_signature(v, |hdsig| {
        let c4 = hd_codomains::<P324_3>(hdsig, &FOUR).expect("four-point chains");
        for h in 0..2 {
            let post = c4.gluing_images[h].clone();
            assert!(
                gluing_images_paired(&post),
                "half {h}: the honest images pair"
            );
            let start = ThetaStructureDim4::new(c4.gluing[h].clone());
            let l = c4.plain[h].len();

            let (sc, fell_back) =
                run_strategy_chain_two_point_checked(&start, &post, l).expect("two-point chain");
            assert!(!fell_back, "half {h}: no fall-back on the honest images");
            assert_eq!(sc.codomains.len(), l);
            for (k, (a, b)) in c4.plain[h].iter().zip(sc.codomains.iter()).enumerate() {
                assert!(a.projective_eq(b), "half {h} step {k}: two-point codomain");
            }

            let mut unpaired = post.clone();
            unpaired[2] = post[3].clone();
            assert!(!gluing_images_paired(&unpaired));
            let (sc_fb, fell_back) = run_strategy_chain_two_point_checked(&start, &unpaired, l)
                .expect("four-point chain on the unpaired images");
            assert!(fell_back, "half {h}: the flag is set");
            let sc_four = run_strategy_chain(&start, &unpaired, l).expect("four-point chain");
            assert_eq!(sc_fb.codomains.len(), sc_four.codomains.len());
            for (k, (a, b)) in sc_four
                .codomains
                .iter()
                .zip(sc_fb.codomains.iter())
                .enumerate()
            {
                assert!(a.projective_eq(b), "half {h} step {k}: fall-back codomain");
            }
            eprintln!(
                "half {h}: fall-back ran {} four-point steps",
                sc_fb.codomains.len()
            );
        }
    });
}

/// The flags of a rejection after the half-chains are those of the halves:
/// on every tampered variant rejected at the middle-codomain or image
/// stage, `hd_verify_with`'s report equals `hd_codomains`'s.
#[test]
fn fall_back_flags_through_rejections() {
    let mut after_halves = 0usize;
    for (i, v) in vectors().iter().enumerate() {
        for (label, sig, pk, msg) in tampered(v) {
            let (verdict, report) = hd_verify_bytes_with::<P324_3>(&sig, &pk, &msg, &TWO_CHECK);
            if !matches!(
                verdict,
                Err(HdReject::MiddleCodomainMismatch) | Err(HdReject::HdImageMismatch)
            ) {
                continue;
            }
            let s = parse_signature::<P324_3>(&sig).unwrap();
            let p = parse_public_key::<P324_3>(&pk).unwrap();
            let chal = challenge_limbs::<P324_3>(&s, &p, &msg).unwrap();
            let hdsig = HdSignature {
                a_pk: p.a_pk.clone(),
                a_com: s.a_com.clone(),
                hint_pk_p: p.hint_pk_p,
                hint_pk_q: p.hint_pk_q,
                hint_com_p: s.hint_com_p,
                hint_com_q: s.hint_com_q,
                chal_limbs: &chal,
                resp_a: s.a,
                resp_b: s.b,
                resp_c_or_d: s.c_or_d,
                q: s.q,
            };
            let diag = hd_codomains::<P324_3>(&hdsig, &TWO_CHECK).expect("half-chains ran");
            assert_eq!(diag.fell_back, report.fell_back, "vector {i}, {label}");
            after_halves += 1;
        }
    }
    eprintln!("flags: {after_halves} tampered variants rejected after the half-chains, flags equal to the diagnostic's");
    assert!(after_halves > 0);
}
