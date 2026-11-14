//! The two-point verifier against the four-point verifier at levels I, III
//! and V: the same verdicts on five vectors per level and on their tampered
//! variants, the same codomain null point at every plain step of both
//! half-chains, no fall-back on an honest signature, and the fall-back
//! branch on perturbed gluing inputs. Feature `compact`.
//!
//! The vectors are five keys with one signature each, from the seeds
//! `"compact level <I|III|V> vector <i>"` (those of `compact_vectors.rs` at
//! level I). The tampered variants of a vector are: the wrong message; one
//! bit of `A_com`, of `q` and of `a`; a truncated signature; the wrong key;
//! and the low bit of nine bytes, every seventh over the last sixty. A flip
//! that lands on the top bit of a transmitted scalar is set apart: such a
//! signature can still verify (the action changes by a map into the
//! 2-torsion), so it is compared across the modes, not required to fail.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::{
    challenge_limbs, gluing_images_paired, hd_codomains, hd_kernel_setup, hd_verify_bytes,
    hd_verify_bytes_with, null_fixed_by_pi, parse_public_key, parse_signature, run_strategy_chain,
    run_strategy_chain_two_point_checked, HdLevel, HdReject, HdSignature, HdVerifyOptions,
    KernelPoints, ThetaStructureDim4,
};
use sqisign_verify::{P324_3, P500_27, P664_17};

const N1: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const N3: usize = sqisign_rs::precomp::p500_27::IBZ_NLIMBS;
const N5: usize = sqisign_rs::precomp::p664_17::IBZ_NLIMBS;
/// Vectors per level.
const COUNT: usize = 5;

const FOUR: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::Four,
};
const TWO_CHECK: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::TwoWithImageCheck,
};
const TWO_NO_CHECK: HdVerifyOptions = HdVerifyOptions {
    kernel_points: KernelPoints::TwoWithoutImageCheck,
};

/// A vector: the public key, the signature and the message.
struct Vector {
    pk: Vec<u8>,
    sig: Vec<u8>,
    msg: Vec<u8>,
    /// A second key of the same seed stream, for the wrong-key negative.
    other_pk: Vec<u8>,
}

/// `count` vectors from the seeds `"<prefix> <i>"`.
fn vectors<L: HdLevel, const N: usize>(
    params: &Params<N>,
    prefix: &str,
    count: usize,
) -> Vec<Vector> {
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let label = format!("{prefix} {i}");
        let mut rng = ShakeRng::new(label.as_bytes(), b"kat");
        let (pk, sk) = compact_keygen::<L, N>(params, &mut rng).expect("keygen");
        let mut pkb = vec![0u8; compact_public_key_bytes::<L>()];
        compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
        let msg = format!("message {i}").into_bytes();
        let mut sig = vec![0u8; compact_signature_bytes::<L>()];
        compact_sign::<L, N>(params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
        let (pk2, _) = compact_keygen::<L, N>(params, &mut rng).expect("keygen 2");
        let mut other_pk = vec![0u8; compact_public_key_bytes::<L>()];
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

/// The three verdicts of one input (Four, then the two two-point modes) and
/// the fall-back flags of the checked two-point run.
#[allow(clippy::type_complexity)]
fn verdicts<L: HdLevel>(
    sig: &[u8],
    pk: &[u8],
    msg: &[u8],
) -> ([Result<(), HdReject>; 3], [bool; 2]) {
    let four = hd_verify_bytes::<L>(sig, pk, msg);
    let (two_check, report) = hd_verify_bytes_with::<L>(sig, pk, msg, &TWO_CHECK);
    let (two_no_check, _) = hd_verify_bytes_with::<L>(sig, pk, msg, &TWO_NO_CHECK);
    ([four, two_check, two_no_check], report.fell_back)
}

/// A variant: its label, the signature bytes, the public key and the
/// message.
type Variant = (String, Vec<u8>, Vec<u8>, Vec<u8>);

/// The bit offsets of the fields of a signature: `q`, `a`, `b`, `c_or_d`.
fn field_bits<L: HdLevel>() -> [usize; 4] {
    let q = 2 * L::P_BITS as usize;
    let a = q + L::E_EMBED as usize - 2;
    let r = L::R as usize;
    [q, a, a + r, a + 2 * r]
}

fn flip(sig: &[u8], bit: usize) -> Vec<u8> {
    let mut t = sig.to_vec();
    t[bit / 8] ^= 1 << (bit % 8);
    t
}

/// The tampered variants of a vector, and the flips that land on the top
/// bit of a transmitted scalar.
fn tampered<L: HdLevel>(v: &Vector) -> (Vec<Variant>, Vec<Variant>) {
    let [q0, a0, b0, cd0] = field_bits::<L>();
    let r = L::R as usize;
    let top_bits = [a0 + r - 1, b0 + r - 1, cd0 + r - 1];
    let n = v.sig.len();
    let mut out = vec![(
        "wrong message".to_string(),
        v.sig.clone(),
        v.pk.clone(),
        b"another message".to_vec(),
    )];
    for (name, bit) in [("A_com bit", 0usize), ("q bit", q0), ("a bit", a0)] {
        out.push((
            name.to_string(),
            flip(&v.sig, bit),
            v.pk.clone(),
            v.msg.clone(),
        ));
    }
    out.push((
        "truncated".to_string(),
        v.sig[..n - 1].to_vec(),
        v.pk.clone(),
        v.msg.clone(),
    ));
    out.push((
        "wrong key".to_string(),
        v.sig.clone(),
        v.other_pk.clone(),
        v.msg.clone(),
    ));
    let mut top = Vec::new();
    for pos in (n - 60..n).step_by(7) {
        let bit = 8 * pos;
        let variant = (
            format!("flip bit {bit}"),
            flip(&v.sig, bit),
            v.pk.clone(),
            v.msg.clone(),
        );
        if top_bits.contains(&bit) {
            top.push(variant);
        } else {
            out.push(variant);
        }
    }
    (out, top)
}

/// (1) The same verdict, and the same reason on a rejection, in the three
/// modes, on every vector and every tampered variant; no fall-back on an
/// honest signature. Returns the counts (accepted, rejected, top-bit flips,
/// of which accepted).
fn check_verdicts<L: HdLevel>(vs: &[Vector]) -> (usize, usize, usize, usize) {
    let (mut accepted, mut rejected, mut top_total, mut top_valid) = (0, 0, 0, 0);
    for (i, v) in vs.iter().enumerate() {
        let (r, fell_back) = verdicts::<L>(&v.sig, &v.pk, &v.msg);
        assert_eq!(r, [Ok(()), Ok(()), Ok(())], "vector {i}");
        assert_eq!(fell_back, [false, false], "vector {i}: no fall-back");
        accepted += 1;
        let (variants, top) = tampered::<L>(v);
        for (label, sig, pk, msg) in variants {
            let (r, _) = verdicts::<L>(&sig, &pk, &msg);
            assert!(r[0].is_err(), "vector {i}, {label}: the four-point verdict");
            assert_eq!(r[1], r[0], "vector {i}, {label}: two-point with check");
            assert_eq!(r[2], r[0], "vector {i}, {label}: two-point without check");
            rejected += 1;
        }
        for (label, sig, pk, msg) in top {
            let (r, _) = verdicts::<L>(&sig, &pk, &msg);
            assert_eq!(r[1], r[0], "vector {i}, {label}: two-point with check");
            assert_eq!(r[2], r[0], "vector {i}, {label}: two-point without check");
            top_total += 1;
            top_valid += r[0].is_ok() as usize;
        }
    }
    (accepted, rejected, top_total, top_valid)
}

/// The parsed signature of a vector as the `hd` layer's input.
fn with_hd_signature<L: HdLevel, T>(
    sig: &[u8],
    pk: &[u8],
    msg: &[u8],
    f: impl FnOnce(&HdSignature<L>) -> T,
) -> T {
    let s = parse_signature::<L>(sig).unwrap();
    let p = parse_public_key::<L>(pk).unwrap();
    let chal = challenge_limbs::<L>(&s, &p, msg).unwrap();
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

/// (2) At every plain step of both half-chains the two-point modes give
/// the codomain null point of the four-point mode, projectively; the gluing
/// images are the same four points, paired by the permutation, and the
/// gluing codomain's null point is fixed by it; no half falls back. Returns
/// the plain steps of the four-point chains, both halves, all vectors.
fn check_codomains<L: HdLevel>(vs: &[Vector]) -> usize {
    let mut steps = 0usize;
    for (i, v) in vs.iter().enumerate() {
        with_hd_signature::<L, _>(&v.sig, &v.pk, &v.msg, |hdsig| {
            let m = hd_kernel_setup::<L>(hdsig).expect("setup").m;
            let c4 = hd_codomains::<L>(hdsig, &FOUR).expect("four-point chains");
            assert_eq!(c4.fell_back, [false, false]);
            // each half-chain as computed has ⌈e/2⌉ steps, of which m in
            // dimension 2 and one gluing
            let e1 = L::E1 as usize;
            assert_eq!(c4.plain[0].len(), e1 - m - 1, "vector {i}: half 1");
            assert_eq!(c4.plain[1].len(), e1 - m - 1, "vector {i}: half 2");
            for (name, opts) in [("with check", &TWO_CHECK), ("without check", &TWO_NO_CHECK)] {
                let c2 = hd_codomains::<L>(hdsig, opts).expect("two-point chains");
                assert_eq!(c2.fell_back, [false, false], "vector {i} {name}");
                for h in 0..2 {
                    assert!(
                        c4.gluing[h].projective_eq(&c2.gluing[h]),
                        "vector {i} {name} half {h}: gluing codomain"
                    );
                    assert!(
                        null_fixed_by_pi(&c2.gluing[h]),
                        "vector {i} {name} half {h}: the gluing null point is fixed"
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
                    assert_eq!(c4.plain[h].len(), c2.plain[h].len());
                    for (k, (a, b)) in c4.plain[h].iter().zip(c2.plain[h].iter()).enumerate() {
                        assert!(a.projective_eq(b), "vector {i} {name} half {h}: step {k}");
                    }
                }
            }
            steps += c4.plain[0].len() + c4.plain[1].len();
        });
    }
    steps
}

/// The conformance run of one level: five vectors, their tampered variants
/// and every plain step. `variants` is the number of tampered variants the
/// construction gives at this level.
fn conformance<L: HdLevel, const N: usize>(params: &Params<N>, roman: &str, variants: usize) {
    let vs = vectors::<L, N>(params, &format!("compact level {roman} vector"), COUNT);
    let (accepted, rejected, top_total, top_valid) = check_verdicts::<L>(&vs);
    let steps = check_codomains::<L>(&vs);
    eprintln!(
        "level {roman}: {accepted} vectors accepted in the three modes with no fall-back; \
         {rejected} tampered variants rejected with the same reason in the three modes; \
         {top_total} flips of a scalar's top bit, {top_valid} of them valid, the same in the \
         three modes; {steps} plain steps with the same codomain in both two-point modes"
    );
    assert_eq!(accepted, COUNT);
    assert_eq!(rejected, variants);
    assert_eq!(rejected + top_total, COUNT * 15);
}

#[test]
fn conformance_level1() {
    conformance::<P324_3, N1>(&level1(), "I", 70);
}

#[test]
fn conformance_level3() {
    conformance::<P500_27, N3>(&level3(), "III", 75);
}

#[test]
fn conformance_level5() {
    conformance::<P664_17, N5>(&level5(), "V", 70);
}

/// Twenty-four signatures at level I: every verdict and every codomain of
/// the four-point verifier, on the honest signatures and on their tampered
/// variants, and no fall-back.
#[test]
fn twenty_four_signatures_level1() {
    let vs = vectors::<P324_3, N1>(&level1(), "compact level I signature", 24);
    let (accepted, rejected, top_total, top_valid) = check_verdicts::<P324_3>(&vs);
    let steps = check_codomains::<P324_3>(&vs);
    eprintln!(
        "level I, 24 signatures: {accepted} accepted, {rejected} tampered variants rejected, \
         {top_total} top-bit flips ({top_valid} valid), {steps} plain steps"
    );
    assert_eq!(accepted, 24);
    assert_eq!(rejected, 24 * 14);
}

/// The fall-back branch: with the gluing images of `g2` and `g3` replaced
/// so that they do not pair with those of `g0` and `g1`, the checked
/// two-point run sets its flag and gives the codomains of the four-point
/// chain on the same inputs; on the honest inputs it does not fall back.
fn fall_back<L: HdLevel, const N: usize>(params: &Params<N>, roman: &str) {
    let v = &vectors::<L, N>(params, &format!("compact level {roman} vector"), 1)[0];
    with_hd_signature::<L, _>(&v.sig, &v.pk, &v.msg, |hdsig| {
        let c4 = hd_codomains::<L>(hdsig, &FOUR).expect("four-point chains");
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
        }
    });
}

#[test]
fn fall_back_on_unpaired_inputs() {
    fall_back::<P324_3, N1>(&level1(), "I");
    fall_back::<P500_27, N3>(&level3(), "III");
    fall_back::<P664_17, N5>(&level5(), "V");
}

/// The flags of a rejection after the half-chains are those of the halves:
/// on every tampered variant rejected at the middle-codomain or image
/// stage, `hd_verify_with`'s report equals `hd_codomains`'s.
#[test]
fn fall_back_flags_through_rejections() {
    type L = P324_3;
    let mut after_halves = 0usize;
    let vs = vectors::<L, N1>(&level1(), "compact level I vector", COUNT);
    for (i, v) in vs.iter().enumerate() {
        for (label, sig, pk, msg) in tampered::<L>(v).0 {
            let (verdict, report) = hd_verify_bytes_with::<L>(&sig, &pk, &msg, &TWO_CHECK);
            if !matches!(
                verdict,
                Err(HdReject::MiddleCodomainMismatch) | Err(HdReject::HdImageMismatch)
            ) {
                continue;
            }
            let diag = with_hd_signature::<L, _>(&sig, &pk, &msg, |hdsig| {
                hd_codomains::<L>(hdsig, &TWO_CHECK).expect("half-chains ran")
            });
            assert_eq!(diag.fell_back, report.fell_back, "vector {i}, {label}");
            after_halves += 1;
        }
    }
    eprintln!("flags: {after_halves} tampered variants rejected after the half-chains, flags equal to the diagnostic's");
    assert!(after_halves > 0);
}
