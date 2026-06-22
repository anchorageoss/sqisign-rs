//! Phase 5: the SQIsignHD SHAKE256 hash-to-challenge, the stage-3 response
//! recovery, and the end-to-end FastVerify orchestration.
//!
//! Ground truth is `test_vectors_l1.json` (curves, challenge, response,
//! recorded `chal`) plus `chain_vectors.json` (the plain-step kernels). The
//! signed message for every recorded signature is 32 zero bytes (the C test
//! `test_sqisignhd.c` uses `unsigned char msg[32] = {0}`), so reproducing the
//! recorded `chal` is exactly a `fips202.c` match for `sha3::Shake256`.

mod hd_common;
use hd_common::{le32, load, parse_coords, parse_fp2, Pt, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::hd::{
    hd_challenge_from_curves, hd_challenge_len, hd_verify, hd_verify_checked, recover_response_cd,
    HdReject, HdVerifyInputs, ThetaPointDim4,
};
use sqisign_verify::Level1;
use std::hint::black_box;
use std::time::Instant;

const CHAIN_VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/chain_vectors.json");

/// Level-1 response modulus exponent `r = ceil(f/2) + 2 = 70`.
const R_LVL1: u32 = 70;
/// The fixed signed message used by the reference test harness.
const MSG: [u8; 32] = [0u8; 32];

fn parse_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}
fn parse_u128(v: &Value) -> u128 {
    v.as_str().unwrap().parse::<u128>().unwrap()
}

/// `q` can exceed 128 bits (q < 2^f = 2^136 at Level 1), but only `q mod 2^r`
/// is needed; read its low 128 bits from the verbatim hex line (`raw_lines[4]`).
fn q_low128(sig: &Value) -> u128 {
    let line = sig["raw_lines"][4].as_str().unwrap();
    let hex = line.rsplit('=').next().unwrap().trim();
    let le = le32(hex);
    let mut b = [0u8; 16];
    b.copy_from_slice(&le[0..16]);
    u128::from_le_bytes(b)
}

/// Extract the recorded challenge as little-endian bytes from the signature's
/// verbatim data line (`raw_lines[7]` = "chal = <hex>").
fn chal_le(sig: &Value) -> [u8; 32] {
    let line = sig["raw_lines"][7].as_str().unwrap();
    let hex = line.rsplit('=').next().unwrap().trim();
    le32(hex)
}

fn parse_k8(v: &Value) -> [Pt; 4] {
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 4);
    core::array::from_fn(|k| parse_coords(&arr[k]))
}
fn parse_kernels(v: &Value) -> Vec<[Pt; 4]> {
    v.as_array().unwrap().iter().map(parse_k8).collect()
}

/// Deliverable 2: the SHAKE256 hash-to-challenge reproduces the recorded `chal`
/// for all 5 signatures - i.e. `sha3::Shake256` matches `fips202.c`.
#[test]
fn challenge_reproduces_recorded_chal() {
    let doc = load(PHASE0_VECTORS);
    let n = hd_challenge_len::<Level1>();
    assert_eq!(n, 32, "Level 1 challenge is 32 bytes");

    let mut count = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let a_com = parse_fp2(&v["signature"]["A_com"]);
        let a_pk = parse_fp2(&v["public_key"]["A_pk"]);

        let mut out = [0u8; 64];
        assert!(
            hd_challenge_from_curves(&a_com, &a_pk, &MSG, &mut out[..n]),
            "curves must be valid"
        );
        let expected = chal_le(&v["signature"]);
        assert_eq!(
            &out[..n],
            &expected[..n],
            "vector {}: recomputed challenge != recorded chal",
            v["index"]
        );
        count += 1;
    }
    assert_eq!(count, 5);
    println!("SHAKE256 hash-to-challenge reproduced recorded chal for all {count} vectors");
}

/// Deliverable 1 (stage 3): the response recovery reproduces the oracle's c,d.
#[test]
fn response_recovery_matches_oracle() {
    let doc = load(PHASE0_VECTORS);
    let mask = (1u128 << R_LVL1) - 1;
    for v in doc["test_vectors"].as_array().unwrap() {
        let sig = &v["signature"];
        let a = parse_i128(&sig["a"]);
        let b = parse_i128(&sig["b"]);
        let c_or_d = parse_i128(&sig["c_or_d"]);
        let q = q_low128(sig);
        let s3 = &v["stage3_image_response"];
        let k = parse_u128(&s3["k"]);
        // One of c,d equals the raw signed c_or_d (possibly negative); compare
        // mod 2^r.
        let c_exp = (parse_i128(&s3["c"]) as u128) & mask;
        let d_exp = (parse_i128(&s3["d"]) as u128) & mask;

        let (c, d) = recover_response_cd(a, b, c_or_d, q, k, R_LVL1);
        assert_eq!(c, c_exp, "vector {}: c mismatch", v["index"]);
        assert_eq!(d, d_exp, "vector {}: d mismatch", v["index"]);

        // The determinant relation a*d - b*c == k*q (mod 2^r).
        let ad = ((a as u128) & mask).wrapping_mul(d) & mask;
        let bc = ((b as u128) & mask).wrapping_mul(c) & mask;
        let kq = (k & mask).wrapping_mul(q & mask) & mask;
        assert_eq!(ad.wrapping_sub(bc) & mask, kq, "determinant relation");
    }
    println!("stage-3 response recovery (c,d) matches the oracle for all 5 vectors");
}

/// Build the per-vector verify inputs from the oracle data.
struct Loaded {
    a_pk: sqisign_verify::Fp2<Level1>,
    a_com: sqisign_verify::Fp2<Level1>,
    chal: [u8; 32],
    f1: Vec<[Pt; 4]>,
    f2: Vec<[Pt; 4]>,
}
fn load_all() -> Vec<Loaded> {
    let main = load(PHASE0_VECTORS);
    let chain = load(CHAIN_VECTORS);
    let mut out = Vec::new();
    for v in main["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let cv = chain["vectors"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["index"].as_u64().unwrap() == vi)
            .unwrap();
        out.push(Loaded {
            a_pk: parse_fp2(&v["public_key"]["A_pk"]),
            a_com: parse_fp2(&v["signature"]["A_com"]),
            chal: chal_le(&v["signature"]),
            f1: parse_kernels(&cv["F1_kernels"]),
            f2: parse_kernels(&cv["F2_dual_kernels"]),
        });
    }
    out
}

/// Deliverable 3: the full verify accepts all 5 valid signatures and rejects
/// tampered inputs.
#[test]
fn end_to_end_accept_and_reject() {
    let n = hd_challenge_len::<Level1>();
    let loaded = load_all();
    assert_eq!(loaded.len(), 5);

    for (vi, l) in loaded.iter().enumerate() {
        let inp = HdVerifyInputs {
            a_pk: l.a_pk.clone(),
            a_com: l.a_com.clone(),
            message: &MSG,
            claimed_chal: &l.chal[..n],
            f1_kernels: &l.f1,
            f2_dual_kernels: &l.f2,
        };
        assert!(hd_verify(&inp), "vector {vi}: valid signature must verify");

        // Tamper the message -> the recomputed challenge no longer matches.
        let mut bad_msg = MSG;
        bad_msg[0] ^= 1;
        let inp_msg = HdVerifyInputs {
            a_pk: l.a_pk.clone(),
            a_com: l.a_com.clone(),
            message: &bad_msg,
            claimed_chal: &l.chal[..n],
            f1_kernels: &l.f1,
            f2_dual_kernels: &l.f2,
        };
        assert_eq!(
            hd_verify_checked(&inp_msg),
            Err(HdReject::ChallengeMismatch),
            "vector {vi}: tampered message must be rejected"
        );

        // Tamper the signature's challenge bytes directly.
        let mut bad_chal = l.chal;
        bad_chal[0] ^= 1;
        let inp_chal = HdVerifyInputs {
            a_pk: l.a_pk.clone(),
            a_com: l.a_com.clone(),
            message: &MSG,
            claimed_chal: &bad_chal[..n],
            f1_kernels: &l.f1,
            f2_dual_kernels: &l.f2,
        };
        assert_eq!(
            hd_verify_checked(&inp_chal),
            Err(HdReject::ChallengeMismatch)
        );

        // Tamper the LAST F1 kernel (the one that determines F1's last
        // codomain) -> the middle-codomain check fails.
        let mut f1_bad = l.f1.clone();
        let last = f1_bad.len() - 1;
        let mut coords: [sqisign_verify::Fp2<Level1>; 16] =
            core::array::from_fn(|i| f1_bad[last][0].coords()[i].clone());
        coords[0] = coords[0].add(&sqisign_verify::Fp2::<Level1>::one());
        f1_bad[last][0] = ThetaPointDim4::new(coords);
        let inp_kt = HdVerifyInputs {
            a_pk: l.a_pk.clone(),
            a_com: l.a_com.clone(),
            message: &MSG,
            claimed_chal: &l.chal[..n],
            f1_kernels: &f1_bad,
            f2_dual_kernels: &l.f2,
        };
        let res = hd_verify_checked(&inp_kt);
        assert!(
            matches!(
                res,
                Err(HdReject::MiddleCodomainMismatch) | Err(HdReject::ChainFailed)
            ),
            "vector {vi}: tampered kernel must be rejected (got {res:?})"
        );
    }
    println!("end-to-end: accept for all 5 valid signatures; reject on tampered message/challenge/kernel");
}

/// Deliverable 4: end-to-end verify timing.
#[test]
fn verify_timing_report() {
    let n = hd_challenge_len::<Level1>();
    let loaded = load_all();

    let iters = 5usize;
    let t0 = Instant::now();
    for _ in 0..iters {
        for l in &loaded {
            let inp = HdVerifyInputs {
                a_pk: l.a_pk.clone(),
                a_com: l.a_com.clone(),
                message: &MSG,
                claimed_chal: &l.chal[..n],
                f1_kernels: &l.f1,
                f2_dual_kernels: &l.f2,
            };
            black_box(hd_verify(&inp));
        }
    }
    let per_sig_ms = t0.elapsed().as_secs_f64() / (iters * loaded.len()) as f64 * 1e3;

    // Isolate the challenge (hash) cost.
    let l = &loaded[0];
    let bn = 2000usize;
    let th = Instant::now();
    let mut out = [0u8; 64];
    for _ in 0..bn {
        black_box(hd_challenge_from_curves(
            &l.a_com,
            &l.a_pk,
            &MSG,
            &mut out[..n],
        ));
    }
    let hash_us = th.elapsed().as_secs_f64() / bn as f64 * 1e6;

    println!("\n=========== PHASE 5 END-TO-END VERIFY TIMING (Level 1, unoptimized) ===========");
    println!(
        "full hd_verify (challenge + dim-4 chain + middle check): {per_sig_ms:.3} ms/signature"
    );
    println!(
        "  of which hash-to-challenge (j-invariants + SHAKE256): {hash_us:.1} us (negligible)"
    );
    println!("  the dim-4 chain dominates (see Phase 3 chain timing).");
    println!("NOTE: stages 1-3 kernel derivation and stage 6 HD-image are not yet computed;");
    println!("      the per-step kernels are supplied from the oracle (see NOTES Phase 5).");
    println!("================================================================================\n");
    assert!(per_sig_ms > 0.0);
}
