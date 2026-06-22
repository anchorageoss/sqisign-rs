//! Phase 7: a stable benchmark for the full raw-bytes verifier
//! ([`hd_verify_bytes_l1`]) at Level 1. Run with:
//!
//! ```text
//! cargo test -p sqisign-verify --release --test hd_bench -- --nocapture
//! ```
//!
//! Reports the mean and min per-signature time over several iterations of all
//! five reference vectors (min is the most stable cross-run number). This is the
//! measurement harness for the before/after optimization table.

mod hd_common;
use hd_common::{load, parse_fp2, PHASE0_VECTORS};

use crypto_bigint::U256;
use serde_json::Value;
use sqisign_verify::hd::{encode_public_key, encode_signature, hd_verify_bytes_l1_bool};
use std::hint::black_box;
use std::time::Instant;

const MSG: [u8; 32] = [0u8; 32];

fn dec_u256(s: &str) -> U256 {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    U256::from_words(limbs)
}
fn dec_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}

fn wire_of(v: &Value) -> (Vec<u8>, Vec<u8>) {
    let sig = &v["signature"];
    let s = encode_signature(
        &parse_fp2(&sig["A_com"]),
        dec_i128(&sig["a"]),
        dec_i128(&sig["b"]),
        dec_i128(&sig["c_or_d"]),
        &dec_u256(sig["q"].as_str().unwrap()),
        sig["hint_com_P"].as_u64().unwrap() as u32,
        sig["hint_com_Q"].as_u64().unwrap() as u32,
    )
    .unwrap();
    let p = encode_public_key(
        &parse_fp2(&v["public_key"]["A_pk"]),
        v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32,
        v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32,
    )
    .unwrap();
    (s.to_vec(), p.to_vec())
}

#[test]
fn verify_bytes_benchmark() {
    let doc = load(PHASE0_VECTORS);
    let wires: Vec<(Vec<u8>, Vec<u8>)> = doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .map(wire_of)
        .collect();
    assert_eq!(wires.len(), 5);

    // Warm-up (also confirms all accept).
    for (s, p) in &wires {
        assert!(
            hd_verify_bytes_l1_bool(s, p, &MSG),
            "valid signature must verify"
        );
    }

    let iters = 8usize;
    let mut min_per_sig = f64::INFINITY;
    let mut total = 0.0;
    for _ in 0..iters {
        let t0 = Instant::now();
        for (s, p) in &wires {
            black_box(hd_verify_bytes_l1_bool(
                black_box(s),
                black_box(p),
                black_box(&MSG),
            ));
        }
        let per_sig = t0.elapsed().as_secs_f64() / wires.len() as f64 * 1e3;
        min_per_sig = min_per_sig.min(per_sig);
        total += per_sig;
    }
    let mean = total / iters as f64;

    println!("\n========= PHASE 7 VERIFY BENCHMARK (Level 1, --release) =========");
    println!(
        "hd_verify_bytes_l1 over {} vectors × {iters} iters:",
        wires.len()
    );
    println!("  mean: {mean:.2} ms/signature");
    println!("  min:  {min_per_sig:.2} ms/signature  (most stable)");
    println!("=================================================================\n");
    assert!(min_per_sig > 0.0);
}
