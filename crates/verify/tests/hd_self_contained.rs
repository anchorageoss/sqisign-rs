//! Phase 5b.6 (front half): the **end-to-end self-contained** Level-1 verifier
//! ([`hd_verify_l1`]). Every quantity the dim-4 chain consumes is derived from
//! the signature alone - no oracle data is read in this file (contrast the
//! Phase-5 `hd_verify.rs`, which fed per-step kernels from `chain_vectors.json`).
//!
//! From `(a_pk, a_com, hints, message, chal, a, b, c_or_d, q)` it runs: the
//! challenge binding, stages 2-3, the Kani norm equation (→ `a1,a2,m`), the
//! canonical bases, the self-derived symplectic/gluing matrices, the two gluing
//! chains, the optimal-strategy chain loop, and the middle-codomain match.
//!
//! Checks: accept all 5 valid signatures; reject tampering of the message, the
//! challenge bytes, the recovered `q`, and a response scalar.

mod hd_common;
use hd_common::{load, parse_fp2, PHASE0_VECTORS};

use crypto_bigint::U256;
use serde_json::Value;
use sqisign_verify::hd::{hd_verify_l1, hd_verify_l1_bool, HdReject, HdSignatureL1};
use std::hint::black_box;
use std::time::Instant;

/// The fixed signed message used by the reference test harness (32 zero bytes).
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
fn chal_limbs_of(sig: &Value) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for ch in sig["chal"].as_str().unwrap().trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    limbs
}
fn limbs_to_le_bytes(limbs: &[u64; 4]) -> [u8; 32] {
    let mut b = [0u8; 32];
    for (i, &l) in limbs.iter().enumerate() {
        b[i * 8..i * 8 + 8].copy_from_slice(&l.to_le_bytes());
    }
    b
}

/// The owned per-vector data the borrowed [`HdSignatureL1`] points into.
struct Owned {
    a_pk: sqisign_verify::Fp2<sqisign_verify::Level1>,
    a_com: sqisign_verify::Fp2<sqisign_verify::Level1>,
    hint_pk_p: u32,
    hint_pk_q: u32,
    hint_com_p: u32,
    hint_com_q: u32,
    chal_limbs: [u64; 4],
    chal_bytes: [u8; 32],
    a: i128,
    b: i128,
    c_or_d: i128,
    q: U256,
}

fn owned_of(v: &Value) -> Owned {
    let sig = &v["signature"];
    let chal_limbs = chal_limbs_of(sig);
    Owned {
        a_pk: parse_fp2(&v["public_key"]["A_pk"]),
        a_com: parse_fp2(&sig["A_com"]),
        hint_pk_p: v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32,
        hint_pk_q: v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32,
        hint_com_p: sig["hint_com_P"].as_u64().unwrap() as u32,
        hint_com_q: sig["hint_com_Q"].as_u64().unwrap() as u32,
        chal_limbs,
        chal_bytes: limbs_to_le_bytes(&chal_limbs),
        a: dec_i128(&sig["a"]),
        b: dec_i128(&sig["b"]),
        c_or_d: dec_i128(&sig["c_or_d"]),
        q: dec_u256(sig["q"].as_str().unwrap()),
    }
}

fn sig_of<'a>(o: &'a Owned, msg: &'a [u8]) -> HdSignatureL1<'a> {
    HdSignatureL1 {
        a_pk: o.a_pk.clone(),
        a_com: o.a_com.clone(),
        hint_pk_p: o.hint_pk_p,
        hint_pk_q: o.hint_pk_q,
        hint_com_p: o.hint_com_p,
        hint_com_q: o.hint_com_q,
        message: msg,
        chal_limbs: &o.chal_limbs,
        claimed_chal: &o.chal_bytes,
        resp_a: o.a,
        resp_b: o.b,
        resp_c_or_d: o.c_or_d,
        q: o.q,
    }
}

#[test]
fn self_contained_accepts_all_valid() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let o = owned_of(v);
        let sig = sig_of(&o, &MSG);
        assert_eq!(
            hd_verify_l1(&sig),
            Ok(()),
            "vec {vi}: self-contained verify must accept the valid signature"
        );
        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "self-contained hd_verify_l1 accepted all {n} valid signatures \
         (everything derived from the signature; no oracle input)"
    );
}

#[test]
fn self_contained_rejects_tampering() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let o = owned_of(v);

        // Tampered message → challenge binding fails.
        let mut bad_msg = MSG;
        bad_msg[0] ^= 1;
        assert_eq!(
            hd_verify_l1(&sig_of(&o, &bad_msg)),
            Err(HdReject::ChallengeMismatch),
            "vec {vi}: tampered message must be rejected"
        );

        // Tampered challenge bytes → binding fails.
        {
            let mut bad = o.clone();
            bad.chal_bytes[0] ^= 1;
            assert_eq!(
                hd_verify_l1(&sig_of(&bad, &MSG)),
                Err(HdReject::ChallengeMismatch),
                "vec {vi}: tampered challenge must be rejected"
            );
        }

        // Tampered recovered q → norm equation / chain / middle match fails.
        {
            let mut bad = o.clone();
            bad.q = bad.q.wrapping_add(&U256::from(2u64)); // keep parity, break the relation
            assert!(
                !hd_verify_l1_bool(&sig_of(&bad, &MSG)),
                "vec {vi}: tampered q must be rejected"
            );
        }

        // Tampered response scalar a → wrong response basis → middle mismatch.
        {
            let mut bad = o.clone();
            bad.a = bad.a.wrapping_add(2); // preserve parity of a (branch unchanged)
            assert!(
                !hd_verify_l1_bool(&sig_of(&bad, &MSG)),
                "vec {vi}: tampered response scalar must be rejected"
            );
        }

        n += 1;
    }
    assert_eq!(n, 5);
    println!("self-contained verify rejects tampered message / challenge / q / response scalar for all {n} vectors");
}

impl Owned {
    fn clone(&self) -> Owned {
        Owned {
            a_pk: self.a_pk.clone(),
            a_com: self.a_com.clone(),
            hint_pk_p: self.hint_pk_p,
            hint_pk_q: self.hint_pk_q,
            hint_com_p: self.hint_com_p,
            hint_com_q: self.hint_com_q,
            chal_limbs: self.chal_limbs,
            chal_bytes: self.chal_bytes,
            a: self.a,
            b: self.b,
            c_or_d: self.c_or_d,
            q: self.q,
        }
    }
}

#[test]
fn self_contained_timing_report() {
    let doc = load(PHASE0_VECTORS);
    let owned: Vec<Owned> = doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .map(owned_of)
        .collect();

    let t0 = Instant::now();
    for o in &owned {
        black_box(hd_verify_l1_bool(&sig_of(o, &MSG)));
    }
    let per_sig_ms = t0.elapsed().as_secs_f64() / owned.len() as f64 * 1e3;

    println!("\n===== PHASE 5b SELF-CONTAINED VERIFY TIMING (Level 1, unoptimized) =====");
    println!("full self-contained hd_verify_l1 per signature: {per_sig_ms:.1} ms");
    println!("  (challenge binding + stages 2-3 + norm equation + canonical bases");
    println!("   + self-derived M1/M2/N_dim2/N_dim4 + 2 gluing chains + 2 strategy");
    println!("   loops + stage 5 middle-codomain match + stage 6 HD-image:");
    println!("   F(T) through F2∘F1 incl. the dual chain). ALL 6 stages.");
    println!("========================================================================\n");
    assert!(per_sig_ms > 0.0);
}
