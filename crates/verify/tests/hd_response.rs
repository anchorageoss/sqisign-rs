//! Phase 5b.3: response image recovery validated against the oracle's
//! recorded `stage3_image_response` for all 5 vectors.
//!
//! This stage consumes the **self-derived** challenge data from 5b.2
//! (`recover_challenge_l1`) - `E_chal`, `w_chal`, and the rescaled basis - and
//! recovers the discrete log `k`, the response scalars `(c, d)`, and the
//! response isogeny images. The central check is that `k` matches the oracle
//! **exactly**: 5b.2's `w_chal` is the oracle's inverse (biextension `weil` vs
//! PARI), and computing `w_com` with the same `weil` makes that convention
//! cancel in the discrete log, as 5b.2 predicted.

mod hd_common;
use hd_common::{fp2_eq, load, parse_fp2, F, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::ec::JacPoint;
use sqisign_verify::hd::{
    jac_to_affine, recover_challenge_l1, recover_response_l1, ResponseScalars,
};
use sqisign_verify::Level1;

const R_LVL1: u32 = 70;

/// Parse a non-negative decimal string into little-endian u64 limbs (the
/// Level-1 challenge is < 2^256).
fn decimal_to_le_limbs(s: &str) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    limbs
}

fn dec_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}
fn dec_u128(v: &Value) -> u128 {
    v.as_str().unwrap().parse::<u128>().unwrap()
}
/// `q` exceeds 128 bits (q < 2^136), but only `q mod 2^r` is used; reduce mod
/// 2^128 by wrapping decimal accumulation (preserves the low 70 bits).
fn dec_u128_wrap(v: &Value) -> u128 {
    let mut acc = 0u128;
    for ch in v.as_str().unwrap().trim().bytes() {
        acc = acc.wrapping_mul(10).wrapping_add((ch - b'0') as u128);
    }
    acc
}

fn check_point(jac: &JacPoint<Level1>, node: &Value, label: &str) {
    let (x, y) = jac_to_affine(jac);
    let ex: F = parse_fp2(&node["x"]);
    let ey: F = parse_fp2(&node["y"]);
    assert!(fp2_eq(&x, &ex), "{label}: x mismatch");
    assert!(fp2_eq(&y, &ey), "{label}: y mismatch");
}

#[test]
fn response_images_match_oracle() {
    let doc = load(PHASE0_VECTORS);
    let mask = (1u128 << R_LVL1) - 1;
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let s1 = &v["stage1_recover_basis"];
        let s3 = &v["stage3_image_response"];
        let sig = &v["signature"];

        // Self-derived challenge data (5b.2).
        let a_pk = parse_fp2(&s1["E_pk_A"]);
        let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        let chal_limbs = decimal_to_le_limbs(sig["chal"].as_str().unwrap());
        let chal = recover_challenge_l1(&a_pk, hp, hq, &chal_limbs).expect("challenge recovery");

        // Stage 3, consuming the self-derived `chal`.
        let a_com = parse_fp2(&s1["E_com_A"]);
        let hcp = sig["hint_com_P"].as_u64().unwrap() as u32;
        let hcq = sig["hint_com_Q"].as_u64().unwrap() as u32;
        let a = dec_i128(&sig["a"]);
        let b = dec_i128(&sig["b"]);
        let c_or_d = dec_i128(&sig["c_or_d"]);
        let q = dec_u128_wrap(&sig["q"]);

        let s = ResponseScalars { a, b, c_or_d, q };
        let rsp = recover_response_l1(&chal, &a_com, hcp, hcq, s).expect("response recovery");

        // (2) The KEY check: the discrete log matches the oracle EXACTLY,
        // confirming the w_chal inverse convention cancels.
        let k_oracle = dec_u128(&s3["k"]);
        assert_eq!(
            rsp.k, k_oracle,
            "vec {vi}: k mismatch (convention did NOT cancel)"
        );

        // w_com is the oracle's inverse (same convention as w_chal).
        let w_com_oracle: F = parse_fp2(&s3["w_com"]);
        assert!(
            fp2_eq(&rsp.w_com.inv(), &w_com_oracle),
            "vec {vi}: w_com is not the oracle's inverse"
        );

        // (1) Response scalars (c, d) mod 2^r.
        let c_exp = (dec_i128(&s3["c"]) as u128) & mask;
        let d_exp = (dec_i128(&s3["d"]) as u128) & mask;
        assert_eq!(rsp.c, c_exp, "vec {vi}: c mismatch");
        assert_eq!(rsp.d, d_exp, "vec {vi}: d mismatch");

        // (3) The determinant relation a·d - b·c ≡ k·q (mod 2^r).
        let ad = ((a as u128) & mask).wrapping_mul(rsp.d) & mask;
        let bc = ((b as u128) & mask).wrapping_mul(rsp.c) & mask;
        let kq = (rsp.k & mask).wrapping_mul(q & mask) & mask;
        assert_eq!(
            ad.wrapping_sub(bc) & mask,
            kq,
            "vec {vi}: determinant relation"
        );

        // Response images and rescaled commitment basis (full points).
        check_point(&rsp.r_com, &s3["R_com"], &format!("vec {vi} R_com"));
        check_point(&rsp.s_com, &s3["S_com"], &format!("vec {vi} S_com"));
        check_point(
            &rsp.phi_rsp_r_com,
            &s3["phi_rsp_R_com"],
            &format!("vec {vi} phi_rsp_R_com"),
        );
        check_point(
            &rsp.phi_rsp_s_com,
            &s3["phi_rsp_S_com"],
            &format!("vec {vi} phi_rsp_S_com"),
        );

        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "stage-3 response image recovery matches the oracle for all {n} vectors; \
         the w_chal inverse convention cancels in the dlog (k exact)"
    );
}
