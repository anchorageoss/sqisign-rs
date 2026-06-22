//! Phase 5b.7: stage 6, the HD-image check, validated against the oracle's
//! recorded `stage6_image_check` (`test_vectors_l1.json`) for all 5 vectors.
//!
//! The oracle records `F(T)` for `T = (P_com, 0, 0, 0)`, the expected
//! `a₁·P_com` / `a₂·P_com`, and the boolean `correct`. The check
//! `F(T) = (±a₁·P_com, ±a₂·P_com, *, 0)` is on **elliptic points**, hence
//! completion-independent - so the self-derived result matches the oracle.
//!
//! The Kani decomposition `(a1, a2)` is unordered up to the parity swap, so the
//! image's first two components are compared to `{a₁·P_com, a₂·P_com}` as a set
//! (by `x`-coordinate, i.e. up to `±`). `F(T)[3]` must be the identity.

mod hd_common;
use hd_common::{load, parse_fp2, F, PHASE0_VECTORS};

use crypto_bigint::U256;
use serde_json::Value;
use sqisign_verify::ec::JacPoint;
use sqisign_verify::hd::{hd_image_l1, hd_verify_l1, jac_to_affine, HdSignatureL1};
use sqisign_verify::Level1;

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

struct Owned {
    a_pk: F,
    a_com: F,
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
fn sig_of<'a>(o: &'a Owned) -> HdSignatureL1<'a> {
    HdSignatureL1 {
        a_pk: o.a_pk.clone(),
        a_com: o.a_com.clone(),
        hint_pk_p: o.hint_pk_p,
        hint_pk_q: o.hint_pk_q,
        hint_com_p: o.hint_com_p,
        hint_com_q: o.hint_com_q,
        message: &MSG,
        chal_limbs: &o.chal_limbs,
        claimed_chal: &o.chal_bytes,
        resp_a: o.a,
        resp_b: o.b,
        resp_c_or_d: o.c_or_d,
        q: o.q,
    }
}

/// Affine `x` of a non-identity Jacobian point.
fn ax(p: &JacPoint<Level1>) -> F {
    jac_to_affine(p).0
}
/// Does `x(p)` equal the oracle affine `x` (i.e. `p = ±(oracle point)`)?
fn x_is(p: &JacPoint<Level1>, oracle_x: &F) -> bool {
    bool::from(ax(p).ct_equal(oracle_x))
}

#[test]
fn stage6_matches_oracle() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let s6 = &v["stage6_image_check"];
        assert!(
            s6["correct"].as_bool().unwrap(),
            "oracle vec {vi}: correct must be true"
        );

        let o = owned_of(v);
        let (ft, a1p, a2p) = hd_image_l1(&sig_of(&o)).expect("HD image computable");

        // Oracle's expected image components (by x, up to ±), as a set.
        let oa1 = parse_fp2(&s6["a1_P_com"]["x"]);
        let oa2 = parse_fp2(&s6["a2_P_com"]["x"]);

        // (a) The self-derived a₁·P_com, a₂·P_com match the oracle's set.
        assert!(
            (x_is(&a1p, &oa1) && x_is(&a2p, &oa2)) || (x_is(&a1p, &oa2) && x_is(&a2p, &oa1)),
            "vec {vi}: self-derived a_i·P_com != oracle (set)"
        );

        // (b) F(T)[0], F(T)[1] match {a₁·P_com, a₂·P_com} as a set (parity swap).
        assert!(
            (x_is(&ft.c[0], &oa1) && x_is(&ft.c[1], &oa2))
                || (x_is(&ft.c[0], &oa2) && x_is(&ft.c[1], &oa1)),
            "vec {vi}: F(T)[0..2] != ±(a1·P_com, a2·P_com) (set)"
        );

        // (c) F(T)[3] is the identity on E_chal (oracle records {inf: true}).
        assert!(
            s6["FT"][3]["inf"].as_bool().unwrap_or(false),
            "oracle FT[3] should be inf"
        );
        assert!(
            bool::from(ft.c[3].z.ct_is_zero()),
            "vec {vi}: F(T)[3] is not the identity"
        );

        // (d) Cross-check against the oracle's recorded FT[0], FT[1] x-coords.
        let oft0 = parse_fp2(&s6["FT"][0]["x"]);
        let oft1 = parse_fp2(&s6["FT"][1]["x"]);
        assert!(
            (x_is(&ft.c[0], &oft0) && x_is(&ft.c[1], &oft1))
                || (x_is(&ft.c[0], &oft1) && x_is(&ft.c[1], &oft0)),
            "vec {vi}: F(T)[0..2] != oracle FT[0..2] (set)"
        );

        // (e) The full verifier (with stage 6) accepts.
        assert_eq!(
            hd_verify_l1(&sig_of(&o)),
            Ok(()),
            "vec {vi}: full verify must accept"
        );
        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "stage 6 HD-image check reproduced the oracle's F(T) (= ±a_i·P_com, FT[3]=0) \
         and `correct` for all {n} vectors; full self-contained verify accepts"
    );
}
