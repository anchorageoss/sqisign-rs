//! Phase 5b.2: challenge isogeny recovery validated against the oracle's
//! recorded `stage2_recover_chal` for all 5 vectors.
//!
//! Stage 2 derives the challenge isogeny `φ_chal : E_pk → E_chal` from the
//! public-key curve and the challenge scalar, the image Weil pairing
//! `w_chal`, and the rescaled image basis `(P_chal_resc, Q_chal_resc)`. We
//! check each against the oracle:
//!
//! * `E_chal_A` and the full points `P_chal_resc`, `Q_chal_resc` are exact
//!   (the latter are pinned to the FESTA square-root sign convention).
//! * `w_chal` is checked up to the global Weil-pairing convention (the dim-2
//!   biextension `weil` vs PARI's `weil_pairing`): an exact match or its
//!   inverse. The downstream discrete log (stage 3) is convention-independent,
//!   so either is correct; the test records which holds.

mod hd_common;
use hd_common::{fp2_eq, load, parse_fp2, F, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::Level1;
use sqisign_verify::hd::{jac_to_affine, recover_challenge_l1};

/// Parse a non-negative decimal string into little-endian u64 limbs (4 limbs
/// cover the Level-1 challenge, which is < 2^256).
fn decimal_to_le_limbs(s: &str) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        assert!(ch.is_ascii_digit(), "non-digit in decimal scalar");
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
        assert_eq!(carry, 0, "decimal scalar exceeds 256 bits");
    }
    limbs
}

/// Assert that a Jacobian point matches the oracle's affine `{x, y}` node.
fn check_point(jac: &sqisign_verify::ec::JacPoint<Level1>, node: &Value, label: &str) {
    let (x, y) = jac_to_affine(jac);
    let ex: F = parse_fp2(&node["x"]);
    let ey: F = parse_fp2(&node["y"]);
    assert!(fp2_eq(&x, &ex), "{label}: x mismatch");
    assert!(fp2_eq(&y, &ey), "{label}: y mismatch");
}

#[test]
fn challenge_recovery_matches_oracle() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    let mut w_inverse_convention = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let s1 = &v["stage1_recover_basis"];
        let s2 = &v["stage2_recover_chal"];

        let a_pk = parse_fp2(&s1["E_pk_A"]);
        let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        let chal = decimal_to_le_limbs(v["signature"]["chal"].as_str().unwrap());

        let rec = recover_challenge_l1(&a_pk, hp, hq, &chal).expect("challenge recovery");

        // Codomain Montgomery coefficient (e_chal is normalised, so a == A/C).
        let e_chal_a: F = parse_fp2(&s2["E_chal_A"]);
        assert!(fp2_eq(&rec.e_chal.a, &e_chal_a), "vec {vi}: E_chal_A mismatch");

        // Rescaled image basis (full points; FESTA-sign sensitive).
        check_point(&rec.p_chal_resc, &s2["P_chal_resc"], &format!("vec {vi} P_chal_resc"));
        check_point(&rec.q_chal_resc, &s2["Q_chal_resc"], &format!("vec {vi} Q_chal_resc"));

        // Image pairing, up to the global Weil convention (value or inverse).
        let w_oracle: F = parse_fp2(&s2["w_chal"]);
        if fp2_eq(&rec.w_chal, &w_oracle) {
            // exact match
        } else if fp2_eq(&rec.w_chal.inv(), &w_oracle) {
            w_inverse_convention += 1;
        } else {
            panic!("vec {vi}: w_chal matches neither the oracle value nor its inverse");
        }

        n += 1;
    }
    assert_eq!(n, 5);
    assert!(
        w_inverse_convention == 0 || w_inverse_convention == n,
        "w_chal convention must be consistent across vectors (got {w_inverse_convention}/{n} inverted)"
    );
    let conv = if w_inverse_convention == 0 {
        "exact"
    } else {
        "inverse (biextension weil = PARI weil^-1; convention-independent downstream)"
    };
    println!("challenge isogeny recovery reproduced the oracle (E_chal, P/Q_chal_resc) for all {n} vectors; w_chal convention: {conv}");
}
