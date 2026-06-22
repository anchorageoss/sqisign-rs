//! Phase 5b.1: HD torsion-basis-from-hint recovery validated against the
//! oracle's recorded `stage1_recover_basis` points for all 5 vectors.
//!
//! The hints come from the signature/public key (`hint_pk_P/Q`,
//! `hint_com_P/Q`); the recovered `(P, Q)` must equal the oracle's affine
//! basis points exactly (the basis is canonical, so this is an exact match,
//! not merely projective).

mod hd_common;
use hd_common::{fp2_eq, load, parse_fp2, F, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::hd::{hd_torsion_basis_l1, jac_to_affine};
use sqisign_verify::Level1;

/// Assert that the lifted point matches the oracle's affine `{x, y}` node.
fn check(jac: &sqisign_verify::ec::JacPoint<Level1>, node: &Value, label: &str) {
    assert!(
        node.get("inf").is_none(),
        "{label}: oracle point is identity?"
    );
    let (x, y) = jac_to_affine(jac);
    let ex: F = parse_fp2(&node["x"]);
    let ey: F = parse_fp2(&node["y"]);
    assert!(fp2_eq(&x, &ex), "{label}: x mismatch");
    assert!(fp2_eq(&y, &ey), "{label}: y mismatch");
}

#[test]
fn torsion_bases_match_oracle() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let s1 = &v["stage1_recover_basis"];

        // Public-key curve basis.
        let a_pk = parse_fp2(&s1["E_pk_A"]);
        let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        let (p_pk, q_pk) = hd_torsion_basis_l1(&a_pk, hp, hq).expect("E_pk basis computable");
        check(&p_pk, &s1["P_pk"], &format!("vec {vi} P_pk"));
        check(&q_pk, &s1["Q_pk"], &format!("vec {vi} Q_pk"));

        // Commitment curve basis.
        let a_com = parse_fp2(&s1["E_com_A"]);
        let hcp = v["signature"]["hint_com_P"].as_u64().unwrap() as u32;
        let hcq = v["signature"]["hint_com_Q"].as_u64().unwrap() as u32;
        let (p_com, q_com) = hd_torsion_basis_l1(&a_com, hcp, hcq).expect("E_com basis computable");
        check(&p_com, &s1["P_com"], &format!("vec {vi} P_com"));
        check(&q_com, &s1["Q_com"], &format!("vec {vi} Q_com"));

        n += 1;
    }
    assert_eq!(n, 5);
    println!("HD torsion-basis-from-hint reproduced the oracle bases (P,Q on E_pk and E_com) for all {n} vectors");
}
