//! Phase 5b.6 (front half): `make_canonical` validated against the oracle's
//! recorded canonical 4-torsion (`product_theta_l1.json`'s `dim1_null_1` /
//! `dim1_null_2`) for all 5 vectors, both half-chains.
//!
//! The canonical basis is completion-independent (the biextension-`weil`
//! inverse cancels in `d1 = inv(a1)·b1` and in the pairing fixup), so it matches
//! the oracle exactly. This replays stages 2-3 to self-derive the response
//! basis, reduces it to the 4-torsion, runs `make_canonical`, and compares the
//! induced dim-1 theta null `(X+Z, X-Z)` to the oracle (projectively).

mod hd_common;
use hd_common::{load, parse_fp2, F};

use serde_json::Value;
use sqisign_verify::ec::jacobian::{jac_add, jac_dbl};
use sqisign_verify::ec::{EcCurve, JacPoint};
use sqisign_verify::Level1;
use sqisign_verify::hd::{
    jac_to_affine, make_canonical, recover_challenge_l1, recover_response_l1, ThetaStructureDim1,
    ResponseScalars,
};

const PRODUCT_THETA: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../sqisignhd-harness/product_theta_l1.json");

/// f - 2 = 70 - 2 doublings take an order-2^70 point to the 4-torsion.
const TO_4TORSION: u32 = 68;

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
fn dec_u128_wrap(v: &Value) -> u128 {
    let mut acc = 0u128;
    for ch in v.as_str().unwrap().trim().bytes() {
        acc = acc.wrapping_mul(10).wrapping_add((ch - b'0') as u128);
    }
    acc
}

fn jac_dbl_n(p: &JacPoint<Level1>, n: u32, curve: &EcCurve<Level1>) -> JacPoint<Level1> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, curve);
    }
    acc
}

fn proj_eq2(a: &[F; 2], b: &[F]) -> bool {
    // (a0:a1) == (b0:b1) projectively
    bool::from(a[0].mul(&b[1]).ct_equal(&b[0].mul(&a[1])))
}

/// The dim-1 null `(x+1, x-1)` induced by the canonical U1 (a 4-torsion point,
/// affine so Z = 1), compared to the oracle's recorded null.
fn check_canonical(u1: &JacPoint<Level1>, oracle_null: &Value, label: &str) {
    let (x, _) = jac_to_affine(u1);
    let s = ThetaStructureDim1::from_torsion(&x, &F::one());
    let want = [parse_fp2(&oracle_null[0]), parse_fp2(&oracle_null[1])];
    assert!(proj_eq2(s.null(), &want), "{label}: canonical dim-1 null mismatch");
}

#[test]
fn make_canonical_matches_oracle() {
    let doc = load(hd_common::PHASE0_VECTORS);
    let pt = load(PRODUCT_THETA);
    let mut n = 0;
    for (v, vp) in doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .zip(pt["vectors"].as_array().unwrap().iter())
    {
        let vi = v["index"].as_u64().unwrap();
        assert_eq!(v["index"], vp["index"]);
        let s1 = &v["stage1_recover_basis"];
        let sig = &v["signature"];

        // Self-derive the response basis (stages 2-3).
        let a_pk = parse_fp2(&s1["E_pk_A"]);
        let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        let chal_limbs = decimal_to_le_limbs(sig["chal"].as_str().unwrap());
        let chal = recover_challenge_l1(&a_pk, hp, hq, &chal_limbs).unwrap();

        let a_com = parse_fp2(&s1["E_com_A"]);
        let hcp = sig["hint_com_P"].as_u64().unwrap() as u32;
        let hcq = sig["hint_com_Q"].as_u64().unwrap() as u32;
        let q = dec_u128_wrap(&sig["q"]);
        let rsp = recover_response_l1(
            &chal,
            &a_com,
            hcp,
            hcq,
            ResponseScalars {
                a: dec_i128(&sig["a"]),
                b: dec_i128(&sig["b"]),
                c_or_d: dec_i128(&sig["c_or_d"]),
                q,
            },
        )
        .unwrap();

        // The canonical basis is the same for both half-chains; check vs F1's.
        let hc = &vp["half_chains"][0];

        // E_com canonical basis from (R_com, S_com) reduced to the 4-torsion.
        let mut e_com = EcCurve::<Level1>::from_a(&a_com).unwrap();
        e_com.normalize_a24();
        let p1 = jac_dbl_n(&rsp.r_com, TO_4TORSION, &e_com);
        let q1 = jac_dbl_n(&rsp.s_com, TO_4TORSION, &e_com);
        let (u1_com, _u2_com, _m_com) = make_canonical(&p1, &q1, &mut e_com).unwrap();
        check_canonical(&u1_com, &hc["dim1_null_1"], &format!("vec {vi} E_com"));

        // E_chal canonical basis from (phi_rsp_R_com, lamb*phi_rsp_S_com).
        let mut e_chal = chal.e_chal.clone();
        let r2 = jac_dbl_n(&rsp.phi_rsp_r_com, TO_4TORSION, &e_chal);
        let mut s2 = jac_dbl_n(&rsp.phi_rsp_s_com, TO_4TORSION, &e_chal);
        // lamb = q^{-1} mod 4 (q odd ⇒ lamb = q mod 4); lamb ∈ {1,3}.
        if (q & 3) == 3 {
            s2 = jac_add(&jac_dbl(&s2, &e_chal), &s2, &e_chal); // 3*s2
        }
        let (u1_chal, _u2_chal, _m_chal) = make_canonical(&r2, &s2, &mut e_chal).unwrap();
        check_canonical(&u1_chal, &hc["dim1_null_2"], &format!("vec {vi} E_chal"));

        n += 1;
    }
    assert_eq!(n, 5);
    println!("make_canonical reproduced the oracle's canonical dim-1 nulls (E_com, E_chal) for all {n} vectors");
}
