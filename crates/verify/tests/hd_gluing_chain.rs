//! Phase 5b.6 (front half): the dim-4 gluing chain
//! (`KaniGluingIsogenyChainDim4Half`) validated against the oracle.
//!
//! Stage 1 (this file) drives the assembly with the **oracle's** completed
//! starting matrices `M1/M2` and gluing matrices `M_gluing_1/2`
//! (`kani_matrices_l1.json`), isolating the dim-4 theta machinery
//! (`base_change_theta_dim4`, the product dim2→4, the Phase-4 `GluingIsogenyDim4`,
//! and the chain's `evaluate`) from the symplectic-completion choice. With the
//! oracle's matrices, the gluing output is completion-matched, so it must equal
//! the oracle's recorded `glue_codomain_null` and `post_glue_basis`
//! (`strategy_vectors.json`) - exactly, projectively - for all 5 vectors, both
//! half-chains. The self-derived-completion variant (the invariant middle match)
//! is exercised in `hd_self_contained.rs`.

mod hd_common;
use hd_common::{load, parse_coords, parse_fp2, Pt, F};

use serde_json::Value;
use sqisign_verify::ec::jacobian::{jac_add, jac_dbl};
use sqisign_verify::ec::pairing::weil;
use sqisign_verify::ec::{EcCurve, JacPoint};
use sqisign_verify::hd::{
    jac_mul_u128, make_canonical, point_matrix_product_k, product_theta_dim2, recover_challenge_l1,
    recover_response_l1, KaniGluingChainHalf, ResponseScalars, ThetaStructureDim1, TuplePoint4,
};
use sqisign_verify::Level1;

const PRODUCT_THETA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/product_theta_l1.json"
);
const KANI: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/kani_matrices_l1.json"
);
const STRATEGY: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/strategy_vectors.json"
);

const TO_4TORSION: u32 = 68;
const F_EXP: u32 = 70; // 2^f-torsion exponent; e1 = e2 = 68 ⇒ modulus 2^(e_i+2) = 2^70.

// parsing

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
    v.as_str().unwrap().trim().parse::<u128>().unwrap()
}
fn dec_u128_wrap(v: &Value) -> u128 {
    let mut acc = 0u128;
    for ch in v.as_str().unwrap().trim().bytes() {
        acc = acc.wrapping_mul(10).wrapping_add((ch - b'0') as u128);
    }
    acc
}
fn parse_m8_u128(v: &Value) -> [[u128; 8]; 8] {
    let rows = v.as_array().unwrap();
    core::array::from_fn(|i| {
        let row = rows[i].as_array().unwrap();
        core::array::from_fn(|j| dec_u128_wrap(&row[j]))
    })
}
fn parse_m8_i64(v: &Value) -> [[i64; 8]; 8] {
    let rows = v.as_array().unwrap();
    core::array::from_fn(|i| {
        let row = rows[i].as_array().unwrap();
        core::array::from_fn(|j| row[j].as_str().unwrap().parse::<i64>().unwrap())
    })
}

// EC helpers

fn jac_dbl_n(p: &JacPoint<Level1>, n: u32, curve: &EcCurve<Level1>) -> JacPoint<Level1> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, curve);
    }
    acc
}
fn weil4(u: &JacPoint<Level1>, v: &JacPoint<Level1>, curve: &mut EcCurve<Level1>) -> F {
    let uv = jac_add(u, &v.neg(), curve);
    weil(2, &u.to_xz(), &v.to_xz(), &uv.to_xz(), curve)
}
fn m0_from_canon(mt: &[[u8; 2]; 2], mu: &[[u8; 2]; 2]) -> [[u8; 4]; 4] {
    [
        [mt[0][0], 0, mt[1][0], 0],
        [0, mu[0][0], 0, mu[1][0]],
        [mt[0][1], 0, mt[1][1], 0],
        [0, mu[0][1], 0, mu[1][1]],
    ]
}
fn dim1_null(u1: &JacPoint<Level1>) -> [F; 2] {
    let (x, _) = sqisign_verify::hd::jac_to_affine(u1);
    ThetaStructureDim1::<Level1>::from_torsion(&x, &F::one())
        .null()
        .clone()
}

/// Everything the gluing chain consumes, self-derived from stages 2-3 plus the
/// canonical bases (all completion-independent). `m1`/`m2`/`mg1`/`mg2` are read
/// from the oracle (the completion-dependent matrices) for stage 1.
struct Setup {
    e_com: EcCurve<Level1>,
    e_chal: EcCurve<Level1>,
    points_m: [JacPoint<Level1>; 4],
    r_com: JacPoint<Level1>,
    s_com: JacPoint<Level1>,
    phi_r: JacPoint<Level1>,
    phi_s: JacPoint<Level1>,
    zero12: [F; 4],
    m0: [[u8; 4]; 4],
    e4: F,
    a1: u128,
    a2: u128,
    q4: u128,
    m: usize,
}

fn build_setup(v: &Value, vk: &Value) -> Setup {
    let s1 = &v["stage1_recover_basis"];
    let sig = &v["signature"];
    let a_pk = parse_fp2(&s1["E_pk_A"]);
    let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
    let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
    let chal_limbs = decimal_to_le_limbs(sig["chal"].as_str().unwrap());
    let chal = recover_challenge_l1(&a_pk, hp, hq, &chal_limbs).unwrap();

    let a_com = parse_fp2(&s1["E_com_A"]);
    let hcp = sig["hint_com_P"].as_u64().unwrap() as u32;
    let hcq = sig["hint_com_Q"].as_u64().unwrap() as u32;
    let q_wrap = dec_u128_wrap(&sig["q"]);
    let rsp = recover_response_l1(
        &chal,
        &a_com,
        hcp,
        hcq,
        ResponseScalars {
            a: dec_i128(&sig["a"]),
            b: dec_i128(&sig["b"]),
            c_or_d: dec_i128(&sig["c_or_d"]),
            q: q_wrap,
        },
    )
    .unwrap();

    let a1 = dec_u128(&vk["a1"]);
    let a2 = dec_u128(&vk["a2"]);
    let m = vk["m"].as_u64().unwrap() as usize;
    let q4 = dec_u128_wrap(&vk["q"]);
    let lamb = q4 & 3;

    let mut e_com = EcCurve::<Level1>::from_a(&a_com).unwrap();
    e_com.normalize_a24();
    let mut e_chal = chal.e_chal.clone();
    e_chal.normalize_a24();

    let p1_4 = jac_dbl_n(&rsp.r_com, TO_4TORSION, &e_com);
    let q1_4 = jac_dbl_n(&rsp.s_com, TO_4TORSION, &e_com);
    let (t1, t2, mt) = make_canonical(&p1_4, &q1_4, &mut e_com).unwrap();
    let r2_4 = jac_dbl_n(&rsp.phi_rsp_r_com, TO_4TORSION, &e_chal);
    let mut s2_4 = jac_dbl_n(&rsp.phi_rsp_s_com, TO_4TORSION, &e_chal);
    if lamb == 3 {
        s2_4 = jac_add(&jac_dbl(&s2_4, &e_chal), &s2_4, &e_chal);
    }
    let (u1, _u2, mu) = make_canonical(&r2_4, &s2_4, &mut e_chal).unwrap();

    let e4 = weil4(&t1, &t2, &mut e_com).inv();
    let zero12 = product_theta_dim2(&dim1_null(&t1), &dim1_null(&u1));
    let m0 = m0_from_canon(&mt, &mu);

    let k = F_EXP - 3 - m as u32; // 67 - m
    let points_m = [
        jac_dbl_n(&rsp.r_com, k, &e_com),
        jac_dbl_n(&rsp.s_com, k, &e_com),
        jac_dbl_n(&rsp.phi_rsp_r_com, k, &e_chal),
        jac_dbl_n(&rsp.phi_rsp_s_com, k, &e_chal),
    ];

    Setup {
        e_com,
        e_chal,
        points_m,
        r_com: rsp.r_com.clone(),
        s_com: rsp.s_com.clone(),
        phi_r: rsp.phi_rsp_r_com.clone(),
        phi_s: rsp.phi_rsp_s_com.clone(),
        zero12,
        m0,
        e4,
        a1,
        a2,
        q4,
        m,
    }
}

/// The full-order kernel basis `B_Kpp = kernel_basis(M, e_i, P1, Q1, R2, lamb·S2)`
/// with `e_i = 68` ⇒ modulus `2^70`.
fn b_kpp(setup: &Setup, m_full: &[[u128; 8]; 8]) -> [TuplePoint4<Level1>; 4] {
    let mask = (1u128 << F_EXP) - 1;
    let lamb = sqisign_verify::hd::inverse_mod_pow2(setup.q4, mask);
    let lamb_s2 = jac_mul_u128(&setup.phi_s, lamb, &setup.e_chal);
    point_matrix_product_k(
        m_full,
        &setup.r_com,
        &setup.s_com,
        &setup.phi_r,
        &lamb_s2,
        mask,
        &setup.e_com,
        &setup.e_chal,
    )
}

#[test]
fn gluing_chain_matches_oracle_with_oracle_matrices() {
    let doc = load(hd_common::PHASE0_VECTORS);
    let pt = load(PRODUCT_THETA);
    let kani = load(KANI);
    let sv = load(STRATEGY);
    let mut n = 0;

    for (((v, _vp), vk), vs) in doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .zip(pt["vectors"].as_array().unwrap().iter())
        .zip(kani["vectors"].as_array().unwrap().iter())
        .zip(sv["vectors"].as_array().unwrap().iter())
    {
        let vi = v["index"].as_u64().unwrap();
        assert_eq!(v["index"], vk["index"]);
        assert_eq!(v["index"], vs["index"]);

        let setup = build_setup(v, vk);

        let m1 = parse_m8_u128(&vk["M1"]);
        let m2 = parse_m8_u128(&vk["M2"]);
        let mg1 = parse_m8_i64(&vk["M_gluing_1"]);
        let mg2 = parse_m8_i64(&vk["M_gluing_2"]);

        for (hc_idx, (dual, m_full, m_glue)) in [(false, &m1, &mg1), (true, &m2, &mg2)]
            .into_iter()
            .enumerate()
        {
            let chain = KaniGluingChainHalf::new(
                &setup.points_m,
                &setup.zero12,
                &setup.m0,
                &setup.e4,
                setup.a1,
                setup.a2,
                setup.q4,
                setup.m,
                m_full,
                m_glue,
                dual,
                &setup.e_com,
                &setup.e_chal,
            )
            .expect("gluing chain computable");

            let hc = &vs["half_chains"][hc_idx];
            let tag = if dual { "F2_dual" } else { "F1" };

            // (a) gluing codomain null.
            let want_null = parse_coords(&hc["glue_codomain_null"]);
            assert!(
                chain.codomain_null().projective_eq(&want_null),
                "vec {vi} {tag}: glue_codomain_null mismatch"
            );

            // (b) post-gluing kernel basis = [chain(T) for T in B_Kpp].
            let basis = b_kpp(&setup, m_full);
            let want_basis = hc["post_glue_basis"].as_array().unwrap();
            assert_eq!(want_basis.len(), 4);
            for (k, t) in basis.iter().enumerate() {
                let got: Pt = chain.evaluate(t);
                let want = parse_coords(&want_basis[k]);
                assert!(
                    got.projective_eq(&want),
                    "vec {vi} {tag}: post_glue_basis[{k}] mismatch"
                );
            }
        }
        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "gluing chain (oracle matrices) reproduced glue_codomain_null + post_glue_basis \
         for all {n} vectors, both half-chains"
    );
}
