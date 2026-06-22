//! Phase 5b.5: dim-1 and product theta structures, validated against the
//! additive oracle `product_theta_l1.json` (from
//! `sqisignhd-harness/extract_product_theta.py`) for all 5 vectors, both
//! half-chains. Theta structures are projective, so comparison is projective
//! (Phase-1 cross-multiply).
//!
//! Checks, per half-chain:
//! * the dim-1 theta null `(X+Z, X-Z)` of `E_com` and `E_chal` from their
//!   canonical 4-torsion points;
//! * `montgomery_to_theta` on a generic point of each curve;
//! * the dim-2 product null `Theta12` = `dim1 ⊗ dim1`, and the product of the
//!   two generic-point images;
//! * the dim-4 product null `domain_product` = `dim2_codomain ⊗ dim2_codomain`
//!   (the gluing's pre-base-change domain, == the oracle's
//!   `theta_null_domain_product`).

mod hd_common;
use hd_common::{fp2_eq, load, parse_fp2, F, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::hd::{product_theta_dim2, product_theta_dim2to4, ThetaStructureDim1};

const PRODUCT_THETA: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../sqisignhd-harness/product_theta_l1.json");

/// Parse a list of `[re, im]` pairs into `Vec<Fp2>`.
fn parse_vec(node: &Value) -> Vec<F> {
    node.as_array().unwrap().iter().map(parse_fp2).collect()
}

/// Projective equality of two equal-length coordinate vectors.
fn proj_eq(a: &[F], b: &[F]) -> bool {
    assert_eq!(a.len(), b.len());
    match a.iter().position(|x| !bool::from(x.ct_is_zero())) {
        None => b.iter().all(|x| bool::from(x.ct_is_zero())),
        Some(p) => a
            .iter()
            .zip(b.iter())
            .all(|(ai, bi)| fp2_eq(&ai.mul(&b[p]), &bi.mul(&a[p]))),
    }
}

#[test]
fn product_theta_structures_match_oracle() {
    let doc = load(PRODUCT_THETA);
    let mut n = 0;
    for v in doc["vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        for hc in v["half_chains"].as_array().unwrap() {
            let chain = hc["chain"].as_str().unwrap();
            let tag = format!("vec {vi} {chain}");

            // dim-1 theta structures on E_com (T1) and E_chal (U1).
            let t1 = parse_vec(&hc["T1_xz"]);
            let u1 = parse_vec(&hc["U1_xz"]);
            let s1 = ThetaStructureDim1::from_torsion(&t1[0], &t1[1]);
            let s2 = ThetaStructureDim1::from_torsion(&u1[0], &u1[1]);
            let n1 = parse_vec(&hc["dim1_null_1"]);
            let n2 = parse_vec(&hc["dim1_null_2"]);
            assert!(proj_eq(s1.null(), &n1), "{tag}: dim1_null_1");
            assert!(proj_eq(s2.null(), &n2), "{tag}: dim1_null_2");

            // montgomery_to_theta on a generic point of each curve.
            let r1 = parse_vec(&hc["conv_pt_1_xz"]);
            let r2 = parse_vec(&hc["conv_pt_2_xz"]);
            let ct1 = s1.montgomery_to_theta(&r1[0], &r1[1]);
            let ct2 = s2.montgomery_to_theta(&r2[0], &r2[1]);
            assert!(proj_eq(&ct1, &parse_vec(&hc["conv_theta_1"])), "{tag}: conv_theta_1");
            assert!(proj_eq(&ct2, &parse_vec(&hc["conv_theta_2"])), "{tag}: conv_theta_2");

            // dim-2 product: Theta12 null and the product of the two images.
            let th12 = product_theta_dim2(s1.null(), s2.null());
            assert!(proj_eq(&th12, &parse_vec(&hc["theta12_null"])), "{tag}: theta12_null");
            let prod_conv = product_theta_dim2(&ct1, &ct2);
            assert!(proj_eq(&prod_conv, &parse_vec(&hc["conv_prod"])), "{tag}: conv_prod");

            // dim-4 product: domain_product = dim2_codomain ⊗ dim2_codomain.
            let cod = parse_vec(&hc["dim2_codomain_null"]);
            let cod4: [F; 4] = core::array::from_fn(|k| cod[k].clone());
            let dp = product_theta_dim2to4(&cod4, &cod4);
            assert!(
                proj_eq(&dp, &parse_vec(&hc["domain_product_null"])),
                "{tag}: domain_product_null"
            );

            n += 1;
        }
    }
    assert_eq!(n, 10, "5 vectors x 2 half-chains");
    println!("dim-1 + product theta structures match the oracle for all {n} half-chains");
}

/// Deliverable 2: the self-built product domain is exactly the
/// `theta_null_domain_product` that the validated Phase-3/4 dim-4 chain
/// consumes (recorded in `test_vectors_l1.json` stage 4). Building it from the
/// dim-2 codomain with `product_theta_dim2to4` reproduces that value, so the
/// product-theta layer correctly forms the gluing's (pre-base-change) domain -
/// which the already-validated Phase 4 gluing then maps to the recorded
/// codomain.
#[test]
fn product_domain_feeds_validated_gluing() {
    let pt = load(PRODUCT_THETA);
    let main = load(PHASE0_VECTORS);
    let mut n = 0;
    for (v, mv) in pt["vectors"]
        .as_array()
        .unwrap()
        .iter()
        .zip(main["test_vectors"].as_array().unwrap().iter())
    {
        assert_eq!(v["index"], mv["index"]);
        let s4 = &mv["stage4_compute_hd"];
        for hc in v["half_chains"].as_array().unwrap() {
            let chain = hc["chain"].as_str().unwrap(); // "F1" | "F2_dual"
            let cod = parse_vec(&hc["dim2_codomain_null"]);
            let cod4: [F; 4] = core::array::from_fn(|k| cod[k].clone());
            let built = product_theta_dim2to4(&cod4, &cod4);
            let recorded = parse_vec(&s4[chain]["theta_null_domain_product"]["coords"]);
            assert!(
                proj_eq(&built, &recorded),
                "vec {} {chain}: product domain != chain's theta_null_domain_product",
                mv["index"]
            );
            n += 1;
        }
    }
    assert_eq!(n, 10);
    println!("self-built product domain == the dim-4 chain's theta_null_domain_product for all {n} half-chains");
}
