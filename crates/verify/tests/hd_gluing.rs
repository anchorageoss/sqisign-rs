//! Phase 4: the dim-4 gluing isogeny + special_image, and the end-to-end
//! chain (gluing + plain steps) for all 5 vectors.
//!
//! Ground truth:
//!   * `tests/gluing_vectors.json` - gluing internals (5-point kernel,
//!     codomain with zero dual null) and `special_image` triples, from
//!     `sqisignhd-harness/extract_gluing.py`.
//!   * `tests/chain_vectors.json` - the plain-step kernels (Phase 3).
//!   * `test_vectors_l1.json` - per-step codomains + middle-codomain match.

mod hd_common;
use hd_common::{load, parse_coords, parse_node, Pt, PHASE0_VECTORS};

use serde_json::Value;
use sqisign_verify::hd::{
    middle_codomain_matches, GluingIsogenyDim4, IsogenyDim4, ThetaPointDim4, GLUING_KERNEL_DIRS,
};

const GLUING_VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/gluing_vectors.json");
const CHAIN_VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/chain_vectors.json");

/// Parse a fixed array of `M` theta points from a JSON list of coord-arrays.
fn parse_points<const M: usize>(v: &Value) -> [Pt; M] {
    let arr = v.as_array().expect("points array");
    assert_eq!(arr.len(), M, "expected {M} points");
    core::array::from_fn(|k| parse_coords(&arr[k]))
}

fn parse_dirs(v: &Value) -> [usize; 5] {
    let arr = v.as_array().expect("dirs array");
    assert_eq!(arr.len(), 5);
    core::array::from_fn(|k| arr[k].as_u64().unwrap() as usize)
}

/// Build the gluing isogeny from an oracle `gluings[..]` entry.
fn gluing_from(entry: &Value) -> GluingIsogenyDim4<sqisign_verify::Level1> {
    let k8: [Pt; 5] = parse_points(&entry["L_K_8"]);
    let dirs = parse_dirs(&entry["L_K_8_ind"]);
    assert_eq!(
        dirs, GLUING_KERNEL_DIRS,
        "kernel directions must be [1,2,4,8,3]"
    );
    GluingIsogenyDim4::from_kernel(&k8, &dirs).expect("gluing codomain computable")
}

/// Deliverable 4a: the gluing codomain matches the oracle and the Phase 0
/// `chain_codomains[chain][0]`, with the expected zero-dual structure.
#[test]
fn gluing_codomain_matches_oracle() {
    let glu = load(GLUING_VECTORS);
    let main = load(PHASE0_VECTORS);

    let mut n = 0;
    for vec in glu["vectors"].as_array().unwrap() {
        let vi = vec["index"].as_u64().unwrap() as usize;
        for entry in vec["gluings"].as_array().unwrap() {
            let chain = entry["chain"].as_str().unwrap();
            let iso = gluing_from(entry);

            // Codomain has zero dual theta-null coordinates (6 at Level 1).
            assert_eq!(
                iso.dual_zero_count(),
                entry["dual_zero_count"].as_u64().unwrap() as usize,
                "dual zero count mismatch"
            );
            assert!(
                iso.dual_zero_count() > 0,
                "gluing codomain must have zero dual coords"
            );

            // Matches the oracle's recorded codomain null point.
            let expected = parse_coords(&entry["codomain_null"]);
            assert!(
                iso.codomain_null().projective_eq(&expected),
                "vector {vi} {chain}: gluing codomain != oracle"
            );

            // Matches the Phase 0 chain_codomains[chain][0] (the gluing is the
            // first codomain of the half-chain).
            let c0 = &main["test_vectors"][vi]["stage4_compute_hd"][chain]["chain_codomains"][0];
            let from_main = parse_node(c0);
            assert!(
                iso.codomain_null().projective_eq(&from_main),
                "vector {vi} {chain}: gluing codomain != test_vectors_l1.json chain_codomains[0]"
            );
            n += 1;
        }
    }
    assert_eq!(n, 10, "expected 2 gluing objects x 5 vectors");
    println!("validated {n} gluing codomains (zero-dual, vs oracle + test_vectors_l1.json)");
}

/// Deliverable 4b: special_image reproduces the reference for captured triples.
#[test]
fn special_image_matches_oracle() {
    let glu = load(GLUING_VECTORS);

    let mut n = 0;
    for vec in glu["vectors"].as_array().unwrap() {
        let gluings = vec["gluings"].as_array().unwrap();
        // Reconstruct both gluing objects once.
        let isos: Vec<_> = gluings.iter().map(gluing_from).collect();

        for pair in vec["special_image_pairs"].as_array().unwrap() {
            let gi = pair["gluing_index"].as_u64().unwrap() as usize;
            let p = parse_coords(&pair["P"]);
            let l_trans: [Pt; 2] = parse_points(&pair["L_trans"]);
            let l_trans_ind: Vec<usize> = pair["L_trans_ind"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();
            let expected = parse_coords(&pair["out"]);

            let got = isos[gi].special_image(&p, &l_trans, &l_trans_ind);
            assert!(got.projective_eq(&expected), "special_image mismatch");
            n += 1;
        }
    }
    assert!(n >= 5, "expected several special_image pairs");
    println!("validated {n} special_image evaluations across 5 vectors");
}

/// Deliverable 4c: the FULL chain (gluing + plain steps) reproduces every
/// codomain and the middle-codomain match, for all 5 vectors.
#[test]
fn full_chain_with_gluing() {
    let glu = load(GLUING_VECTORS);
    let chain_doc = load(CHAIN_VECTORS);
    let main = load(PHASE0_VECTORS);

    // Index the gluing entries and plain kernels by (vector, chain).
    let chain_vecs = chain_doc["vectors"].as_array().unwrap();

    let mut total_codomains = 0usize;
    for gvec in glu["vectors"].as_array().unwrap() {
        let vi = gvec["index"].as_u64().unwrap() as usize;
        let cvec = chain_vecs
            .iter()
            .find(|c| c["index"].as_u64().unwrap() as usize == vi)
            .unwrap();

        let mut last: [Option<Pt>; 2] = [None, None]; // [F1, F2_dual]
        for (ci, (chain, plain_key)) in [("F1", "F1_kernels"), ("F2_dual", "F2_dual_kernels")]
            .iter()
            .enumerate()
        {
            let gentry = gvec["gluings"]
                .as_array()
                .unwrap()
                .iter()
                .find(|e| e["chain"].as_str().unwrap() == *chain)
                .unwrap();

            // Full codomain sequence: index 0 from the gluing, 1.. from plain.
            let mut full: Vec<Pt> = Vec::new();
            full.push(gluing_from(gentry).codomain_null().clone());
            for k8j in cvec[plain_key].as_array().unwrap() {
                let k8: [Pt; 4] = parse_points(k8j);
                let iso = IsogenyDim4::from_kernel(&k8).expect("plain step computable");
                full.push(iso.codomain_null().clone());
            }

            // Reproduces every recorded chain codomain.
            let recorded = main["test_vectors"][vi]["stage4_compute_hd"][chain]["chain_codomains"]
                .as_array()
                .unwrap();
            assert_eq!(full.len(), recorded.len(), "chain length mismatch {chain}");
            for (i, cod) in full.iter().enumerate() {
                let expected = parse_node(&recorded[i]);
                assert!(
                    cod.projective_eq(&expected),
                    "vector {vi} {chain}: full-chain codomain {i} mismatch"
                );
            }
            total_codomains += full.len();
            last[ci] = full.last().cloned();
        }

        // Middle-codomain match from the fully computed chains (incl. gluing).
        let f1_last = last[0].as_ref().unwrap();
        let f2_last = last[1].as_ref().unwrap();
        assert!(
            middle_codomain_matches(f1_last, f2_last),
            "vector {vi}: middle-codomain check failed"
        );
        assert_eq!(
            main["test_vectors"][vi]["stage5_codomain_check"]["match"].as_bool(),
            Some(true)
        );
        // Tie to the recorded middle codomains.
        let c1 = parse_node(&main["test_vectors"][vi]["stage5_codomain_check"]["C1_zero"]);
        assert!(
            f1_last.projective_eq(&c1),
            "vector {vi}: F1 last != C1.zero()"
        );
        let hc2 = parse_node(&main["test_vectors"][vi]["stage5_codomain_check"]["HC2_zero"]);
        let my_hc2 = ThetaPointDim4::new(sqisign_verify::hd::hadamard(f2_last.coords()));
        assert!(
            my_hc2.projective_eq(&hc2),
            "vector {vi}: H(F2 last) != HC2.zero()"
        );
    }
    println!(
        "FULL chain (gluing + plain) reproduced {total_codomains} codomains and the \
         middle-codomain match for all 5 vectors"
    );
}
