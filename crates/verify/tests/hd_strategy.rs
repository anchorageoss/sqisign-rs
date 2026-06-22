//! Phase 5b.6: the dim-4 optimal-strategy chain loop derives each half-chain's
//! per-step kernels from a single post-gluing kernel basis, and reproduces the
//! oracle chain + the middle-codomain match with those SELF-DERIVED kernels.
//!
//! Ground truth:
//! * `strategy_vectors.json` - the post-gluing kernel basis + gluing codomain
//!   per half-chain (the loop's starting point);
//! * `chain_vectors.json` - the oracle per-step plain kernels `K_8`;
//! * `test_vectors_l1.json` - `chain_codomains` and thus the middle match.
//!
//! Comparison is projective (theta points up to scaling). The derived kernels
//! and codomains are completion-independent given the starting basis (the
//! strategy only picks the doubling/pushforward path), so they match the oracle
//! exactly even though the symplectic completion is self-consistent (5b.4).

mod hd_common;
use hd_common::{load, parse_coords, parse_node, Pt, PHASE0_VECTORS};

use serde_json::Value;
use std::hint::black_box;
use std::time::Instant;
use sqisign_verify::hd::{middle_codomain_matches, run_strategy_chain, ThetaStructureDim4};

const STRATEGY_VECTORS: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../sqisignhd-harness/strategy_vectors.json");
const CHAIN_VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/chain_vectors.json");

/// Parse `[K0, K1, K2, K3]` (four 16-coord theta points) into a `[Pt; 4]`.
fn parse_basis(node: &Value) -> [Pt; 4] {
    let arr = node.as_array().unwrap();
    assert_eq!(arr.len(), 4);
    core::array::from_fn(|i| parse_coords(&arr[i]))
}

/// Run one half-chain's strategy loop and validate the derived kernels and
/// codomains against the oracle. Returns the last codomain (middle null point).
fn run_and_check(hc: &Value, oracle_kernels: &Value, chain_codomains: &Value, tag: &str) -> Pt {
    let start = ThetaStructureDim4::new(parse_coords(&hc["glue_codomain_null"]));
    let basis = parse_basis(&hc["post_glue_basis"]);
    let l = hc["n_plain"].as_u64().unwrap() as usize;

    // Gluing codomain captured here must equal chain_codomains[0].
    let cc = chain_codomains.as_array().unwrap();
    assert!(
        start.null_point().projective_eq(&parse_node(&cc[0])),
        "{tag}: glue codomain != chain_codomains[0]"
    );

    let sc = run_strategy_chain(&start, &basis, l).expect("strategy chain computable");

    // Derived per-step kernels match the oracle's recorded kernels (projectively).
    let ok = oracle_kernels.as_array().unwrap();
    assert_eq!(sc.kernels.len(), l, "{tag}: kernel count");
    assert_eq!(ok.len(), l, "{tag}: oracle kernel count");
    for (i, k8) in sc.kernels.iter().enumerate() {
        let want = ok[i].as_array().unwrap();
        for (j, kj) in k8.iter().enumerate() {
            assert!(
                kj.projective_eq(&parse_coords(&want[j])),
                "{tag}: derived kernel {i} point {j} mismatch"
            );
        }
    }

    // Derived codomains match chain_codomains[1..] (chain_codomains[0] = gluing).
    assert_eq!(sc.codomains.len(), l, "{tag}: codomain count");
    for (i, cod) in sc.codomains.iter().enumerate() {
        assert!(
            cod.projective_eq(&parse_node(&cc[i + 1])),
            "{tag}: derived codomain {i} mismatch"
        );
    }

    sc.last_codomain().expect("non-empty chain").clone()
}

#[test]
fn strategy_loop_derives_kernels_and_middle_match() {
    let sv = load(STRATEGY_VECTORS);
    let cv = load(CHAIN_VECTORS);
    let main = load(PHASE0_VECTORS);
    let mut n = 0;
    for (vs, (vc, vm)) in sv["vectors"].as_array().unwrap().iter().zip(
        cv["vectors"]
            .as_array()
            .unwrap()
            .iter()
            .zip(main["test_vectors"].as_array().unwrap().iter()),
    ) {
        let vi = vs["index"].as_u64().unwrap();
        assert_eq!(vs["index"], vc["index"]);
        assert_eq!(vs["index"], vm["index"]);
        let hcs = vs["half_chains"].as_array().unwrap();
        let s4 = &vm["stage4_compute_hd"];

        let f1_last = run_and_check(
            &hcs[0],
            &vc["F1_kernels"],
            &s4["F1"]["chain_codomains"],
            &format!("vec {vi} F1"),
        );
        let f2_last = run_and_check(
            &hcs[1],
            &vc["F2_dual_kernels"],
            &s4["F2_dual"]["chain_codomains"],
            &format!("vec {vi} F2_dual"),
        );

        // The headline invariant: middle-codomain match with SELF-DERIVED kernels.
        assert!(
            middle_codomain_matches(&f1_last, &f2_last),
            "vec {vi}: middle-codomain match failed with self-derived kernels"
        );
        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "strategy loop derived all per-step kernels + codomains (== oracle) and the \
         middle-codomain match held for all {n} vectors (self-derived kernels)"
    );
}

/// Inter-step integrity: tampering the derived-chain input (the post-gluing
/// kernel basis) must break the middle-codomain match. This is the integrity
/// the Phase 5 skeleton (oracle per-step kernels) could not test: with the
/// kernels *derived* from a single basis, corrupting that basis propagates
/// through every step and fails the match (or makes a step uncomputable).
#[test]
fn tampered_basis_breaks_middle_match() {
    let sv = load(STRATEGY_VECTORS);
    let cv = load(CHAIN_VECTORS);
    let main = load(PHASE0_VECTORS);
    let mut n = 0;
    for (vs, (vc, vm)) in sv["vectors"].as_array().unwrap().iter().zip(
        cv["vectors"]
            .as_array()
            .unwrap()
            .iter()
            .zip(main["test_vectors"].as_array().unwrap().iter()),
    ) {
        let hcs = vs["half_chains"].as_array().unwrap();
        let s4 = &vm["stage4_compute_hd"];

        // Honest F2_dual last codomain (untampered).
        let f2_last = run_and_check(
            &hcs[1],
            &vc["F2_dual_kernels"],
            &s4["F2_dual"]["chain_codomains"],
            "tamper-baseline F2_dual",
        );

        // Tamper one coordinate of F1's post-gluing basis, then derive F1.
        let start = ThetaStructureDim4::new(parse_coords(&hcs[0]["glue_codomain_null"]));
        let mut basis = parse_basis(&hcs[0]["post_glue_basis"]);
        let mut coords: [sqisign_verify::Fp2<sqisign_verify::Level1>; 16] =
            core::array::from_fn(|i| basis[0].coords()[i].clone());
        coords[0] = coords[0].add(&sqisign_verify::Fp2::<sqisign_verify::Level1>::one());
        basis[0] = sqisign_verify::hd::ThetaPointDim4::new(coords);
        let l = hcs[0]["n_plain"].as_u64().unwrap() as usize;

        let bad_matches = match run_strategy_chain(&start, &basis, l) {
            Some(sc) => match sc.last_codomain() {
                Some(f1_last) => middle_codomain_matches(f1_last, &f2_last),
                None => false,
            },
            None => false, // a step became uncomputable - also a rejection
        };
        assert!(!bad_matches, "vec {}: tampered basis still matched", vm["index"]);
        n += 1;
    }
    assert_eq!(n, 5);
    println!("tampering the derived-chain input breaks the middle-codomain match for all {n} vectors");
}

#[test]
fn strategy_chain_timing() {
    let sv = load(STRATEGY_VECTORS);
    // Pre-parse the starting bases/structures so timing covers the derivation only.
    struct HC {
        start: ThetaStructureDim4<sqisign_verify::Level1>,
        basis: [Pt; 4],
        l: usize,
    }
    let mut loaded: Vec<[HC; 2]> = Vec::new();
    for vs in sv["vectors"].as_array().unwrap() {
        let hcs = vs["half_chains"].as_array().unwrap();
        let mk = |hc: &Value| HC {
            start: ThetaStructureDim4::new(parse_coords(&hc["glue_codomain_null"])),
            basis: parse_basis(&hc["post_glue_basis"]),
            l: hc["n_plain"].as_u64().unwrap() as usize,
        };
        loaded.push([mk(&hcs[0]), mk(&hcs[1])]);
    }

    let t0 = Instant::now();
    for v in &loaded {
        for hc in v {
            let sc = run_strategy_chain(&hc.start, &hc.basis, hc.l).expect("chain");
            black_box(sc.last_codomain().unwrap().coords()[0].clone());
        }
    }
    let per_sig_ms = t0.elapsed().as_secs_f64() / loaded.len() as f64 * 1e3;

    println!("\n=========== PHASE 5b.6 STRATEGY-LOOP TIMING (Level 1, unoptimized) ===========");
    println!("per-signature kernel derivation (both half-chains, strategy walk +");
    println!("  ~938 doublings + ~1882 pushforwards + the plain (2,2,2,2) steps): {per_sig_ms:.1} ms");
    println!("NOTE: this is the chain-derivation cost from the gluing output. The full");
    println!("      self-contained verify additionally needs stages 1-3 (done) and the");
    println!("      gluing-chain self-derivation (the remaining 5b.6 piece; see NOTES).");
    println!("==============================================================================\n");
    assert!(per_sig_ms > 0.0);
}
