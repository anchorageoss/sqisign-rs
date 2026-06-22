//! Phase 3: drive the plain-step half-chains F1 and F2_dual, validate every
//! per-step codomain against the oracle, confirm the middle-codomain match for
//! all 5 vectors, and report the headline chain timing.
//!
//! Ground truth: `tests/chain_vectors.json` (ordered per-step kernels + the
//! reference operation counts, from `sqisignhd-harness/extract_chain.py`) and
//! the Phase 0 `test_vectors_l1.json` (per-step codomains + the recorded
//! middle-codomain match).

mod hd_common;
use hd_common::{load, parse_coords, parse_node, zero_point, Pt, PHASE0_VECTORS};

use serde_json::Value;
use std::hint::black_box;
use std::time::Instant;
use sqisign_verify::hd::{
    hadamard, middle_codomain_matches, run_half_chain, run_half_chain_collect, IsogenyDim4,
    ThetaPointDim4, ThetaStructureDim4,
};

const CHAIN_VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/chain_vectors.json");

fn parse_k8(v: &Value) -> [Pt; 4] {
    let arr = v.as_array().expect("K_8 array");
    assert_eq!(arr.len(), 4, "K_8 must have 4 points");
    core::array::from_fn(|k| parse_coords(&arr[k]))
}

fn parse_kernels(v: &Value) -> Vec<[Pt; 4]> {
    v.as_array().expect("kernels array").iter().map(parse_k8).collect()
}

/// Deliverable 2: both half-chains reproduce the oracle per-step codomains,
/// each codomain feeding the next.
#[test]
fn half_chains_match_oracle_per_step() {
    let chain_doc = load(CHAIN_VECTORS);
    let main = load(PHASE0_VECTORS);

    let mut total_steps = 0usize;
    for vec in chain_doc["vectors"].as_array().unwrap() {
        let vi = vec["index"].as_u64().unwrap() as usize;
        for (chain, key) in [("F1", "F1_kernels"), ("F2_dual", "F2_dual_kernels")] {
            let kernels = parse_kernels(&vec[key]);
            let mut buf = vec![zero_point(); kernels.len()];
            let n = run_half_chain_collect(&kernels, &mut buf).expect("chain computable");
            assert_eq!(n, kernels.len());

            // chain_codomains index 0 is the gluing codomain; plain step i sits
            // at index i + 1.
            let cods = main["test_vectors"][vi]["stage4_compute_hd"][chain]["chain_codomains"]
                .as_array()
                .unwrap();
            for (i, cod) in buf.iter().enumerate().take(n) {
                let expected = parse_node(&cods[i + 1]);
                assert!(
                    cod.projective_eq(&expected),
                    "vector {vi} {chain} plain step {} codomain mismatch",
                    i + 1
                );
            }
            total_steps += n;
        }
    }
    println!("validated {total_steps} plain-step codomains across 5 vectors (both half-chains)");
}

/// Deliverable 3: the middle-codomain check returns the oracle's accept result
/// for all 5 vectors, computed from the chained codomains.
#[test]
fn middle_codomain_matches_for_all_vectors() {
    let chain_doc = load(CHAIN_VECTORS);
    let main = load(PHASE0_VECTORS);

    let mut n = 0;
    for vec in chain_doc["vectors"].as_array().unwrap() {
        let vi = vec["index"].as_u64().unwrap() as usize;
        let f1 = parse_kernels(&vec["F1_kernels"]);
        let f2 = parse_kernels(&vec["F2_dual_kernels"]);

        let f1_last = run_half_chain(&f1).expect("F1 chain");
        let f2_last = run_half_chain(&f2).expect("F2_dual chain");

        // The heart of verification.
        assert!(
            middle_codomain_matches(&f1_last, &f2_last),
            "vector {vi}: middle-codomain check failed"
        );

        // Matches the oracle's recorded accept result.
        let s5 = &main["test_vectors"][vi]["stage5_codomain_check"];
        assert_eq!(s5["match"].as_bool(), Some(true), "oracle match must be accept");

        // Tie the computed codomains to the oracle's recorded C1/HC2.
        let c1 = parse_node(&s5["C1_zero"]);
        assert!(f1_last.projective_eq(&c1), "vector {vi}: F1 last != recorded C1.zero()");
        let hc2 = parse_node(&s5["HC2_zero"]);
        let my_hc2 = ThetaPointDim4::new(hadamard(f2_last.coords()));
        assert!(
            my_hc2.projective_eq(&hc2),
            "vector {vi}: Hadamard(F2_dual last) != recorded HC2.zero()"
        );
        n += 1;
    }
    assert_eq!(n, 5, "expected all 5 vectors");
    println!("middle-codomain match confirmed (accept) for all {n} vectors");
}

/// Deliverable 4: the headline timing report.
#[test]
fn chain_timing_report() {
    // Parse all kernels OUTSIDE the timed region.
    let chain_doc = load(CHAIN_VECTORS);
    struct Sig {
        f1: Vec<[Pt; 4]>,
        f2: Vec<[Pt; 4]>,
        doublings: u64,
        images: u64,
        steps: usize,
    }
    let parsed: Vec<Sig> = chain_doc["vectors"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            let f1 = parse_kernels(&v["F1_kernels"]);
            let f2 = parse_kernels(&v["F2_dual_kernels"]);
            let steps = f1.len() + f2.len();
            Sig {
                f1,
                f2,
                doublings: v["ref_dim4_doublings"].as_u64().unwrap(),
                images: v["ref_plain_image_calls"].as_u64().unwrap(),
                steps,
            }
        })
        .collect();

    // Headline: codomain-reconstruction chain time per signature (both halves).
    let iters = 5usize;
    let t0 = Instant::now();
    for _ in 0..iters {
        for s in &parsed {
            let a = run_half_chain(&s.f1).unwrap();
            let b = run_half_chain(&s.f2).unwrap();
            black_box((a, b));
        }
    }
    let elapsed = t0.elapsed();
    let n_sig = (iters * parsed.len()) as f64;
    let per_sig_ms = elapsed.as_secs_f64() / n_sig * 1e3;
    let avg_steps = parsed.iter().map(|s| s.steps).sum::<usize>() as f64 / parsed.len() as f64;

    // Micro-bench the primitives for an honest full-chain estimate.
    let sample_k = &parsed[0].f1[0];
    let bn = 300usize;
    let tk = Instant::now();
    for _ in 0..bn {
        black_box(IsogenyDim4::from_kernel(sample_k));
    }
    let t_from_kernel = tk.elapsed().as_secs_f64() / bn as f64;

    let iso = IsogenyDim4::from_kernel(sample_k).unwrap();
    let sample_pt = &parsed[0].f1[0][0];
    let ti = Instant::now();
    for _ in 0..bn {
        black_box(iso.image(sample_pt));
    }
    let t_image = ti.elapsed().as_secs_f64() / bn as f64;

    // Doubling cost is a fixed-op-count formula, independent of the point/struct
    // values, so any suitable structure gives a representative timing.
    let mut st = ThetaStructureDim4::new(iso.codomain_null().clone());
    let t_double = if st.precompute() {
        let td = Instant::now();
        for _ in 0..bn {
            black_box(st.double(sample_pt));
        }
        td.elapsed().as_secs_f64() / bn as f64
    } else {
        0.0
    };

    let avg_doublings = parsed.iter().map(|s| s.doublings).sum::<u64>() as f64 / parsed.len() as f64;
    let avg_images = parsed.iter().map(|s| s.images).sum::<u64>() as f64 / parsed.len() as f64;
    let est_full_ms =
        (avg_steps * t_from_kernel + avg_doublings * t_double + avg_images * t_image) * 1e3;

    println!("\n========= PHASE 3 CHAIN TIMING (Level 1, unoptimized, this machine) =========");
    println!("plain steps per signature (both half-chains): {avg_steps:.0}");
    println!("HEADLINE -- codomain-reconstruction chain time: {per_sig_ms:.3} ms/signature");
    println!("            (= {avg_steps:.0} from_kernel calls; EXCLUDES the strategy");
    println!("             doublings + pushforwards that derive the kernels, and the");
    println!("             gluing step -- both deferred)");
    println!(
        "measured primitives: from_kernel {:.1} us | image {:.1} us | double {:.1} us",
        t_from_kernel * 1e6,
        t_image * 1e6,
        t_double * 1e6
    );
    println!("reference per-signature counts: ~{avg_doublings:.0} dim-4 doublings, ~{avg_images:.0} plain image-evals");
    println!("MODELED full plain-chain time (codomain + doublings + pushforwards):");
    println!("            ~{est_full_ms:.1} ms/signature");
    println!("ORDER OF MAGNITUDE: TENS of milliseconds per signature (plain steps; excludes gluing).");
    println!("=============================================================================\n");

    assert!(per_sig_ms > 0.0);
    assert!(est_full_ms > per_sig_ms);
}
