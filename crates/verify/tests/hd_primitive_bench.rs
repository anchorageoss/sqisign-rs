//! Phase 8a - microbenchmark of the dim-4 strategy-loop primitives.
//!
//! The optimal-strategy generator ([`sqisign_verify::hd`] `run_strategy_chain`)
//! balances point doublings against isogeny image-evaluations using a cost
//! ratio `mul_c = (doubling cost) / (image cost)`. This measures the
//! post-Phase-7 costs of one doubling, one `from_kernel` (codomain
//! reconstruction), and one `image`, and prints the ratio that should feed the
//! strategy generator. Run with:
//!
//! ```text
//! cargo test -p sqisign-verify --release --test hd_primitive_bench -- --nocapture
//! ```

mod hd_common;
use hd_common::{load, parse_coords, Pt};

use serde_json::Value;
use sqisign_verify::hd::{IsogenyDim4, ThetaStructureDim4};
use std::hint::black_box;
use std::time::Instant;

const ISOGENY_STEP_VECTORS: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/isogeny_step_vectors.json");

fn parse_k8(v: &Value) -> [Pt; 4] {
    let arr = v.as_array().expect("K_8 array");
    core::array::from_fn(|k| parse_coords(&arr[k]))
}

fn time_ns<F: FnMut()>(iters: u32, mut f: F) -> f64 {
    // warm-up
    for _ in 0..(iters / 10).max(1) {
        f();
    }
    let t = Instant::now();
    for _ in 0..iters {
        f();
    }
    t.elapsed().as_secs_f64() / iters as f64 * 1e9
}

#[test]
fn primitive_cost_microbench() {
    let doc = load(ISOGENY_STEP_VECTORS);
    let cases = doc["cases"].as_array().expect("cases");
    // A generic (non-zero-dual) step: the generic `image` formula applies.
    let case = cases
        .iter()
        .find(|c| !c["has_zero_dual"].as_bool().unwrap())
        .expect("a generic (non-zero-dual) case");

    let k8 = parse_k8(&case["K_8"]);
    let mut domain = ThetaStructureDim4::new(parse_coords(&case["domain_null"]));
    assert!(domain.precompute(), "domain null must be suitable");
    let iso = IsogenyDim4::from_kernel(&k8).expect("codomain computable");
    let input = parse_coords(&case["image_pairs"][0]["in"]);

    // One doubling on the domain structure (what the strategy descends with).
    let dbl_ns = time_ns(5000, || {
        black_box(domain.double_iter(black_box(&k8[0]), 1));
    });
    // One codomain reconstruction (once per step, strategy-independent).
    let fk_ns = time_ns(2000, || {
        black_box(IsogenyDim4::from_kernel(black_box(&k8)));
    });
    // One isogeny image-evaluation (the strategy's per-point pushforward).
    let img_ns = time_ns(5000, || {
        black_box(iso.image(black_box(&input)));
    });

    let mul_c = dbl_ns / img_ns;
    println!("\n========= PHASE 8a PRIMITIVE COSTS (Level 1, --release) =========");
    println!("  doubling      : {dbl_ns:8.0} ns");
    println!("  from_kernel   : {fk_ns:8.0} ns");
    println!("  image         : {img_ns:8.0} ns");
    println!("  mul_c = doubling/image = {mul_c:.3}   (strategy currently uses 1.000)");
    println!("  (from_kernel/image = {:.3}, fixed per step, not in the trade-off)", fk_ns / img_ns);
    println!("=================================================================\n");
    assert!(dbl_ns > 0.0 && img_ns > 0.0 && fk_ns > 0.0);

    // 8a strategy comparison: does mul_c=1.605 change the strategy?
    // Model the per-half-chain cost (doublings * dbl + images * img) of the
    // current (mul_c=1.0) strategy vs the measured-ratio one, over the chain
    // lengths a Level-1 half-chain actually uses (n_plain = 67 - v2(a2)).
    println!("===== 8a STRATEGY COMPARISON (cost via measured primitives) =====");
    let mut any_diff = false;
    for l in 60..=67usize {
        let s0 = sqisign_verify::hd::optimised_strategy(l, 1.0);
        let s1 = sqisign_verify::hd::optimised_strategy(l, mul_c);
        let (d0, i0) = walk_counts(&s0, l);
        let (d1, i1) = walk_counts(&s1, l);
        let c0 = d0 as f64 * dbl_ns + i0 as f64 * img_ns;
        let c1 = d1 as f64 * dbl_ns + i1 as f64 * img_ns;
        let diff = s0 != s1;
        any_diff |= diff;
        // ×4 points/basis, ×2 half-chains => per-verify modeled delta.
        let per_verify_us = (c0 - c1) * 4.0 * 2.0 / 1000.0;
        println!(
            "  l={l}: changed={diff}  dbls {d0}->{d1}  imgs {i0}->{i1}  modeled Δ/verify = {per_verify_us:+.1} µs"
        );
    }
    println!("  strategy differs for some l: {any_diff}");
    println!("=================================================================\n");
}

/// Count per-basis doublings and image-evaluations a strategy walk performs for
/// an `l`-step chain (mirrors `run_strategy_chain`'s bookkeeping; the per-step
/// `from_kernel` is fixed and excluded from the doubling/image trade-off).
fn walk_counts(strategy: &[u32], l: usize) -> (u64, u64) {
    let mut strat_idx = 0usize;
    let mut level: Vec<u32> = vec![0];
    let mut stack_depth: u64 = 1;
    let (mut doublings, mut images) = (0u64, 0u64);
    for k in 0..l {
        let target = (l - 1 - k) as u32;
        let mut prev: u32 = level.iter().sum();
        while prev != target {
            let s = strategy[strat_idx];
            level.push(s);
            prev += s;
            doublings += s as u64;
            stack_depth += 1;
            strat_idx += 1;
        }
        level.pop();
        stack_depth -= 1; // pop the order-8 kernel of this step
        images += stack_depth; // push the remaining saved bases through
    }
    (doublings, images)
}

