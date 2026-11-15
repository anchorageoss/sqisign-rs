//! mike-step <out.csv> [runs]
//!
//! Times MIKE-I's key generation and shared secret (each computes one chain
//! of e = 228 dimension-4 2-isogenies from a product of elliptic curves),
//! `runs` times each (default 200) after five warm-up calls, with `rdtsc`
//! and the wall clock around the call, and writes one CSV line per call:
//! operation, iteration, cycles, ns. Pin it to a core with taskset.

use mike_rs::mike::MIKE_I;
use std::io::Write;
use std::time::Instant;

fn tsc() -> u64 {
    // SAFETY: rdtsc is available on every x86-64 processor.
    unsafe { core::arch::x86_64::_rdtsc() }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut out = std::fs::File::create(&argv[1]).expect("out.csv");
    let runs: usize = argv.get(2).map(|s| s.parse().unwrap()).unwrap_or(200);
    writeln!(
        out,
        "# MIKE-I (mike_rs 0d79121, e = 228 steps): rdtsc cycles and nanoseconds per call"
    )
    .unwrap();
    writeln!(out, "operation,iteration,cycles,ns").unwrap();
    let (_, sk_a) = MIKE_I.keygen_seeded(&[0u8; 48]).unwrap();
    let (pk_b, _) = MIKE_I.keygen_seeded(&[1u8; 48]).unwrap();
    for i in 0..runs + 5 {
        let mut seed = [0u8; 48];
        seed[..8].copy_from_slice(&(i as u64).to_le_bytes());
        let (t, c) = (Instant::now(), tsc());
        let kp = MIKE_I.keygen_seeded(&seed).unwrap();
        let (c, ns) = (tsc() - c, t.elapsed().as_nanos());
        std::hint::black_box(kp);
        if i >= 5 {
            writeln!(out, "keygen,{},{c},{ns}", i - 5).unwrap();
        }
    }
    for i in 0..runs + 5 {
        let (t, c) = (Instant::now(), tsc());
        let ss = sk_a.shared_secret(std::hint::black_box(&pk_b)).unwrap();
        let (c, ns) = (tsc() - c, t.elapsed().as_nanos());
        std::hint::black_box(ss);
        if i >= 5 {
            writeln!(out, "shared_secret,{},{c},{ns}", i - 5).unwrap();
        }
    }
}
