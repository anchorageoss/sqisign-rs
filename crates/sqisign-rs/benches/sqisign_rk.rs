//! Criterion benchmarks for the SQIsign-RK key-randomization API.
//!
//! Benchmarked at Level 1 (the NIST target). Uses only the public API:
//! `rand_pk`, `rand_sk`, `ver_key` (rng-free) plus the base keygen/sign/verify.
//! Run with `cargo bench --features sqisign-rk`.

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use sqisign_rs::keygen::keypair;
use sqisign_rs::sign::{generate_compact, sign};
use sqisign_rs::sqisign_rk::{
    rand_pk, rand_pk_compact, rand_sk, rand_sk_compact, ver_key, ver_key_compact,
};
use sqisign_rs::{Level1, Verifier};

type L = Level1;

fn bench_rand_pk(c: &mut Criterion) {
    let mut group = c.benchmark_group("rand_pk");
    let mut rng = rand::thread_rng();
    let (pk, _sk) = keypair::<L>(&mut rng);
    let rr = b"bench-randomness-0";

    group.bench_function("derive_public_key", |b| {
        b.iter(|| black_box(rand_pk(black_box(&pk), black_box(rr))))
    });
    group.finish();
}

fn bench_rand_sk(c: &mut Criterion) {
    let mut group = c.benchmark_group("rand_sk");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let rr = b"bench-randomness-0";

    group.bench_function("derive_secret_key", |b| {
        b.iter(|| black_box(rand_sk(black_box(&sk), black_box(&pk), black_box(rr))))
    });
    group.finish();
}

fn bench_ver_key(c: &mut Criterion) {
    let mut group = c.benchmark_group("ver_key");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);

    group.bench_function("verify_keypair", |b| {
        b.iter(|| black_box(ver_key(black_box(&pk), black_box(&sk))))
    });
    group.finish();
}

fn bench_derived_sign_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("derived_sign_verify");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let rr = b"bench-randomness-0";
    let pk_child = rand_pk(&pk, rr);
    let sk_child = rand_sk(&sk, &pk, rr);
    let msg = b"benchmark message";

    group.bench_function("sign_derived", |b| {
        b.iter(|| black_box(sign::<L>(&sk_child, &pk_child, msg, &mut rng).expect("sign")))
    });

    let sig = sign::<L>(&sk_child, &pk_child, msg, &mut rng).expect("sign");
    group.bench_function("verify_derived", |b| {
        b.iter(|| black_box(pk_child.verify(msg, &sig).is_ok()))
    });
    group.finish();
}

fn bench_full_derivation_round_trip(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_round_trip");
    group.sample_size(10); // RandSK is slow, reduce samples.
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let msg = b"round trip message";

    group.bench_function("derive_and_sign", |b| {
        let mut counter: u32 = 0;
        b.iter(|| {
            let rr = counter.to_be_bytes();
            counter = counter.wrapping_add(1);
            let pk_child = rand_pk(&pk, &rr);
            let sk_child = rand_sk(&sk, &pk, &rr);
            black_box(sign::<L>(&sk_child, &pk_child, msg, &mut rng).expect("sign"))
        })
    });
    group.finish();
}

fn bench_sequential_derivation(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequential_derivation");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);

    // Nested public derivation: pk -> pk' -> pk''.
    group.bench_function("rand_pk_depth_2", |b| {
        b.iter(|| {
            let child = rand_pk(&pk, b"level-1");
            black_box(rand_pk(&child, b"level-2"))
        })
    });

    // Nested secret derivation: sk -> sk' -> sk''.
    group.bench_function("rand_sk_depth_2", |b| {
        b.iter(|| {
            let child_pk = rand_pk(&pk, b"level-1");
            let child_sk = rand_sk(&sk, &pk, b"level-1");
            let grandchild_pk = rand_pk(&child_pk, b"level-2");
            black_box(rand_sk(&child_sk, &child_pk, b"level-2"));
            black_box(grandchild_pk)
        })
    });
    group.finish();
}

// ---- Compact (dim-4, 108-byte) variants, Level 1 ----

fn bench_rand_pk_compact(c: &mut Criterion) {
    let mut group = c.benchmark_group("rand_pk_compact");
    let mut rng = rand::thread_rng();
    let (pk, _sk) = generate_compact(&mut rng);
    let rr = b"bench-randomness-0";

    group.bench_function("derive_public_key", |b| {
        b.iter(|| black_box(rand_pk_compact(black_box(&pk), black_box(rr))))
    });
    group.finish();
}

fn bench_rand_sk_compact(c: &mut Criterion) {
    let mut group = c.benchmark_group("rand_sk_compact");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = generate_compact(&mut rng);
    let rr = b"bench-randomness-0";

    group.bench_function("derive_secret_key", |b| {
        b.iter(|| {
            black_box(rand_sk_compact(
                black_box(&sk),
                black_box(&pk),
                black_box(rr),
            ))
        })
    });
    group.finish();
}

fn bench_ver_key_compact(c: &mut Criterion) {
    let mut group = c.benchmark_group("ver_key_compact");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = generate_compact(&mut rng);

    group.bench_function("verify_keypair", |b| {
        b.iter(|| black_box(ver_key_compact(black_box(&pk), black_box(&sk))))
    });
    group.finish();
}

fn bench_compact_derived_sign_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("compact_derived_sign_verify");
    group.sample_size(10);
    let mut rng = rand::thread_rng();
    let (pk, sk) = generate_compact(&mut rng);
    let rr = b"bench-randomness-0";
    let pk_child = rand_pk_compact(&pk, rr);
    let sk_child = rand_sk_compact(&sk, &pk, rr);
    let msg = b"benchmark message";

    group.bench_function("sign_derived", |b| {
        b.iter(|| black_box(sk_child.sign(msg, &mut rng).expect("compact sign")))
    });

    let sig = sk_child.sign(msg, &mut rng).expect("compact sign");
    group.bench_function("verify_derived", |b| {
        b.iter(|| black_box(pk_child.verify(msg, &sig).is_ok()))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_rand_pk,
    bench_rand_sk,
    bench_ver_key,
    bench_derived_sign_verify,
    bench_full_derivation_round_trip,
    bench_sequential_derivation,
    bench_rand_pk_compact,
    bench_rand_sk_compact,
    bench_ver_key_compact,
    bench_compact_derived_sign_verify,
);
criterion_main!(benches);
