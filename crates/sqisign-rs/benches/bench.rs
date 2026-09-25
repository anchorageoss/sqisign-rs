//! Round-3 key generation, signing and verification at levels I, III and V
//! in `rdtsc` cycles (criterion means; the unit the C reference's own
//! benchmark reports), the default build. `SQISIGN_FORCE_PORTABLE=1` runs
//! the portable field kernels of the x86-64 backends instead.
//!
//!     cargo bench -p sqisign-rs --bench bench -- "level I$"

mod cycles;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cycles::Cycles;
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{keygen, level1, level3, level5, sign, Params};
use sqisign_verify::compressed::{
    compress, compressed_from_bytes, compressed_to_bytes, verify_compressed_prepared,
    MAX_COMPRESSED_BYTES,
};
use sqisign_verify::fp::FpBackend;
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::sqisign::{
    public_key_from_bytes, public_key_to_bytes, signature_from_bytes, signature_to_bytes, verify,
    verify_batch, verify_prepared, BatchItem, PreparedPublicKey, MAX_PUBLICKEY_BYTES,
    MAX_SIGNATURE_BYTES,
};
use sqisign_verify::{P324_3, P500_27, P664_17};
use std::time::Duration;

fn kernels_from_env() -> &'static str {
    if std::env::var("SQISIGN_FORCE_PORTABLE").is_ok_and(|v| v == "1") {
        sqisign_verify::fp::force_portable_arithmetic();
    }
    sqisign_verify::fp::arithmetic_backend()
}

fn level<L: FpBackend + PrimePrecomp, const N: usize>(
    c: &mut Criterion<Cycles>,
    params: &Params<N>,
    name: &str,
) {
    let vp = &params.verify;
    let mut g = c.benchmark_group(format!("SQIsign (Rust) {name}"));
    g.sample_size(20).measurement_time(Duration::from_secs(20));
    let mut r = ShakeRng::new(name.as_bytes(), b"bnc");
    g.bench_function(BenchmarkId::new("keygen", name), |b| {
        b.iter(|| keygen::<L, N>(params, &mut r).unwrap())
    });
    let (pk, sk) = keygen::<L, N>(params, &mut r).unwrap();
    let mut i = 0u64;
    g.bench_function(BenchmarkId::new("sign", name), |b| {
        b.iter(|| {
            i += 1;
            sign::<L, N>(params, &pk, &sk, &i.to_le_bytes(), &mut r).unwrap()
        })
    });
    let msg = b"criterion";
    let sig = sign::<L, N>(params, &pk, &sk, msg, &mut r).unwrap();
    let mut pkb = [0u8; MAX_PUBLICKEY_BYTES];
    let pklen = public_key_to_bytes(&pk, &mut pkb);
    let mut sigb = [0u8; MAX_SIGNATURE_BYTES];
    let siglen = signature_to_bytes(vp, &sig, &mut sigb).unwrap();
    g.sample_size(100).measurement_time(Duration::from_secs(10));
    // verification from bytes: decode the key and the signature, verify
    g.bench_function(
        BenchmarkId::new(format!("verify from bytes ({siglen} B)"), name),
        |b| {
            b.iter(|| {
                let pk = public_key_from_bytes::<L>(vp, &pkb[..pklen]).unwrap();
                let sig = signature_from_bytes::<L>(vp, &sigb[..siglen]).unwrap();
                assert!(verify(vp, &pk, &sig, msg))
            })
        },
    );
    // the key's canonical basis prepared once, the signature decoded
    let key = PreparedPublicKey::new(vp, &pk).unwrap();
    g.bench_function(
        BenchmarkId::new(
            format!("verify from bytes, prepared key ({siglen} B)"),
            name,
        ),
        |b| {
            b.iter(|| {
                let sig = signature_from_bytes::<L>(vp, &sigb[..siglen]).unwrap();
                assert!(verify_prepared(vp, &key, &sig, msg))
            })
        },
    );
    // throughput: a batch of 16 signatures under the prepared key, the
    // per-signature cost (criterion reports the batch; divide by 16)
    let batch: Vec<([u8; 8], Vec<u8>)> = (0..16u64)
        .map(|i| {
            let m = i.to_le_bytes();
            let s = sign::<L, N>(params, &pk, &sk, &m, &mut r).unwrap();
            let mut sb = [0u8; MAX_SIGNATURE_BYTES];
            let n = signature_to_bytes(vp, &s, &mut sb).unwrap();
            (m, sb[..n].to_vec())
        })
        .collect();
    let items: Vec<BatchItem<'_, L>> = batch
        .iter()
        .map(|(m, s)| BatchItem {
            key: &key,
            msg: &m[..],
            signature: &s[..],
        })
        .collect();
    g.bench_function(
        BenchmarkId::new("verify_batch of 16, prepared key", name),
        |b| {
            b.iter(|| {
                let mut out = [false; 16];
                assert_eq!(verify_batch(vp, &items, &mut out), 16)
            })
        },
    );
    // the compressed signature (three matrix entries; the fourth recovered
    // from two Weil pairings and a discrete logarithm)
    let c = compress(vp, &sig).unwrap();
    let mut cb = [0u8; MAX_COMPRESSED_BYTES];
    let clen = compressed_to_bytes(vp, &c, &mut cb).unwrap();
    g.bench_function(
        BenchmarkId::new(
            format!("verify compressed from bytes, prepared key ({clen} B)"),
            name,
        ),
        |b| {
            b.iter(|| {
                let c = compressed_from_bytes::<L>(vp, &cb[..clen]).unwrap();
                assert!(verify_compressed_prepared(vp, &key, &c, msg))
            })
        },
    );
    g.finish();
}

/// Dependent chains of field operations, the latency of one operation at
/// each level (comparable with prism-rs's `field` bench and the C harness).
fn field<L: FpBackend>(c: &mut Criterion<Cycles>, label: &str, enc: usize) {
    use sqisign_verify::fp::{Fp, Fp2};
    use sqisign_verify::rng::Rng;
    let mut r = ShakeRng::new(b"field bench", b"fld");
    let mut bytes = [0u8; 4 * 84];
    assert!(r.fill(&mut bytes[..4 * enc]));
    let a = Fp2::<L> {
        re: Fp::<L>::decode_reduce(&bytes[..enc]),
        im: Fp::<L>::decode_reduce(&bytes[enc..2 * enc]),
    };
    let b = Fp2::<L> {
        re: Fp::<L>::decode_reduce(&bytes[2 * enc..3 * enc]),
        im: Fp::<L>::decode_reduce(&bytes[3 * enc..4 * enc]),
    };
    let mut g = c.benchmark_group(label);
    g.sample_size(200);
    g.bench_function("Fp mul, chain of 1000", |bn| {
        bn.iter(|| {
            let mut x = criterion::black_box(&a).re.clone();
            for _ in 0..1000 {
                x = x.mul(&b.re);
            }
            x
        })
    });
    g.bench_function("Fp2 mul, chain of 1000", |bn| {
        bn.iter(|| {
            let mut x = criterion::black_box(&a).clone();
            for _ in 0..1000 {
                x = x.mul(&b);
            }
            x
        })
    });
    g.bench_function("Fp2 sqr, chain of 1000", |bn| {
        bn.iter(|| {
            let mut x = criterion::black_box(&a).clone();
            for _ in 0..1000 {
                x = x.sqr();
            }
            x
        })
    });
    g.bench_function("Fp2 add, chain of 1000", |bn| {
        bn.iter(|| {
            let mut x = criterion::black_box(&a).clone();
            for _ in 0..1000 {
                x = x.add(&b);
            }
            x
        })
    });
    g.finish();
}

fn all(c: &mut Criterion<Cycles>) {
    eprintln!("field arithmetic: {}", kernels_from_env());
    field::<P324_3>(c, "field level I (p324_3)", 41);
    field::<P500_27>(c, "field level III (p500_27)", 64);
    field::<P664_17>(c, "field level V (p664_17)", 84);
    level::<P324_3, { sqisign_rs::precomp::p324_3::IBZ_NLIMBS }>(c, &level1(), "level I");
    level::<P500_27, { sqisign_rs::precomp::p500_27::IBZ_NLIMBS }>(c, &level3(), "level III");
    level::<P664_17, { sqisign_rs::precomp::p664_17::IBZ_NLIMBS }>(c, &level5(), "level V");
}

/// The compact (dimension-4) format at level I (feature `compact`): key
/// generation, signing (with the good-degree rejection loop) and
/// verification from bytes. No C reference exists at these parameters.
#[cfg(feature = "compact")]
fn compact_level1(c: &mut Criterion<Cycles>) {
    use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign};
    use sqisign_rs::sqisign::level1;
    use sqisign_verify::hd::hd_verify_bytes;
    const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
    let params = level1();
    let mut g = c.benchmark_group("SQIsign compact (Rust) level I");
    g.sample_size(20).measurement_time(Duration::from_secs(20));
    let mut r = ShakeRng::new(b"compact bench", b"bnc");
    g.bench_function(BenchmarkId::new("keygen", "level I"), |b| {
        b.iter(|| compact_keygen::<P324_3, N>(&params, &mut r).unwrap())
    });
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
    let mut pkb = [0u8; 84];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let mut i = 0u64;
    let mut sig = [0u8; 142];
    g.bench_function(BenchmarkId::new("sign", "level I"), |b| {
        b.iter(|| {
            i += 1;
            compact_sign::<P324_3, N>(&params, &pk, &sk, &i.to_le_bytes(), &mut r, &mut sig)
                .unwrap()
        })
    });
    let msg = b"criterion";
    compact_sign::<P324_3, N>(&params, &pk, &sk, msg, &mut r, &mut sig).unwrap();
    g.sample_size(30).measurement_time(Duration::from_secs(20));
    g.bench_function(
        BenchmarkId::new("verify from bytes (142 B)", "level I"),
        |b| b.iter(|| assert!(hd_verify_bytes::<P324_3>(&sig, &pkb, msg).is_ok())),
    );
    g.finish();
}

#[cfg(feature = "compact")]
criterion_group! {
    name = compact_benches;
    config = Criterion::default().with_measurement(Cycles);
    targets = compact_level1
}

criterion_group! {
    name = benches;
    config = Criterion::default().with_measurement(Cycles);
    targets = all
}
#[cfg(feature = "compact")]
criterion_main!(benches, compact_benches);
#[cfg(not(feature = "compact"))]
criterion_main!(benches);
