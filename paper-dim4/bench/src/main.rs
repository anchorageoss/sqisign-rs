//! One session for the dimension-4 paper: key generation, signing and
//! verification of SQIsign round 2 and round 3, dimension 2 and dimension 4,
//! measured by one binary with one method (`rdtsc` cycles and wall time
//! around each call, medians over a fixed number of runs, seeded
//! randomness). The round-3 code is the crate in this repository; the
//! round-2 code is the last round-2 release, `v0.4.31`, pulled from GitHub
//! as a dependency. The output is Markdown in the format paper-dim4/BENCH.md's session
//! section uses; the paper's figure scripts read that section.
//!
//! Also reported: the compact signer's rejection loop (samples and primality
//! tests per signature, from `compact_sign_stats`) and the peak heap of one
//! compact verification (a counting allocator around the call).

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::time::Instant;

// ---- a counting allocator, for the peak heap of one verification ---------

struct Counting;
static CUR: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            let c = CUR.fetch_add(layout.size(), Relaxed) + layout.size();
            PEAK.fetch_max(c, Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        CUR.fetch_sub(layout.size(), Relaxed);
        System.dealloc(ptr, layout)
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let p = System.realloc(ptr, layout, new_size);
        if !p.is_null() {
            if new_size >= layout.size() {
                let c = CUR.fetch_add(new_size - layout.size(), Relaxed) + new_size - layout.size();
                PEAK.fetch_max(c, Relaxed);
            } else {
                CUR.fetch_sub(layout.size() - new_size, Relaxed);
            }
        }
        p
    }
}

#[global_allocator]
static ALLOC: Counting = Counting;

fn heap_mark() {
    PEAK.store(CUR.load(Relaxed), Relaxed);
}
fn heap_peak_since_mark() -> usize {
    PEAK.load(Relaxed)
}
fn heap_now() -> usize {
    CUR.load(Relaxed)
}

fn vm_hwm_kib() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = s.lines().find(|l| l.starts_with("VmHWM:"))?;
    line.split_whitespace().nth(1)?.parse().ok()
}

// ---- timing ---------------------------------------------------------------

#[cfg(target_arch = "x86_64")]
fn cycles() -> u64 {
    // SAFETY: rdtsc is available on every x86-64 CPU.
    unsafe { core::arch::x86_64::_rdtsc() }
}
#[cfg(not(target_arch = "x86_64"))]
fn cycles() -> u64 {
    0
}

struct Stats {
    runs: usize,
    median_mcyc: f64,
    mean_mcyc: f64,
    min_mcyc: f64,
    max_mcyc: f64,
    median_ms: f64,
}

fn measure<F: FnMut()>(runs: usize, mut f: F) -> Stats {
    let mut cyc = Vec::with_capacity(runs);
    let mut ns = Vec::with_capacity(runs);
    for _ in 0..runs {
        let t0 = Instant::now();
        let c0 = cycles();
        f();
        let c1 = cycles();
        let t1 = t0.elapsed();
        cyc.push((c1 - c0) as f64 / 1e6);
        ns.push(t1.as_secs_f64() * 1e3);
    }
    let mut sc = cyc.clone();
    sc.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut sn = ns.clone();
    sn.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = |v: &Vec<f64>| {
        let n = v.len();
        if n % 2 == 1 {
            v[n / 2]
        } else {
            (v[n / 2 - 1] + v[n / 2]) / 2.0
        }
    };
    Stats {
        runs,
        median_mcyc: med(&sc),
        mean_mcyc: cyc.iter().sum::<f64>() / runs as f64,
        min_mcyc: sc[0],
        max_mcyc: sc[runs - 1],
        median_ms: med(&sn),
    }
}

struct Row {
    impl_: &'static str,
    round: u8,
    dim: u8,
    level: &'static str,
    op: &'static str,
    st: Stats,
}

struct Size {
    impl_: &'static str,
    round: u8,
    dim: u8,
    level: &'static str,
    pk: usize,
    sig: usize,
}

const RUNS_KEYGEN: usize = 50;
const RUNS_SIGN: usize = 50;
const RUNS_VERIFY: usize = 200;
const RUNS_COMPACT_SIGN: usize = 100;
const RUNS_COMPACT_VERIFY: usize = 100;
const MSG: &[u8] = b"dimension-4 paper session";

const R3: &str = "SQIsign (Rust, this crate)";
const R2: &str = "SQIsign (Rust, v0.4.31)";

// ---- round 3, dimension 2 (the current crate) -------------------------------

fn r3_level<L: sqisign_rs::SigningLevel>(
    level: &'static str,
    rng: &mut StdRng,
    rows: &mut Vec<Row>,
    sizes: &mut Vec<Size>,
) {
    use sqisign_rs::{generate, CompressedSignature, PublicKey, Signature};
    let st = measure(RUNS_KEYGEN, || {
        let (pk, sk) = generate::<L>(rng);
        black_box((pk, sk));
    });
    rows.push(Row { impl_: R3, round: 3, dim: 2, level, op: "keygen", st });
    let (pk, sk) = generate::<L>(rng);
    let mut i = 0u64;
    let st = measure(RUNS_SIGN, || {
        i += 1;
        let sig = sk.sign(&i.to_le_bytes(), rng).unwrap();
        black_box(sig);
    });
    rows.push(Row { impl_: R3, round: 3, dim: 2, level, op: "sign", st });
    let sig = sk.sign(MSG, rng).unwrap();
    let pkb = pk.to_bytes();
    let sigb = sig.to_bytes();
    sizes.push(Size { impl_: R3, round: 3, dim: 2, level, pk: pkb.len(), sig: sigb.len() });
    let st = measure(RUNS_VERIFY, || {
        let pk = PublicKey::<L>::from_bytes(&pkb).unwrap();
        let sig = Signature::<L>::from_bytes(&sigb).unwrap();
        assert!(pk.verify(MSG, &sig).is_ok());
    });
    rows.push(Row { impl_: R3, round: 3, dim: 2, level, op: "verify from bytes", st });
    let cb = sig.compress().to_bytes();
    let st = measure(RUNS_VERIFY, || {
        let pk = PublicKey::<L>::from_bytes(&pkb).unwrap();
        let c = CompressedSignature::<L>::from_bytes(&cb).unwrap();
        assert!(pk.verify_compressed(MSG, &c).is_ok());
    });
    rows.push(Row { impl_: R3, round: 3, dim: 2, level, op: "verify compressed from bytes", st });
    sizes.push(Size { impl_: "SQIsign compressed (Rust, this crate)", round: 3, dim: 2, level, pk: pkb.len(), sig: cb.len() });
}

// ---- round 3, dimension 4 (the compact format, level I) ---------------------

struct Loop {
    signatures: usize,
    samples_mean: f64,
    samples_median: f64,
    samples_min: u32,
    samples_max: u32,
    tests_mean: f64,
}

struct Mem {
    peak_heap_bytes: usize,
    heap_before: usize,
    vm_hwm_kib: Option<u64>,
}

fn r3_compact(rows: &mut Vec<Row>, sizes: &mut Vec<Size>) -> (Loop, Mem) {
    use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign_stats};
    use sqisign_rs::mp::ShakeRng;
    use sqisign_rs::sqisign::level1;
    use sqisign_rs::{CompactPublicKey, CompactSignature, Level1 as P324_3};
    const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
    let params = level1();
    let mut r = ShakeRng::new(b"dimension-4 paper session", b"bnc");
    let level = "I";
    let st = measure(RUNS_KEYGEN, || {
        let kp = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
        black_box(kp);
    });
    rows.push(Row { impl_: R3, round: 3, dim: 4, level, op: "keygen", st });
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
    let mut pkb = [0u8; 84];
    let pklen = compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let mut sig = [0u8; 142];
    let mut i = 0u64;
    let mut samples: Vec<u32> = Vec::new();
    let mut tests: Vec<u32> = Vec::new();
    let st = measure(RUNS_COMPACT_SIGN, || {
        i += 1;
        let (_, s) =
            compact_sign_stats::<P324_3, N>(&params, &pk, &sk, &i.to_le_bytes(), &mut r, &mut sig).unwrap();
        samples.push(s.response_samples);
        tests.push(s.primality_tests);
    });
    rows.push(Row { impl_: R3, round: 3, dim: 4, level, op: "sign", st });
    let (siglen, _) = compact_sign_stats::<P324_3, N>(&params, &pk, &sk, MSG, &mut r, &mut sig).unwrap();
    sizes.push(Size { impl_: R3, round: 3, dim: 4, level, pk: pklen, sig: siglen });
    let pkb = &pkb[..pklen];
    let sigb = &sig[..siglen];
    let st = measure(RUNS_COMPACT_VERIFY, || {
        let pk = CompactPublicKey::<P324_3>::from_bytes(pkb).unwrap();
        let s = CompactSignature::<P324_3>::from_bytes(sigb).unwrap();
        assert!(pk.verify_compact(MSG, &s).is_ok());
    });
    rows.push(Row { impl_: R3, round: 3, dim: 4, level, op: "verify from bytes", st });
    // peak heap of one verification
    let heap_before = heap_now();
    heap_mark();
    {
        let pk = CompactPublicKey::<P324_3>::from_bytes(pkb).unwrap();
        let s = CompactSignature::<P324_3>::from_bytes(sigb).unwrap();
        assert!(pk.verify_compact(MSG, &s).is_ok());
    }
    let peak = heap_peak_since_mark() - heap_before;
    let mem = Mem { peak_heap_bytes: peak, heap_before, vm_hwm_kib: vm_hwm_kib() };
    let mut ss = samples.clone();
    ss.sort();
    let n = ss.len();
    let lp = Loop {
        signatures: n,
        samples_mean: samples.iter().map(|&x| x as f64).sum::<f64>() / n as f64,
        samples_median: if n % 2 == 1 { ss[n / 2] as f64 } else { (ss[n / 2 - 1] + ss[n / 2]) as f64 / 2.0 },
        samples_min: ss[0],
        samples_max: ss[n - 1],
        tests_mean: tests.iter().map(|&x| x as f64).sum::<f64>() / n as f64,
    };
    (lp, mem)
}

// ---- round 2 (v0.4.31): dimension 2 at three levels, dimension 4 at level I --

macro_rules! r2_level {
    ($lvl:ty, $name:expr, $rng:expr, $rows:expr, $sizes:expr) => {{
        use sqisign_rs_r2::{generate, PublicKey, Signature, Verifier};
        let level: &'static str = $name;
        let rng: &mut StdRng = $rng;
        let st = measure(RUNS_KEYGEN, || {
            let (pk, sk) = generate::<$lvl>(rng);
            black_box((pk, sk));
        });
        $rows.push(Row { impl_: R2, round: 2, dim: 2, level, op: "keygen", st });
        let (pk, sk) = generate::<$lvl>(rng);
        let mut i = 0u64;
        let st = measure(RUNS_SIGN, || {
            i += 1;
            let sig = sk.sign(&i.to_le_bytes(), rng).unwrap();
            black_box(sig);
        });
        $rows.push(Row { impl_: R2, round: 2, dim: 2, level, op: "sign", st });
        let sig = sk.sign(MSG, rng).unwrap();
        let pkb = pk.to_bytes();
        let sigb = sig.to_bytes();
        $sizes.push(Size { impl_: R2, round: 2, dim: 2, level, pk: pkb.len(), sig: sigb.len() });
        let st = measure(RUNS_VERIFY, || {
            let pk = PublicKey::<$lvl>::from_bytes(&pkb).unwrap();
            let sig = Signature::<$lvl>::from_bytes(&sigb).unwrap();
            assert!(pk.verify(MSG, &sig).is_ok());
        });
        $rows.push(Row { impl_: R2, round: 2, dim: 2, level, op: "verify from bytes", st });
    }};
}

fn r2_compact(rng: &mut StdRng, rows: &mut Vec<Row>, sizes: &mut Vec<Size>) {
    use sqisign_rs_r2::{generate_compact, CompactPublicKey, CompactSignature, Level1, Verifier};
    let level = "I";
    let st = measure(RUNS_KEYGEN, || {
        let kp = generate_compact(rng);
        black_box(kp);
    });
    rows.push(Row { impl_: R2, round: 2, dim: 4, level, op: "keygen", st });
    let (pk, sk) = generate_compact(rng);
    let mut i = 0u64;
    let st = measure(RUNS_COMPACT_SIGN, || {
        i += 1;
        let sig = sk.sign(&i.to_le_bytes(), rng).unwrap();
        black_box(sig);
    });
    rows.push(Row { impl_: R2, round: 2, dim: 4, level, op: "sign", st });
    let sig = sk.sign(MSG, rng).unwrap();
    let pkb = pk.to_bytes();
    let sigb = sig.to_bytes();
    sizes.push(Size { impl_: R2, round: 2, dim: 4, level, pk: pkb.len(), sig: sigb.len() });
    let st = measure(RUNS_COMPACT_VERIFY, || {
        let pk = CompactPublicKey::<Level1>::from_bytes(&pkb).unwrap();
        let s = CompactSignature::<Level1>::from_bytes(&sigb).unwrap();
        assert!(pk.verify(MSG, &s).is_ok());
    });
    rows.push(Row { impl_: R2, round: 2, dim: 4, level, op: "verify from bytes", st });
}

fn main() {
    let only_r3 = std::env::args().any(|a| a == "--r3-only");
    let mut rows = Vec::new();
    let mut sizes = Vec::new();
    let mut rng = StdRng::seed_from_u64(20260926);

    eprintln!("round 3, dimension 2 ...");
    r3_level::<sqisign_rs::Level1>("I", &mut rng, &mut rows, &mut sizes);
    r3_level::<sqisign_rs::Level3>("III", &mut rng, &mut rows, &mut sizes);
    r3_level::<sqisign_rs::Level5>("V", &mut rng, &mut rows, &mut sizes);
    eprintln!("round 3, dimension 4 ...");
    let (lp, mem) = r3_compact(&mut rows, &mut sizes);
    if !only_r3 {
        eprintln!("round 2, dimension 2 ...");
        r2_level!(sqisign_rs_r2::Level1, "I", &mut rng, rows, sizes);
        r2_level!(sqisign_rs_r2::Level3, "III", &mut rng, rows, sizes);
        r2_level!(sqisign_rs_r2::Level5, "V", &mut rng, rows, sizes);
        eprintln!("round 2, dimension 4 ...");
        r2_compact(&mut rng, &mut rows, &mut sizes);
    }

    println!("#### Sizes observed (bytes)\n");
    println!("| implementation | round | dimension | level | public key | signature |");
    println!("|---|---|---|---|---|---|");
    for s in &sizes {
        println!("| {} | {} | {} | {} | {} | {} |", s.impl_, s.round, s.dim, s.level, s.pk, s.sig);
    }
    println!("\n#### Key generation, signing, verification (Rust; megacycles over `runs`, and wall-clock milliseconds)\n");
    println!("| implementation | round | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |");
    println!("|---|---|---|---|---|---|---|---|---|---|---|");
    for r in &rows {
        let s = &r.st;
        println!(
            "| {} | {} | {} | {} | {} | {} | {:.1} | {:.1} | {:.1} | {:.1} | {:.2} |",
            r.impl_, r.round, r.dim, r.level, r.op, s.runs, s.median_mcyc, s.mean_mcyc, s.min_mcyc, s.max_mcyc, s.median_ms
        );
    }
    println!("\n#### The compact signer's rejection loop (round 3, level I, the `sign` rows' signatures)\n");
    println!("| signatures | response samples per signature (mean) | median | min | max | primality tests per signature (mean) |");
    println!("|---|---|---|---|---|---|");
    println!(
        "| {} | {:.1} | {:.1} | {} | {} | {:.1} |",
        lp.signatures, lp.samples_mean, lp.samples_median, lp.samples_min, lp.samples_max, lp.tests_mean
    );
    println!("\n#### Memory of one compact verification (round 3, level I)\n");
    println!("| peak heap during the call (bytes) | heap live before the call (bytes) | process VmHWM after (KiB) |");
    println!("|---|---|---|");
    println!(
        "| {} | {} | {} |",
        mem.peak_heap_bytes,
        mem.heap_before,
        mem.vm_hwm_kib.map(|k| k.to_string()).unwrap_or_else(|| "n/a".into())
    );
}
