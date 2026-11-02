//! The dimension-4 paper's session binary (SECENG-2423, SECENG-2424): one
//! block per invocation, every iteration written as one CSV line, so that a
//! session script can run the blocks round-robin on one pinned core and a
//! report script can recompute every statistic from the raw files.
//!
//!     dim4-session <block> --out <dir> --round <k> [--runs <n>] [--runs-keygen <n>]
//!                  [--runs-sign <n>] [--runs-verify <n>] [--sigs <n>] [--seed <s>]
//!
//! Blocks: `dim2-I`, `dim2-III`, `dim2-V` (round-3 standard format: keygen,
//! sign, verify from bytes, verify compressed from bytes), `dim4-keygen`,
//! `dim4-sign` (with the signer's rejection-loop counts), `dim4-verify-four`,
//! `dim4-verify-two-check`, `dim4-verify-two-nocheck` (the compact format,
//! from bytes, at `KernelPoints::Four`, `TwoWithImageCheck`,
//! `TwoWithoutImageCheck`), `profile-two-check` and `profile-four` (the
//! verifier's phases through `hd_verify_profiled`, decoding timed here; with
//! the cargo feature `op-counters` also the `Fp2` operation counts per
//! phase, SECENG-2434), and `fp2-ops` (the cycles per `Fp2` operation in
//! this binary, for the floor table of SECENG-2434).
//! The three dimension-4 verify blocks and the two profile blocks of one
//! round verify the same signatures: the signature set is seeded by the
//! round number alone. Timing is `rdtsc` and the wall clock around each
//! call, as in the 2026-09-26 session; the old all-in-one Markdown run is
//! `legacy`.
//!
//! Round-2 rows (the `v0.4.31` crate) are behind the cargo feature `round2`
//! and are not part of the session.

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fs::{self, File, OpenOptions};
use std::hint::black_box;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use sqisign_verify::fp::opcount::{self, OPS, OP_NAMES};
use sqisign_verify::fp::Fp2;
use sqisign_verify::hd::{
    challenge_limbs, hd_verify_bytes_with, hd_verify_profiled_record, parse_public_key,
    parse_signature, HdSignature, HdVerifyOptions, KernelPoints, SLOTS, SLOT_NAMES, SUB_NAMES,
    SUB_SLOTS,
};
use sqisign_verify::P324_3;

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

/// One raw file: `round,iteration,cycles,ns,extra` per line, appended.
struct Raw {
    file: File,
    round: u32,
    iteration: u64,
}

impl Raw {
    fn open(dir: &Path, name: &str, header: &str, round: u32) -> Raw {
        let path = dir.join(format!("{name}.csv"));
        let new = !path.exists();
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .expect("raw file");
        if new {
            writeln!(file, "# {header}").unwrap();
            writeln!(file, "round,iteration,cycles,ns,extra").unwrap();
        }
        Raw {
            file,
            round,
            iteration: 0,
        }
    }
    fn line(&mut self, cyc: u64, ns: u128, extra: &str) {
        writeln!(
            self.file,
            "{},{},{},{},{}",
            self.round, self.iteration, cyc, ns, extra
        )
        .unwrap();
        self.iteration += 1;
    }
}

/// Time `f` once: cycles and nanoseconds.
fn timed<T, F: FnOnce() -> T>(f: F) -> (T, u64, u128) {
    let t0 = Instant::now();
    let c0 = cycles();
    let r = f();
    let c1 = cycles();
    let ns = t0.elapsed().as_nanos();
    (r, c1 - c0, ns)
}

/// Repeat `f` `runs` times, one raw line each; `f` returns the extra column.
fn repeat<F: FnMut() -> String>(raw: &mut Raw, runs: usize, mut f: F) {
    for _ in 0..runs {
        let (extra, cyc, ns) = timed(&mut f);
        raw.line(cyc, ns, &extra);
    }
}

// ---- arguments ------------------------------------------------------------

struct Args {
    block: String,
    out: PathBuf,
    round: u32,
    runs: usize,
    runs_keygen: usize,
    runs_sign: usize,
    runs_verify: usize,
    sigs: usize,
    seed: u64,
}

fn parse_args() -> Args {
    let mut a = Args {
        block: String::new(),
        out: PathBuf::from("session-out"),
        round: 1,
        runs: 50,
        runs_keygen: 13,
        runs_sign: 13,
        runs_verify: 50,
        sigs: 5,
        seed: 20261101,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let v = |i: &mut usize| -> String {
            *i += 1;
            argv.get(*i).cloned().expect("value")
        };
        match argv[i].as_str() {
            "--out" => a.out = PathBuf::from(v(&mut i)),
            "--round" => a.round = v(&mut i).parse().unwrap(),
            "--runs" => a.runs = v(&mut i).parse().unwrap(),
            "--runs-keygen" => a.runs_keygen = v(&mut i).parse().unwrap(),
            "--runs-sign" => a.runs_sign = v(&mut i).parse().unwrap(),
            "--runs-verify" => a.runs_verify = v(&mut i).parse().unwrap(),
            "--sigs" => a.sigs = v(&mut i).parse().unwrap(),
            "--seed" => a.seed = v(&mut i).parse().unwrap(),
            s if !s.starts_with("--") && a.block.is_empty() => a.block = s.to_string(),
            s => panic!("unknown argument {s}"),
        }
        i += 1;
    }
    a
}

// ---- dimension 2 (the round-3 standard format, this crate) ------------------

const MSG: &[u8] = b"dimension-4 paper session";

fn dim2_level<L: sqisign_rs::SigningLevel>(level: &str, a: &Args, raw_dir: &Path) {
    use sqisign_rs::{generate, CompressedSignature, PublicKey, Signature};
    let mut rng = StdRng::seed_from_u64(a.seed ^ (a.round as u64) ^ (level.len() as u64) << 32);
    let tag = format!("dim2-{level}");
    let mut raw = Raw::open(
        raw_dir,
        &format!("{tag}-keygen"),
        &format!("dimension 2 | {level} | keygen (generate)"),
        a.round,
    );
    repeat(&mut raw, a.runs_keygen, || {
        let kp = generate::<L>(&mut rng);
        black_box(kp);
        String::new()
    });
    let (pk, sk) = generate::<L>(&mut rng);
    let mut raw = Raw::open(
        raw_dir,
        &format!("{tag}-sign"),
        &format!("dimension 2 | {level} | sign (a fresh message each iteration)"),
        a.round,
    );
    let mut i = 0u64;
    repeat(&mut raw, a.runs_sign, || {
        i += 1;
        let sig = sk.sign(&i.to_le_bytes(), &mut rng).unwrap();
        black_box(sig);
        String::new()
    });
    let sig = sk.sign(MSG, &mut rng).unwrap();
    let pkb = pk.to_bytes();
    let sigb = sig.to_bytes();
    let mut raw = Raw::open(
        raw_dir,
        &format!("{tag}-verify"),
        &format!("dimension 2 | {level} | verify from bytes (standard format)"),
        a.round,
    );
    repeat(&mut raw, a.runs_verify, || {
        let pk = PublicKey::<L>::from_bytes(&pkb).unwrap();
        let sig = Signature::<L>::from_bytes(&sigb).unwrap();
        assert!(pk.verify(MSG, &sig).is_ok());
        String::new()
    });
    let cb = sig.compress().to_bytes();
    let mut raw = Raw::open(
        raw_dir,
        &format!("{tag}-verify-compressed"),
        &format!("dimension 2 | {level} | verify compressed from bytes"),
        a.round,
    );
    repeat(&mut raw, a.runs_verify, || {
        let pk = PublicKey::<L>::from_bytes(&pkb).unwrap();
        let c = CompressedSignature::<L>::from_bytes(&cb).unwrap();
        assert!(pk.verify_compressed(MSG, &c).is_ok());
        String::new()
    });
}

// ---- dimension 4 (the compact format, level I) ------------------------------

const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const PK_BYTES: usize = 84;
const SIG_BYTES: usize = 142;

/// A key and `n` signatures of this round, the same in every block of the
/// round that asks for them.
struct Set {
    pk: [u8; PK_BYTES],
    sigs: Vec<([u8; SIG_BYTES], Vec<u8>)>,
}

fn signature_set(a: &Args) -> Set {
    use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign};
    use sqisign_rs::mp::ShakeRng;
    use sqisign_rs::sqisign::level1;
    let params = level1();
    let label = format!("dimension-4 session round {}", a.round);
    let mut r = ShakeRng::new(label.as_bytes(), b"set");
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
    let mut pkb = [0u8; PK_BYTES];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let mut sigs = Vec::with_capacity(a.sigs);
    for j in 0..a.sigs {
        let msg = format!("round {} signature {j}", a.round).into_bytes();
        let mut sig = [0u8; SIG_BYTES];
        compact_sign::<P324_3, N>(&params, &pk, &sk, &msg, &mut r, &mut sig).unwrap();
        sigs.push((sig, msg));
    }
    Set { pk: pkb, sigs }
}

fn dim4_keygen(a: &Args, raw_dir: &Path) {
    use sqisign_rs::compact::compact_keygen;
    use sqisign_rs::mp::ShakeRng;
    use sqisign_rs::sqisign::level1;
    let params = level1();
    let label = format!("dimension-4 session keygen round {}", a.round);
    let mut r = ShakeRng::new(label.as_bytes(), b"kgn");
    let mut raw = Raw::open(
        raw_dir,
        "dim4-keygen",
        "dimension 4 | I | keygen (compact_keygen)",
        a.round,
    );
    repeat(&mut raw, a.runs_keygen, || {
        let kp = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
        black_box(kp);
        String::new()
    });
}

fn dim4_sign(a: &Args, raw_dir: &Path) {
    use sqisign_rs::compact::{compact_keygen, compact_sign_stats};
    use sqisign_rs::mp::ShakeRng;
    use sqisign_rs::sqisign::level1;
    let params = level1();
    let label = format!("dimension-4 session sign round {}", a.round);
    let mut r = ShakeRng::new(label.as_bytes(), b"sgn");
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
    let mut sig = [0u8; SIG_BYTES];
    let mut raw = Raw::open(
        raw_dir,
        "dim4-sign",
        "dimension 4 | I | sign (compact_sign_stats, a fresh message each iteration); extra = response samples;primality tests",
        a.round,
    );
    let mut i = 0u64;
    repeat(&mut raw, a.runs_sign, || {
        i += 1;
        let (_, s) =
            compact_sign_stats::<P324_3, N>(&params, &pk, &sk, &i.to_le_bytes(), &mut r, &mut sig)
                .unwrap();
        format!("{};{}", s.response_samples, s.primality_tests)
    });
}

fn dim4_verify(a: &Args, raw_dir: &Path, name: &str, mode: KernelPoints, label: &str) {
    let set = signature_set(a);
    let opts = HdVerifyOptions {
        kernel_points: mode,
    };
    let mut raw = Raw::open(
        raw_dir,
        name,
        &format!("dimension 4 | I | verify from bytes, {label} (hd_verify_bytes_with); extra = signature index;fell_back F1;fell_back F2"),
        a.round,
    );
    let mut j = 0usize;
    repeat(&mut raw, a.runs_verify, || {
        let (sig, msg) = &set.sigs[j % set.sigs.len()];
        let (verdict, report) = hd_verify_bytes_with::<P324_3>(sig, &set.pk, msg, &opts);
        assert_eq!(verdict, Ok(()));
        let extra = format!(
            "{};{};{}",
            j % set.sigs.len(),
            report.fell_back[0] as u8,
            report.fell_back[1] as u8
        );
        j += 1;
        extra
    });
}

/// The phase profile: decoding timed here, the verifier's slots from
/// `hd_verify_profiled`, the whole call timed around both.
fn profile(a: &Args, raw_dir: &Path, name: &str, mode: KernelPoints, label: &str) {
    let set = signature_set(a);
    let opts = HdVerifyOptions {
        kernel_points: mode,
    };
    let path = raw_dir.join(format!("{name}.csv"));
    let new = !path.exists();
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .unwrap();
    if new {
        writeln!(file, "# dimension 4 | I | phase profile of verify from bytes, {label}; cycles per phase, both half-chains added; total = rdtsc around decode and verify").unwrap();
        let slots: Vec<String> = (0..SLOTS).map(|k| format!("slot{k}")).collect();
        writeln!(
            file,
            "round,iteration,signature,total,decode,{}",
            slots.join(",")
        )
        .unwrap();
        for (k, n) in SLOT_NAMES.iter().enumerate() {
            writeln!(file, "# slot{k} = {n}").unwrap();
        }
    }
    // The sub-slots (SECENG-2428), working evidence, in a second file.
    let sub_path = raw_dir.join(format!("{name}-sub.csv"));
    let sub_new = !sub_path.exists();
    let mut subf = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&sub_path)
        .unwrap();
    if sub_new {
        writeln!(subf, "# dimension 4 | I | the sub-slots of the setup phase of verify from bytes, {label}; cycles per sub-slot, both half-chains added where a sub-slot runs per half").unwrap();
        let subs: Vec<String> = (0..SUB_SLOTS).map(|k| format!("sub{k}")).collect();
        writeln!(subf, "round,iteration,signature,{}", subs.join(",")).unwrap();
        for (k, n) in SUB_NAMES.iter().enumerate() {
            writeln!(subf, "# sub{k} = {n}").unwrap();
        }
    }
    // The Fp2 operation counts per phase (SECENG-2434): all zero unless the
    // binary has the cargo feature op-counters; the remainder is the whole
    // call minus the decode and the eight slots.
    let regions: Vec<String> = ["decode".to_string()]
        .into_iter()
        .chain((0..SLOTS).map(|k| format!("slot{k}")))
        .chain(["total".to_string(), "remainder".to_string()])
        .collect();
    let ops_path = raw_dir.join(format!("{name}-ops.csv"));
    let ops_new = !ops_path.exists();
    let mut opsf = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&ops_path)
        .unwrap();
    if ops_new {
        writeln!(opsf, "# dimension 4 | I | the Fp2 operation counts per phase of verify from bytes, {label}; counters {}; unmatched = the spans the profiler could not charge", if opcount::ENABLED { "compiled in" } else { "not compiled in, every count zero" }).unwrap();
        let cols: Vec<String> = regions
            .iter()
            .flat_map(|r| OP_NAMES.iter().map(move |o| format!("{r}_{o}")))
            .collect();
        writeln!(
            opsf,
            "round,iteration,signature,unmatched,{}",
            cols.join(",")
        )
        .unwrap();
    }
    let subops_path = raw_dir.join(format!("{name}-subops.csv"));
    let subops_new = !subops_path.exists();
    let mut subopsf = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&subops_path)
        .unwrap();
    if subops_new {
        writeln!(subopsf, "# dimension 4 | I | the Fp2 operation counts per sub-slot of verify from bytes, {label}; counters {}", if opcount::ENABLED { "compiled in" } else { "not compiled in, every count zero" }).unwrap();
        let cols: Vec<String> = (0..SUB_SLOTS)
            .flat_map(|k| OP_NAMES.iter().map(move |o| format!("sub{k}_{o}")))
            .collect();
        writeln!(subopsf, "round,iteration,signature,{}", cols.join(",")).unwrap();
    }
    let clock: &dyn Fn() -> u64 = &cycles;
    for it in 0..a.runs {
        let j = it % set.sigs.len();
        let (sig, msg) = &set.sigs[j];
        let o0 = opcount::snapshot();
        let c0 = cycles();
        let s = parse_signature::<P324_3>(sig).unwrap();
        let p = parse_public_key::<P324_3>(&set.pk).unwrap();
        let chal = challenge_limbs::<P324_3>(&s, &p, msg).unwrap();
        let hdsig = HdSignature {
            a_pk: p.a_pk.clone(),
            a_com: s.a_com.clone(),
            hint_pk_p: p.hint_pk_p,
            hint_pk_q: p.hint_pk_q,
            hint_com_p: s.hint_com_p,
            hint_com_q: s.hint_com_q,
            chal_limbs: &chal,
            resp_a: s.a,
            resp_b: s.b,
            resp_c_or_d: s.c_or_d,
            q: s.q,
        };
        let c1 = cycles();
        let o1 = opcount::snapshot();
        let (verdict, _report, rec) = hd_verify_profiled_record::<P324_3>(&hdsig, &opts, clock);
        let c2 = cycles();
        let o2 = opcount::snapshot();
        assert_eq!(verdict, Ok(()));
        let (slots, subs) = (rec.slots, rec.subs);
        let subcols: Vec<String> = subs.iter().map(|v| v.to_string()).collect();
        writeln!(subf, "{},{},{},{}", a.round, it, j, subcols.join(",")).unwrap();
        let decode: [u64; OPS] = core::array::from_fn(|k| o1[k] - o0[k]);
        let total: [u64; OPS] = core::array::from_fn(|k| o2[k] - o0[k]);
        let remainder: [u64; OPS] = core::array::from_fn(|k| {
            total[k] - decode[k] - rec.ops.iter().map(|s| s[k]).sum::<u64>()
        });
        let mut opcols: Vec<String> = Vec::new();
        for region in [decode]
            .iter()
            .chain(rec.ops.iter())
            .chain([total, remainder].iter())
        {
            opcols.extend(region.iter().map(|v| v.to_string()));
        }
        writeln!(
            opsf,
            "{},{},{},{},{}",
            a.round,
            it,
            j,
            rec.unmatched,
            opcols.join(",")
        )
        .unwrap();
        let subopcols: Vec<String> = rec
            .sub_ops
            .iter()
            .flat_map(|s| s.iter().map(|v| v.to_string()))
            .collect();
        writeln!(subopsf, "{},{},{},{}", a.round, it, j, subopcols.join(",")).unwrap();
        let cols: Vec<String> = slots.iter().map(|v| v.to_string()).collect();
        writeln!(
            file,
            "{},{},{},{},{},{}",
            a.round,
            it,
            j,
            c2 - c0,
            c1 - c0,
            cols.join(",")
        )
        .unwrap();
    }
}

// ---- the cost of the field operations in this binary (SECENG-2434) ---------

type F = Fp2<P324_3>;

/// A field element from a seed, by squarings and additions.
fn elt(seed: u64) -> F {
    let mut x = F::from_small(seed | 1);
    let s = F::from_small(seed ^ 0x5555);
    for _ in 0..40 {
        x.sqr_assign();
        x.add_assign(&s);
    }
    x
}

/// A timed loop repeated `reps` times: the median and the minimum of its cycles.
fn med_min<G: FnMut() -> u64>(reps: usize, mut g: G) -> (u64, u64) {
    let mut v: Vec<u64> = (0..reps).map(|_| g()).collect();
    v.sort_unstable();
    (v[reps / 2], v[0])
}

/// The cycles per `Fp2` operation in this binary, on the kernels in use: a
/// loop of `n` operations, each on the previous result (`dep`), the median
/// and the minimum over the repetitions; `lane16` runs sixteen independent
/// elements per iteration (the access pattern of the theta coordinates) and
/// reports cycles per element operation; the composites are the dimension-4
/// kernels on one point, cycles per call.
fn fp2_ops(a: &Args, raw_dir: &Path) {
    use sqisign_verify::hd::arith::{
        hadamard_assign, pointwise_mul_assign, pointwise_square_assign, to_squared_theta_assign,
    };
    use sqisign_verify::hd::field::products_of_others;
    use sqisign_verify::hd::isogeny::apply_plain_image;
    use sqisign_verify::hd::{IsogenyDim4, ThetaPointDim4, ThetaStructureDim4};
    let path = raw_dir.join("fp2-ops.csv");
    let new = !path.exists();
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .unwrap();
    if new {
        writeln!(file, "# dimension 4 | I | cycles per Fp2 operation in this binary (fp2-ops): a loop of n operations, median and minimum over reps repetitions, rdtsc; dep = each operation on the previous result; lane16 = sixteen independent elements per iteration, cycles per element operation; composite = the dimension-4 kernels on one point, cycles per call").unwrap();
        writeln!(
            file,
            "round,op,variant,n,reps,median_per_op,min_per_op,backend"
        )
        .unwrap();
    }
    let backend = sqisign_verify::fp::arithmetic_backend();
    let reps = 11usize;
    let mut row = |op: &str, variant: &str, n: u64, per: u64, (med, min): (u64, u64)| {
        writeln!(
            file,
            "{},{},{},{},{},{:.2},{:.2},{}",
            a.round,
            op,
            variant,
            n,
            reps,
            med as f64 / per as f64,
            min as f64 / per as f64,
            backend
        )
        .unwrap();
    };
    let a0 = elt(3);
    let b = elt(7);
    macro_rules! dep {
        ($op:expr, $n:expr, |$x:ident| $body:expr) => {{
            let n: u64 = $n;
            let r = med_min(reps, || {
                let mut $x = a0.clone();
                let t = cycles();
                for _ in 0..n {
                    $body;
                }
                black_box(&$x);
                cycles() - t
            });
            row($op, "dep", n, n, r);
        }};
    }
    dep!("add", 2_000_000, |x| x.add_assign(black_box(&b)));
    dep!("add_out", 2_000_000, |x| x = x.add(black_box(&b)));
    dep!("sub", 2_000_000, |x| x.sub_assign(black_box(&b)));
    dep!("neg", 2_000_000, |x| x = x.neg());
    dep!("mul", 1_000_000, |x| x.mul_assign(black_box(&b)));
    dep!("mul_out", 1_000_000, |x| x = x.mul(black_box(&b)));
    dep!("sqr", 1_000_000, |x| x.sqr_assign());
    dep!("sqr_out", 1_000_000, |x| x = x.sqr());
    dep!("mul_small", 1_000_000, |x| x = x.mul_small(3));
    dep!("half", 1_000_000, |x| x = x.half());
    dep!("inv", 20_000, |x| x = x.inv());
    dep!("sqrt", 5_000, |x| x = x.sqrt());
    // The Legendre symbol returns a bit; one addition keeps the chain dependent.
    dep!("legendre_plus_add", 5_000, |x| {
        black_box(x.is_square());
        x.add_assign(black_box(&b))
    });
    let p0: [F; 16] = core::array::from_fn(|i| elt(11 + i as u64));
    let w: [F; 16] = core::array::from_fn(|i| elt(101 + i as u64));
    macro_rules! lanes {
        ($op:expr, $variant:expr, $n:expr, $per:expr, |$p:ident| $body:expr) => {{
            let n: u64 = $n;
            let r = med_min(reps, || {
                let mut $p = p0.clone();
                let t = cycles();
                for _ in 0..n {
                    $body;
                    black_box(&mut $p);
                }
                cycles() - t
            });
            row($op, $variant, n, n * $per, r);
        }};
    }
    lanes!("add", "lane16", 100_000, 16, |p| for k in 0..16 {
        p[k].add_assign(&w[k])
    });
    lanes!("sub", "lane16", 100_000, 16, |p| for k in 0..16 {
        p[k].sub_assign(&w[k])
    });
    lanes!("mul", "lane16", 100_000, 16, |p| for k in 0..16 {
        p[k].mul_assign(&w[k])
    });
    lanes!("sqr", "lane16", 100_000, 16, |p| for k in 0..16 {
        p[k].sqr_assign()
    });
    lanes!("hadamard16", "composite", 200_000, 1, |p| hadamard_assign(
        &mut p
    ));
    lanes!("pointwise_mul16", "composite", 200_000, 1, |p| {
        pointwise_mul_assign(&mut p, &w)
    });
    lanes!("pointwise_sqr16", "composite", 200_000, 1, |p| {
        pointwise_square_assign(&mut p)
    });
    lanes!("to_squared_theta16", "composite", 200_000, 1, |p| {
        to_squared_theta_assign(&mut p)
    });
    lanes!("products_of_others16", "composite", 50_000, 1, |p| {
        p = products_of_others(&p)
    });
    let mut st = ThetaStructureDim4::<P324_3>::new(ThetaPointDim4::new(w.clone()));
    assert!(st.precompute());
    {
        let n = 50_000u64;
        let r = med_min(reps, || {
            let mut pt = ThetaPointDim4::new(p0.clone());
            let t = cycles();
            for _ in 0..n {
                pt = st.double(black_box(&pt));
            }
            black_box(&pt);
            cycles() - t
        });
        row("double", "composite", n, n, r);
        let r = med_min(reps, || {
            let mut pt = ThetaPointDim4::new(p0.clone());
            let t = cycles();
            for _ in 0..n {
                pt = apply_plain_image(&w, black_box(&pt));
            }
            black_box(&pt);
            cycles() - t
        });
        row("image", "composite", n, n, r);
    }
    {
        let n = 5_000u64;
        let k8: [ThetaPointDim4<P324_3>; 4] = core::array::from_fn(|i| {
            ThetaPointDim4::new(core::array::from_fn(|j| {
                elt(1000 + 16 * i as u64 + j as u64)
            }))
        });
        let r = med_min(reps, || {
            let t = cycles();
            for _ in 0..n {
                let iso = IsogenyDim4::from_kernel(black_box(&k8)).unwrap();
                let mut c = iso.codomain().clone();
                assert!(c.precompute());
                black_box(&c);
            }
            cycles() - t
        });
        row("from_kernel_precompute", "composite", n, n, r);
    }
    eprintln!("fp2-ops done, backend {backend}");
}

// ---- the legacy all-in-one run (the 2026-09-26 format) ----------------------

#[cfg(feature = "round2")]
mod legacy;

fn main() {
    let a = parse_args();
    if a.block == "legacy" {
        #[cfg(feature = "round2")]
        {
            legacy::run();
            return;
        }
        #[cfg(not(feature = "round2"))]
        panic!("the legacy run needs the cargo feature round2");
    }
    let raw_dir = a.out.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    eprintln!(
        "block {} round {} out {}",
        a.block,
        a.round,
        a.out.display()
    );
    match a.block.as_str() {
        "dim2-I" => dim2_level::<sqisign_rs::Level1>("I", &a, &raw_dir),
        "dim2-III" => dim2_level::<sqisign_rs::Level3>("III", &a, &raw_dir),
        "dim2-V" => dim2_level::<sqisign_rs::Level5>("V", &a, &raw_dir),
        "dim4-keygen" => dim4_keygen(&a, &raw_dir),
        "dim4-sign" => dim4_sign(&a, &raw_dir),
        "dim4-verify-four" => dim4_verify(
            &a,
            &raw_dir,
            "dim4-verify-four",
            KernelPoints::Four,
            "four points (KernelPoints::Four)",
        ),
        "dim4-verify-two-check" => dim4_verify(
            &a,
            &raw_dir,
            "dim4-verify-two-check",
            KernelPoints::TwoWithImageCheck,
            "two points with the image check (KernelPoints::TwoWithImageCheck)",
        ),
        "dim4-verify-two-nocheck" => dim4_verify(
            &a,
            &raw_dir,
            "dim4-verify-two-nocheck",
            KernelPoints::TwoWithoutImageCheck,
            "two points without the image check (KernelPoints::TwoWithoutImageCheck)",
        ),
        "profile-two-check" => profile(
            &a,
            &raw_dir,
            "profile-two-check",
            KernelPoints::TwoWithImageCheck,
            "two points with the image check",
        ),
        "profile-four" => profile(
            &a,
            &raw_dir,
            "profile-four",
            KernelPoints::Four,
            "four points",
        ),
        "fp2-ops" => fp2_ops(&a, &raw_dir),
        "basis-dim2" => basis_dim2(&a, &raw_dir),
        other => panic!("unknown block {other}"),
    }
}

/// The torsion basis from hints per curve, ours against the dimension-2
/// verifier's (basis-dim2): the dimension-2 verifier's
/// ec_curve_to_basis_2f_from_hint on its own level-I keys at its exponent
/// (challenge bits plus the extra torsion), and on the curves of this
/// round's set (E_pk and each E_com) at that exponent and at the full
/// exponent 324 with a hint computed outside the timing; ours
/// (hd_torsion_basis) on the same curves with the wire hints. Cycles per
/// basis, median and minimum over reps repetitions of n bases, rdtsc.
fn basis_dim2(a: &Args, raw_dir: &Path) {
    use sqisign_rs::generate;
    use sqisign_verify::ec::basis::{ec_curve_to_basis_2f_from_hint, ec_curve_to_basis_2f_to_hint};
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::hd::{hd_torsion_basis, parse_public_key, parse_signature};
    use sqisign_verify::precomp::{p324_3::TWO_ADIC_EXPONENT, EXTRA_TORSION};
    use sqisign_verify::{Level, Level1, PublicKey as PublicKeyDim2};
    let path = raw_dir.join("basis-dim2.csv");
    let new = !path.exists();
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .unwrap();
    if new {
        writeln!(file, "# the torsion basis from hints per curve (basis-dim2): cycles per basis, median and minimum over reps repetitions of n bases, rdtsc; kind = dim2-native (the dimension-2 verifier's function on its own level-I keys, f = challenge bits + 2), dim2-on-ours-f130 and dim2-on-ours-f324 (the same function on this round's E_pk and E_com with a hint computed outside the timing), ours (hd_torsion_basis with the wire hints, the full exponent)").unwrap();
        writeln!(file, "round,kind,curve,f,n,reps,median_per_basis,min_per_basis").unwrap();
    }
    let reps = 11usize;
    let n = 100u64;
    let mut row = |kind: &str, curve: &str, f: u32, (med, min): (u64, u64)| {
        writeln!(
            file,
            "{},{},{},{},{},{},{:.1},{:.1}",
            a.round,
            kind,
            curve,
            f,
            n,
            reps,
            med as f64 / n as f64,
            min as f64 / n as f64
        )
        .unwrap();
    };
    let measure = |body: &mut dyn FnMut()| -> (u64, u64) {
        let mut v: Vec<u64> = (0..reps)
            .map(|_| {
                let t = cycles();
                for _ in 0..n {
                    body();
                }
                cycles() - t
            })
            .collect();
        v.sort_unstable();
        (v[reps / 2], v[0])
    };
    // the dimension-2 verifier on its own keys
    let f2 = Level1::PARAMS.challenge_bits + EXTRA_TORSION;
    let mut rng = StdRng::seed_from_u64(a.seed ^ (a.round as u64) ^ 0xba515);
    for k in 0..5 {
        let (pk, _sk) = generate::<Level1>(&mut rng);
        let pkb = pk.to_bytes();
        let pk2 = PublicKeyDim2::<Level1>::from_bytes(&pkb).unwrap();
        let inner = pk2.inner();
        let curve0 = inner.curve.clone();
        let hint = inner.hint_pk;
        let r = measure(&mut || {
            let mut c = curve0.clone();
            black_box(ec_curve_to_basis_2f_from_hint(&mut c, f2, hint));
        });
        row("dim2-native", &format!("E_pk key {k}"), f2, r);
    }
    // this round's curves: E_pk and each E_com
    let set = signature_set(a);
    let p = parse_public_key::<P324_3>(&set.pk).unwrap();
    let mut curves: Vec<(String, Fp2<P324_3>, u32, u32)> =
        vec![("E_pk".to_string(), p.a_pk.clone(), p.hint_pk_p, p.hint_pk_q)];
    for (j, (sig, _msg)) in set.sigs.iter().enumerate() {
        let s = parse_signature::<P324_3>(sig).unwrap();
        curves.push((format!("E_com sig {j}"), s.a_com.clone(), s.hint_com_p, s.hint_com_q));
    }
    for (name, a_coeff, hp, hq) in &curves {
        let r = measure(&mut || {
            black_box(hd_torsion_basis::<P324_3>(a_coeff, *hp, *hq));
        });
        row("ours", name, TWO_ADIC_EXPONENT, r);
        for f in [f2, TWO_ADIC_EXPONENT] {
            let mut c = EcCurve::<P324_3>::from_a(a_coeff).unwrap();
            let (_basis, hint) = ec_curve_to_basis_2f_to_hint(&mut c, f, 0).unwrap();
            let r = measure(&mut || {
                let mut cc = c.clone();
                black_box(ec_curve_to_basis_2f_from_hint(&mut cc, f, hint));
            });
            row(&format!("dim2-on-ours-f{f}"), name, f, r);
        }
    }
}
