//! The paper's benchmark protocol (Section "Setup" of
//! paper-dim4/main_standalone.tex), one pass per invocation:
//!
//!     paper --out <dir> [--only <block>]
//!
//! Every call is timed with `rdtsc` and the wall clock around it and
//! written as one CSV line `iteration,cycles,ns,extra` under `<dir>/raw/`,
//! so that the report script recomputes every statistic from the raw
//! files. Run counts, per the paper:
//!
//! * dimension 2, each level: 52 key generations, 52 signatures, 200
//!   verifications (standard and compressed), from bytes;
//! * compact, level I: 52 key generations, 100 signatures (with the
//!   good-degree loop's counts), 200 verifications of five signatures in
//!   each mode (four points; two points with the check at the gluing; two
//!   points without it, as an extra row);
//! * compact, levels III and V: 20 key generations, 52 signatures, and for
//!   verification six signatures, 15 runs each, in each mode;
//! * the phase profile of the four-point and two-point verifiers at each
//!   level, with the pre-chain checks on their own: the decode, the hint
//!   and pairing-order checks, the goodness test.
//!
//! Verification is from bytes: decode the key, decode the signature,
//! verify. The compact signatures are in the bit-packed layout. The heap
//! figures come from the `paper-heap` binary, which has the counting
//! allocator; this one uses the system allocator.

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fs::{self, File};
use std::hint::black_box;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign_stats,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::{
    challenge_limbs, hd_verify_bytes_with, hd_verify_profiled_full, parse_public_key,
    parse_signature,
    HdLevel, HdSignature, HdVerifyOptions, KernelPoints, SLOTS, SLOT_NAMES,
};
use sqisign_verify::{P324_3, P500_27, P664_17};

fn cycles() -> u64 {
    // SAFETY: rdtsc is available on every x86-64 CPU.
    unsafe { core::arch::x86_64::_rdtsc() }
}

struct Raw {
    file: File,
    iteration: u64,
}

impl Raw {
    fn open(dir: &Path, name: &str, header: &str) -> Raw {
        let mut file = File::create(dir.join(format!("{name}.csv"))).expect("raw file");
        writeln!(file, "# {header}").unwrap();
        writeln!(file, "iteration,cycles,ns,extra").unwrap();
        Raw { file, iteration: 0 }
    }
    /// Time `f` once and write its line; `f` returns the extra column.
    fn time<F: FnOnce() -> String>(&mut self, f: F) {
        let t0 = Instant::now();
        let c0 = cycles();
        let extra = f();
        let c1 = cycles();
        let ns = t0.elapsed().as_nanos();
        writeln!(self.file, "{},{},{},{}", self.iteration, c1 - c0, ns, extra).unwrap();
        self.iteration += 1;
    }
}

const MSG: &[u8] = b"dimension-4 paper session";

fn dim2<L: sqisign_rs::SigningLevel>(level: &str, raw: &Path) {
    use sqisign_rs::{generate, CompressedSignature, PublicKey, Signature};
    let mut rng = StdRng::seed_from_u64(0x2026_1114 ^ ((level.len() as u64) << 32));
    let mut f = Raw::open(
        raw,
        &format!("dim2-{level}-keygen"),
        &format!("dimension 2 | {level} | keygen, 52 runs"),
    );
    for _ in 0..52 {
        f.time(|| {
            black_box(generate::<L>(&mut rng));
            String::new()
        });
    }
    let (pk, sk) = generate::<L>(&mut rng);
    let mut f = Raw::open(
        raw,
        &format!("dim2-{level}-sign"),
        &format!("dimension 2 | {level} | sign, a fresh message each run, 52 runs"),
    );
    for i in 0..52u64 {
        f.time(|| {
            black_box(sk.sign(&i.to_le_bytes(), &mut rng).unwrap());
            String::new()
        });
    }
    let sig = sk.sign(MSG, &mut rng).unwrap();
    let (pkb, sigb) = (pk.to_bytes(), sig.to_bytes());
    let mut f = Raw::open(
        raw,
        &format!("dim2-{level}-verify"),
        &format!(
            "dimension 2 | {level} | verify from bytes, standard format ({} + {} bytes), 200 runs",
            pkb.len(),
            sigb.len()
        ),
    );
    for _ in 0..200 {
        f.time(|| {
            let pk = PublicKey::<L>::from_bytes(black_box(&pkb)).unwrap();
            let sig = Signature::<L>::from_bytes(black_box(&sigb)).unwrap();
            assert!(pk.verify(MSG, &sig).is_ok());
            String::new()
        });
    }
    let cb = sig.compress().to_bytes();
    let mut f = Raw::open(
        raw,
        &format!("dim2-{level}-verify-compressed"),
        &format!(
            "dimension 2 | {level} | verify compressed from bytes ({} bytes), 200 runs",
            cb.len()
        ),
    );
    for _ in 0..200 {
        f.time(|| {
            let pk = PublicKey::<L>::from_bytes(black_box(&pkb)).unwrap();
            let c = CompressedSignature::<L>::from_bytes(black_box(&cb)).unwrap();
            assert!(pk.verify_compressed(MSG, &c).is_ok());
            String::new()
        });
    }
}

const MODES: [(&str, KernelPoints); 3] = [
    ("four", KernelPoints::Four),
    ("two", KernelPoints::TwoWithImageCheck),
    ("two-nocheck", KernelPoints::TwoWithoutImageCheck),
];

/// The compact format at one level: `keygens` key generations, `signs`
/// signatures with the loop's counts, and `sigs` signatures verified
/// `runs_per_sig` times each in each mode, the signatures taken in turn.
fn compact<L: HdLevel, const N: usize>(
    level: &str,
    params: &Params<N>,
    keygens: usize,
    signs: usize,
    sigs: usize,
    runs_per_sig: usize,
    raw: &Path,
) {
    let (pk_len, sig_len) = (
        compact_public_key_bytes::<L>(),
        compact_signature_bytes::<L>(),
    );
    let mut r = ShakeRng::new(format!("paper session {level} keygen").as_bytes(), b"kgn");
    let mut f = Raw::open(
        raw,
        &format!("dim4-{level}-keygen"),
        &format!("compact | {level} | keygen, {keygens} runs"),
    );
    for _ in 0..keygens {
        f.time(|| {
            black_box(compact_keygen::<L, N>(params, &mut r).unwrap());
            String::new()
        });
    }

    let mut r = ShakeRng::new(format!("paper session {level} sign").as_bytes(), b"sgn");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut r).unwrap();
    let mut sig = vec![0u8; sig_len];
    let mut f = Raw::open(
        raw,
        &format!("dim4-{level}-sign"),
        &format!("compact | {level} | sign, a fresh message each run, {signs} runs; extra = response samples;primality tests;restarts"),
    );
    for i in 0..signs as u64 {
        f.time(|| {
            let (_, s) =
                compact_sign_stats::<L, N>(params, &pk, &sk, &i.to_le_bytes(), &mut r, &mut sig)
                    .unwrap();
            format!("{};{};{}", s.response_samples, s.primality_tests, s.restarts)
        });
    }

    // the signatures every mode verifies
    let mut r = ShakeRng::new(format!("paper session {level} set").as_bytes(), b"set");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut r).unwrap();
    let mut pkb = vec![0u8; pk_len];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let set: Vec<(Vec<u8>, Vec<u8>)> = (0..sigs)
        .map(|j| {
            let msg = format!("level {level} signature {j}").into_bytes();
            let mut s = vec![0u8; sig_len];
            compact_sign_stats::<L, N>(params, &pk, &sk, &msg, &mut r, &mut s).unwrap();
            (s, msg)
        })
        .collect();
    for (name, mode) in MODES {
        let opts = HdVerifyOptions {
            kernel_points: mode,
        };
        let mut f = Raw::open(
            raw,
            &format!("dim4-{level}-verify-{name}"),
            &format!("compact | {level} | verify from bytes ({pk_len} + {sig_len} bytes, packed layout), {mode:?}, {sigs} signatures x {runs_per_sig} runs in turn; extra = signature;fell back F1;fell back F2"),
        );
        for it in 0..sigs * runs_per_sig {
            let j = it % sigs;
            let (s, msg) = &set[j];
            f.time(|| {
                let (verdict, report) =
                    hd_verify_bytes_with::<L>(black_box(s), black_box(&pkb), msg, &opts);
                assert_eq!(verdict, Ok(()));
                format!(
                    "{j};{};{}",
                    report.fell_back[0] as u8, report.fell_back[1] as u8
                )
            });
        }
    }

    // the phase profile, four points and two
    for (name, mode) in [MODES[0], MODES[1]] {
        let opts = HdVerifyOptions {
            kernel_points: mode,
        };
        let mut file = File::create(raw.join(format!("dim4-{level}-profile-{name}.csv"))).unwrap();
        writeln!(file, "# compact | {level} | phase profile of verify from bytes, {mode:?}: cycles per phase, both half-chains added; total = rdtsc around decode and verify").unwrap();
        for (k, n) in SLOT_NAMES.iter().enumerate() {
            writeln!(file, "# slot{k} = {n}").unwrap();
        }
        let slots: Vec<String> = (0..SLOTS).map(|k| format!("slot{k}")).collect();
        writeln!(file, "# checks_hints_pairing = the canonical hints of both curves and the order of e(P, Q) (sub-slot 16); checks_goodness = the goodness test on q (sub-slot 17); both lie inside the slots, decode does not").unwrap();
        writeln!(
            file,
            "iteration,signature,total,decode,{},checks_hints_pairing,checks_goodness",
            slots.join(",")
        )
        .unwrap();
        let clock: &dyn Fn() -> u64 = &cycles;
        for it in 0..sigs * 10 {
            let j = it % sigs;
            let (s, msg) = &set[j];
            let c0 = cycles();
            let ps = parse_signature::<L>(s).unwrap();
            let pp = parse_public_key::<L>(&pkb).unwrap();
            let chal = challenge_limbs::<L>(&ps, &pp, msg).unwrap();
            let hdsig = HdSignature {
                a_pk: pp.a_pk.clone(),
                a_com: ps.a_com.clone(),
                hint_pk_p: pp.hint_pk_p,
                hint_pk_q: pp.hint_pk_q,
                hint_com_p: ps.hint_com_p,
                hint_com_q: ps.hint_com_q,
                chal_limbs: &chal,
                resp_a: ps.a,
                resp_b: ps.b,
                resp_c_or_d: ps.c_or_d,
                q: ps.q,
            };
            let c1 = cycles();
            let (verdict, _, sl, subs) = hd_verify_profiled_full::<L>(&hdsig, &opts, clock);
            let c2 = cycles();
            assert_eq!(verdict, Ok(()));
            let cols: Vec<String> = sl.iter().map(|v| v.to_string()).collect();
            writeln!(
                file,
                "{it},{j},{},{},{},{},{}",
                c2 - c0,
                c1 - c0,
                cols.join(","),
                subs[16],
                subs[17]
            )
            .unwrap();
        }
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let (mut out, mut only) = (PathBuf::from("paper-session-out"), None::<String>);
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out" => {
                i += 1;
                out = PathBuf::from(&argv[i]);
            }
            "--only" => {
                i += 1;
                only = Some(argv[i].clone());
            }
            s => panic!("unknown argument {s}"),
        }
        i += 1;
    }
    let raw = out.join("raw");
    fs::create_dir_all(&raw).unwrap();
    let want = |b: &str| only.as_deref().is_none_or(|o| o == b);
    eprintln!(
        "field backend: {}",
        sqisign_verify::fp::arithmetic_backend()
    );
    fs::write(
        out.join("backend.txt"),
        format!("{}\n", sqisign_verify::fp::arithmetic_backend()),
    )
    .unwrap();
    const N1: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
    const N3: usize = sqisign_rs::precomp::p500_27::IBZ_NLIMBS;
    const N5: usize = sqisign_rs::precomp::p664_17::IBZ_NLIMBS;
    let blocks: [(&str, Box<dyn Fn()>); 6] = [
        ("dim2-I", Box::new(|| dim2::<sqisign_rs::Level1>("I", &raw))),
        ("dim2-III", Box::new(|| dim2::<sqisign_rs::Level3>("III", &raw))),
        ("dim2-V", Box::new(|| dim2::<sqisign_rs::Level5>("V", &raw))),
        (
            "dim4-I",
            Box::new(|| compact::<P324_3, N1>("I", &level1(), 52, 100, 5, 40, &raw)),
        ),
        (
            "dim4-III",
            Box::new(|| compact::<P500_27, N3>("III", &level3(), 20, 52, 6, 15, &raw)),
        ),
        (
            "dim4-V",
            Box::new(|| compact::<P664_17, N5>("V", &level5(), 20, 52, 6, 15, &raw)),
        ),
    ];
    for (name, run) in blocks.iter() {
        if want(name) {
            let t = Instant::now();
            run();
            eprintln!("{name}: {:.1} s", t.elapsed().as_secs_f64());
        }
    }
}
