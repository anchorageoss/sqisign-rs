//! Generates the per-prime field backends of `sqisign-verify`.
//!
//! Every SQIsign / PRISM prime has the form `p = c * 2^e - 1` with a small
//! odd cofactor `c`. Field elements are held in `n` unsaturated limbs of
//! `r` bits each (`n * r` a few bits above `log2 p`), and arithmetic runs
//! in Montgomery form with `R = 2^(n*r)`. Because `2^(r*(n-1))` divides
//! `p + 1`, the Montgomery reduction step is a single multiply-add of
//! each freshly emitted low limb into the top limb position, which is
//! what makes the unrolled schoolbook `modmul` below fast.
//!
//! The emitted code is straight-line and constant-time in the field
//! elements. Only the choice of `(c, e, r, n)` and the derived constants
//! are prime-specific; everything else is the same text with `n`
//! substituted. This binary owns that text so a new parameter set is a
//! one-line table entry, and CI regenerates every backend and diffs it
//! against the checked-in copy.
//!
//! Usage:
//!   gen-fp --check DIR    regenerate all backends and diff against DIR
//!   gen-fp --write DIR    regenerate all backends into DIR
//!   gen-fp NAME           print one backend to stdout

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use std::fmt::Write as _;
use std::path::Path;

/// One parameter set. `radix` and `nlimbs` are the layout the SQIsign
/// reference implementation uses for the same prime, so that both
/// codebases make identical headroom assumptions.
struct PrimeSpec {
    /// Module and file name, e.g. `p324_3`.
    name: &'static str,
    /// Odd cofactor `c` in `p = c * 2^e - 1`.
    cofactor: u64,
    /// Two-adic exponent `e`.
    e: u32,
    /// Bits per unsaturated limb.
    radix: u32,
    /// Number of limbs.
    nlimbs: usize,
    /// Canonical encoding length `ceil(log2(p) / 8)`.
    encoded_bytes: usize,
    /// Where the prime comes from, for the module docs.
    origin: &'static str,
    /// Whether an x86-64 inline-assembly backend exists for this prime
    /// (`fp/<name>_asm.rs`), the implementation on x86-64; the generated
    /// trait implementation is then compiled out on that architecture.
    asm: bool,
}

const PRIMES: &[PrimeSpec] = &[
    PrimeSpec {
        name: "p324_3",
        cofactor: 3,
        e: 324,
        radix: 55,
        nlimbs: 6,
        encoded_bytes: 41,
        origin: "SQIsign round 3, NIST level I",
        asm: true,
    },
    PrimeSpec {
        name: "p500_27",
        cofactor: 27,
        e: 500,
        radix: 57,
        nlimbs: 9,
        encoded_bytes: 64,
        origin: "SQIsign round 3, NIST level III (also round 2, level V)",
        asm: true,
    },
    PrimeSpec {
        name: "p664_17",
        cofactor: 17,
        e: 664,
        radix: 61,
        nlimbs: 11,
        encoded_bytes: 84,
        origin: "SQIsign round 3, NIST level V",
        asm: true,
    },
];

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.as_slice() {
        [flag, dir] if flag == "--check" => {
            let mut bad = 0;
            for spec in PRIMES {
                let path = Path::new(dir).join(format!("{}.rs", spec.name));
                let want = rustfmt(generate(spec));
                let have = std::fs::read_to_string(&path).unwrap_or_default();
                if want != have {
                    eprintln!("stale: {}", path.display());
                    bad += 1;
                }
            }
            if bad > 0 {
                std::process::exit(1);
            }
            eprintln!("all {} backends up to date", PRIMES.len());
        }
        [flag, dir] if flag == "--write" => {
            for spec in PRIMES {
                let path = Path::new(dir).join(format!("{}.rs", spec.name));
                std::fs::write(&path, rustfmt(generate(spec))).expect("write backend");
                eprintln!("wrote {}", path.display());
            }
        }
        [name] => match PRIMES.iter().find(|s| s.name == name) {
            Some(spec) => print!("{}", rustfmt(generate(spec))),
            None => {
                eprintln!("unknown prime {name}");
                std::process::exit(2);
            }
        },
        _ => {
            eprintln!("usage: gen-fp (--check DIR | --write DIR | <name>)");
            std::process::exit(2);
        }
    }
}

/// Pipe generated source through `rustfmt` so the checked-in backends are
/// formatting-clean and `--check` compares like with like.
fn rustfmt(src: String) -> String {
    use std::io::Write as _;
    use std::process::{Command, Stdio};
    let mut child = Command::new("rustfmt")
        .args(["--edition", "2021", "--emit", "stdout"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("rustfmt must be installed (rustup component add rustfmt)");
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(src.as_bytes())
        .expect("feed rustfmt");
    let out = child.wait_with_output().expect("rustfmt exit");
    assert!(out.status.success(), "rustfmt rejected generated code");
    String::from_utf8(out.stdout).expect("utf8")
}

// ---------------------------------------------------------------------
// Arithmetic on the spec (build time only)
// ---------------------------------------------------------------------

fn prime(s: &PrimeSpec) -> BigUint {
    (BigUint::from(s.cofactor) << s.e) - BigUint::one()
}

/// Split `x < 2^(radix * nlimbs)` into little-endian radix-`2^radix`
/// limbs. The top limb may hold more than `radix` bits.
fn to_limbs(x: &BigUint, s: &PrimeSpec) -> Vec<u64> {
    let mask = (BigUint::one() << s.radix) - BigUint::one();
    let mut v = Vec::with_capacity(s.nlimbs);
    let mut t = x.clone();
    for i in 0..s.nlimbs {
        let limb = if i + 1 == s.nlimbs {
            t.clone()
        } else {
            &t & &mask
        };
        v.push(limb.to_u64().expect("limb fits in u64"));
        t >>= s.radix;
    }
    assert!(t.is_zero(), "value does not fit in {} limbs", s.nlimbs);
    v
}

fn fmt_limb_const(out: &mut String, doc: &str, name: &str, limbs: &[u64]) {
    writeln!(out, "/// {doc}").unwrap();
    writeln!(out, "pub const {name}: [u64; NLIMBS] = [").unwrap();
    for l in limbs {
        writeln!(out, "    0x{l:016x},").unwrap();
    }
    writeln!(out, "];").unwrap();
    writeln!(out).unwrap();
}

/// Inverse of a small integer `k` modulo `p`: the unique `x` with
/// `k * x = 1 mod p`, found as `(j * p + 1) / k` for the `j` that makes
/// the division exact.
fn small_inverse(k: u64, p: &BigUint) -> BigUint {
    for j in 1..k {
        let num = p * BigUint::from(j) + BigUint::one();
        let (q, r) = num.div_rem(&BigUint::from(k));
        if r.is_zero() {
            return q;
        }
    }
    unreachable!("k coprime to p")
}

// ---------------------------------------------------------------------
// Addition chain for the square-root progenitor exponent (p - 3) / 4
// ---------------------------------------------------------------------

/// Window width for the fixed sliding-window exponentiation. 16 odd
/// powers are precomputed; the chain then needs one squaring per bit of
/// the exponent and roughly one multiplication per `WINDOW + 1` bits.
const WINDOW: usize = 5;

enum ChainOp {
    /// Square the accumulator this many times.
    Sqr(usize),
    /// Multiply the accumulator by table entry `j` (`w^(2j+1)`).
    Mul(usize),
}

/// Left-to-right sliding-window decomposition of `exp`. The first op is
/// always a `Mul` that initialises the accumulator from the table.
fn sliding_window(exp: &BigUint) -> Vec<ChainOp> {
    let bits: Vec<bool> = (0..exp.bits()).map(|i| exp.bit(i)).collect();
    let mut ops = Vec::new();
    let mut i = bits.len() as isize - 1;
    let mut pending_sqr = 0usize;
    let mut first = true;
    while i >= 0 {
        if !bits[i as usize] {
            pending_sqr += 1;
            i -= 1;
            continue;
        }
        // Window [j, i] with bit j set, at most WINDOW bits wide.
        let mut j = (i - WINDOW as isize + 1).max(0);
        while !bits[j as usize] {
            j += 1;
        }
        let width = (i - j + 1) as usize;
        let mut w = 0usize;
        for k in (j..=i).rev() {
            w = (w << 1) | usize::from(bits[k as usize]);
        }
        debug_assert!(w % 2 == 1);
        if first {
            ops.push(ChainOp::Mul(w / 2));
            first = false;
        } else {
            ops.push(ChainOp::Sqr(pending_sqr + width));
            ops.push(ChainOp::Mul(w / 2));
        }
        pending_sqr = 0;
        i = j - 1;
    }
    if pending_sqr > 0 {
        ops.push(ChainOp::Sqr(pending_sqr));
    }
    ops
}

// ---------------------------------------------------------------------
// Code emission
// ---------------------------------------------------------------------

fn generate(s: &PrimeSpec) -> String {
    let n = s.nlimbs;
    let r = s.radix;
    let p = prime(s);
    let bits = p.bits();
    assert_eq!(
        s.encoded_bytes as u64,
        bits.div_ceil(8),
        "encoded_bytes must be ceil(log2 p / 8)"
    );
    assert!(
        (n as u32) * r > bits as u32,
        "limb layout too small for the prime"
    );
    let low = r * (n as u32 - 1);
    assert!(
        s.e >= low,
        "2^(radix*(nlimbs-1)) must divide p+1 for the folding reduction"
    );
    let ptop: u64 = s.cofactor << (s.e - low);
    let two_ptop = ptop.checked_mul(2).expect("2p top limb fits");

    let big_r = BigUint::one() << (r * n as u32);
    let one = &big_r % &p;
    let two_inv = (small_inverse(2, &p) * &big_r) % &p;
    let three_inv = (small_inverse(3, &p) * &big_r) % &p;
    let r2_nres = (&big_r * &big_r) % &p;
    let chunk_bytes = s.encoded_bytes - 1;
    assert!(
        ((8 * chunk_bytes) as u64) < bits,
        "decode_reduce chunks must be smaller than p"
    );
    let chunk_r = ((BigUint::one() << (8 * chunk_bytes)) * &big_r) % &p;
    let exp = (&p - BigUint::from(3u32)) >> 2;
    let chain = sliding_window(&exp);

    let mut o = String::new();
    let top = n - 1;
    let top_m1 = top - 1;
    let upper = s.name.to_uppercase();
    // an inline-assembly backend replaces this one's trait implementation
    // on x86-64; the functions here then go unused
    let asm_note = if s.asm {
        "\n//!\n//! On x86-64 the trait implementation lives in `<name>_asm.rs` (saturated\n//! limbs; `mulx`/`adcx`/`adox` or portable kernels chosen at run time), and\n//! this module is compiled but unused; elsewhere it is the implementation.\n#![cfg_attr(target_arch = \"x86_64\", allow(dead_code, unused_imports))]"
            .replace("<name>", s.name)
    } else {
        String::new()
    };

    // ---- header and constants ------------------------------------------
    writeln!(
        o,
        "//! Field arithmetic modulo `p = {c} * 2^{e} - 1` ({bits} bits).\n\
         //!\n\
         //! Origin: {origin}.\n\
         //!\n\
         //! Elements are {n} limbs of {r} bits (unsaturated, {cap} bits of\n\
         //! storage) in Montgomery form with `R = 2^{cap}`. Arithmetic results\n\
         //! are lazily reduced to `[0, 2p)`; encoding and comparison reduce\n\
         //! fully. `2^{low}` divides `p + 1`, so the Montgomery reduction folds\n\
         //! each emitted low limb `v` back in as `v * {c} * 2^{shift}` at limb\n\
         //! position {top}.\n\
         //!\n\
         //! GENERATED by `tools/gen-fp`. Do not edit; change the generator.{asm_note}\n",
        c = s.cofactor,
        e = s.e,
        origin = s.origin,
        cap = r * n as u32,
        shift = s.e - low,
    )
    .unwrap();
    writeln!(
        o,
        "use super::FpBackend;\n\
         use crate::params::{upper};\n\
         use hybrid_array::Array;\n\
         use subtle::Choice;\n"
    )
    .unwrap();
    writeln!(o, "/// Number of 64-bit limbs in an `Fp` element.").unwrap();
    writeln!(o, "pub const NLIMBS: usize = {n};\n").unwrap();
    writeln!(o, "/// Bit width of each unsaturated limb.").unwrap();
    writeln!(o, "pub const RADIX: u32 = {r};\n").unwrap();
    writeln!(o, "/// `2^RADIX - 1`.").unwrap();
    writeln!(o, "pub const MASK: u64 = (1u64 << RADIX) - 1;\n").unwrap();
    writeln!(o, "/// Canonical encoding length in bytes.").unwrap();
    writeln!(o, "pub const ENCODED_BYTES: usize = {};\n", s.encoded_bytes).unwrap();
    writeln!(
        o,
        "/// `(p + 1) / 2^{low}`: the contribution of `p + 1` to the top limb.\n\
         /// Montgomery folding constant in `modmul` and `modsqr`."
    )
    .unwrap();
    writeln!(o, "pub const PTOP: u64 = 0x{ptop:x};\n").unwrap();
    writeln!(
        o,
        "/// `2 * PTOP`, used to add or subtract `2p` in the lazy reductions."
    )
    .unwrap();
    writeln!(o, "pub const TWO_PTOP: u64 = 0x{two_ptop:x};\n").unwrap();
    writeln!(o, "/// The prime `p` as canonical little-endian bytes.").unwrap();
    writeln!(o, "pub const PRIME_LE_BYTES: [u8; ENCODED_BYTES] = [").unwrap();
    let pbytes = p.to_bytes_le();
    for chunk in pbytes.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|b| format!("0x{b:02x}")).collect();
        writeln!(o, "    {},", row.join(", ")).unwrap();
    }
    writeln!(o, "];\n").unwrap();
    writeln!(o, "/// `0`.").unwrap();
    writeln!(o, "pub const ZERO: [u64; NLIMBS] = [0; NLIMBS];\n").unwrap();
    fmt_limb_const(
        &mut o,
        "Montgomery form of `1` (`R mod p`).",
        "ONE",
        &to_limbs(&one, s),
    );
    fmt_limb_const(
        &mut o,
        "Montgomery form of `2^-1 mod p`.",
        "TWO_INV",
        &to_limbs(&two_inv, s),
    );
    fmt_limb_const(
        &mut o,
        "Montgomery form of `3^-1 mod p`.",
        "THREE_INV",
        &to_limbs(&three_inv, s),
    );
    fmt_limb_const(
        &mut o,
        "`R^2 mod p`: multiplying a canonical integer by this in Montgomery\n/// multiplication converts it to Montgomery form.",
        "R2_NRES",
        &to_limbs(&r2_nres, s),
    );
    fmt_limb_const(
        &mut o,
        &format!(
            "Montgomery form of `2^{}`: shifts the accumulator by one\n/// {chunk_bytes}-byte chunk in `fp_decode_reduce`.",
            8 * chunk_bytes
        ),
        "CHUNK_R",
        &to_limbs(&chunk_r, s),
    );

    // ---- carry handling ------------------------------------------------
    writeln!(
        o,
        "/// Propagate carries so limbs `0..{top}` are in `[0, 2^RADIX)` and the top\n\
         /// limb holds the overflow. Returns all-ones if the value is negative\n\
         /// (top bit of the top limb set), else zero.\n\
         #[inline]\n\
         pub(crate) fn prop(n: &mut [u64; NLIMBS]) -> u64 {{\n\
         \x20   let mut carry: i64 = n[0] as i64;\n\
         \x20   carry >>= RADIX;\n\
         \x20   n[0] &= MASK;\n\
         \x20   for limb in n.iter_mut().take({top}).skip(1) {{\n\
         \x20       carry = carry.wrapping_add(*limb as i64);\n\
         \x20       *limb = (carry as u64) & MASK;\n\
         \x20       carry >>= RADIX;\n\
         \x20   }}\n\
         \x20   n[{top}] = n[{top}].wrapping_add(carry as u64);\n\
         \x20   0u64.wrapping_sub(n[{top}] >> 63)\n\
         }}\n"
    )
    .unwrap();
    writeln!(
        o,
        "/// [`prop`], then add `p` back if the value went negative, then [`prop`]\n\
         /// again. Returns `1` if the fix-up happened.\n\
         #[inline]\n\
         pub(crate) fn flatten(n: &mut [u64; NLIMBS]) -> u32 {{\n\
         \x20   let carry = prop(n);\n\
         \x20   n[0] = n[0].wrapping_sub(1 & carry);\n\
         \x20   n[{top}] = n[{top}].wrapping_add(PTOP & carry);\n\
         \x20   let _ = prop(n);\n\
         \x20   (carry & 1) as u32\n\
         }}\n\n\
         /// Final subtraction of `p`: `n -= p`, undone if that went negative.\n\
         /// Returns `1` if the input was already `< p`.\n\
         #[inline]\n\
         pub(crate) fn modfsb(n: &mut [u64; NLIMBS]) -> u32 {{\n\
         \x20   n[0] = n[0].wrapping_add(1);\n\
         \x20   n[{top}] = n[{top}].wrapping_sub(PTOP);\n\
         \x20   flatten(n)\n\
         }}\n"
    )
    .unwrap();

    // ---- add / sub / neg -----------------------------------------------
    writeln!(o, "/// `n <- a + b mod 2p` (lazily reduced).\n#[inline]").unwrap();
    writeln!(
        o,
        "pub(crate) fn modadd(n: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{"
    )
    .unwrap();
    for i in 0..n {
        writeln!(o, "    n[{i}] = a[{i}].wrapping_add(b[{i}]);").unwrap();
    }
    writeln!(
        o,
        "    n[0] = n[0].wrapping_add(2);\n\
         \x20   n[{top}] = n[{top}].wrapping_sub(TWO_PTOP);\n\
         \x20   let carry = prop(n);\n\
         \x20   n[0] = n[0].wrapping_sub(2 & carry);\n\
         \x20   n[{top}] = n[{top}].wrapping_add(TWO_PTOP & carry);\n\
         \x20   let _ = prop(n);\n}}\n"
    )
    .unwrap();
    writeln!(o, "/// `n <- a - b mod 2p` (lazily reduced).\n#[inline]").unwrap();
    writeln!(
        o,
        "pub(crate) fn modsub(n: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{"
    )
    .unwrap();
    for i in 0..n {
        writeln!(o, "    n[{i}] = a[{i}].wrapping_sub(b[{i}]);").unwrap();
    }
    writeln!(
        o,
        "    let carry = prop(n);\n\
         \x20   n[0] = n[0].wrapping_sub(2 & carry);\n\
         \x20   n[{top}] = n[{top}].wrapping_add(TWO_PTOP & carry);\n\
         \x20   let _ = prop(n);\n}}\n"
    )
    .unwrap();
    writeln!(o, "/// `n <- -b mod 2p` (lazily reduced).\n#[inline]").unwrap();
    writeln!(
        o,
        "pub(crate) fn modneg(n: &mut [u64; NLIMBS], b: &[u64; NLIMBS]) {{"
    )
    .unwrap();
    for i in 0..n {
        writeln!(o, "    n[{i}] = 0u64.wrapping_sub(b[{i}]);").unwrap();
    }
    writeln!(
        o,
        "    let carry = prop(n);\n\
         \x20   n[0] = n[0].wrapping_sub(2 & carry);\n\
         \x20   n[{top}] = n[{top}].wrapping_add(TWO_PTOP & carry);\n\
         \x20   let _ = prop(n);\n}}\n"
    )
    .unwrap();

    // ---- modmul ----------------------------------------------------------
    writeln!(
        o,
        "/// `c <- a * b * R^-1 mod 2p`. Schoolbook {n}x{n} product in radix `2^{r}`\n\
         /// with interleaved Montgomery folding: each emitted low limb `v_i` is\n\
         /// added back at position `i + {top}` as `v_i * PTOP`, which is\n\
         /// `v_i * (p + 1) = v_i (mod p)` and divides the result by `R`.\n\
         #[inline]\n\
         pub(crate) fn modmul(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], b: &[u64; NLIMBS]) {{\n\
         \x20   let ptop = PTOP as u128;\n\
         \x20   let mut t: u128 = 0;"
    )
    .unwrap();
    for k in 0..(2 * n - 1) {
        let lo = k.saturating_sub(n - 1);
        let hi = k.min(n - 1);
        for i in lo..=hi {
            let j = k - i;
            writeln!(
                o,
                "    t = t.wrapping_add((a[{i}] as u128) * (b[{j}] as u128));"
            )
            .unwrap();
        }
        if k >= n - 1 {
            writeln!(
                o,
                "    t = t.wrapping_add((v{} as u128) * ptop);",
                k - (n - 1)
            )
            .unwrap();
        }
        if k < n {
            writeln!(o, "    let v{k} = (t as u64) & MASK;").unwrap();
        } else {
            writeln!(o, "    c[{}] = (t as u64) & MASK;", k - n).unwrap();
        }
        writeln!(o, "    t >>= RADIX;").unwrap();
    }
    writeln!(o, "    c[{top}] = t as u64;\n}}\n").unwrap();

    // ---- modsqr ----------------------------------------------------------
    writeln!(
        o,
        "/// `c <- a * a * R^-1 mod 2p`. Each cross term `a_i * a_j`, `i != j`, is\n\
         /// computed once and doubled.\n\
         #[inline]\n\
         pub(crate) fn modsqr(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS]) {{\n\
         \x20   let ptop = PTOP as u128;\n\
         \x20   let mut t: u128 = 0;\n\
         \x20   let mut tot: u128;"
    )
    .unwrap();
    for k in 0..(2 * n - 1) {
        let lo = k.saturating_sub(n - 1);
        let hi = k.min(n - 1);
        // cross terms i < j
        let mut cross = Vec::new();
        let mut square = None;
        for i in lo..=hi {
            let j = k - i;
            if i < j {
                cross.push((i, j));
            } else if i == j {
                square = Some(i);
            }
        }
        writeln!(o, "    tot = 0;").unwrap();
        for (i, j) in &cross {
            writeln!(
                o,
                "    tot = tot.wrapping_add((a[{i}] as u128) * (a[{j}] as u128));"
            )
            .unwrap();
        }
        if !cross.is_empty() {
            writeln!(o, "    tot = tot.wrapping_mul(2);").unwrap();
        }
        if let Some(i) = square {
            writeln!(
                o,
                "    tot = tot.wrapping_add((a[{i}] as u128) * (a[{i}] as u128));"
            )
            .unwrap();
        }
        writeln!(o, "    t = t.wrapping_add(tot);").unwrap();
        if k >= n - 1 {
            writeln!(
                o,
                "    t = t.wrapping_add((v{} as u128) * ptop);",
                k - (n - 1)
            )
            .unwrap();
        }
        if k < n {
            writeln!(o, "    let v{k} = (t as u64) & MASK;").unwrap();
        } else {
            writeln!(o, "    c[{}] = (t as u64) & MASK;", k - n).unwrap();
        }
        writeln!(o, "    t >>= RADIX;").unwrap();
    }
    writeln!(o, "    c[{top}] = t as u64;\n}}\n").unwrap();

    // ---- helpers -----------------------------------------------------------
    writeln!(
        o,
        "/// Square `a` in place `n` times.\n\
         #[inline]\n\
         pub(crate) fn modnsqr(a: &mut [u64; NLIMBS], n: u32) {{\n\
         \x20   for _ in 0..n {{\n\
         \x20       let mut tmp = [0u64; NLIMBS];\n\
         \x20       modsqr(&mut tmp, a);\n\
         \x20       *a = tmp;\n\
         \x20   }}\n\
         }}\n\n\
         /// `a <- a * b`.\n\
         #[inline]\n\
         fn mul_assign(a: &mut [u64; NLIMBS], b: &[u64; NLIMBS]) {{\n\
         \x20   let mut tmp = [0u64; NLIMBS];\n\
         \x20   modmul(&mut tmp, a, b);\n\
         \x20   *a = tmp;\n\
         }}\n"
    )
    .unwrap();

    // ---- modpro ------------------------------------------------------------
    let n_sqr: usize = chain
        .iter()
        .map(|op| if let ChainOp::Sqr(k) = op { *k } else { 0 })
        .sum();
    let n_mul = chain
        .iter()
        .filter(|op| matches!(op, ChainOp::Mul(_)))
        .count();
    writeln!(
        o,
        "/// Square-root progenitor: `z <- w^((p-3)/4) mod p`, computed with a\n\
         /// fixed {WINDOW}-bit sliding-window chain over the public exponent\n\
         /// ({n_sqr} squarings, {n_mul} multiplications after a 16-entry table).\n\
         /// Since `p = 3 mod 4`, `sqrt(w) = w * z` and `1/w = w * z^4`.\n\
         #[inline]\n\
         pub(crate) fn modpro(z: &mut [u64; NLIMBS], w: &[u64; NLIMBS]) {{\n\
         \x20   // tab[j] = w^(2j + 1)\n\
         \x20   let mut tab = [[0u64; NLIMBS]; {tab}];\n\
         \x20   let mut w2 = [0u64; NLIMBS];\n\
         \x20   modsqr(&mut w2, w);\n\
         \x20   tab[0] = *w;\n\
         \x20   for j in 1..{tab} {{\n\
         \x20       let prev = tab[j - 1];\n\
         \x20       modmul(&mut tab[j], &prev, &w2);\n\
         \x20   }}",
        tab = 1 << (WINDOW - 1)
    )
    .unwrap();
    let mut first = true;
    for op in &chain {
        match op {
            ChainOp::Mul(j) if first => {
                writeln!(o, "    *z = tab[{j}];").unwrap();
                first = false;
            }
            ChainOp::Mul(j) => writeln!(o, "    mul_assign(z, &tab[{j}]);").unwrap(),
            ChainOp::Sqr(k) => writeln!(o, "    modnsqr(z, {k});").unwrap(),
        }
    }
    writeln!(o, "}}\n").unwrap();

    // ---- inverse, conversions, predicates ------------------------------------
    writeln!(
        o,
        "/// `z <- 1 / x mod p` (`0` for `x = 0`). `h`, if given, is the precomputed\n\
         /// progenitor `x^((p-3)/4)`.\n\
         #[inline]\n\
         pub(crate) fn modinv(z: &mut [u64; NLIMBS], x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) {{\n\
         \x20   let mut t = [0u64; NLIMBS];\n\
         \x20   match h {{\n\
         \x20       None => modpro(&mut t, x),\n\
         \x20       Some(h) => t = *h,\n\
         \x20   }}\n\
         \x20   modnsqr(&mut t, 2);\n\
         \x20   modmul(z, x, &t);\n\
         }}\n\n\
         /// Canonical integer to Montgomery form.\n\
         #[inline]\n\
         pub(crate) fn nres(n: &mut [u64; NLIMBS], m: &[u64; NLIMBS]) {{\n\
         \x20   modmul(n, m, &R2_NRES);\n\
         }}\n\n\
         /// Montgomery form to canonical integer in `[0, p)`.\n\
         #[inline]\n\
         pub(crate) fn redc(m: &mut [u64; NLIMBS], n: &[u64; NLIMBS]) {{\n\
         \x20   let mut c = [0u64; NLIMBS];\n\
         \x20   c[0] = 1;\n\
         \x20   modmul(m, n, &c);\n\
         \x20   let _ = modfsb(m);\n\
         }}\n\n\
         /// `1` if `a == 0 mod p`, else `0`.\n\
         #[inline]\n\
         pub(crate) fn modis0(a: &[u64; NLIMBS]) -> u32 {{\n\
         \x20   let mut c = [0u64; NLIMBS];\n\
         \x20   redc(&mut c, a);\n\
         \x20   let mut d: u64 = 0;\n\
         \x20   for limb in c.iter() {{\n\
         \x20       d |= *limb;\n\
         \x20   }}\n\
         \x20   (1 & (d.wrapping_sub(1) >> RADIX)) as u32\n\
         }}\n\n\
         /// `1` if `a == 1 mod p`, else `0`.\n\
         #[inline]\n\
         pub(crate) fn modis1(a: &[u64; NLIMBS]) -> u32 {{\n\
         \x20   let mut c = [0u64; NLIMBS];\n\
         \x20   redc(&mut c, a);\n\
         \x20   let mut d: u64 = 0;\n\
         \x20   for limb in c.iter().skip(1) {{\n\
         \x20       d |= *limb;\n\
         \x20   }}\n\
         \x20   (1 & (d.wrapping_sub(1) >> RADIX) & ((c[0] ^ 1).wrapping_sub(1) >> RADIX)) as u32\n\
         }}\n\n\
         /// `1` if `a == b mod p`, else `0`.\n\
         #[inline]\n\
         pub(crate) fn modcmp(a: &[u64; NLIMBS], b: &[u64; NLIMBS]) -> u32 {{\n\
         \x20   let mut c = [0u64; NLIMBS];\n\
         \x20   let mut d = [0u64; NLIMBS];\n\
         \x20   redc(&mut c, a);\n\
         \x20   redc(&mut d, b);\n\
         \x20   let mut eq: u64 = 1;\n\
         \x20   for i in 0..NLIMBS {{\n\
         \x20       eq &= (c[i] ^ d[i]).wrapping_sub(1) >> RADIX;\n\
         \x20   }}\n\
         \x20   (eq & 1) as u32\n\
         }}\n\n\
         /// `1` if `x` is a square (or zero), else `0`. `h` is the optional\n\
         /// precomputed progenitor.\n\
         #[inline]\n\
         pub(crate) fn modqr(x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) -> u32 {{\n\
         \x20   let mut r = [0u64; NLIMBS];\n\
         \x20   match h {{\n\
         \x20       None => {{\n\
         \x20           let mut pro = [0u64; NLIMBS];\n\
         \x20           modpro(&mut pro, x);\n\
         \x20           modsqr(&mut r, &pro);\n\
         \x20       }}\n\
         \x20       Some(h) => modsqr(&mut r, h),\n\
         \x20   }}\n\
         \x20   mul_assign(&mut r, x);\n\
         \x20   modis1(&r) | modis0(x)\n\
         }}\n\n\
         /// `r <- sqrt(x)` for square `x` (`x^((p+1)/4)`), defined up to sign.\n\
         #[inline]\n\
         pub(crate) fn modsqrt(r: &mut [u64; NLIMBS], x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>) {{\n\
         \x20   let mut y = [0u64; NLIMBS];\n\
         \x20   match h {{\n\
         \x20       None => modpro(&mut y, x),\n\
         \x20       Some(h) => y = *h,\n\
         \x20   }}\n\
         \x20   modmul(r, &y, x);\n\
         }}\n\n\
         /// `a <- x` in Montgomery form.\n\
         #[inline]\n\
         pub(crate) fn modint(a: &mut [u64; NLIMBS], x: u64) {{\n\
         \x20   let mut m = [0u64; NLIMBS];\n\
         \x20   m[0] = x;\n\
         \x20   nres(a, &m);\n\
         }}\n\n\
         /// `c <- a * x mod 2p` for a small integer `x`.\n\
         #[inline]\n\
         pub(crate) fn modmli(c: &mut [u64; NLIMBS], a: &[u64; NLIMBS], x: u64) {{\n\
         \x20   let mut t = [0u64; NLIMBS];\n\
         \x20   modint(&mut t, x);\n\
         \x20   modmul(c, a, &t);\n\
         }}\n\n\
         /// Shift left by `n < RADIX` bits, carrying into the next limb.\n\
         #[inline]\n\
         pub(crate) fn modshl(a: &mut [u64; NLIMBS], n: u32) {{\n\
         \x20   a[{top}] = (a[{top}] << n) | (a[{top_m1}] >> (RADIX - n));\n\
         \x20   for i in (1..{top}).rev() {{\n\
         \x20       a[i] = ((a[i] << n) & MASK) | (a[i - 1] >> (RADIX - n));\n\
         \x20   }}\n\
         \x20   a[0] = (a[0] << n) & MASK;\n\
         }}\n\n\
         /// Shift right by `n < RADIX` bits, returning the bits shifted out.\n\
         #[inline]\n\
         pub(crate) fn modshr(a: &mut [u64; NLIMBS], n: u32) -> u64 {{\n\
         \x20   let r = a[0] & ((1u64 << n) - 1);\n\
         \x20   for i in 0..{top} {{\n\
         \x20       a[i] = (a[i] >> n) | ((a[i + 1] << (RADIX - n)) & MASK);\n\
         \x20   }}\n\
         \x20   a[{top}] >>= n;\n\
         \x20   r\n\
         }}\n\n\
         /// Constant-time conditional swap of `f` and `g` when `b == 1`.\n\
         #[inline]\n\
         pub(crate) fn modcsw(b: u64, g: &mut [u64; NLIMBS], f: &mut [u64; NLIMBS]) {{\n\
         \x20   let r: u64 = 0x3cc3_c33c_5aa5_a55a;\n\
         \x20   let c0 = (1u64.wrapping_sub(b)).wrapping_add(r);\n\
         \x20   let c1 = b.wrapping_add(r);\n\
         \x20   for i in 0..NLIMBS {{\n\
         \x20       let s = g[i];\n\
         \x20       let t = f[i];\n\
         \x20       let w = r.wrapping_mul(t.wrapping_add(s));\n\
         \x20       f[i] = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s)).wrapping_sub(w);\n\
         \x20       g[i] = c0.wrapping_mul(s).wrapping_add(c1.wrapping_mul(t)).wrapping_sub(w);\n\
         \x20   }}\n\
         }}\n"
    )
    .unwrap();

    // ---- encoding ------------------------------------------------------------
    writeln!(
        o,
        "/// Canonical little-endian encoding of `a` into `ENCODED_BYTES` bytes.\n\
         #[inline]\n\
         pub(crate) fn fp_encode(out: &mut [u8], a: &[u64; NLIMBS]) {{\n\
         \x20   debug_assert!(out.len() >= ENCODED_BYTES);\n\
         \x20   let mut c = [0u64; NLIMBS];\n\
         \x20   redc(&mut c, a);\n\
         \x20   for byte in out.iter_mut().take(ENCODED_BYTES) {{\n\
         \x20       *byte = (c[0] & 0xff) as u8;\n\
         \x20       let _ = modshr(&mut c, 8);\n\
         \x20   }}\n\
         }}\n\n\
         /// Decode `ENCODED_BYTES` canonical little-endian bytes. Returns all-ones\n\
         /// if the value was in `[0, p)`, else zero with `out` zeroed.\n\
         #[inline]\n\
         pub(crate) fn fp_decode(out: &mut [u64; NLIMBS], bytes: &[u8]) -> u32 {{\n\
         \x20   *out = [0; NLIMBS];\n\
         \x20   if bytes.len() < ENCODED_BYTES {{\n\
         \x20       return 0;\n\
         \x20   }}\n\
         \x20   for i in (0..ENCODED_BYTES).rev() {{\n\
         \x20       modshl(out, 8);\n\
         \x20       out[0] = out[0].wrapping_add(bytes[i] as u64);\n\
         \x20   }}\n\
         \x20   let ok = 0u64.wrapping_sub(modfsb(out) as u64);\n\
         \x20   let canonical = *out;\n\
         \x20   nres(out, &canonical);\n\
         \x20   for limb in out.iter_mut() {{\n\
         \x20       *limb &= ok;\n\
         \x20   }}\n\
         \x20   ok as u32\n\
         }}\n\n\
         /// Bytes consumed per step of [`fp_decode_reduce`]. One byte short of the\n\
         /// encoding length, so every chunk is below `p` and decodes without\n\
         /// reduction.\n\
         const CHUNK_BYTES: usize = ENCODED_BYTES - 1;\n\n\
         /// Decode a little-endian byte string of any length modulo `p`. Never\n\
         /// fails; used to map hash output into the field.\n\
         #[inline]\n\
         pub(crate) fn fp_decode_reduce(out: &mut [u64; NLIMBS], bytes: &[u8]) {{\n\
         \x20   *out = [0; NLIMBS];\n\
         \x20   let mut len = bytes.len();\n\
         \x20   let rem = len % CHUNK_BYTES;\n\
         \x20   if rem != 0 {{\n\
         \x20       let mut tmp = [0u8; ENCODED_BYTES];\n\
         \x20       tmp[..rem].copy_from_slice(&bytes[len - rem..]);\n\
         \x20       let _ = fp_decode(out, &tmp);\n\
         \x20       len -= rem;\n\
         \x20   }}\n\
         \x20   while len > 0 {{\n\
         \x20       len -= CHUNK_BYTES;\n\
         \x20       mul_assign(out, &CHUNK_R);\n\
         \x20       let mut tmp = [0u8; ENCODED_BYTES];\n\
         \x20       tmp[..CHUNK_BYTES].copy_from_slice(&bytes[len..len + CHUNK_BYTES]);\n\
         \x20       let mut chunk = [0u64; NLIMBS];\n\
         \x20       let _ = fp_decode(&mut chunk, &tmp);\n\
         \x20       let acc = *out;\n\
         \x20       modadd(out, &acc, &chunk);\n\
         \x20   }}\n\
         }}\n"
    )
    .unwrap();

    // ---- trait impl ------------------------------------------------------------
    let asm_gate = if s.asm {
        "#[cfg(not(target_arch = \"x86_64\"))]\n"
    } else {
        ""
    };
    writeln!(
        o,
        "type Limbs = Array<u64, <{upper} as crate::params::Prime>::FpLimbs>;\n\n\
         #[inline]\n\
         fn as_arr(a: &Limbs) -> &[u64; NLIMBS] {{\n\
         \x20   <&[u64; NLIMBS]>::try_from(&a[..]).expect(\"invariant: FpLimbs == NLIMBS\")\n\
         }}\n\n\
         #[inline]\n\
         fn as_arr_mut(a: &mut Limbs) -> &mut [u64; NLIMBS] {{\n\
         \x20   <&mut [u64; NLIMBS]>::try_from(&mut a[..]).expect(\"invariant: FpLimbs == NLIMBS\")\n\
         }}\n\n\
         {asm_gate}impl FpBackend for {upper} {{\n\
         \x20   #[inline]\n\
         \x20   fn set_zero(out: &mut Limbs) {{\n\
         \x20       *as_arr_mut(out) = ZERO;\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn set_one(out: &mut Limbs) {{\n\
         \x20       *as_arr_mut(out) = ONE;\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn set_small(out: &mut Limbs, val: u64) {{\n\
         \x20       modint(as_arr_mut(out), val);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn is_equal(a: &Limbs, b: &Limbs) -> Choice {{\n\
         \x20       Choice::from(modcmp(as_arr(a), as_arr(b)) as u8)\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn is_zero(a: &Limbs) -> Choice {{\n\
         \x20       Choice::from(modis0(as_arr(a)) as u8)\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn copy(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       *as_arr_mut(out) = *as_arr(a);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn add(out: &mut Limbs, a: &Limbs, b: &Limbs) {{\n\
         \x20       modadd(as_arr_mut(out), as_arr(a), as_arr(b));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn sub(out: &mut Limbs, a: &Limbs, b: &Limbs) {{\n\
         \x20       modsub(as_arr_mut(out), as_arr(a), as_arr(b));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn neg(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       modneg(as_arr_mut(out), as_arr(a));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn mul(out: &mut Limbs, a: &Limbs, b: &Limbs) {{\n\
         \x20       modmul(as_arr_mut(out), as_arr(a), as_arr(b));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn sqr(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       modsqr(as_arr_mut(out), as_arr(a));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn inv(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       let a = *as_arr(a);\n\
         \x20       modinv(as_arr_mut(out), &a, None);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn sqrt(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       let a = *as_arr(a);\n\
         \x20       modsqrt(as_arr_mut(out), &a, None);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn is_square(a: &Limbs) -> Choice {{\n\
         \x20       Choice::from(modqr(as_arr(a), None) as u8)\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn half(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       modmul(as_arr_mut(out), &TWO_INV, as_arr(a));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn div3(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       modmul(as_arr_mut(out), &THREE_INV, as_arr(a));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn exp3div4(out: &mut Limbs, a: &Limbs) {{\n\
         \x20       let a = *as_arr(a);\n\
         \x20       modpro(as_arr_mut(out), &a);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn mul_small(out: &mut Limbs, a: &Limbs, val: u32) {{\n\
         \x20       modmli(as_arr_mut(out), as_arr(a), val as u64);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn encode(out: &mut [u8], a: &Limbs) {{\n\
         \x20       fp_encode(out, as_arr(a));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn decode(out: &mut Limbs, bytes: &[u8]) -> Choice {{\n\
         \x20       Choice::from((fp_decode(as_arr_mut(out), bytes) & 1) as u8)\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn decode_reduce(out: &mut Limbs, bytes: &[u8]) {{\n\
         \x20       fp_decode_reduce(as_arr_mut(out), bytes);\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn cswap(a: &mut Limbs, b: &mut Limbs, ctl: Choice) {{\n\
         \x20       modcsw(ctl.unwrap_u8() as u64, as_arr_mut(a), as_arr_mut(b));\n\
         \x20   }}\n\n\
         \x20   #[inline]\n\
         \x20   fn select(out: &mut Limbs, a0: &Limbs, a1: &Limbs, ctl: Choice) {{\n\
         \x20       let cw = 0u64.wrapping_sub(ctl.unwrap_u8() as u64);\n\
         \x20       let a0 = as_arr(a0);\n\
         \x20       let a1 = as_arr(a1);\n\
         \x20       let out = as_arr_mut(out);\n\
         \x20       for i in 0..NLIMBS {{\n\
         \x20           out[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]));\n\
         \x20       }}\n\
         \x20   }}\n\
         }}"
    )
    .unwrap();

    o
}
