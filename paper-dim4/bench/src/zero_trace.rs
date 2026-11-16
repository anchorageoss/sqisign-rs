//! All-zero vectors in the compact verifier under tampering (feature
//! `instrument`; never part of a timing session):
//!
//!     zero_trace --level 1|3|5 --sigs <n> --random <n> --out <file>
//!
//! For `n` signatures of one key at the level, every tampered input of the
//! classes below is given to the verifier in its three modes (four points,
//! two points with the check at the gluing, two points without). For each
//! class the record has the verdicts, the number of runs in which some
//! vector of `sqisign_verify::hd::zero_trace::SITES` was all zero, the
//! totals per site, and every accepted input. The classes `top-bit` (with
//! or without the third scalar) and `negation` are the two malleabilities
//! the paper states; an acceptance in any other class is reported as
//! UNEXPECTED and the exit status is 2.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::self_contained::hd_verify_without_hint_check;
use sqisign_verify::hd::zero_trace::{reset, snapshot, SITES};
use sqisign_verify::hd::{
    canonical_hints, challenge_limbs, hd_verify_parsed_with, is_good_degree, parse_public_key,
    parse_signature, HdLevel, HdReject, HdSignature, HdVerifyOptions, KernelPoints,
    ParsedPublicKey, ParsedSignature, U6, W,
};
use sqisign_verify::{Fp2, P324_3, P500_27, P664_17};

const MODES: [(KernelPoints, &str); 3] = [
    (KernelPoints::Four, "four"),
    (KernelPoints::TwoWithImageCheck, "two+check"),
    (KernelPoints::TwoWithoutImageCheck, "two"),
];

/// The sites that are not an all-zero vector: counted apart.
const NOT_ZERO_VECTOR: [usize; 4] = [
    sqisign_verify::hd::zero_trace::XEQ_LEFT_Z0,
    sqisign_verify::hd::zero_trace::XEQ_RIGHT_Z0,
    sqisign_verify::hd::zero_trace::IMAGE_C3_NOT_ID,
    sqisign_verify::hd::zero_trace::D2_DUAL_NOT_PRODUCT,
];

#[derive(Default)]
struct Class {
    inputs: u64,
    runs: u64,
    verdicts: BTreeMap<String, u64>,
    accepted: u64,
    runs_with_zero: u64,
    accepted_with_zero: u64,
    site_hits: Vec<u64>,
    site_visits: Vec<u64>,
    accepted_list: Vec<String>,
}

struct Tally {
    classes: BTreeMap<&'static str, Class>,
    order: Vec<&'static str>,
}

impl Tally {
    fn class(&mut self, name: &'static str) -> &mut Class {
        if !self.classes.contains_key(name) {
            self.order.push(name);
            self.classes.insert(
                name,
                Class {
                    site_hits: vec![0; SITES.len()],
                    site_visits: vec![0; SITES.len()],
                    ..Class::default()
                },
            );
        }
        self.classes.get_mut(name).unwrap()
    }

    /// One tampered input, in the three modes.
    fn input<L: HdLevel>(
        &mut self,
        class: &'static str,
        what: &str,
        sig: &ParsedSignature<L>,
        pk: &ParsedPublicKey<L>,
        msg: &[u8],
        hint_check: bool,
    ) {
        self.class(class).inputs += 1;
        for (mode, mode_name) in MODES {
            let opts = HdVerifyOptions {
                kernel_points: mode,
            };
            reset();
            let verdict = if hint_check {
                hd_verify_parsed_with::<L>(sig, pk, msg, &opts).0
            } else {
                match challenge_limbs::<L>(sig, pk, msg) {
                    None => Err(HdReject::BadCurve),
                    Some(chal) => hd_verify_without_hint_check::<L>(
                        &HdSignature {
                            a_pk: pk.a_pk.clone(),
                            a_com: sig.a_com.clone(),
                            hint_pk_p: pk.hint_pk_p,
                            hint_pk_q: pk.hint_pk_q,
                            hint_com_p: sig.hint_com_p,
                            hint_com_q: sig.hint_com_q,
                            chal_limbs: &chal,
                            resp_a: sig.a,
                            resp_b: sig.b,
                            resp_c_or_d: sig.c_or_d,
                            q: sig.q,
                        },
                        &opts,
                    ),
                }
            };
            let snap = snapshot();
            let c = self.class(class);
            c.runs += 1;
            let mut zero = false;
            for (i, (visits, hits)) in snap.iter().enumerate() {
                c.site_visits[i] += visits;
                c.site_hits[i] += hits;
                if *hits > 0 && !NOT_ZERO_VECTOR.contains(&i) {
                    zero = true;
                }
            }
            c.runs_with_zero += zero as u64;
            let key = match verdict {
                Ok(()) => "ACCEPTED".to_string(),
                Err(e) => format!("{e:?}"),
            };
            *c.verdicts.entry(key).or_insert(0) += 1;
            if verdict.is_ok() {
                c.accepted += 1;
                c.accepted_with_zero += zero as u64;
                c.accepted_list
                    .push(format!("{what} [{mode_name}] zero vector: {zero}"));
            }
        }
    }

    /// A tampered byte string: the decoder first.
    fn bytes<L: HdLevel>(
        &mut self,
        class: &'static str,
        what: &str,
        sig: &[u8],
        pk: &[u8],
        msg: &[u8],
    ) {
        match (parse_signature::<L>(sig), parse_public_key::<L>(pk)) {
            (Ok(s), Ok(p)) => self.input::<L>(class, what, &s, &p, msg, true),
            _ => {
                let c = self.class(class);
                c.inputs += 1;
                c.runs += 1;
                *c.verdicts.entry("MalformedInput".into()).or_insert(0) += 1;
            }
        }
    }
}

struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn w(&mut self, r: u32) -> W {
        let l: Vec<u64> = (0..4).map(|_| self.next()).collect();
        W::from_limbs(&l) & W::mask(r)
    }
}

fn flip(bytes: &[u8], bit: usize) -> Vec<u8> {
    let mut t = bytes.to_vec();
    t[bit / 8] ^= 1 << (bit % 8);
    t
}

/// The good degrees `2^e - p` for primes `p = a1^2 + a2^2` with a chosen
/// `v2(a2) = m`: the first `a1` odd for which the verifier's own test
/// passes.
fn degree_with_m(e: u32, m: u32) -> Option<(U6, u64)> {
    let a2sq = U6::pow2(2 * m);
    (1u64..4000).step_by(2).find_map(|a1| {
        let p = a2sq.wrapping_add(&U6::from_u64(a1 * a1));
        let q = U6::pow2(e).wrapping_sub(&p);
        is_good_degree(e, &q).then_some((q, a1))
    })
}

fn with_canonical_com<L: HdLevel>(s: &ParsedSignature<L>) -> Option<ParsedSignature<L>> {
    let (hp, hq) = canonical_hints::<L>(&s.a_com)?;
    let mut t = s.clone();
    t.hint_com_p = hp;
    t.hint_com_q = hq;
    Some(t)
}

fn with_canonical_pk<L: HdLevel>(p: &ParsedPublicKey<L>) -> Option<ParsedPublicKey<L>> {
    let (hp, hq) = canonical_hints::<L>(&p.a_pk)?;
    let mut t = p.clone();
    t.hint_pk_p = hp;
    t.hint_pk_q = hq;
    Some(t)
}

fn run<L: HdLevel, const N: usize>(
    params: &Params<N>,
    label: &[u8],
    sigs: usize,
    random: usize,
) -> Tally {
    let mut t = Tally {
        classes: BTreeMap::new(),
        order: Vec::new(),
    };
    let r = L::R;
    let e = L::E_EMBED;
    let top = W::pow2(r - 1);
    let mask = W::mask(r);
    let neg = |x: W| x.wrapping_neg() & mask;
    let mut rng = ShakeRng::new(label, b"zero");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut rng).expect("keygen");
    let mut pkb = vec![0u8; compact_public_key_bytes::<L>()];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let ppk = parse_public_key::<L>(&pkb).unwrap();
    let mut xs = XorShift(0x9e3779b97f4a7c15 ^ (e as u64));
    let small = |v: u64| Fp2::<L>::from_small(v);
    let i = Fp2::<L>::i_element();

    for n in 0..sigs {
        let msg = [n as u8, 0x5a, 0xc3, L::R as u8];
        let mut sig = vec![0u8; compact_signature_bytes::<L>()];
        compact_sign::<L, N>(params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
        let ps = parse_signature::<L>(&sig).unwrap();
        eprintln!("signature {n}");

        // the signature itself: the reference for the counters
        t.input::<L>("honest", "the signature", &ps, &ppk, &msg, true);
        assert_eq!(t.class("honest").accepted, 3 * (n as u64 + 1));

        // every bit of the signature, of the key, and the message
        for bit in 0..8 * sig.len() {
            let what = format!("sig {n} bit {bit}");
            t.bytes::<L>("sig-bit", &what, &flip(&sig, bit), &pkb, &msg);
        }
        for bit in 0..8 * pkb.len() {
            let what = format!("sig {n} pk bit {bit}");
            t.bytes::<L>("pk-bit", &what, &sig, &flip(&pkb, bit), &msg);
        }
        for bit in 0..8 * msg.len() {
            let what = format!("sig {n} msg bit {bit}");
            t.bytes::<L>("msg-bit", &what, &sig, &pkb, &flip(&msg, bit));
        }
        t.bytes::<L>("msg-bit", "empty message", &sig, &pkb, &[]);

        // every bit of the commitment curve, with the canonical hints of
        // the new curve (the real verifier), and with the hint check off
        for bit in 0..2 * L::P_BITS as usize {
            let Ok(s) = parse_signature::<L>(&flip(&sig, bit)) else {
                continue;
            };
            let what = format!("sig {n} A_com bit {bit}");
            if let Some(c) = with_canonical_com(&s) {
                t.input::<L>("acom-bit+canonical-hints", &what, &c, &ppk, &msg, true);
            }
            t.input::<L>("acom-bit, hint check off", &what, &s, &ppk, &msg, false);
        }
        // every bit of the key's curve with its canonical hints
        for bit in 0..8 * L::FP2_BYTES {
            let Ok(p) = parse_public_key::<L>(&flip(&pkb, bit)) else {
                continue;
            };
            if let Some(c) = with_canonical_pk(&p) {
                let what = format!("sig {n} A_pk bit {bit}");
                t.input::<L>("apk-bit+canonical-hints", &what, &ps, &c, &msg, true);
            }
        }

        // the top bit of the scalars: the seven subsets
        for subset in 1..8u8 {
            let mut s = ps.clone();
            if subset & 1 != 0 {
                s.a = s.a ^ top;
            }
            if subset & 2 != 0 {
                s.b = s.b ^ top;
            }
            if subset & 4 != 0 {
                s.c_or_d = s.c_or_d ^ top;
            }
            let class = if subset & 4 != 0 {
                "top-bit, with the third scalar"
            } else {
                "top-bit"
            };
            let what = format!("sig {n} top bits {subset:03b} (c_or_d, b, a)");
            t.input::<L>(class, &what, &s, &ppk, &msg, true);
        }
        // the negations
        for subset in 1..8u8 {
            let mut s = ps.clone();
            if subset & 1 != 0 {
                s.a = neg(s.a);
            }
            if subset & 2 != 0 {
                s.b = neg(s.b);
            }
            if subset & 4 != 0 {
                s.c_or_d = neg(s.c_or_d);
            }
            let class = if subset == 7 {
                "negation"
            } else {
                "negation of some scalars"
            };
            let what = format!("sig {n} negated {subset:03b} (c_or_d, b, a)");
            t.input::<L>(class, &what, &s, &ppk, &msg, true);
        }

        // a and b both even
        let even = |x: W| x & (mask ^ W::ONE);
        let both_even: [(W, W, &str); 6] = [
            (even(ps.a), even(ps.b), "low bits cleared"),
            (W::from_u64(0), W::from_u64(0), "a = b = 0"),
            (W::from_u64(0), W::from_u64(2), "a = 0, b = 2"),
            (W::from_u64(2), W::from_u64(0), "a = 2, b = 0"),
            (top, top, "a = b = 2^(r-1)"),
            (even(ps.a), W::from_u64(0), "a even, b = 0"),
        ];
        for (a, b, name) in both_even {
            for (cd, cd_name) in [
                (ps.c_or_d, "kept"),
                (W::from_u64(0), "0"),
                (W::from_u64(1), "1"),
            ] {
                let mut s = ps.clone();
                (s.a, s.b, s.c_or_d) = (a, b, cd);
                let what = format!("sig {n} {name}, c_or_d {cd_name}");
                t.input::<L>("both-even", &what, &s, &ppk, &msg, true);
            }
        }

        // scalars with zeros and ones
        let z = W::from_u64(0);
        let o = W::from_u64(1);
        let simple: [(W, W, W, &str); 14] = [
            (z, ps.b, ps.c_or_d, "a = 0"),
            (ps.a, z, ps.c_or_d, "b = 0"),
            (ps.a, ps.b, z, "c_or_d = 0"),
            (o, z, z, "(1, 0, 0)"),
            (z, o, z, "(0, 1, 0)"),
            (o, z, o, "(1, 0, 1)"),
            (z, o, o, "(0, 1, 1)"),
            (o, o, z, "(1, 1, 0)"),
            (o, o, o, "(1, 1, 1)"),
            (mask, mask, mask, "all ones"),
            (o, z, ps.c_or_d, "(1, 0, c_or_d)"),
            (z, o, ps.c_or_d, "(0, 1, c_or_d)"),
            (ps.b, ps.a, ps.c_or_d, "a and b swapped"),
            (ps.a, ps.b, ps.a, "c_or_d = a"),
        ];
        for (a, b, cd, name) in simple {
            let mut s = ps.clone();
            (s.a, s.b, s.c_or_d) = (a, b, cd);
            t.input::<L>(
                "zero-and-one scalars",
                &format!("sig {n} {name}"),
                &s,
                &ppk,
                &msg,
                true,
            );
        }

        // other good degrees: the neighbours, and 2^e - p for the small
        // primes and for p = a1^2 + 2^(2m) over the range of m and past it
        let mut degrees: Vec<(U6, String)> = Vec::new();
        for dir in [true, false] {
            let mut q = ps.q;
            let mut found = 0;
            while found < 3 {
                q = if dir {
                    q.wrapping_add(&U6::from_u64(4))
                } else {
                    q.wrapping_sub(&U6::from_u64(4))
                };
                if is_good_degree(e, &q) {
                    found += 1;
                    degrees.push((q, format!("good neighbour {found} {dir}")));
                }
            }
        }
        for p in [5u64, 13, 17, 29, 37, 41, 257, 65537] {
            let q = U6::pow2(e).wrapping_sub(&U6::from_u64(p));
            degrees.push((q, format!("2^e - {p}")));
        }
        for m in [1u32, 2, 3, 10, 20, 30, 38, 39, 40, 41, 42, 50, 60] {
            if 2 * m + 1 < e {
                if let Some((q, a1)) = degree_with_m(e, m) {
                    degrees.push((q, format!("2^e - ({a1}^2 + 2^{})", 2 * m)));
                }
            }
        }
        degrees.push((U6::from_u64(3), "q = 3".into()));
        degrees.push((U6::from_u64(7), "q = 7".into()));
        for (q, name) in &degrees {
            let mut s = ps.clone();
            s.q = *q;
            let what = format!("sig {n} q: {name}");
            t.input::<L>("other degree", &what, &s, &ppk, &msg, true);
            // and with scalars that are not the signature's
            s.a = o;
            s.b = z;
            s.c_or_d = z;
            t.input::<L>(
                "other degree, scalars (1, 0, 0)",
                &what,
                &s,
                &ppk,
                &msg,
                true,
            );
            s.a = xs.w(r) | W::ONE;
            s.b = xs.w(r);
            s.c_or_d = xs.w(r);
            t.input::<L>("other degree, random scalars", &what, &s, &ppk, &msg, true);
        }

        // random scalars under the signature's degree
        for k in 0..random {
            let mut s = ps.clone();
            s.a = xs.w(r);
            s.b = xs.w(r);
            s.c_or_d = xs.w(r);
            if k % 3 == 1 {
                s.a = ps.a;
            }
            if k % 3 == 2 {
                s.c_or_d = ps.c_or_d;
            }
            let what = format!("sig {n} random scalars {k}");
            t.input::<L>("random scalars", &what, &s, &ppk, &msg, true);
        }

        // crafted curves, each with its canonical hints
        let curves: Vec<(Fp2<L>, &str)> = vec![
            (ppk.a_pk.clone(), "A = A_pk (as the commitment only)"),
            (small(0), "A = 0"),
            (small(6), "A = 6"),
            (small(6).neg(), "A = -6"),
            (small(2), "A = 2"),
            (small(2).neg(), "A = -2"),
            (small(1), "A = 1"),
            (small(3), "A = 3"),
            (i.clone(), "A = i"),
            (i.mul(&small(6)), "A = 6 i"),
            (ppk.a_pk.neg(), "A = -A_pk"),
            (ppk.a_pk.conjugate(), "A = conj(A_pk)"),
            (ps.a_com.neg(), "A = -A_com"),
            (ps.a_com.conjugate(), "A = conj(A_com)"),
        ];
        for (a, name) in &curves {
            let mut s = ps.clone();
            s.a_com = a.clone();
            let what = format!("sig {n} commitment {name}");
            match with_canonical_com(&s) {
                Some(c) => {
                    t.input::<L>("crafted commitment curve", &what, &c, &ppk, &msg, true);
                    let mut c1 = c.clone();
                    (c1.a, c1.b, c1.c_or_d) = (o, z, z);
                    t.input::<L>("crafted commitment curve", &what, &c1, &ppk, &msg, true);
                    let mut c2 = c.clone();
                    (c2.a, c2.b, c2.c_or_d) = (xs.w(r) | W::ONE, xs.w(r), xs.w(r));
                    t.input::<L>("crafted commitment curve", &what, &c2, &ppk, &msg, true);
                }
                None => t.input::<L>("crafted commitment curve", &what, &s, &ppk, &msg, true),
            }
            let mut p = ppk.clone();
            p.a_pk = a.clone();
            if *name == "A = -A_com" || *name == "A = conj(A_com)" {
                // as a key: the commitment curve itself and its conjugate
                p.a_pk = if *name == "A = -A_com" {
                    ps.a_com.clone()
                } else {
                    ps.a_com.conjugate()
                };
            }
            let what = format!("sig {n} key {name}");
            if bool::from(p.a_pk.ct_equal(&ppk.a_pk)) {
                continue;
            }
            match with_canonical_pk(&p) {
                Some(c) => t.input::<L>("crafted key curve", &what, &ps, &c, &msg, true),
                None => t.input::<L>("crafted key curve", &what, &ps, &p, &msg, true),
            }
        }
        // the commitment curve with a zero coordinate
        for (a, name) in [
            (
                Fp2::<L>::decode(&{
                    let mut enc = ps.a_com.encode().to_vec();
                    let half = enc.len() / 2;
                    enc[half..].iter_mut().for_each(|b| *b = 0);
                    enc
                }),
                "Im(A_com) = 0",
            ),
            (
                Fp2::<L>::decode(&{
                    let mut enc = ps.a_com.encode().to_vec();
                    let half = enc.len() / 2;
                    enc[..half].iter_mut().for_each(|b| *b = 0);
                    enc
                }),
                "Re(A_com) = 0",
            ),
        ] {
            let mut s = ps.clone();
            s.a_com = a.expect("a field element");
            let what = format!("sig {n} commitment {name}");
            let c = with_canonical_com(&s).unwrap_or(s);
            t.input::<L>("crafted commitment curve", &what, &c, &ppk, &msg, true);
        }
    }
    t
}

fn report(t: &Tally, level: u32, sigs: usize) -> (String, u64) {
    let mut s = String::new();
    let mut unexpected = 0;
    writeln!(s, "level {level}: {sigs} signatures of one key; each tampered input in the three verifier modes (a run is one input in one mode; an input the decoder rejects is one run)").unwrap();
    for name in &t.order {
        let c = &t.classes[name];
        writeln!(s, "\n[{name}]").unwrap();
        writeln!(
            s,
            "  inputs {}, runs {}, accepted runs {}, runs with an all-zero vector {}, accepted runs with an all-zero vector {}",
            c.inputs, c.runs, c.accepted, c.runs_with_zero, c.accepted_with_zero
        )
        .unwrap();
        let v: Vec<String> = c.verdicts.iter().map(|(k, n)| format!("{k} {n}")).collect();
        writeln!(s, "  verdicts: {}", v.join(", ")).unwrap();
        for (i, site) in SITES.iter().enumerate() {
            if c.site_hits[i] > 0 {
                writeln!(
                    s,
                    "  site hit: {site}: {} of {} visits",
                    c.site_hits[i], c.site_visits[i]
                )
                .unwrap();
            }
        }
        let expected = matches!(
            *name,
            "honest" | "top-bit" | "top-bit, with the third scalar" | "negation"
        );
        if c.accepted > 0 && !expected {
            unexpected += c.accepted;
            writeln!(s, "  UNEXPECTED ACCEPTANCE").unwrap();
        }
        if *name != "honest" {
            for a in &c.accepted_list {
                writeln!(s, "  accepted: {a}").unwrap();
            }
        }
    }
    let total = |f: &dyn Fn(&Class) -> u64| -> u64 {
        t.order
            .iter()
            .filter(|n| **n != "honest")
            .map(|n| f(&t.classes[n]))
            .sum()
    };
    writeln!(
        s,
        "\ntotal over the tampered classes: inputs {}, runs {}, accepted runs {} (unexpected {unexpected}), runs with an all-zero vector {}, accepted runs with an all-zero vector {}",
        total(&|c| c.inputs),
        total(&|c| c.runs),
        total(&|c| c.accepted),
        total(&|c| c.runs_with_zero),
        total(&|c| c.accepted_with_zero),
    )
    .unwrap();
    (s, unexpected)
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let (mut level, mut sigs, mut random, mut out) = (1u32, 1usize, 20usize, String::new());
    let mut i = 0;
    while i < argv.len() {
        let v = argv.get(i + 1).cloned().unwrap_or_default();
        match argv[i].as_str() {
            "--level" => level = v.parse().unwrap(),
            "--sigs" => sigs = v.parse().unwrap(),
            "--random" => random = v.parse().unwrap(),
            "--out" => out = v,
            s => panic!("unknown argument {s}"),
        }
        i += 2;
    }
    let t = match level {
        1 => run::<P324_3, { sqisign_rs::precomp::p324_3::IBZ_NLIMBS }>(
            &level1(),
            b"zero trace I",
            sigs,
            random,
        ),
        3 => run::<P500_27, { sqisign_rs::precomp::p500_27::IBZ_NLIMBS }>(
            &level3(),
            b"zero trace III",
            sigs,
            random,
        ),
        5 => run::<P664_17, { sqisign_rs::precomp::p664_17::IBZ_NLIMBS }>(
            &level5(),
            b"zero trace V",
            sigs,
            random,
        ),
        _ => panic!("level 1, 3 or 5"),
    };
    let (text, unexpected) = report(&t, level, sigs);
    print!("{text}");
    if !out.is_empty() {
        fs::write(&out, &text).unwrap();
    }
    if unexpected > 0 {
        std::process::exit(2);
    }
}
