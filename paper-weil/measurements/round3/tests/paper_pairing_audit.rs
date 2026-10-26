//! Numerical experiments for the Weil pairing compression note on the 300
//! official NIST round-3 KAT vectors (`kat/*.rsp` of the-sqisign at tag
//! `nist-v3`, 100 per level), run on sqisign-rs 0.7.
//!
//! Task 1 (the theorem checks of the paper's "What was checked"): the
//! verifier's canonical bases of `E_chl[2^N]` and `E_aux[2^N]`, `N = e_rsp
//! + 2`, the Weil pairings `omega = e(A, B)` and `eta = e(C, D)`, the
//! discrete logarithm `y = log_omega(eta^-1)` in `mu_{2^N}`, and against
//! `det M` from the four transmitted entries: `y + 2^(g-1) == det M mod
//! 2^g` with `g = N - 1`, the first bit at which the two differ,
//! `omega^det eta` a primitive fourth root of unity and which of the two
//! (`u = 1` or `3`), `v_2(det M)`, the parities of the entries, and that the
//! basis hints recomputed from the curves equal the transmitted ones.
//!
//! Task 2 (the Tate sweep of the paper's "Pairing selection"): for every
//! `m` in `[N, F]`, one reduced Tate ladder of `2^m` steps per curve, (a) on
//! the verifier's order-`2^N` bases and (c) on the full-order deterministic
//! bases (`TorsionBasisToHint(E, F, N)`, whose `[2^(F-N)]` multiple is the
//! verifier's basis) scaled to order `2^m`; the 2-adic orders of the values,
//! the discrete logarithm `y_m` in the subgroup the challenge value
//! generates, and whether `y_m + 2^(g-1) == det M mod 2^g` and whether `y_m
//! == det M mod 2^(g-1) = 2^e_rsp`; then the smallest `m` at which every
//! vector of a level passes, against the bound `ceil((F + g) / 2)`.
//!
//! Task 3 (end to end): the decompressor run with `det M` taken from one
//! Tate ladder per curve on the scaled full-order bases instead of the Weil
//! pairings, at `m = bound - 1` and `m = bound`, the reconstructed
//! signature handed to the standard verifier; the pass counts, and at
//! `bound - 1` whether the recovered dropped entry differs from the
//! signer's in bit `e_rsp` only.
//!
//! Timing (level I, KAT entry 0): the verifier's pairing stage, the bases
//! from the hints at order `2^N` and two Weil pairings, against the Tate
//! alternative, the bases at order `2^m*` and one Tate ladder per curve,
//! `m* = PAPER_TATE_M` (default the bound). Medians of 200 iterations after
//! 5 warm-ups, each on a copy of its input state.
//!
//! Run in release from a copy in `crates/kat/tests` (see the README):
//!   cargo test -p sqisign-kat --release --test paper_pairing_audit -- --nocapture --test-threads=1

#![allow(clippy::too_many_arguments)]

use rayon::prelude::*;
use sqisign_verify::ec::basis::{
    difference_point, ec_curve_to_basis_2f_from_hint, ec_curve_to_basis_2f_to_hint,
};
use sqisign_verify::ec::isogeny::iso_isogeny_2chain_from;
use sqisign_verify::ec::pairing::{fp2_dlog_2e, reduced_tate, weil};
use sqisign_verify::ec::point::{ec_dbl_iter, ec_dbl_iter_basis, ec_ladder3pt};
use sqisign_verify::ec::{EcBasis, EcCurve, MAX_ORDER_WORDS};
use sqisign_verify::fp::Fp2;
use sqisign_verify::precomp::EXTRA_TORSION;
use sqisign_verify::compressed::compress;
use sqisign_verify::sqisign::{
    public_key_from_bytes, signature_from_bytes, verify, Scalar, Signature, HD_EXTRA_TORSION,
};
use sqisign_verify::{Level, Level1, Level3, Level5, Prime};
use std::collections::BTreeMap;
use std::hint::black_box;
use std::time::Instant;

const KAT_L1: &str = include_str!("../kat/PQCsignKAT_270_SQIsign_p324_3.rsp");
const KAT_L3: &str = include_str!("../kat/PQCsignKAT_417_SQIsign_p500_27.rsp");
const KAT_L5: &str = include_str!("../kat/PQCsignKAT_549_SQIsign_p664_17.rsp");

// ---------------------------------------------------------------------------
// KAT parsing
// ---------------------------------------------------------------------------

struct Entry {
    count: usize,
    msg: Vec<u8>,
    pk: Vec<u8>,
    sm: Vec<u8>,
}

fn parse_rsp(text: &str) -> Vec<Entry> {
    let mut out = Vec::new();
    let mut cur: Option<Entry> = None;
    for line in text.lines() {
        let Some((k, v)) = line.split_once('=') else {
            continue;
        };
        let (k, v) = (k.trim(), v.trim());
        if k == "count" {
            out.extend(cur.take());
            cur = Some(Entry {
                count: v.parse().unwrap(),
                msg: vec![],
                pk: vec![],
                sm: vec![],
            });
            continue;
        }
        let e = cur.as_mut().expect("count first");
        match k {
            "msg" => e.msg = hex::decode(v).unwrap(),
            "pk" => e.pk = hex::decode(v).unwrap(),
            "sm" => e.sm = hex::decode(v).unwrap(),
            _ => {}
        }
    }
    out.extend(cur.take());
    out
}

// ---------------------------------------------------------------------------
// Multiprecision on 16-limb little-endian two's complement integers
// ---------------------------------------------------------------------------

const W: usize = 16;
type Wide = [u64; W];

fn wide(s: &[u64]) -> Wide {
    let mut out = [0u64; W];
    out[..s.len()].copy_from_slice(s);
    out
}

fn mul_wide(a: &[u64], b: &[u64]) -> Wide {
    let mut out = [0u64; W];
    for i in 0..a.len() {
        let mut carry: u128 = 0;
        for j in 0..b.len() {
            if i + j >= W {
                break;
            }
            let p = (a[i] as u128) * (b[j] as u128) + (out[i + j] as u128) + carry;
            out[i + j] = p as u64;
            carry = p >> 64;
        }
        if i + b.len() < W {
            out[i + b.len()] = carry as u64;
        }
    }
    out
}

fn sub_wide(a: &Wide, b: &Wide) -> Wide {
    let mut out = [0u64; W];
    let mut borrow = 0u64;
    for i in 0..W {
        let (d1, b1) = a[i].overflowing_sub(b[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        out[i] = d2;
        borrow = (b1 | b2) as u64;
    }
    out
}

fn add_wide(a: &Wide, b: &Wide) -> Wide {
    let mut out = [0u64; W];
    let mut carry = 0u128;
    for i in 0..W {
        let s = a[i] as u128 + b[i] as u128 + carry;
        out[i] = s as u64;
        carry = s >> 64;
    }
    out
}

fn mask(a: &mut Wide, nbits: usize) {
    let q = nbits / 64;
    let r = nbits % 64;
    if q < W {
        if r != 0 {
            a[q] &= (1u64 << r) - 1;
        } else {
            a[q] = 0;
        }
        for limb in a[q + 1..].iter_mut() {
            *limb = 0;
        }
    }
}

fn masked(a: &Wide, nbits: usize) -> Wide {
    let mut out = *a;
    mask(&mut out, nbits);
    out
}

fn pow2(k: usize) -> Wide {
    let mut out = [0u64; W];
    out[k / 64] = 1u64 << (k % 64);
    out
}

fn lowest_set_bit(a: &[u64]) -> Option<usize> {
    a.iter()
        .enumerate()
        .find(|(_, l)| **l != 0)
        .map(|(i, l)| 64 * i + l.trailing_zeros() as usize)
}

/// The first bit below `nbits` at which `a` and `b` differ.
fn first_mismatch(a: &Wide, b: &Wide, nbits: usize) -> Option<usize> {
    let mut d = [0u64; W];
    for i in 0..W {
        d[i] = a[i] ^ b[i];
    }
    mask(&mut d, nbits);
    lowest_set_bit(&d)
}

fn agree(a: &Wide, b: &Wide, nbits: usize) -> bool {
    first_mismatch(a, b, nbits).is_none()
}

/// `det M = M00 M11 - M01 M10` from the four entries, not reduced.
fn ground_truth_det<L: Level>(sig: &Signature<L>) -> Wide {
    let m = &sig.mat_bchall_can_to_bchall;
    let w = L::ORDER_WORDS;
    let p1 = mul_wide(&m[0][0][..w], &m[1][1][..w]);
    let p2 = mul_wide(&m[0][1][..w], &m[1][0][..w]);
    sub_wide(&p1, &p2)
}

// ---------------------------------------------------------------------------
// Fp2 helpers
// ---------------------------------------------------------------------------

fn fp2_eq<L: Level>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

fn fp2_is_one<L: Level>(a: &Fp2<L>) -> bool {
    bool::from(a.ct_is_one())
}

/// Smallest `a` with `x^(2^a) == 1`, or `None` within `cap` squarings.
fn ord_2adic<L: Level>(x: &Fp2<L>, cap: usize) -> Option<usize> {
    if fp2_is_one(x) {
        return Some(0);
    }
    let mut t = x.clone();
    for k in 1..=cap {
        t = t.sqr();
        if fp2_is_one(&t) {
            return Some(k);
        }
    }
    None
}

fn pow2k<L: Level>(x: &Fp2<L>, k: usize) -> Fp2<L> {
    let mut r = x.clone();
    for _ in 0..k {
        r = r.sqr();
    }
    r
}

/// `y` with `arg = base^y` in the subgroup of order `2^e`, as the
/// decompressor calls it (`fp2_dlog_2e(y, arg, base^-1, e)`).
fn dlog<L: Level>(arg: &Fp2<L>, base: &Fp2<L>, e: usize) -> Option<Wide> {
    let mut y = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e(&mut y[..L::ORDER_WORDS], arg, &base.inv(), e as u32)?;
    Some(masked(&wide(&y), e))
}

fn ord_str(o: Option<usize>) -> String {
    match o {
        Some(k) => format!("2^{k}"),
        None => "not-2-power".to_string(),
    }
}

// ---------------------------------------------------------------------------
// Per-vector setup: what the verifier computes before the pairings
// ---------------------------------------------------------------------------

struct Setup<L: Level> {
    idx: usize,
    sig: Signature<L>,
    det: Wide,
    e_chall: EcCurve<L>,
    /// The verifier's basis of `E_chl[2^N]` (`TorsionBasisFromHint`).
    b_chall: EcBasis<L>,
    /// `TorsionBasisToHint(E_chl, F, N)`: order `2^F`, `[2^(F-N)]` of it is `b_chall`.
    b_chall_full: EcBasis<L>,
    e_aux: EcCurve<L>,
    b_aux: EcBasis<L>,
    b_aux_full: EcBasis<L>,
    hint_chall_ok: bool,
    hint_aux_ok: bool,
    scaled_chall_ok: bool,
    scaled_aux_ok: bool,
}

fn basis_eq<L: Level>(a: &EcBasis<L>, b: &EcBasis<L>) -> bool {
    bool::from(a.p.ct_equal(&b.p) & a.q.ct_equal(&b.q) & a.pmq.ct_equal(&b.pmq))
}

/// A curve with its verifier basis at order `2^n` and the full-order basis
/// adjusted to it (`ec_curve_to_basis_2f_to_hint(curve, F, n)`).
fn bases<L: Level>(curve: &mut EcCurve<L>, n: u32, hint: u8) -> (EcBasis<L>, EcBasis<L>, bool, bool) {
    let fc = L::TWO_ADIC_EXPONENT;
    let b = ec_curve_to_basis_2f_from_hint(curve, n, hint).expect("basis from hint");
    let (full, h) = ec_curve_to_basis_2f_to_hint(curve, fc, n).expect("basis to hint");
    let scaled = ec_dbl_iter_basis(&full, (fc - n) as usize, curve);
    let scaled_ok = basis_eq(&scaled, &b);
    (b, full, h == hint, scaled_ok)
}

fn setup<L: Level>(idx: usize, e: &Entry) -> Setup<L> {
    assert_eq!(e.count, idx, "KAT entries in order");
    let params = &L::PARAMS;
    let pk = public_key_from_bytes::<L>(params, &e.pk).expect("pk decode");
    let sig = signature_from_bytes::<L>(params, &e.sm[..params.signature_bytes]).expect("sig decode");
    let n = params.response_bits + HD_EXTRA_TORSION;

    // E_chl exactly as the verifier derives it (sqisign::challenge_state)
    let f_chall = params.challenge_bits + EXTRA_TORSION;
    let mut e_pk = pk.curve.clone();
    let bas = ec_curve_to_basis_2f_from_hint(&mut e_pk, f_chall, pk.hint_pk).expect("pk basis");
    let kernel = ec_ladder3pt(
        &sig.chall_coeff[..L::ORDER_WORDS],
        &bas.p,
        &bas.q,
        &bas.pmq,
        &e_pk,
    )
    .expect("challenge kernel");
    let mut e_chall = iso_isogeny_2chain_from(&e_pk, &kernel, f_chall).expect("challenge chain");
    let (b_chall, b_chall_full, hint_chall_ok, scaled_chall_ok) = bases(&mut e_chall, n, sig.hint_chall);

    let mut e_aux = EcCurve::<L>::from_a(&sig.e_aux_a).expect("aux curve");
    let (b_aux, b_aux_full, hint_aux_ok, scaled_aux_ok) = bases(&mut e_aux, n, sig.hint_aux);

    let det = ground_truth_det(&sig);
    Setup {
        idx,
        sig,
        det,
        e_chall,
        b_chall,
        b_chall_full,
        e_aux,
        b_aux,
        b_aux_full,
        hint_chall_ok,
        hint_aux_ok,
        scaled_chall_ok,
        scaled_aux_ok,
    }
}

fn header<L: Level>(label: &str, task: &str) {
    let params = &L::PARAMS;
    let n = params.response_bits + HD_EXTRA_TORSION;
    println!("==================== {task}: {label} ====================");
    println!(
        "p = {} 2^{} - 1, lambda = {}, e_rsp = RESPONSE_BITS = {}, N = e_rsp + 2 = {}, g = N - 1 = {}, F = {}, field kernels: {}",
        L::COFACTOR,
        L::TWO_ADIC_EXPONENT,
        params.challenge_bits,
        params.response_bits,
        n,
        n - 1,
        L::TWO_ADIC_EXPONENT,
        sqisign_kat::arithmetic_from_env()
    );
}

// ---------------------------------------------------------------------------
// Task 1: the theorem checks
// ---------------------------------------------------------------------------

struct T1 {
    idx: usize,
    hints_ok: bool,
    scaled_ok: bool,
    ord_omega: Option<usize>,
    ord_eta: Option<usize>,
    dlog_ok: bool,
    /// `y + 2^(g-1) == det M mod 2^g`
    agree_g: bool,
    /// first bit below N at which `y + 2^(g-1)` and `det M` differ; N if none
    k: usize,
    /// 0: `z == 1`; 1: `z == -1`; 2: `z^2 == -1`; 3: none
    z_class: u8,
    /// `zeta = omega^(2^(N-2) u)`: 1, 3, or 0 if neither
    u: u8,
    v2_det: Option<usize>,
    m00_even: bool,
    m01_even: bool,
    m10_even: bool,
    m11_even: bool,
}

fn task1_vector<L: Level>(s: &Setup<L>) -> T1 {
    let params = &L::PARAMS;
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let g = n - 1;
    let mut e_chall = s.e_chall.clone();
    let mut e_aux = s.e_aux.clone();

    // the decompressor's pairings and discrete logarithm (compressed::pairing_dlog)
    let omega = weil::<L>(n as u32, &s.b_chall.p, &s.b_chall.q, &s.b_chall.pmq, &mut e_chall);
    let eta = weil::<L>(n as u32, &s.b_aux.p, &s.b_aux.q, &s.b_aux.pmq, &mut e_aux);
    let y = dlog::<L>(&eta.inv(), &omega, n);
    let det_n = masked(&s.det, n);
    let (dlog_ok, agree_g, k) = match y {
        Some(y) => {
            let yp = masked(&add_wide(&y, &pow2(g - 1)), n);
            (true, agree(&yp, &det_n, g), first_mismatch(&yp, &det_n, n).unwrap_or(n))
        }
        None => (false, false, 0),
    };

    // z = omega^(det M mod 2^N) eta; zeta = omega^(2^(N-2) u)
    let z = omega.pow_vartime(&det_n[..]).mul(&eta);
    let minus_one = Fp2::<L>::one().neg();
    let z_class = if fp2_is_one(&z) {
        0
    } else if fp2_eq(&z, &minus_one) {
        1
    } else if fp2_eq(&z.sqr(), &minus_one) {
        2
    } else {
        3
    };
    let i0 = pow2k(&omega, n - 2);
    let u = if fp2_eq(&z, &i0) {
        1
    } else if fp2_eq(&z, &i0.inv()) {
        3
    } else {
        0
    };

    let m = &s.sig.mat_bchall_can_to_bchall;
    let even = |x: &[u64]| x[0] & 1 == 0;
    T1 {
        idx: s.idx,
        hints_ok: s.hint_chall_ok && s.hint_aux_ok,
        scaled_ok: s.scaled_chall_ok && s.scaled_aux_ok,
        ord_omega: ord_2adic(&omega, n + 2),
        ord_eta: ord_2adic(&eta, n + 2),
        dlog_ok,
        agree_g,
        k,
        z_class,
        u,
        v2_det: lowest_set_bit(&s.det),
        m00_even: even(&m[0][0]),
        m01_even: even(&m[0][1]),
        m10_even: even(&m[1][0]),
        m11_even: even(&m[1][1]),
    }
}

fn run_task1<L: Level>(label: &str, kat: &str) {
    let entries = parse_rsp(kat);
    assert_eq!(entries.len(), 100, "KAT entries");
    let params = &L::PARAMS;
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let g = n - 1;
    let results: Vec<T1> = entries
        .par_iter()
        .enumerate()
        .map(|(i, e)| task1_vector::<L>(&setup::<L>(i, e)))
        .collect();
    header::<L>(label, "TASK 1");
    let cnt = |f: &dyn Fn(&T1) -> bool| results.iter().filter(|r| f(r)).count();
    println!(
        "recomputed basis hints equal the transmitted ones: {}/100; [2^(F-N)] of the full-order basis equals the verifier's basis, both curves: {}/100",
        cnt(&|r| r.hints_ok),
        cnt(&|r| r.scaled_ok)
    );
    println!(
        "order(omega) == 2^N on {}/100, order(eta) == 2^N on {}/100",
        cnt(&|r| r.ord_omega == Some(n)),
        cnt(&|r| r.ord_eta == Some(n))
    );
    println!("--- Theorem: y = log_omega(eta^-1) in mu_{{2^N}}, g = N - 1 = {g} ---");
    println!("dlog succeeded on {}/100", cnt(&|r| r.dlog_ok));
    println!(
        "y + 2^(g-1) == det M mod 2^g on {}/100; first bit at which y + 2^(g-1) and det M differ: == g on {}/100 (u = 3), > g on {}/100 (u = 1, they agree mod 2^N), < g on {}/100",
        cnt(&|r| r.agree_g),
        cnt(&|r| r.k == g),
        cnt(&|r| r.k > g),
        cnt(&|r| r.k < g)
    );
    let mut k_hist: BTreeMap<usize, usize> = BTreeMap::new();
    for r in &results {
        *k_hist.entry(r.k).or_default() += 1;
    }
    println!("histogram of that first differing bit (N = {n} means none below N): {k_hist:?}");
    println!("--- omega^(det M) eta ---");
    println!(
        "== 1: {}/100, == -1: {}/100, primitive fourth root of unity: {}/100, none: {}/100; u = 1 on {}/100, u = 3 on {}/100",
        cnt(&|r| r.z_class == 0),
        cnt(&|r| r.z_class == 1),
        cnt(&|r| r.z_class == 2),
        cnt(&|r| r.z_class == 3),
        cnt(&|r| r.u == 1),
        cnt(&|r| r.u == 3)
    );
    println!("--- Corollary: v_2(det M) = 0; the pivot ---");
    let mut v2_hist: BTreeMap<String, usize> = BTreeMap::new();
    for r in &results {
        let v = r.v2_det.map_or("det=0".to_string(), |v| v.to_string());
        *v2_hist.entry(v).or_default() += 1;
    }
    println!(
        "v_2(det M) histogram {v2_hist:?}; det M odd on {}/100",
        cnt(&|r| r.v2_det == Some(0))
    );
    println!(
        "M00 and M01 both even: {}/100; M00 even (pivot path, M11 transmitted): {}/100; M01 even: {}/100, M10 even: {}/100, M11 even: {}/100",
        cnt(&|r| r.m00_even && r.m01_even),
        cnt(&|r| r.m00_even),
        cnt(&|r| r.m01_even),
        cnt(&|r| r.m10_even),
        cnt(&|r| r.m11_even)
    );
    println!("--- per vector ---");
    for r in &results {
        println!(
            "KAT {:3}: hints_ok={} scaled_ok={} ord(omega)={} ord(eta)={} dlog_ok={} agree_g={} k={} zclass={} u={} v2det={:?} even(M00,M01,M10,M11)=({},{},{},{})",
            r.idx, r.hints_ok, r.scaled_ok, ord_str(r.ord_omega), ord_str(r.ord_eta), r.dlog_ok,
            r.agree_g, r.k, r.z_class, r.u, r.v2_det, r.m00_even, r.m01_even, r.m10_even, r.m11_even
        );
    }
    assert!(results.iter().all(|r| r.dlog_ok && r.agree_g && r.z_class == 2 && r.u != 0));
    assert!(results.iter().all(|r| r.v2_det == Some(0) && !(r.m00_even && r.m01_even)));
    assert!(results.iter().all(|r| r.hints_ok && r.scaled_ok));
}

#[test]
fn task1_formula_audit_l1() {
    run_task1::<Level1>("level I (p324_3)", KAT_L1);
}

#[test]
fn task1_formula_audit_l3() {
    run_task1::<Level3>("level III (p500_27)", KAT_L3);
}

#[test]
fn task1_formula_audit_l5() {
    run_task1::<Level5>("level V (p664_17)", KAT_L5);
}

// ---------------------------------------------------------------------------
// Task 2: the Tate sweep
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Cell {
    /// 2-adic orders of the challenge and auxiliary values
    ord_c: Option<usize>,
    ord_a: Option<usize>,
    /// the dlog in the subgroup of order ord_c succeeded
    dlog_ok: bool,
    /// first bit at which y_m and det M differ, capped at min(ord_c, N)
    k: usize,
    /// `y_m + 2^(g-1) == det M mod 2^g`
    ok_g: bool,
    /// `y_m == det M mod 2^(g-1)`
    ok_e: bool,
}

const NV: usize = 2;
const VARIANTS: [&str; NV] = ["a: verifier's order-2^N bases", "c: full-order bases scaled to 2^m"];

struct Sweep {
    rows: Vec<(usize, [Cell; NV])>,
}

fn tate_cell<L: Level>(
    m: usize,
    bc: &EcBasis<L>,
    e_chall: &mut EcCurve<L>,
    ba: &EcBasis<L>,
    e_aux: &mut EcCurve<L>,
    det: &Wide,
    n: usize,
) -> Cell {
    let g = n - 1;
    let tc = reduced_tate::<L>(m as u32, &bc.p, &bc.q, &bc.pmq, e_chall);
    let ta = reduced_tate::<L>(m as u32, &ba.p, &ba.q, &ba.pmq, e_aux);
    let ord_c = ord_2adic(&tc, m);
    let ord_a = ord_2adic(&ta, m);
    let (dlog_ok, k, ok_g, ok_e) = match (ord_c, ord_a) {
        (Some(a_c), Some(a_a)) if a_a <= a_c && a_c > 0 => match dlog::<L>(&ta.inv(), &tc, a_c) {
            Some(y) => {
                let cap = a_c.min(n);
                let det_n = masked(det, n);
                let k = first_mismatch(&y, &det_n, cap).unwrap_or(cap);
                let yp = masked(&add_wide(&y, &pow2(g - 1)), n);
                let ok_g = a_c >= g && agree(&yp, &det_n, g);
                let ok_e = a_c >= g - 1 && agree(&y, &det_n, g - 1);
                (true, k, ok_g, ok_e)
            }
            None => (false, 0, false, false),
        },
        _ => (false, 0, false, false),
    };
    Cell {
        ord_c,
        ord_a,
        dlog_ok,
        k,
        ok_g,
        ok_e,
    }
}

fn sweep_vector<L: Level>(s: &Setup<L>) -> Sweep {
    let params = &L::PARAMS;
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let fc = L::TWO_ADIC_EXPONENT as usize;
    let mut e_chall = s.e_chall.clone();
    let mut e_aux = s.e_aux.clone();
    let mut rows = Vec::new();
    for m in n..=fc {
        let a = tate_cell::<L>(m, &s.b_chall, &mut e_chall, &s.b_aux, &mut e_aux, &s.det, n);
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m, &mut e_chall);
        let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - m, &mut e_aux);
        let c = tate_cell::<L>(m, &bc, &mut e_chall, &ba, &mut e_aux, &s.det, n);
        rows.push((m, [a, c]));
    }
    Sweep { rows }
}

#[derive(Clone, PartialEq, Eq)]
struct Agg {
    n: usize,
    /// per variant: (ord_c min, max), (ord_a min, max), both orders == 2^max(2m - F, 0), dlog failures, min k, ok_g, ok_e
    v: [(usize, usize, usize, usize, usize, usize, usize, usize, usize); NV],
}

const NOT_POW2: usize = usize::MAX - 1;

fn ord_range(mn: usize, mx: usize) -> String {
    let one = |v: usize| {
        if v == NOT_POW2 {
            "x".to_string()
        } else {
            format!("2^{v}")
        }
    };
    if mn == mx {
        one(mn)
    } else {
        format!("{}..{}", one(mn), one(mx))
    }
}

fn agg_body(a: &Agg) -> String {
    let mut s = String::new();
    for (i, name) in VARIANTS.iter().enumerate() {
        let v = a.v[i];
        s.push_str(&format!(
            " [{}] ordC={} ordA={} both==2^max(2m-F,0): {} dlogfail={} min_k={} ok_g={} ok_e={} |",
            name,
            ord_range(v.0, v.1),
            ord_range(v.2, v.3),
            v.4,
            v.5,
            v.6,
            v.7,
            v.8
        ));
    }
    s
}

fn run_task2<L: Level>(label: &str, kat: &str) {
    let entries = parse_rsp(kat);
    assert_eq!(entries.len(), 100);
    let params = &L::PARAMS;
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let g = n - 1;
    let fc = L::TWO_ADIC_EXPONENT as usize;
    let t0 = Instant::now();
    let results: Vec<Sweep> = entries
        .par_iter()
        .enumerate()
        .map(|(i, e)| sweep_vector::<L>(&setup::<L>(i, e)))
        .collect();
    let elapsed = t0.elapsed();
    header::<L>(label, "TASK 2");
    println!(
        "reduced_tate(m, X, Y, X - Y, E): one cubical ladder of 2^m steps, final exponent (p^2 - 1)/2^m; y_m: T_aux^-1 = T_chl^(y_m) in the subgroup generated by T_chl; wall time {elapsed:?}"
    );
    println!("variant a: the verifier's bases of order 2^N = 2^{n}; variant c: TorsionBasisToHint(E, F, N) scaled by [2^(F-m)]; ok_g: y_m + 2^(g-1) == det M mod 2^g (g = {g}); ok_e: y_m == det M mod 2^(g-1) = 2^e_rsp");

    let mut table: Vec<(usize, Agg)> = Vec::new();
    for m in n..=fc {
        let pred = (2 * m as i64 - fc as i64).max(0) as usize;
        let mut agg = Agg {
            n: results.len(),
            v: [(usize::MAX, 0, usize::MAX, 0, 0, 0, usize::MAX, 0, 0); NV],
        };
        for r in &results {
            let (_, cells) = r.rows.iter().find(|(mm, _)| *mm == m).unwrap();
            for i in 0..NV {
                let c = cells[i];
                let v = &mut agg.v[i];
                let oc = c.ord_c.unwrap_or(NOT_POW2);
                let oa = c.ord_a.unwrap_or(NOT_POW2);
                v.0 = v.0.min(oc);
                v.1 = v.1.max(oc);
                v.2 = v.2.min(oa);
                v.3 = v.3.max(oa);
                if c.ord_c == Some(pred) && c.ord_a == Some(pred) {
                    v.4 += 1;
                }
                if !c.dlog_ok {
                    v.5 += 1;
                }
                v.6 = v.6.min(c.k);
                v.7 += c.ok_g as usize;
                v.8 += c.ok_e as usize;
            }
        }
        table.push((m, agg));
    }
    println!("--- full table, one row per m ---");
    for (m, a) in &table {
        println!("m={m:3} 2m-F={:4} n={} |{}", 2 * *m as i64 - fc as i64, a.n, agg_body(a));
    }
    println!("--- compressed table (runs of m with identical cells) ---");
    let mut i = 0;
    while i < table.len() {
        let mut j = i;
        while j + 1 < table.len() && table[j + 1].1 == table[i].1 {
            j += 1;
        }
        let (m0, a) = &table[i];
        let m1 = table[j].0;
        let mr = if *m0 == m1 {
            format!("m={m0}")
        } else {
            format!("m={m0}..{m1}")
        };
        println!("{mr:>12} 2m-F=[{},{}] |{}", 2 * *m0 as i64 - fc as i64, 2 * m1 as i64 - fc as i64, agg_body(a));
        i = j + 1;
    }

    println!("--- thresholds ---");
    let bound_g = (fc + g).div_ceil(2);
    let bound_e = (fc + g - 1).div_ceil(2);
    for (vi, name) in VARIANTS.iter().enumerate() {
        let first_g = table.iter().find(|(_, a)| a.v[vi].7 == a.n).map(|(m, _)| *m);
        let first_e = table.iter().find(|(_, a)| a.v[vi].8 == a.n).map(|(m, _)| *m);
        let max_g = table.iter().map(|(_, a)| a.v[vi].7).max().unwrap();
        let max_e = table.iter().map(|(_, a)| a.v[vi].8).max().unwrap();
        println!(
            "  [{name}] mod 2^g (g = {g}): smallest m with 100/100 = {first_g:?}, bound ceil((F + g)/2) = {bound_g}, most vectors at any m = {max_g}; mod 2^e_rsp: smallest m with 100/100 = {first_e:?}, bound ceil((F + e_rsp)/2) = {bound_e}, most at any m = {max_e}"
        );
        // per-vector smallest working m
        let mut hist_g: BTreeMap<Option<usize>, usize> = BTreeMap::new();
        let mut hist_e: BTreeMap<Option<usize>, usize> = BTreeMap::new();
        for r in &results {
            let fg = r.rows.iter().find(|(_, c)| c[vi].ok_g).map(|(m, _)| *m);
            let fe = r.rows.iter().find(|(_, c)| c[vi].ok_e).map(|(m, _)| *m);
            *hist_g.entry(fg).or_default() += 1;
            *hist_e.entry(fe).or_default() += 1;
        }
        println!("    per-vector smallest m, mod 2^g: {hist_g:?}; mod 2^e_rsp: {hist_e:?}");
        // the orders against 2^max(2m-F, 0) at every m
        let ord_ok = results
            .iter()
            .filter(|r| {
                r.rows.iter().all(|(m, c)| {
                    let a = (2 * *m as i64 - fc as i64).max(0) as usize;
                    c[vi].ord_c == Some(a) && c[vi].ord_a == Some(a)
                })
            })
            .count();
        println!("    order(T_chl) == order(T_aux) == 2^max(2m-F,0) at every m in [N, F]: {ord_ok}/100");
    }
    // the sanity row for the paper
    let ok_c = table.iter().find(|(m, _)| *m == bound_g).map(|(_, a)| a.v[1].7).unwrap();
    let below = table.iter().find(|(m, _)| *m + 1 == bound_g).map(|(_, a)| a.v[1].7).unwrap_or(0);
    println!("variant c at m = {bound_g}: ok_g on {ok_c}/100; at m = {}: {below}/100", bound_g - 1);
}

#[test]
fn task2_tate_sweep_l1() {
    run_task2::<Level1>("level I (p324_3)", KAT_L1);
}

#[test]
fn task2_tate_sweep_l3() {
    run_task2::<Level3>("level III (p500_27)", KAT_L3);
}

#[test]
fn task2_tate_sweep_l5() {
    run_task2::<Level5>("level V (p664_17)", KAT_L5);
}

// ---------------------------------------------------------------------------
// Task 3: the decompressor on one Tate ladder per curve, end to end
// ---------------------------------------------------------------------------

fn words(t: u32) -> usize {
    (t as usize).div_ceil(64)
}

fn mask_2exp(x: &mut Scalar, t: u32) {
    for (i, w) in x.iter_mut().enumerate() {
        let lo = 64 * i as u32;
        if lo >= t {
            *w = 0;
        } else if lo + 64 > t {
            *w &= (1u64 << (t - lo)) - 1;
        }
    }
}

fn reduced(x: &Scalar, t: u32) -> Scalar {
    let mut y = *x;
    mask_2exp(&mut y, t);
    y
}

fn mul_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let n = words(t);
    let mut out = [0u64; MAX_ORDER_WORDS];
    for i in 0..n {
        let mut carry = 0u128;
        for j in 0..n - i {
            let acc = out[i + j] as u128 + (a[i] as u128) * (b[j] as u128) + carry;
            out[i + j] = acc as u64;
            carry = acc >> 64;
        }
    }
    mask_2exp(&mut out, t);
    out
}

fn add_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let mut out = [0u64; MAX_ORDER_WORDS];
    let mut carry = 0u128;
    for i in 0..words(t) {
        let s = a[i] as u128 + b[i] as u128 + carry;
        out[i] = s as u64;
        carry = s >> 64;
    }
    mask_2exp(&mut out, t);
    out
}

fn sub_mod_2exp(a: &Scalar, b: &Scalar, t: u32) -> Scalar {
    let mut out = [0u64; MAX_ORDER_WORDS];
    let mut borrow = 0u64;
    for i in 0..words(t) {
        let (d, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = d.overflowing_sub(borrow);
        out[i] = d;
        borrow = (b1 | b2) as u64;
    }
    mask_2exp(&mut out, t);
    out
}

fn inv_mod_2exp(a: &Scalar, t: u32) -> Option<Scalar> {
    if a[0] & 1 == 0 {
        return None;
    }
    let mut two = [0u64; MAX_ORDER_WORDS];
    two[0] = 2;
    let mut x = reduced(a, t);
    let mut bits = 3u32;
    while bits < t {
        let ax = mul_mod_2exp(a, &x, t);
        let f = sub_mod_2exp(&two, &ax, t);
        x = mul_mod_2exp(&x, &f, t);
        bits *= 2;
    }
    Some(x)
}

/// The decompressor's reconstruction (compressed::reconstruct) from the
/// three transmitted entries and `det M` modulo `2^(e_rsp + 1)`.
fn reconstruct<L: Level>(
    r: u32,
    c: &sqisign_verify::compressed::CompressedSignature<L>,
    det: &Scalar,
) -> Option<[[Scalar; 2]; 2]> {
    let t = r + 1;
    let (m00, m01, m_var) = (c.m00, c.m01, c.m_var);
    let (m10, m11) = if m00[0] & 1 == 1 {
        let num = add_mod_2exp(det, &mul_mod_2exp(&m01, &m_var, t), t);
        (m_var, mul_mod_2exp(&num, &inv_mod_2exp(&m00, t)?, t))
    } else {
        let num = sub_mod_2exp(&mul_mod_2exp(&m00, &m_var, t), det, t);
        (mul_mod_2exp(&num, &inv_mod_2exp(&m01, t)?, t), m_var)
    };
    Some([[m00, m01], [m10, m11]])
}

struct T3 {
    /// the standard verifier accepts the reconstructed signature
    pass: bool,
    /// the recovered dropped entry equals the signer's modulo 2^(e_rsp + 1)
    entry_exact: bool,
    /// it differs from the signer's in bit e_rsp only
    entry_bit_only: bool,
    /// the Tate value had the expected order and the dlog succeeded
    dlog_ok: bool,
}

fn task3_vector<L: Level>(s: &Setup<L>, pk: &sqisign_verify::sqisign::PublicKey<L>, msg: &[u8], m: usize) -> T3 {
    let params = &L::PARAMS;
    let r = params.response_bits;
    let n = (r + HD_EXTRA_TORSION) as usize;
    let g = n - 1;
    let fc = L::TWO_ADIC_EXPONENT as usize;
    let mut e_chall = s.e_chall.clone();
    let mut e_aux = s.e_aux.clone();
    let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m, &mut e_chall);
    let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - m, &mut e_aux);
    let tc = reduced_tate::<L>(m as u32, &bc.p, &bc.q, &bc.pmq, &mut e_chall);
    let ta = reduced_tate::<L>(m as u32, &ba.p, &ba.q, &ba.pmq, &mut e_aux);
    let a = 2 * m - fc;
    let y = match ord_2adic(&tc, m) {
        Some(ord) if ord == a => dlog::<L>(&ta.inv(), &tc, a),
        _ => None,
    };
    let Some(y) = y else {
        return T3 { pass: false, entry_exact: false, entry_bit_only: false, dlog_ok: false };
    };
    // det M = y_m + 2^(g-1) mod 2^g, with y_m as the ladder gives it, below 2^(2m-F)
    let mut y_s = [0u64; MAX_ORDER_WORDS];
    y_s.copy_from_slice(&y[..MAX_ORDER_WORDS]);
    let mut two = [0u64; MAX_ORDER_WORDS];
    two[((g - 1) / 64) as usize] |= 1u64 << ((g - 1) % 64);
    let det = add_mod_2exp(&y_s, &two, g as u32);
    let c = compress(params, &s.sig).expect("compress");
    let Some(mat) = reconstruct::<L>(r, &c, &det) else {
        return T3 { pass: false, entry_exact: false, entry_bit_only: false, dlog_ok: true };
    };
    let rec = Signature {
        e_aux_a: s.sig.e_aux_a.clone(),
        mat_bchall_can_to_bchall: mat,
        chall_coeff: s.sig.chall_coeff,
        hint_aux: s.sig.hint_aux,
        hint_chall: s.sig.hint_chall,
    };
    let pass = verify(params, pk, &rec, msg);
    let dropped = if s.sig.mat_bchall_can_to_bchall[0][0][0] & 1 == 1 { (1, 1) } else { (1, 0) };
    let orig = reduced(&s.sig.mat_bchall_can_to_bchall[dropped.0][dropped.1], g as u32);
    let got = reduced(&mat[dropped.0][dropped.1], g as u32);
    let mut diff = [0u64; MAX_ORDER_WORDS];
    for i in 0..MAX_ORDER_WORDS {
        diff[i] = orig[i] ^ got[i];
    }
    let entry_exact = diff.iter().all(|w| *w == 0);
    let mut bit = [0u64; MAX_ORDER_WORDS];
    bit[(r / 64) as usize] |= 1u64 << (r % 64);
    let entry_bit_only = diff == bit;
    T3 { pass, entry_exact, entry_bit_only, dlog_ok: true }
}

fn run_task3<L: Level>(label: &str, kat: &str) {
    let entries = parse_rsp(kat);
    assert_eq!(entries.len(), 100);
    let params = &L::PARAMS;
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let g = n - 1;
    let fc = L::TWO_ADIC_EXPONENT as usize;
    let bound = (fc + g).div_ceil(2);
    header::<L>(label, "TASK 3");
    println!("the decompressor with det M from one Tate ladder per curve on the scaled full-order bases (y_m + 2^(g-1) mod 2^g, y_m below 2^(2m-F)), the reconstructed signature given to the standard verifier");
    for m in [bound - 1, bound] {
        let results: Vec<T3> = entries
            .par_iter()
            .enumerate()
            .map(|(i, e)| {
                let s = setup::<L>(i, e);
                let pk = public_key_from_bytes::<L>(params, &e.pk).expect("pk");
                task3_vector::<L>(&s, &pk, &e.msg, m)
            })
            .collect();
        let cnt = |f: &dyn Fn(&T3) -> bool| results.iter().filter(|r| f(r)).count();
        println!(
            "m = {m} (2m - F = {}): dlog ok {}/100; standard verifier accepts {}/100; recovered dropped entry equals the signer's mod 2^g {}/100, differs in bit e_rsp only {}/100, differs otherwise {}/100; accepted with a wrong entry {}/100, rejected with the right entry {}/100",
            2 * m - fc,
            cnt(&|r| r.dlog_ok),
            cnt(&|r| r.pass),
            cnt(&|r| r.entry_exact),
            cnt(&|r| r.entry_bit_only),
            cnt(&|r| r.dlog_ok && !r.entry_exact && !r.entry_bit_only),
            cnt(&|r| r.pass && !r.entry_exact),
            cnt(&|r| !r.pass && r.entry_exact)
        );
        if m == bound {
            assert!(results.iter().all(|r| r.pass && r.entry_exact));
        } else {
            assert!(results.iter().all(|r| r.entry_exact || r.entry_bit_only));
        }
    }
}

#[test]
fn task3_end_to_end_tate_l1() {
    run_task3::<Level1>("level I (p324_3)", KAT_L1);
}

#[test]
fn task3_end_to_end_tate_l3() {
    run_task3::<Level3>("level III (p500_27)", KAT_L3);
}

#[test]
fn task3_end_to_end_tate_l5() {
    run_task3::<Level5>("level V (p664_17)", KAT_L5);
}

// ---------------------------------------------------------------------------
// Timing, level I, KAT entry 0
// ---------------------------------------------------------------------------

const ITERS: usize = 200;

/// Median of `ITERS` timings of `f` after 5 warm-ups, each on a fresh input
/// from `setup` (the routines mutate their curves).
fn measure<I, T>(mut setup: impl FnMut() -> I, mut f: impl FnMut(I) -> T) -> f64 {
    for _ in 0..5 {
        black_box(f(setup()));
    }
    let mut v = Vec::with_capacity(ITERS);
    for _ in 0..ITERS {
        let input = setup();
        let t = Instant::now();
        let out = f(input);
        let dt = t.elapsed();
        black_box(out);
        v.push(dt.as_secs_f64() * 1e6);
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (v[ITERS / 2 - 1] + v[ITERS / 2]) / 2.0
}

#[test]
fn timing_l1_kat0() {
    type L = Level1;
    let params = &L::PARAMS;
    let entries = parse_rsp(KAT_L1);
    let s = setup::<L>(0, &entries[0]);
    let n = params.response_bits + HD_EXTRA_TORSION;
    let g = n - 1;
    let fc = L::TWO_ADIC_EXPONENT;
    let m_star: u32 = std::env::var("PAPER_TATE_M")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or((fc + g).div_ceil(2));
    let (hc, ha) = (s.sig.hint_chall, s.sig.hint_aux);
    let e_chall0 = s.e_chall.clone();
    let e_aux0 = s.e_aux.clone();

    // the verifier's stage: bases at 2^N from the hints, two Weil pairings
    let bases_n = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            (
                ec_curve_to_basis_2f_from_hint(&mut ec, n, hc).unwrap(),
                ec_curve_to_basis_2f_from_hint(&mut ea, n, ha).unwrap(),
            )
        },
    );
    let weil_one = measure(
        || e_chall0.clone(),
        |mut ec| weil::<L>(n, &s.b_chall.p, &s.b_chall.q, &s.b_chall.pmq, &mut ec),
    );
    let weil_both = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            (
                weil::<L>(n, &s.b_chall.p, &s.b_chall.q, &s.b_chall.pmq, &mut ec),
                weil::<L>(n, &s.b_aux.p, &s.b_aux.q, &s.b_aux.pmq, &mut ea),
            )
        },
    );
    let stage_weil = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            let bc = ec_curve_to_basis_2f_from_hint(&mut ec, n, hc).unwrap();
            let ba = ec_curve_to_basis_2f_from_hint(&mut ea, n, ha).unwrap();
            (
                weil::<L>(n, &bc.p, &bc.q, &bc.pmq, &mut ec),
                weil::<L>(n, &ba.p, &ba.q, &ba.pmq, &mut ea),
            )
        },
    );
    let mut e1 = e_chall0.clone();
    let mut e2 = e_aux0.clone();
    let omega = weil::<L>(n, &s.b_chall.p, &s.b_chall.q, &s.b_chall.pmq, &mut e1);
    let eta = weil::<L>(n, &s.b_aux.p, &s.b_aux.q, &s.b_aux.pmq, &mut e2);
    let (eta_inv, omega_inv) = (eta.inv(), omega.inv());
    let dlog_n = measure(
        || [0u64; MAX_ORDER_WORDS],
        |mut y| {
            fp2_dlog_2e::<L>(&mut y[..L::ORDER_WORDS], &eta_inv, &omega_inv, n).unwrap();
            y
        },
    );

    // the Tate alternative: bases at 2^m*, one ladder per curve
    let bases_m = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            (
                ec_curve_to_basis_2f_from_hint(&mut ec, m_star, hc).unwrap(),
                ec_curve_to_basis_2f_from_hint(&mut ea, m_star, ha).unwrap(),
            )
        },
    );
    let mut e1 = e_chall0.clone();
    let mut e2 = e_aux0.clone();
    let bc_m = ec_dbl_iter_basis(&s.b_chall_full, (fc - m_star) as usize, &mut e1);
    let ba_m = ec_dbl_iter_basis(&s.b_aux_full, (fc - m_star) as usize, &mut e2);
    let tate_one = measure(
        || e_chall0.clone(),
        |mut ec| reduced_tate::<L>(m_star, &bc_m.p, &bc_m.q, &bc_m.pmq, &mut ec),
    );
    let tate_both = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            (
                reduced_tate::<L>(m_star, &bc_m.p, &bc_m.q, &bc_m.pmq, &mut ec),
                reduced_tate::<L>(m_star, &ba_m.p, &ba_m.q, &ba_m.pmq, &mut ea),
            )
        },
    );
    let stage_tate = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            let bc = ec_curve_to_basis_2f_from_hint(&mut ec, m_star, hc).unwrap();
            let ba = ec_curve_to_basis_2f_from_hint(&mut ea, m_star, ha).unwrap();
            (
                reduced_tate::<L>(m_star, &bc.p, &bc.q, &bc.pmq, &mut ec),
                reduced_tate::<L>(m_star, &ba.p, &ba.q, &ba.pmq, &mut ea),
            )
        },
    );
    // what the integrated verifier still needs after the Tate ladders: the
    // two points of each basis brought down to order 2^N for the basis
    // change and the chain, and the difference point at that order (the
    // verifier's basis, TorsionBasisFromHint at N, is exactly that)
    let b_from_m = |bm: &EcBasis<L>, e: &mut EcCurve<L>| {
        let p = ec_dbl_iter(&bm.p, (m_star - n) as usize, e);
        let q = ec_dbl_iter(&bm.pmq, (m_star - n) as usize, e);
        let r = difference_point(&p, &q, e);
        EcBasis::new(p, r, q)
    };
    let (mut e3, mut e4) = (e_chall0.clone(), e_aux0.clone());
    let b_chall_m = ec_curve_to_basis_2f_from_hint(&mut e3, m_star, hc).unwrap();
    let b_aux_m = ec_curve_to_basis_2f_from_hint(&mut e4, m_star, ha).unwrap();
    let reduced_chall = b_from_m(&b_chall_m, &mut e3);
    let reduced_aux = b_from_m(&b_aux_m, &mut e4);
    let same_basis = |a: &EcBasis<L>, b: &EcBasis<L>| {
        bool::from(a.p.ct_equal(&b.p) & a.q.ct_equal(&b.q) & a.pmq.ct_equal(&b.pmq))
    };
    let reduction_ok = same_basis(&reduced_chall, &s.b_chall) && same_basis(&reduced_aux, &s.b_aux);
    let reduce = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| (b_from_m(&b_chall_m, &mut ec), b_from_m(&b_aux_m, &mut ea)),
    );
    let stage_tate_full = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            let bc = ec_curve_to_basis_2f_from_hint(&mut ec, m_star, hc).unwrap();
            let ba = ec_curve_to_basis_2f_from_hint(&mut ea, m_star, ha).unwrap();
            let t1 = reduced_tate::<L>(m_star, &bc.p, &bc.q, &bc.pmq, &mut ec);
            let t2 = reduced_tate::<L>(m_star, &ba.p, &ba.q, &ba.pmq, &mut ea);
            (t1, t2, b_from_m(&bc, &mut ec), b_from_m(&ba, &mut ea))
        },
    );
    // the Tate recovery at m* on this vector, and its discrete logarithm
    let tc = reduced_tate::<L>(m_star, &bc_m.p, &bc_m.q, &bc_m.pmq, &mut e1);
    let ta = reduced_tate::<L>(m_star, &ba_m.p, &ba_m.q, &ba_m.pmq, &mut e2);
    let a = ord_2adic(&tc, m_star as usize).expect("tate value order");
    let y = dlog::<L>(&ta.inv(), &tc, a).expect("tate dlog");
    let det_n = masked(&s.det, n as usize);
    let yp = masked(&add_wide(&y, &pow2(g as usize - 1)), n as usize);
    let recovered = a >= g as usize && agree(&yp, &det_n, g as usize);
    let (ta_inv, tc_inv) = (ta.inv(), tc.inv());
    let dlog_a = measure(
        || [0u64; MAX_ORDER_WORDS],
        |mut y| {
            fp2_dlog_2e::<L>(&mut y[..L::ORDER_WORDS], &ta_inv, &tc_inv, a as u32).unwrap();
            y
        },
    );

    println!(
        "==================== TIMING: level I (p324_3), KAT entry 0, N = {n}, g = {g}, F = {fc}, m* = {m_star}, medians of {ITERS} iterations after 5 warm-ups, {} ====================",
        sqisign_kat::arithmetic_from_env()
    );
    println!("Weil path (the verifier's): bases from the hints at 2^N on both curves ({} doublings per point after the cofactor): {bases_n:.1} us; one Weil pairing at 2^N (two ladders of {} steps): {weil_one:.1} us; both Weil pairings: {weil_both:.1} us; stage (bases + both pairings): {stage_weil:.1} us; discrete logarithm at 2^N: {dlog_n:.1} us", fc - n, n - 1);
    println!("Tate path: bases from the hints at 2^m* on both curves ({} doublings per point after the cofactor): {bases_m:.1} us; one Tate ladder of 2^m* steps ({} steps, final exponentiation with {} squarings): {tate_one:.1} us; both ladders: {tate_both:.1} us; stage (bases + both ladders): {stage_tate:.1} us; value order 2^{a} (2m* - F = {}), discrete logarithm at that order: {dlog_a:.1} us; recovers det M mod 2^g on this vector: {recovered}", fc - m_star, m_star - 1, fc - m_star, 2 * m_star as i64 - fc as i64);
    println!("Tate path, continued: bringing the two points of each basis from order 2^m* down to 2^N ({} doublings per point) and the difference point at 2^N, both curves: {reduce:.1} us (gives the verifier's bases: {reduction_ok}); stage with that reduction (bases at 2^m* + both ladders + reduction): {stage_tate_full:.1} us", m_star - n);
    println!("stage ratio Tate / Weil: {:.3} without the reduction, {:.3} with it; pairings alone Tate / Weil: {:.3}", stage_tate / stage_weil, stage_tate_full / stage_weil, tate_both / weil_both);
    assert!(recovered && reduction_ok);
}
