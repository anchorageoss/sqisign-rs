//! Numerical experiments for the SQIsign compression paper, run on the 300
//! official NIST round-2 KAT vectors (100 per level) at tag v0.4.31.
//!
//! Task 1 (formula audit): canonical bases exactly as the decompressor builds
//! them, Weil pairings w_f and w_g, the production discrete logarithm, the
//! congruence with det M, the candidate relations between w_f and w_g, the
//! direct isotropy check on P' = [M00]P + [M10]Q, Q' = [M01]P + [M11]Q, plus
//! Check A (primitivity of w_f^det * w_g) and Check B (v_2(det M) against r').
//!
//! Task 2 (Tate order sweep): reduced Tate pairings whose cubical ladder runs
//! 2^m steps and whose final exponentiation is (p^2 - 1)/2^m, for every m in
//! [f, F], on (a) the decompressor's order-2^f / order-2^g points and (b) the
//! order-2^F canonical basis points, in both argument orders.
//!
//! Timing (Level 1, KAT entry 0): production Weil at 2^f, one Tate ladder at
//! m = PAPER_TATE_M (environment variable, default f) and one at m = F.
//!
//! Run (release only; `--test-threads=1` keeps the printed output in order,
//! the tests parallelise over vectors internally):
//!   cargo test -p sqisign-verify --release --features kat-compat \
//!       --test paper_pairing_audit -- --nocapture --test-threads=1

#![cfg(feature = "kat-compat")]

use hybrid_array::typenum::Unsigned;
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::isogeny::ec_eval_small_chain;
use sqisign_verify::ec::pairing::{fp2_dlog_2e_pub, reduced_tate, weil};
use sqisign_verify::ec::point::{ec_biscalar_mul, ec_dbl_iter, ec_dbl_iter_basis, xadd};
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::params::{Level1, Level3, Level5, SecurityLevel};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::theta::HD_EXTRA_TORSION;
use sqisign_verify::verify::compute_challenge_curve;
use sqisign_verify::{PublicKey, Scalar, Signature};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const KAT_LVL1: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp"
));
const KAT_LVL3: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_529_SQIsign_lvl3.rsp"
));
const KAT_LVL5: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp"
));

// ---------------------------------------------------------------------------
// KAT parsing
// ---------------------------------------------------------------------------

fn hex_to_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

struct KatEntry {
    pk: Vec<u8>,
    sm: Vec<u8>,
    mlen: usize,
}

fn parse_kat_entries(content: &str) -> Vec<KatEntry> {
    let mut entries = Vec::new();
    let mut pk = None;
    let mut sm = None;
    let mut mlen = None;
    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("mlen = ") {
            mlen = Some(rest.parse::<usize>().unwrap());
        } else if let Some(rest) = line.strip_prefix("pk = ") {
            pk = Some(hex_to_bytes(rest));
        } else if let Some(rest) = line.strip_prefix("sm = ") {
            sm = Some(hex_to_bytes(rest));
        }
        if pk.is_some() && sm.is_some() && mlen.is_some() {
            entries.push(KatEntry {
                pk: pk.take().unwrap(),
                sm: sm.take().unwrap(),
                mlen: mlen.take().unwrap(),
            });
        }
    }
    entries
}

// ---------------------------------------------------------------------------
// Multiprecision helpers on 16-limb little-endian arrays
// ---------------------------------------------------------------------------

const W: usize = 16;
type Wide = [u64; W];

fn mp_mul_wide(a: &[u64], b: &[u64]) -> Wide {
    let mut out = [0u64; W];
    for i in 0..a.len() {
        let mut carry: u64 = 0;
        for j in 0..b.len() {
            if i + j >= W {
                break;
            }
            let p = (a[i] as u128) * (b[j] as u128) + (out[i + j] as u128) + (carry as u128);
            out[i + j] = p as u64;
            carry = (p >> 64) as u64;
        }
        if i + b.len() < W {
            out[i + b.len()] = carry;
        }
    }
    out
}

fn mp_sub_wide(a: &Wide, b: &Wide) -> Wide {
    let mut out = [0u64; W];
    let mut borrow: u64 = 0;
    for i in 0..W {
        let (d1, b1) = a[i].overflowing_sub(b[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        out[i] = d2;
        borrow = (b1 as u64) + (b2 as u64);
    }
    out
}

fn mp_mul_small(a: &Wide, c: u64) -> Wide {
    mp_mul_wide(a, &[c])
}

fn mp_shl(a: &Wide, shift: usize) -> Wide {
    let mut out = [0u64; W];
    let q = shift / 64;
    let r = shift % 64;
    for i in (0..W).rev() {
        if i < q {
            break;
        }
        let lo = a[i - q] << r;
        let hi = if r > 0 && i >= q + 1 {
            a[i - q - 1] >> (64 - r)
        } else {
            0
        };
        out[i] = lo | hi;
    }
    out
}

fn mask_bits(a: &mut [u64], nbits: usize) {
    let q = nbits / 64;
    let r = nbits % 64;
    if q < a.len() {
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
    mask_bits(&mut out, nbits);
    out
}

/// (2^nbits - a) mod 2^nbits.
fn mp_neg_mod(a: &Wide, nbits: usize) -> Wide {
    let zero = [0u64; W];
    masked(&mp_sub_wide(&zero, a), nbits)
}

fn pow2(k: usize) -> Wide {
    let mut out = [0u64; W];
    out[k / 64] = 1u64 << (k % 64);
    out
}

fn lowest_set_bit(a: &[u64]) -> Option<usize> {
    for (i, &limb) in a.iter().enumerate() {
        if limb != 0 {
            return Some(i * 64 + limb.trailing_zeros() as usize);
        }
    }
    None
}

/// First differing bit below `nbits`, or None if the low `nbits` bits agree.
fn first_mismatch(a: &[u64], b: &[u64], nbits: usize) -> Option<usize> {
    let mut d = [0u64; W];
    for i in 0..a.len().min(b.len()).min(W) {
        d[i] = a[i] ^ b[i];
    }
    mask_bits(&mut d, nbits);
    lowest_set_bit(&d)
}

fn mp_eq(a: &Wide, b: &Wide, nbits: usize) -> bool {
    first_mismatch(a, b, nbits).is_none()
}

fn wide_from(s: &[u64]) -> Wide {
    let mut out = [0u64; W];
    out[..s.len()].copy_from_slice(s);
    out
}

fn scalar_wide<L: FpBackend>(s: &Scalar<L>) -> Wide {
    wide_from(s.digits())
}

/// det(M) = M00*M11 - M01*M10 as a 1024-bit two's complement integer from the
/// four transmitted entries (not reduced; callers mask).
fn ground_truth_det<L: FpBackend>(sig: &Signature<L>) -> Wide {
    let m = sig.mat();
    let p1 = mp_mul_wide(m[0][0].digits(), m[1][1].digits());
    let p2 = mp_mul_wide(m[0][1].digits(), m[1][0].digits());
    mp_sub_wide(&p1, &p2)
}

// ---------------------------------------------------------------------------
// Fp2 helpers
// ---------------------------------------------------------------------------

fn fp2_eq<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

fn fp2_is_one<L: FpBackend>(a: &Fp2<L>) -> bool {
    bool::from(a.ct_is_one())
}

fn fp2_pow2k<L: FpBackend>(x: &Fp2<L>, k: usize) -> Fp2<L> {
    let mut r = x.clone();
    for _ in 0..k {
        r = r.sqr();
    }
    r
}

/// Smallest `a` with `x^(2^a) == 1`, or None if not reached within `cap`
/// squarings (then the order is not a power of two dividing 2^cap).
fn fp2_ord_2adic<L: FpBackend>(x: &Fp2<L>, cap: usize) -> Option<usize> {
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

fn ord_str(o: Option<usize>) -> String {
    match o {
        Some(k) => format!("2^{}", k),
        None => "NOT-2-POWER".to_string(),
    }
}

// ---------------------------------------------------------------------------
// Simple scoped thread pool (no extra dependency)
// ---------------------------------------------------------------------------

fn par_map<T: Send, F: Fn(usize) -> T + Sync>(n: usize, f: F) -> Vec<T> {
    let threads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1)
        .min(n)
        .max(1);
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<T>>> = (0..n).map(|_| Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n {
                    break;
                }
                let r = f(i);
                *results[i].lock().unwrap() = Some(r);
            });
        }
    });
    results
        .into_iter()
        .map(|m| m.into_inner().unwrap().unwrap())
        .collect()
}

// ---------------------------------------------------------------------------
// Per-vector setup: everything the decompressor computes before the dlog
// ---------------------------------------------------------------------------

struct VecSetup<L: SecurityLevel> {
    idx: usize,
    bt: usize,
    trl: usize,
    pow_dim2: usize,
    /// f = pow_dim2 + HD_EXTRA_TORSION + trl = E_RSP - bt + 2 (E_chall torsion exponent, = N)
    f: usize,
    /// g = pow_dim2 + HD_EXTRA_TORSION = E_RSP - bt - trl + 2 (E_aux torsion exponent)
    g: usize,
    /// d = det_precision = pow_dim2 + trl = E_RSP - bt (bits the production code keeps)
    d: usize,
    sig: Signature<L>,
    det: Wide,
    hint_chall_ok: bool,
    hint_aux_ok: bool,
    e_chall: EcCurve<L>,
    /// Canonical basis of E_chall[2^F] (after the odd cofactor, before the 2^(F-f) doublings).
    b_chall_full: EcBasis<L>,
    /// = [2^(F-f)] b_chall_full: the decompressor's points.
    b_chall: EcBasis<L>,
    ppq_chall: EcPoint<L>,
    ppq_chall_full: EcPoint<L>,
    e_aux: EcCurve<L>,
    b_aux_full: EcBasis<L>,
    /// = [2^(F-g)] b_aux_full
    b_aux: EcBasis<L>,
    ppq_aux: EcPoint<L>,
    ppq_aux_full: EcPoint<L>,
}

fn canonical_basis<L: FpBackend + LevelPrecomp>(curve: &mut EcCurve<L>) -> (EcBasis<L>, u8) {
    ec_curve_to_basis_2f_to_hint(
        curve,
        L::F_CHR,
        L::basis_e0_px_bytes(),
        L::basis_e0_qx_bytes(),
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
        L::torsion_even_power(),
    )
    .expect("canonical basis")
}

fn setup<L: FpBackend + LevelPrecomp>(idx: usize, entry: &KatEntry) -> VecSetup<L> {
    let pk = PublicKey::<L>::from_bytes(&entry.pk).expect("pk decode");
    let sig_len = entry.sm.len() - entry.mlen;
    let sig = Signature::<L>::from_bytes(&entry.sm[..sig_len]).expect("sig decode");

    let bt = sig.backtracking() as usize;
    let trl = sig.two_resp_length() as usize;
    let pow_dim2 = L::E_RSP as usize - trl - bt;
    let f = pow_dim2 + HD_EXTRA_TORSION as usize + trl;
    let g = pow_dim2 + HD_EXTRA_TORSION as usize;
    let d = pow_dim2 + trl;
    let fc = L::F_CHR as usize;

    let det = ground_truth_det(&sig);

    // E_chall and its canonical basis, exactly as CompressedSignature::decompress.
    let mut e_chall =
        compute_challenge_curve::<L>(sig.chall_coeff(), sig.backtracking(), pk.curve(), pk.hint_pk())
            .expect("challenge curve");
    let (b_chall_full, hint_chall) = canonical_basis::<L>(&mut e_chall);
    let b_chall = ec_dbl_iter_basis(&b_chall_full, fc - f, &mut e_chall);
    let ppq_chall = xadd(&b_chall.p, &b_chall.q, &b_chall.pmq);
    let ppq_chall_full = xadd(&b_chall_full.p, &b_chall_full.q, &b_chall_full.pmq);

    let mut e_aux = EcCurve::<L>::from_a(sig.e_aux_a()).expect("aux curve");
    let (b_aux_full, hint_aux) = canonical_basis::<L>(&mut e_aux);
    let b_aux = ec_dbl_iter_basis(&b_aux_full, fc - g, &mut e_aux);
    let ppq_aux = xadd(&b_aux.p, &b_aux.q, &b_aux.pmq);
    let ppq_aux_full = xadd(&b_aux_full.p, &b_aux_full.q, &b_aux_full.pmq);

    e_chall.normalize_a24();
    e_aux.normalize_a24();

    VecSetup {
        idx,
        bt,
        trl,
        pow_dim2,
        f,
        g,
        d,
        hint_chall_ok: hint_chall == sig.hint_chall(),
        hint_aux_ok: hint_aux == sig.hint_aux(),
        sig,
        det,
        e_chall,
        b_chall_full,
        b_chall,
        ppq_chall,
        ppq_chall_full,
        e_aux,
        b_aux_full,
        b_aux,
        ppq_aux,
        ppq_aux_full,
    }
}

/// Production dlog convention: find y with `arg = base^y` in the 2^e subgroup,
/// calling `fp2_dlog_2e_pub(y, arg, base^-1, e)`.
fn dlog<L: FpBackend>(arg: &Fp2<L>, base: &Fp2<L>, e: usize) -> Option<Wide> {
    let nw = L::MpLimbs::USIZE;
    let mut y = [0u64; W];
    fp2_dlog_2e_pub::<L>(&mut y[..nw], arg, &base.inv(), e as u32)?;
    mask_bits(&mut y, e);
    Some(y)
}

// ---------------------------------------------------------------------------
// Task 1
// ---------------------------------------------------------------------------

struct T1 {
    idx: usize,
    bt: usize,
    trl: usize,
    f: usize,
    g: usize,
    d: usize,
    hints_ok: bool,
    ord_wf: Option<usize>,
    ord_wg: Option<usize>,
    dlog_ok: bool,
    /// first mismatch bit of y vs det M (mod 2^f); f when they agree on all f bits
    k: usize,
    /// (label, w_g == w_f^exp, w_g == w_f^-exp) after raising both sides to 2^raise
    rels: Vec<(String, bool, bool)>,
    /// e_{2^f}(P', Q') == w_f^{det M mod 2^f}
    e_pq_is_wf_det: bool,
    /// s = 1, 2, 4: (e^s * w_g^s == 1, e^s * w_g^-s == 1), pre-chain at order 2^f
    c_pre: [(bool, bool); 3],
    /// same after pushing P', Q' through the 2^trl chain (order 2^g); None if trl == 0
    c_post: Option<[(bool, bool); 3]>,
    /// Check A: 0: z == 1, 1: z == -1, 2: z^2 == -1, 3: none
    z_class: u8,
    /// Check A: (det - y) mod 2^N == k * 2^(N-2) -> k in 0..=3; 4: other
    diff_class: u8,
    /// Check B
    v2_det: Option<usize>,
    m00_even: bool,
    m01_even: bool,
    m10_even: bool,
    m11_even: bool,
}

fn task1_vector<L: FpBackend + LevelPrecomp>(s: &VecSetup<L>) -> T1 {
    let fc = L::F_CHR as usize;
    let cof = L::p_cofactor_for_2f()[0];
    let mut e_chall = s.e_chall.clone();
    let mut e_aux = s.e_aux.clone();

    // Production pairings (formats.rs decompress): weil(f, P, Q, P+Q) and weil(g, Pa, Qa, Pa+Qa).
    let w_f = weil::<L>(s.f as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall);
    let w_g = weil::<L>(s.g as u32, &s.b_aux.p, &s.b_aux.q, &s.ppq_aux, &mut e_aux);

    // Production dlog: fp2_dlog_2e_pub(det, omega_aux^-1, omega_f^-1, f):
    // arg = w_g^-1, base = w_f  =>  w_g^-1 = w_f^y.
    let y = dlog::<L>(&w_g.inv(), &w_f, s.f);
    let det_f = masked(&s.det, s.f);
    let (dlog_ok, k, yv) = match y {
        Some(yv) => (true, first_mismatch(&det_f, &yv, s.f).unwrap_or(s.f), yv),
        None => (false, 0, [0u64; W]),
    };

    // (b) candidate relations.
    let mut rels = Vec::new();
    let test_rel = |exp: &Wide, raise: usize| -> (bool, bool) {
        let lhs = fp2_pow2k(&w_g, raise);
        let rhs = fp2_pow2k(&w_f.pow_vartime(&exp[..]), raise);
        (fp2_eq(&lhs, &rhs), fp2_eq(&lhs, &rhs.inv()))
    };
    let neg_det_f = mp_neg_mod(&det_f, s.f);
    for (label, exp) in [
        ("det", det_f),
        ("det*c", masked(&mp_mul_small(&det_f, cof), s.f)),
        ("det*2", masked(&mp_shl(&det_f, 1), s.f)),
        ("det*4", masked(&mp_shl(&det_f, 2), s.f)),
        ("det*8", masked(&mp_shl(&det_f, 3), s.f)),
        ("det*16", masked(&mp_shl(&det_f, 4), s.f)),
    ] {
        for raise in [0usize, s.f - s.g, 2, s.f - s.g + 2] {
            let (pos, neg) = test_rel(&exp, raise);
            rels.push((format!("{} raise=2^{}", label, raise), pos, neg));
        }
    }
    // The literal four candidates (w_g^-1 = w_f^{+-det} is the same statement as
    // w_g = w_f^{-+det}); recorded for completeness with the negated exponent.
    {
        let (pos, neg) = test_rel(&neg_det_f, 0);
        rels.push(("-det raise=2^0 (w_g = w_f^-det | w_g = w_f^det)".to_string(), pos, neg));
    }

    // (c) P' = [M00]P + [M10]Q, Q' = [M01]P + [M11]Q, exactly the verifier's basis change.
    let nw = L::MpLimbs::USIZE;
    let m = s.sig.mat();
    let p1 = ec_biscalar_mul(m[0][0].digits(), m[1][0].digits(), s.f, &s.b_chall, &e_chall)
        .expect("P'");
    let q1 = ec_biscalar_mul(m[0][1].digits(), m[1][1].digits(), s.f, &s.b_chall, &e_chall)
        .expect("Q'");
    let s0 = masked(&mp_sub_wide(&scalar_wide(&m[0][0]), &scalar_wide(&m[0][1])), s.f);
    let s1 = masked(&mp_sub_wide(&scalar_wide(&m[1][0]), &scalar_wide(&m[1][1])), s.f);
    let pmq1 = ec_biscalar_mul(&s0[..nw], &s1[..nw], s.f, &s.b_chall, &e_chall).expect("P'-Q'");
    let ppq1 = xadd(&p1, &q1, &pmq1);
    let e_pq = weil::<L>(s.f as u32, &p1, &q1, &ppq1, &mut e_chall);
    let e_pq_is_wf_det = fp2_eq(&e_pq, &w_f.pow_vartime(&det_f[..]));
    let iso = |e: &Fp2<L>| -> [(bool, bool); 3] {
        let mut out = [(false, false); 3];
        for (i, sexp) in [0usize, 1, 2].iter().enumerate() {
            let es = fp2_pow2k(e, *sexp);
            let ws = fp2_pow2k(&w_g, *sexp);
            out[i] = (fp2_is_one(&es.mul(&ws)), fp2_is_one(&es.mul(&ws.inv())));
        }
        out
    };
    let c_pre = iso(&e_pq);
    let c_post = if s.trl > 0 {
        // two_response_isogeny_verify_inner: kernel = [2^g](Q' if M00, M10 even else P').
        let kernel_is_q = m[0][0].digits()[0] & 1 == 0 && m[1][0].digits()[0] & 1 == 0;
        let ker = if kernel_is_q { q1.clone() } else { p1.clone() };
        let mut e2 = e_chall.clone();
        let ker = ec_dbl_iter(&ker, s.g, &mut e2);
        let mut pts = [p1.clone(), q1.clone(), pmq1.clone()];
        ec_eval_small_chain(&mut e2, &ker, s.trl as i32, &mut pts, false).expect("2^trl chain");
        let ppq2 = xadd(&pts[0], &pts[1], &pts[2]);
        let e_post = weil::<L>(s.g as u32, &pts[0], &pts[1], &ppq2, &mut e2);
        Some(iso(&e_post))
    } else {
        None
    };

    // Check A: z = w_f^{det mod 2^N} * w_g.
    let z = w_f.pow_vartime(&det_f[..]).mul(&w_g);
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
    let diff = masked(&mp_sub_wide(&det_f, &yv), s.f);
    let quarter = pow2(s.f - 2);
    let mut diff_class = 4u8;
    for kq in 0..4usize {
        let cand = masked(&mp_mul_small(&quarter, kq as u64), s.f);
        if mp_eq(&diff, &cand, s.f) {
            diff_class = kq as u8;
        }
    }

    // Check B.
    let v2_det = lowest_set_bit(&s.det);
    let even = |x: &Scalar<L>| x.digits()[0] & 1 == 0;

    let _ = fc;
    T1 {
        idx: s.idx,
        bt: s.bt,
        trl: s.trl,
        f: s.f,
        g: s.g,
        d: s.d,
        hints_ok: s.hint_chall_ok && s.hint_aux_ok,
        ord_wf: fp2_ord_2adic(&w_f, fc + 2),
        ord_wg: fp2_ord_2adic(&w_g, fc + 2),
        dlog_ok,
        k,
        rels,
        e_pq_is_wf_det,
        c_pre,
        c_post,
        z_class,
        diff_class,
        v2_det,
        m00_even: even(&m[0][0]),
        m01_even: even(&m[0][1]),
        m10_even: even(&m[1][0]),
        m11_even: even(&m[1][1]),
    }
}

fn run_task1<L: FpBackend + LevelPrecomp>(label: &str, kat: &str) {
    let entries = parse_kat_entries(kat);
    assert_eq!(entries.len(), 100, "expected 100 KAT entries");
    let n = entries.len();
    let fc = L::F_CHR as usize;
    let e_rsp = L::E_RSP as usize;

    let results = par_map(n, |i| {
        let s = setup::<L>(i, &entries[i]);
        task1_vector::<L>(&s)
    });

    println!("==================== TASK 1: {} ====================", label);
    println!(
        "E_RSP = {} (params/level*.rs const E_RSP), HD_EXTRA_TORSION = {} (theta/mod.rs:19), F = F_CHR = {} (params/level*.rs const F_CHR), c = (p+1)/2^F = {} (precomp/level*/mod.rs P_COFACTOR_FOR_2F), LAMBDA = {}",
        e_rsp,
        HD_EXTRA_TORSION,
        fc,
        L::p_cofactor_for_2f()[0],
        L::LAMBDA
    );
    println!(
        "N = E_RSP + HD_EXTRA_TORSION - n_bt (verify.rs:112 check_canonical_basis_change_matrix; = f below); matrix entry wire bytes = floor((E_RSP + 9)/8) = {} (types.rs:520)",
        (e_rsp + 9) / 8
    );
    println!("per-vector: f = pow_dim2 + 2 + trl (formats.rs:636), g = pow_dim2 + 2 (formats.rs:637), d = det_precision = pow_dim2 + trl (formats.rs:635), pow_dim2 = E_RSP - trl - bt (formats.rs:628)");
    println!(
        "pow_dim2 of vector 0 = {} (sanity: E_RSP - trl - bt)",
        setup::<L>(0, &entries[0]).pow_dim2
    );
    let range = |get: &dyn Fn(&T1) -> usize| {
        let mn = results.iter().map(get).min().unwrap();
        let mx = results.iter().map(get).max().unwrap();
        (mn, mx)
    };
    let (bt_mn, bt_mx) = range(&|r| r.bt);
    let (trl_mn, trl_mx) = range(&|r| r.trl);
    let (f_mn, f_mx) = range(&|r| r.f);
    let (g_mn, g_mx) = range(&|r| r.g);
    let (d_mn, d_mx) = range(&|r| r.d);
    println!(
        "ranges over 100 vectors: n_bt in [{}, {}], r' = two_resp_length in [{}, {}], f = N in [{}, {}], g in [{}, {}], d in [{}, {}]",
        bt_mn, bt_mx, trl_mn, trl_mx, f_mn, f_mx, g_mn, g_mx, d_mn, d_mx
    );
    let mut bt_hist: BTreeMap<usize, usize> = BTreeMap::new();
    let mut trl_hist: BTreeMap<usize, usize> = BTreeMap::new();
    for r in &results {
        *bt_hist.entry(r.bt).or_default() += 1;
        *trl_hist.entry(r.trl).or_default() += 1;
    }
    println!("n_bt histogram: {:?}", bt_hist);
    println!("r'  histogram: {:?}", trl_hist);
    println!(
        "recomputed canonical-basis hints equal the transmitted hints on {}/100 vectors",
        results.iter().filter(|r| r.hints_ok).count()
    );

    // Orders of the pairings.
    let wf_ok = results.iter().filter(|r| r.ord_wf == Some(r.f)).count();
    let wg_ok = results.iter().filter(|r| r.ord_wg == Some(r.g)).count();
    println!("order(w_f) == 2^f on {}/100, order(w_g) == 2^g on {}/100", wf_ok, wg_ok);

    // (a)
    println!("--- (a) y = log_{{w_f}}(w_g^-1) at order 2^f [fp2_dlog_2e_pub(det, w_g^-1, w_f^-1, f)] vs det M ---");
    println!("dlog succeeded on {}/100", results.iter().filter(|r| r.dlog_ok).count());
    let k_min = results.iter().map(|r| r.k).min().unwrap();
    let k_max = results.iter().map(|r| r.k).max().unwrap();
    println!("largest k with y == det M mod 2^k on ALL vectors (absolute): k = {} (per-vector k ranges over [{}, {}])", k_min, k_min, k_max);
    let eq_d = results.iter().filter(|r| r.k == r.d).count();
    let gt_d = results.iter().filter(|r| r.k > r.d).count();
    let lt_d = results.iter().filter(|r| r.k < r.d).count();
    println!(
        "per vector, first mismatch bit k vs d = E_RSP - n_bt: k == d on {}/100, k > d on {}/100, k < d on {}/100",
        eq_d, gt_d, lt_d
    );
    let mut kd_hist: BTreeMap<i64, usize> = BTreeMap::new();
    for r in &results {
        *kd_hist.entry(r.k as i64 - r.d as i64).or_default() += 1;
    }
    println!("histogram of (k - d): {:?}", kd_hist);
    let mut by_bt_trl: BTreeMap<(usize, usize), (usize, usize, usize)> = BTreeMap::new();
    for r in &results {
        let e = by_bt_trl.entry((r.bt, r.trl)).or_insert((usize::MAX, 0, 0));
        e.0 = e.0.min(r.k);
        e.1 = e.1.max(r.k);
        e.2 += 1;
    }
    for ((bt, trl), (kmn, kmx, cnt)) in &by_bt_trl {
        println!(
            "  (n_bt={}, r'={}): {} vectors, k in [{}, {}], d = {}, f = {}",
            bt,
            trl,
            cnt,
            kmn,
            kmx,
            e_rsp - bt,
            e_rsp - bt + 2
        );
    }

    // (b)
    println!("--- (b) relations between w_g and w_f^exp (both sides raised to 2^raise; raise=2^(f-g) is the task's 2^(f-g)) ---");
    let nrel = results[0].rels.len();
    for j in 0..nrel {
        let pos = results.iter().filter(|r| r.rels[j].1).count();
        let neg = results.iter().filter(|r| r.rels[j].2).count();
        // the raise label is per vector (f-g = r'); summarise by name of the first vector and by r' buckets
        let mut per_trl_pos: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for r in &results {
            let e = per_trl_pos.entry(r.trl).or_insert((0, 0));
            if r.rels[j].1 {
                e.0 += 1;
            }
            if r.rels[j].2 {
                e.1 += 1;
            }
        }
        println!(
            "  [{}] {:<44} w_g == w_f^exp: {:>3}/100   w_g == w_f^-exp: {:>3}/100   (by r': {:?})",
            j, results[0].rels[j].0, pos, neg, per_trl_pos
        );
    }
    println!("  (raise labels above are from vector 0; for a vector with r' = t the entries 'raise=2^{{f-g}}' use 2^t and 'raise=2^{{f-g+2}}' use 2^(t+2))");

    // (c)
    println!("--- (c) direct isotropy: e = e_{{2^f}}(P', Q') with P' = [M00]P + [M10]Q, Q' = [M01]P + [M11]Q ---");
    println!(
        "e_{{2^f}}(P', Q') == w_f^(det M mod 2^f) exactly on {}/100",
        results.iter().filter(|r| r.e_pq_is_wf_det).count()
    );
    for (i, sv) in [1usize, 2, 4].iter().enumerate() {
        let pos = results.iter().filter(|r| r.c_pre[i].0).count();
        let neg = results.iter().filter(|r| r.c_pre[i].1).count();
        println!(
            "  s = {}: e^s * w_g^(+s) == 1 on {:>3}/100   e^s * w_g^(-s) == 1 on {:>3}/100   (pre-chain, order 2^f on E_chall)",
            sv, pos, neg
        );
    }
    let npost = results.iter().filter(|r| r.c_post.is_some()).count();
    if npost > 0 {
        for (i, sv) in [1usize, 2, 4].iter().enumerate() {
            let pos = results
                .iter()
                .filter(|r| r.c_post.map_or(false, |c| c[i].0))
                .count();
            let neg = results
                .iter()
                .filter(|r| r.c_post.map_or(false, |c| c[i].1))
                .count();
            println!(
                "  s = {}: e^s * w_g^(+s) == 1 on {:>3}/{}   e^s * w_g^(-s) == 1 on {:>3}/{}   (post 2^r' chain, order 2^g on E_chall', r' > 0 vectors only)",
                sv, pos, npost, neg, npost
            );
        }
    }

    // Check A
    println!("--- Check A: z = w_f^(det M mod 2^N) * w_g ---");
    let cnt = |c: u8| results.iter().filter(|r| r.z_class == c).count();
    println!(
        "z == 1: {}/100, z == -1: {}/100, z^2 == -1 (primitive 4th root): {}/100, none: {}/100",
        cnt(0),
        cnt(1),
        cnt(2),
        cnt(3)
    );
    let dc = |c: u8| results.iter().filter(|r| r.diff_class == c).count();
    println!(
        "(det M - y) mod 2^N: == 0: {}/100, == 2^(N-2): {}/100, == 2*2^(N-2): {}/100, == 3*2^(N-2): {}/100, other: {}/100",
        dc(0),
        dc(1),
        dc(2),
        dc(3),
        dc(4)
    );

    // Check B
    println!("--- Check B: v_2(det M) tabulated against r' ---");
    let mut tab: BTreeMap<usize, BTreeMap<String, usize>> = BTreeMap::new();
    for r in &results {
        let v = match r.v2_det {
            Some(v) => v.to_string(),
            None => "det=0".to_string(),
        };
        *tab.entry(r.trl).or_default().entry(v).or_default() += 1;
    }
    for (trl, hist) in &tab {
        println!("  r' = {}: v_2(det M) histogram {:?}", trl, hist);
    }
    println!(
        "v_2(det M) == r' on {}/100 vectors",
        results.iter().filter(|r| r.v2_det == Some(r.trl)).count()
    );
    println!(
        "both M[0][0] and M[0][1] even: {}/100; both M[0][0] and M[1][0] even: {}/100; M00 even: {}/100, M01 even: {}/100, M10 even: {}/100, M11 even: {}/100",
        results.iter().filter(|r| r.m00_even && r.m01_even).count(),
        results.iter().filter(|r| r.m00_even && r.m10_even).count(),
        results.iter().filter(|r| r.m00_even).count(),
        results.iter().filter(|r| r.m01_even).count(),
        results.iter().filter(|r| r.m10_even).count(),
        results.iter().filter(|r| r.m11_even).count(),
    );
    println!("convention (verify.rs:260-283 matrix_scalar_application_even_basis): P' = [M[0][0]]P + [M[1][0]]Q, Q' = [M[0][1]]P + [M[1][1]]Q (column j of M = coefficients of the j-th new point); 2^r' kernel is Q' iff M[0][0] and M[1][0] are both even (verify.rs:389)");

    // Per-vector dump (compact).
    println!("--- per-vector ---");
    for r in &results {
        println!(
            "KAT {:3}: bt={} trl={} f={} g={} d={} ord(w_f)={} ord(w_g)={} dlog_ok={} k={} k-d={} zclass={} diffclass={} v2det={:?} pre(s=1,2,4)={:?} post={:?}",
            r.idx,
            r.bt,
            r.trl,
            r.f,
            r.g,
            r.d,
            ord_str(r.ord_wf),
            ord_str(r.ord_wg),
            r.dlog_ok,
            r.k,
            r.k as i64 - r.d as i64,
            r.z_class,
            r.diff_class,
            r.v2_det,
            r.c_pre,
            r.c_post
        );
    }
    assert!(results.iter().all(|r| r.dlog_ok));
    assert!(results.iter().all(|r| r.k >= r.d), "production precision claim violated");
}

#[test]
fn task1_formula_audit_l1() {
    run_task1::<Level1>("Level 1 (p = 5*2^248 - 1)", KAT_LVL1);
}

#[test]
fn task1_formula_audit_l3() {
    run_task1::<Level3>("Level 3 (p = 65*2^376 - 1)", KAT_LVL3);
}

#[test]
fn task1_formula_audit_l5() {
    run_task1::<Level5>("Level 5 (p = 27*2^500 - 1)", KAT_LVL5);
}

// ---------------------------------------------------------------------------
// Task 2: Tate order sweep
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Cell {
    ord_c: Option<usize>,
    ord_a: Option<usize>,
    /// value^(2^m) == 1 on both curves
    in_mu_2m: bool,
    /// dlog taken in mu_{2^m} (the sweep's nominal order)
    dlog_ok: bool,
    /// bits of det M recovered by y: [y vs det, y vs -det, 2^r' y vs det, 2^r' y vs -det],
    /// each the first mismatching bit (comparisons capped at min(order, f) bits)
    k: [usize; 4],
    /// dlog taken in mu_{2^a} where 2^a = order(T_chall) (needs order(T_aux) | 2^a)
    dlog_o_ok: bool,
    k_o: [usize; 4],
}

/// Bits of det M recovered by a dlog value `y` known modulo 2^e: first mismatch
/// of y vs det, y vs -det, 2^trl*y vs det, 2^trl*y vs -det. Comparisons are
/// capped at f bits because the matrix entries are only defined modulo 2^f.
fn kbits(y: &Wide, e: usize, det: &Wide, trl: usize, f: usize) -> [usize; 4] {
    let nb = e.min(f);
    let det_nb = masked(det, nb);
    let neg_nb = mp_neg_mod(&det_nb, nb);
    let k_pos = first_mismatch(&det_nb, y, nb).unwrap_or(nb);
    let k_neg = first_mismatch(&neg_nb, y, nb).unwrap_or(nb);
    let nb2 = (e + trl).min(f);
    let det_nb2 = masked(det, nb2);
    let neg_nb2 = mp_neg_mod(&det_nb2, nb2);
    let ysh = masked(&mp_shl(&masked(y, e), trl), nb2);
    let k_pos_sh = first_mismatch(&det_nb2, &ysh, nb2).unwrap_or(nb2);
    let k_neg_sh = first_mismatch(&neg_nb2, &ysh, nb2).unwrap_or(nb2);
    [k_pos, k_neg, k_pos_sh, k_neg_sh]
}

const NV: usize = 6;
const VARIANTS: [&str; NV] = [
    "a:t(P,Q)",
    "a:t(Q,P)",
    "b:t(R,S)",
    "b:t(S,R)",
    "c:t(A,B)",
    "c:t(B,A)",
];

struct SweepVec {
    bt: usize,
    trl: usize,
    f: usize,
    g: usize,
    d: usize,
    rows: Vec<(usize, [Cell; NV])>,
}

fn tate_cell<L: FpBackend + LevelPrecomp>(
    m: usize,
    pc: &EcPoint<L>,
    qc: &EcPoint<L>,
    ppqc: &EcPoint<L>,
    e_chall: &mut EcCurve<L>,
    pa: &EcPoint<L>,
    qa: &EcPoint<L>,
    ppqa: &EcPoint<L>,
    e_aux: &mut EcCurve<L>,
    det: &Wide,
    trl: usize,
    f: usize,
) -> Cell {
    let fc = L::F_CHR;
    let cof = L::p_cofactor_for_2f();
    let tc = reduced_tate::<L>(m as u32, pc, qc, ppqc, e_chall, fc, cof);
    let ta = reduced_tate::<L>(m as u32, pa, qa, ppqa, e_aux, fc, cof);
    let ord_c = fp2_ord_2adic(&tc, m);
    let ord_a = fp2_ord_2adic(&ta, m);
    let in_mu_2m = ord_c.is_some() && ord_a.is_some();
    let (dlog_ok, k) = match dlog::<L>(&ta.inv(), &tc, m) {
        Some(y) => (true, kbits(&y, m, det, trl, f)),
        None => (false, [0; 4]),
    };
    let (dlog_o_ok, k_o) = match (ord_c, ord_a) {
        (Some(a_c), Some(a_a)) if a_a <= a_c && a_c > 0 => match dlog::<L>(&ta.inv(), &tc, a_c) {
            Some(y) => (true, kbits(&y, a_c, det, trl, f)),
            None => (false, [0; 4]),
        },
        _ => (false, [0; 4]),
    };
    Cell {
        ord_c,
        ord_a,
        in_mu_2m,
        dlog_ok,
        k,
        dlog_o_ok,
        k_o,
    }
}

fn sweep_vector<L: FpBackend + LevelPrecomp>(s: &VecSetup<L>) -> SweepVec {
    let fc = L::F_CHR as usize;
    let mut e_chall = s.e_chall.clone();
    let mut e_aux = s.e_aux.clone();
    let mut rows = Vec::new();
    for m in s.f..=fc {
        let a_pq = tate_cell::<L>(
            m,
            &s.b_chall.p,
            &s.b_chall.q,
            &s.ppq_chall,
            &mut e_chall,
            &s.b_aux.p,
            &s.b_aux.q,
            &s.ppq_aux,
            &mut e_aux,
            &s.det,
            s.trl,
            s.f,
        );
        let a_qp = tate_cell::<L>(
            m,
            &s.b_chall.q,
            &s.b_chall.p,
            &s.ppq_chall,
            &mut e_chall,
            &s.b_aux.q,
            &s.b_aux.p,
            &s.ppq_aux,
            &mut e_aux,
            &s.det,
            s.trl,
            s.f,
        );
        let b_rs = tate_cell::<L>(
            m,
            &s.b_chall_full.p,
            &s.b_chall_full.q,
            &s.ppq_chall_full,
            &mut e_chall,
            &s.b_aux_full.p,
            &s.b_aux_full.q,
            &s.ppq_aux_full,
            &mut e_aux,
            &s.det,
            s.trl,
            s.f,
        );
        let b_sr = tate_cell::<L>(
            m,
            &s.b_chall_full.q,
            &s.b_chall_full.p,
            &s.ppq_chall_full,
            &mut e_chall,
            &s.b_aux_full.q,
            &s.b_aux_full.p,
            &s.ppq_aux_full,
            &mut e_aux,
            &s.det,
            s.trl,
            s.f,
        );
        // Variant c: scale the order-2^F bases down to order 2^m first
        // (F - m doublings of P, Q and P-Q on each curve), then pair at order 2^m.
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m, &mut e_chall);
        let ppqc = xadd(&bc.p, &bc.q, &bc.pmq);
        let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - m, &mut e_aux);
        let ppqa = xadd(&ba.p, &ba.q, &ba.pmq);
        let c_ab = tate_cell::<L>(
            m, &bc.p, &bc.q, &ppqc, &mut e_chall, &ba.p, &ba.q, &ppqa, &mut e_aux, &s.det, s.trl, s.f,
        );
        let c_ba = tate_cell::<L>(
            m, &bc.q, &bc.p, &ppqc, &mut e_chall, &ba.q, &ba.p, &ppqa, &mut e_aux, &s.det, s.trl, s.f,
        );
        rows.push((m, [a_pq, a_qp, b_rs, b_sr, c_ab, c_ba]));
    }
    SweepVec {
        bt: s.bt,
        trl: s.trl,
        f: s.f,
        g: s.g,
        d: s.d,
        rows,
    }
}

/// Per-variant aggregate over the participating vectors at one m.
#[derive(Clone, Copy, PartialEq, Eq)]
struct VAgg {
    oc_min: usize,
    oc_max: usize,
    oa_min: usize,
    oa_max: usize,
    notmu: usize,
    /// vectors whose ord_chall and ord_aux both equal 2^max(2m-F, 0)
    pred: usize,
    fail: usize,
    /// min over vectors of k[j] (dlog at order 2^m)
    k: [usize; 4],
    /// number of vectors with k[j] >= d
    reach: [usize; 4],
    fail_o: usize,
    k_o: [usize; 4],
    reach_o: [usize; 4],
}

#[derive(Clone, PartialEq, Eq)]
struct Agg {
    n: usize,
    v: [VAgg; NV],
    /// min/max over participating vectors of 2f - F (clipped at 0) and of d
    tf_min: i64,
    tf_max: i64,
    d_min: usize,
    d_max: usize,
}

const NOT_POW2: usize = usize::MAX - 1;

fn ord_range(mn: usize, mx: usize) -> String {
    let one = |v: usize| {
        if v == NOT_POW2 {
            "x".to_string()
        } else {
            format!("2^{}", v)
        }
    };
    if mn == mx {
        one(mn)
    } else {
        format!("{}..{}", one(mn), one(mx))
    }
}

fn agg_line(m: usize, fc: usize, a: &Agg) -> String {
    let mut line = format!(
        "m={:3} n={:3} 2m-F={:4} 2f-F=[{},{}] d=[{},{}] |",
        m,
        a.n,
        2 * m as i64 - fc as i64,
        a.tf_min,
        a.tf_max,
        a.d_min,
        a.d_max
    );
    for (i, name) in VARIANTS.iter().enumerate() {
        let v = a.v[i];
        line.push_str(&format!(
            " {} ordC={} ordA={} notmu={} pred2m-F={} | dlog@2^m: fail={} k(+,-,+sh,-sh)={}/{}/{}/{} reachd={}/{}/{}/{} | dlog@ord: fail={} k={}/{}/{}/{} reachd={}/{}/{}/{} |",
            name,
            ord_range(v.oc_min, v.oc_max),
            ord_range(v.oa_min, v.oa_max),
            v.notmu,
            v.pred,
            v.fail,
            v.k[0],
            v.k[1],
            v.k[2],
            v.k[3],
            v.reach[0],
            v.reach[1],
            v.reach[2],
            v.reach[3],
            v.fail_o,
            v.k_o[0],
            v.k_o[1],
            v.k_o[2],
            v.k_o[3],
            v.reach_o[0],
            v.reach_o[1],
            v.reach_o[2],
            v.reach_o[3],
        ));
    }
    line
}

fn run_task2<L: FpBackend + LevelPrecomp>(label: &str, kat: &str) {
    let entries = parse_kat_entries(kat);
    assert_eq!(entries.len(), 100);
    let n = entries.len();
    let fc = L::F_CHR as usize;

    let t0 = Instant::now();
    let results = par_map(n, |i| {
        let s = setup::<L>(i, &entries[i]);
        sweep_vector::<L>(&s)
    });
    let elapsed = t0.elapsed();

    println!("==================== TASK 2: {} ====================", label);
    println!(
        "F = v2(p+1) = {}, c = {}; reduced_tate(m, X, Y, X+Y, E, F, c) = ladder of 2^m steps on X, final exponent (p-1)*c*2^(F-m) = (p^2-1)/2^m (pairing.rs reduced_tate); dlog y: T_aux^-1 = T_chall^y in mu_{{2^m}} [fp2_dlog_2e_pub(y, T_aux^-1, T_chall^-1, m)]",
        fc,
        L::p_cofactor_for_2f()[0]
    );
    println!("variant a: decompressor's points, P,Q of order 2^f on E_chall and Pa,Qa of order 2^g on E_aux; variant b: order-2^F canonical bases R,S (basis_to_hint at F_CHR, before the 2^(F-f) resp. 2^(F-g) doublings): P = [2^(F-f)]R, Q = [2^(F-f)]S, Pa = [2^(F-g)]Ra, Qa = [2^(F-g)]Sa");
    println!("k+ = min over vectors of the first bit where y differs from det M (mod 2^m); k- = same against -det M; reachd = number of vectors whose k reaches d = E_RSP - n_bt; notmu = vectors whose value is not in mu_{{2^m}}; wall time for the sweep: {:?}", elapsed);

    let f_min = results.iter().map(|r| r.f).min().unwrap();
    let mut table: Vec<(usize, Agg)> = Vec::new();
    for m in f_min..=fc {
        let part: Vec<(&SweepVec, &[Cell; NV])> = results
            .iter()
            .filter_map(|r| r.rows.iter().find(|(mm, _)| *mm == m).map(|(_, c)| (r, c)))
            .collect();
        if part.is_empty() {
            continue;
        }
        let empty = VAgg {
            oc_min: usize::MAX,
            oc_max: 0,
            oa_min: usize::MAX,
            oa_max: 0,
            notmu: 0,
            pred: 0,
            fail: 0,
            k: [usize::MAX; 4],
            reach: [0; 4],
            fail_o: 0,
            k_o: [usize::MAX; 4],
            reach_o: [0; 4],
        };
        let pred_ord = (2 * m as i64 - fc as i64).max(0) as usize;
        let mut agg = Agg {
            n: part.len(),
            v: [empty; NV],
            tf_min: i64::MAX,
            tf_max: i64::MIN,
            d_min: usize::MAX,
            d_max: 0,
        };
        for (r, cells) in &part {
            let tf = (2 * r.f as i64 - fc as i64).max(0);
            agg.tf_min = agg.tf_min.min(tf);
            agg.tf_max = agg.tf_max.max(tf);
            agg.d_min = agg.d_min.min(r.d);
            agg.d_max = agg.d_max.max(r.d);
            for i in 0..NV {
                let c = cells[i];
                let v = &mut agg.v[i];
                let oc = c.ord_c.unwrap_or(NOT_POW2);
                let oa = c.ord_a.unwrap_or(NOT_POW2);
                if c.ord_c == Some(pred_ord) && c.ord_a == Some(pred_ord) {
                    v.pred += 1;
                }
                v.oc_min = v.oc_min.min(oc);
                v.oc_max = v.oc_max.max(oc);
                v.oa_min = v.oa_min.min(oa);
                v.oa_max = v.oa_max.max(oa);
                if !c.in_mu_2m {
                    v.notmu += 1;
                }
                if !c.dlog_ok {
                    v.fail += 1;
                }
                if !c.dlog_o_ok {
                    v.fail_o += 1;
                }
                for j in 0..4 {
                    v.k[j] = v.k[j].min(c.k[j]);
                    v.k_o[j] = v.k_o[j].min(c.k_o[j]);
                    if c.dlog_ok && c.k[j] >= r.d {
                        v.reach[j] += 1;
                    }
                    if c.dlog_o_ok && c.k_o[j] >= r.d {
                        v.reach_o[j] += 1;
                    }
                }
            }
        }
        table.push((m, agg));
    }

    println!("--- full table, one row per m ---");
    for (m, a) in &table {
        println!("{}", agg_line(*m, fc, a));
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
        let mrange = if *m0 == m1 {
            format!("m={}", m0)
        } else {
            format!("m={}..{}", m0, m1)
        };
        let body = agg_line(*m0, fc, a);
        let body = body.splitn(2, '|').nth(1).unwrap_or("");
        println!(
            "{:>14} n={:3} 2m-F=[{},{}] 2f-F=[{},{}] d=[{},{}] |{}",
            mrange,
            a.n,
            2 * *m0 as i64 - fc as i64,
            2 * m1 as i64 - fc as i64,
            a.tf_min,
            a.tf_max,
            a.d_min,
            a.d_max,
            body
        );
        i = j + 1;
    }

    println!("--- smallest m at which all 100 vectors reach d (per variant; comparisons: y vs det, y vs -det, 2^r' y vs det, 2^r' y vs -det) ---");
    const CMP: [&str; 4] = ["y=+det", "y=-det", "2^r'y=+det", "2^r'y=-det"];
    for (vi, name) in VARIANTS.iter().enumerate() {
        for j in 0..4 {
            let first = table
                .iter()
                .find(|(_, a)| a.v[vi].reach[j] == a.n && a.n == n)
                .map(|(m, _)| *m);
            let best = table.iter().map(|(_, a)| a.v[vi].reach[j]).max().unwrap_or(0);
            let best_m = table
                .iter()
                .filter(|(_, a)| a.v[vi].reach[j] == best && best > 0)
                .map(|(m, _)| *m)
                .collect::<Vec<_>>();
            let first_o = table
                .iter()
                .find(|(_, a)| a.v[vi].reach_o[j] == a.n && a.n == n)
                .map(|(m, _)| *m);
            let best_o = table.iter().map(|(_, a)| a.v[vi].reach_o[j]).max().unwrap_or(0);
            println!(
                "  {} [{}]: dlog@2^m: smallest m with 100/100 reaching d = {:?}; max vectors reaching d over m = {} (at m = {:?}) | dlog@ord: smallest m = {:?}; max = {}",
                name, CMP[j], first, best, best_m, first_o, best_o
            );
        }
    }

    println!("--- per n_bt group: smallest m at which every vector of the group reaches d (dlog@2^m, or dlog@ord when the full-order dlog fails), vs ceil((F + d + 1)/2) ---");
    let mut bts: Vec<usize> = results.iter().map(|r| r.bt).collect();
    bts.sort_unstable();
    bts.dedup();
    for (vi, name) in VARIANTS.iter().enumerate() {
        for j in 0..4 {
            let mut parts = Vec::new();
            for &bt in &bts {
                let d = L::E_RSP as usize - bt;
                let group: Vec<&SweepVec> = results.iter().filter(|r| r.bt == bt).collect();
                let f_bt = L::E_RSP as usize + 2 - bt;
                let mut first: Option<usize> = None;
                for m in f_bt..=fc {
                    let all = group.iter().all(|r| {
                        let (_, c) = r.rows.iter().find(|(mm, _)| *mm == m).unwrap();
                        let c = c[vi];
                        (c.dlog_ok && c.k[j] >= r.d) || (!c.dlog_ok && c.dlog_o_ok && c.k_o[j] >= r.d)
                    });
                    if all {
                        first = Some(m);
                        break;
                    }
                }
                parts.push(format!(
                    "bt={} ({} vec, d={}, ceil((F+d+1)/2)={}): {:?}",
                    bt,
                    group.len(),
                    d,
                    (fc + d + 2) / 2,
                    first
                ));
            }
            println!("  {} [{}]: {}", name, CMP[j], parts.join("; "));
        }
    }

    // Variant c, per vector: smallest m at which 2^r' y == det M mod 2^d (dlog at
    // the actual order when the full-order dlog fails) vs ceil((F + d - r')/2).
    for vi in [4usize, 5] {
        let mut ok = 0usize;
        let mut mism = Vec::new();
        let mut ord_ok = 0usize;
        for (idx, r) in results.iter().enumerate() {
            let first = r
                .rows
                .iter()
                .find(|(_, c)| {
                    let c = c[vi];
                    (c.dlog_ok && c.k[2] >= r.d) || (!c.dlog_ok && c.dlog_o_ok && c.k_o[2] >= r.d)
                })
                .map(|(m, _)| *m);
            let pred = (fc + r.d - r.trl).div_ceil(2);
            if first == Some(pred) {
                ok += 1;
            } else {
                mism.push((idx, r.bt, r.trl, r.d, first, pred));
            }
            // order prediction 2^max(2m-F,0) on both curves for every m
            if r.rows.iter().all(|(m, c)| {
                let a = (2 * *m as i64 - fc as i64).max(0) as usize;
                c[vi].ord_c == Some(a) && c[vi].ord_a == Some(a)
            }) {
                ord_ok += 1;
            }
        }
        println!(
            "variant {}: smallest working m == ceil((F + d - r')/2) on {}/100 vectors; mismatches (idx, bt, r', d, first, pred): {:?}; order(t_chall) == order(t_aux) == 2^max(2m-F,0) for every m in [f, F] on {}/100 vectors",
            VARIANTS[vi], ok, mism, ord_ok
        );
    }

    // All vectors at m = F: order-level dlog bits for the decompressor's points
    // (variant a) and the shifted comparison for the order-2^F points (variant b).
    println!("--- all vectors at m = F: a:t(P,Q) ord_chall, ord_aux, dlog@ord k=[+,-,+sh,-sh]; b:t(R,S) dlog@2^F k=[+,-,+sh,-sh] ---");
    let mut hist_a: BTreeMap<(usize, usize), BTreeMap<(usize, usize, usize), usize>> = BTreeMap::new();
    for (idx, r) in results.iter().enumerate() {
        let (_, cells) = r.rows.iter().find(|(m, _)| *m == fc).expect("m = F row");
        let a = cells[0];
        let b = cells[2];
        println!(
            "KAT {:3}: bt={} trl={} f={} g={} d={} | a@F: ordC={} ordA={} dlog@ord ok={} k_o={:?} | b@F: ordC={} ordA={} dlog ok={} k={:?}",
            idx, r.bt, r.trl, r.f, r.g, r.d,
            ord_str(a.ord_c), ord_str(a.ord_a), a.dlog_o_ok, a.k_o,
            ord_str(b.ord_c), ord_str(b.ord_a), b.dlog_ok, b.k
        );
        let key = (a.ord_c.unwrap_or(NOT_POW2), a.ord_a.unwrap_or(NOT_POW2), if a.dlog_o_ok { a.k_o[0] } else { usize::MAX });
        *hist_a.entry((r.bt, r.trl)).or_default().entry(key).or_default() += 1;
    }
    println!("--- variant a, histogram by (n_bt, r') of (log2 ord_chall, log2 ord_aux, k_o[y=+det]) at m = F (k_o = MAX means dlog@ord failed); 2f-F and 2g-F for reference ---");
    for ((bt, trl), h) in &hist_a {
        let f = L::E_RSP as usize + 2 - bt;
        let g = f - trl;
        println!(
            "  (n_bt={}, r'={}): 2f-F={} 2g-F={} : {:?}",
            bt, trl, 2 * f as i64 - fc as i64, 2 * g as i64 - fc as i64, h
        );
    }
    let all_b_ok = results.iter().all(|r| {
        let (_, cells) = r.rows.iter().find(|(m, _)| *m == fc).unwrap();
        cells[2].dlog_ok && cells[2].k[2] >= r.d && cells[3].dlog_ok && cells[3].k[2] >= r.d
    });
    println!("variant b at m = F: 2^r' * y == det M mod 2^d on all 100 vectors (both argument orders): {}", all_b_ok);

    // Constancy in m of the variant-a cells, per vector, over the whole sweep [f, F].
    let same = |x: &Cell, y: &Cell| {
        x.ord_c == y.ord_c
            && x.ord_a == y.ord_a
            && x.dlog_ok == y.dlog_ok
            && x.k == y.k
            && x.dlog_o_ok == y.dlog_o_ok
            && x.k_o == y.k_o
    };
    let const_a = results
        .iter()
        .filter(|r| {
            let first = &r.rows[0].1;
            r.rows.iter().all(|(_, c)| same(&c[0], &first[0]) && same(&c[1], &first[1]))
        })
        .count();
    println!(
        "variant a: vectors whose t(P,Q) and t(Q,P) cells (orders, dlog results, k, k_o) are identical for every m in [f, F]: {}/100",
        const_a
    );
    let ord_pred = results
        .iter()
        .filter(|r| {
            let c = &r.rows[0].1[0];
            c.ord_c == Some((2 * r.f as i64 - fc as i64).max(0) as usize)
                && c.ord_a == Some((2 * r.g as i64 - fc as i64).max(0) as usize)
        })
        .count();
    println!(
        "variant a: vectors with log2 ord(t_m(P,Q)) == max(2f-F,0) and log2 ord(t_m(Pa,Qa)) == max(2g-F,0): {}/100",
        ord_pred
    );

    // Per-vector detail for the first 3 vectors at m = f, m = f + 1, m = F.
    println!("--- per-vector detail (vectors 0..3) ---");
    for (idx, r) in results.iter().enumerate().take(3) {
        println!("KAT {}: bt={} trl={} f={} g={} d={}", idx, r.bt, r.trl, r.f, r.g, r.d);
        for (m, cells) in r.rows.iter().filter(|(m, _)| *m == r.f || *m == r.f + 1 || *m == fc) {
            for (vi, c) in cells.iter().enumerate() {
                println!(
                    "   m={} {}: ord_chall={} ord_aux={} in_mu={} dlog@2^m ok={} k={:?} | dlog@ord ok={} k={:?}",
                    m,
                    VARIANTS[vi],
                    ord_str(c.ord_c),
                    ord_str(c.ord_a),
                    c.in_mu_2m,
                    c.dlog_ok,
                    c.k,
                    c.dlog_o_ok,
                    c.k_o
                );
            }
        }
    }
}

#[test]
fn task2_tate_sweep_l1() {
    run_task2::<Level1>("Level 1 (p = 5*2^248 - 1, F = 248)", KAT_LVL1);
}

#[test]
fn task2_tate_sweep_l3() {
    run_task2::<Level3>("Level 3 (p = 65*2^376 - 1, F = 376)", KAT_LVL3);
}

#[test]
fn task2_tate_sweep_l5() {
    run_task2::<Level5>("Level 5 (p = 27*2^500 - 1, F = 500)", KAT_LVL5);
}

// ---------------------------------------------------------------------------
// Timing, Level 1, KAT entry 0
// ---------------------------------------------------------------------------

fn median_us<F: FnMut()>(iters: usize, mut f: F) -> f64 {
    let mut v = Vec::with_capacity(iters);
    for _ in 0..iters {
        let t = Instant::now();
        f();
        v.push(t.elapsed().as_nanos() as f64 / 1000.0);
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[iters / 2]
}

#[test]
fn timing_l1_kat0() {
    type L = Level1;
    let entries = parse_kat_entries(KAT_LVL1);
    let s = setup::<L>(0, &entries[0]);
    let fc = L::F_CHR as usize;
    let cof = L::p_cofactor_for_2f();
    let iters = 201;
    let m_star: usize = std::env::var("PAPER_TATE_M")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(s.f);
    let mut e_chall = s.e_chall.clone();

    let w = median_us(iters, || {
        let r = weil::<L>(s.f as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall);
        std::hint::black_box(r);
    });
    let t_f = median_us(iters, || {
        let r = reduced_tate::<L>(s.f as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall, fc as u32, cof);
        std::hint::black_box(r);
    });
    let t_star = median_us(iters, || {
        let r = reduced_tate::<L>(m_star as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall, fc as u32, cof);
        std::hint::black_box(r);
    });
    let t_full_a = median_us(iters, || {
        let r = reduced_tate::<L>(fc as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall, fc as u32, cof);
        std::hint::black_box(r);
    });
    let t_full_b = median_us(iters, || {
        let r = reduced_tate::<L>(fc as u32, &s.b_chall_full.p, &s.b_chall_full.q, &s.ppq_chall_full, &mut e_chall, fc as u32, cof);
        std::hint::black_box(r);
    });
    let m_c: usize = std::env::var("PAPER_TATE_MC")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(fc);
    let mut e_aux = s.e_aux.clone();
    let dbl_only = median_us(iters, || {
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m_c, &mut e_chall);
        std::hint::black_box(bc);
    });
    let t_c_one = median_us(iters, || {
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m_c, &mut e_chall);
        let ppqc = xadd(&bc.p, &bc.q, &bc.pmq);
        let r = reduced_tate::<L>(m_c as u32, &bc.p, &bc.q, &ppqc, &mut e_chall, fc as u32, cof);
        std::hint::black_box(r);
    });
    let t_c_both = median_us(iters, || {
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - m_c, &mut e_chall);
        let ppqc = xadd(&bc.p, &bc.q, &bc.pmq);
        let r1 = reduced_tate::<L>(m_c as u32, &bc.p, &bc.q, &ppqc, &mut e_chall, fc as u32, cof);
        let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - m_c, &mut e_aux);
        let ppqa = xadd(&ba.p, &ba.q, &ba.pmq);
        let r2 = reduced_tate::<L>(m_c as u32, &ba.p, &ba.q, &ppqa, &mut e_aux, fc as u32, cof);
        std::hint::black_box((r1, r2));
    });
    let w_both = median_us(iters, || {
        let r1 = weil::<L>(s.f as u32, &s.b_chall.p, &s.b_chall.q, &s.ppq_chall, &mut e_chall);
        let r2 = weil::<L>(s.g as u32, &s.b_aux.p, &s.b_aux.q, &s.ppq_aux, &mut e_aux);
        std::hint::black_box((r1, r2));
    });
    // Production decompressor stage: F - f doublings of the E_chall basis triple,
    // F - g doublings of the E_aux basis triple, then the two Weil pairings
    // (formats.rs:657-668), measured like for like with the variant c stage.
    let prod_dbl = median_us(iters, || {
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - s.f, &mut e_chall);
        let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - s.g, &mut e_aux);
        std::hint::black_box((bc, ba));
    });
    let prod_stage = median_us(iters, || {
        let bc = ec_dbl_iter_basis(&s.b_chall_full, fc - s.f, &mut e_chall);
        let ppqc = xadd(&bc.p, &bc.q, &bc.pmq);
        let r1 = weil::<L>(s.f as u32, &bc.p, &bc.q, &ppqc, &mut e_chall);
        let ba = ec_dbl_iter_basis(&s.b_aux_full, fc - s.g, &mut e_aux);
        let ppqa = xadd(&ba.p, &ba.q, &ba.pmq);
        let r2 = weil::<L>(s.g as u32, &ba.p, &ba.q, &ppqa, &mut e_aux);
        std::hint::black_box((r1, r2));
    });
    println!("==================== TIMING: Level 1, KAT entry 0 (bt={}, trl={}, f={}, F={}), medians of {} iterations ====================", s.bt, s.trl, s.f, fc, iters);
    println!("production Weil e_{{2^f}}(P,Q) on E_chall (two cubical ladders of 2^{} steps): {:.1} us", s.f, w);
    println!("single Tate ladder, m = f = {} (order-2^f points): {:.1} us", s.f, t_f);
    println!("single Tate ladder, m = PAPER_TATE_M = {} (order-2^f points): {:.1} us", m_star, t_star);
    println!("single Tate ladder, m = F = {} (order-2^f points): {:.1} us", fc, t_full_a);
    println!("single Tate ladder, m = F = {} (order-2^F points R,S): {:.1} us", fc, t_full_b);
    println!("variant c at m = PAPER_TATE_MC = {}: {} doublings of the E_chall basis triple (P, Q, P-Q) alone: {:.1} us; doublings + x(P+Q) + one Tate ladder of 2^{} steps on E_chall: {:.1} us; same on both curves (four basis points R,S,R_aux,S_aux doubled, two ladders): {:.1} us", m_c, fc - m_c, dbl_only, m_c, t_c_one, t_c_both);
    println!("production Weil on both curves (e_{{2^f}} on E_chall + e_{{2^g}} on E_aux, four ladders): {:.1} us", w_both);
    println!("production doublings alone ({} on the E_chall triple + {} on the E_aux triple): {:.1} us; production pairing stage (those doublings + x(P+Q) + both Weil pairings): {:.1} us", fc - s.f, fc - s.g, prod_dbl, prod_stage);
}
