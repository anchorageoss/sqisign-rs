//! Verification cost breakdown for the Weil pairing compression note on
//! round 3: level I, every component of the standard path and of the
//! compressed path timed on a copy of the state its own path feeds it, on
//! each of the 100 KAT vectors (the median of 200 iterations after 5
//! warm-ups per vector), reported as the median over the vectors. The
//! same per-vector method times the Tate alternative to the pairing stage
//! (the bases from the hints at order `2^m*`, `m* = ceil((F + g)/2)`, one
//! Tate ladder per curve, and the reduction of the two points of each
//! basis to order `2^N` with the difference point, which gives the
//! verifier's bases), so that the paper's pairing-stage comparison and its
//! cost table rest on one run.
//!
//! Both paths run under a prepared public key (`PreparedPublicKey`, the
//! per-key basis work done once). The standard path is
//! `sqisign::verify_prepared` on a decoded signature; the compressed path is
//! `compressed::verify_compressed_prepared`, which is integrated: it derives
//! the challenge curve and the two bases once, adds the two Weil pairings,
//! the discrete logarithm and the recovery of the dropped entry, and
//! continues as the standard path on the recovered matrix.
//!
//! Run in release from a copy in `crates/kat/tests` (see the README):
//!   cargo test -p sqisign-kat --release --test paper_timing -- --nocapture --test-threads=1

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::Instant;

use sqisign_verify::compressed::{
    compress, compressed_from_bytes, compressed_to_bytes, verify_compressed_prepared,
    CompressedSignature, MAX_COMPRESSED_BYTES,
};
use sqisign_verify::ec::basis::{difference_point, ec_curve_to_basis_2f_from_hint};
use sqisign_verify::ec::isogeny::iso_isogeny_2chain_from;
use sqisign_verify::ec::pairing::{fp2_dlog_2e, reduced_tate, weil};
use sqisign_verify::ec::point::{ec_biscalar_mul_verif, ec_dbl_iter, ec_ladder3pt};
use sqisign_verify::ec::{EcBasis, EcCurve, MAX_ORDER_WORDS};
use sqisign_verify::fp::Fp2;
use sqisign_verify::precomp::EXTRA_TORSION;
use sqisign_verify::sqisign::{
    hash_to_challenge, public_key_from_bytes, signature_from_bytes, verify_prepared,
    PreparedPublicKey, PublicKey, Scalar, VerifyParams, HD_EXTRA_TORSION,
};
use sqisign_verify::theta::chain::theta_chain_compute_and_eval;
use sqisign_verify::theta::{ChainMode, ThetaCoupleCurve, ThetaKernelCouplePoints};
use sqisign_verify::{Level, Level1};

const KAT_L1: &str = include_str!("../kat/PQCsignKAT_270_SQIsign_p324_3.rsp");

const ITERS: usize = 200;

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
// Timing
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Stat {
    median_us: f64,
}

/// Time `f` over `ITERS` iterations after 5 warm-ups. `setup` runs outside
/// the timed region and gives each iteration a fresh input (the routines
/// mutate their curves).
fn measure<I, T>(mut setup: impl FnMut() -> I, mut f: impl FnMut(I) -> T) -> Stat {
    for _ in 0..5 {
        black_box(f(setup()));
    }
    let mut samples = Vec::with_capacity(ITERS);
    for _ in 0..ITERS {
        let input = setup();
        let t0 = Instant::now();
        let out = f(input);
        let dt = t0.elapsed();
        black_box(out);
        samples.push(dt.as_secs_f64() * 1e6);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Stat {
        median_us: (samples[ITERS / 2 - 1] + samples[ITERS / 2]) / 2.0,
    }
}

// ---------------------------------------------------------------------------
// The verifier's internals that are not public, copied: the matrix check
// (sqisign::check_canonical_basis_change_matrix) and the recovery of the
// dropped entry (compressed::reconstruct and its arithmetic modulo 2^t)
// ---------------------------------------------------------------------------

fn bit_length(x: &Scalar) -> u32 {
    for i in (0..MAX_ORDER_WORDS).rev() {
        if x[i] != 0 {
            return 64 * i as u32 + 64 - x[i].leading_zeros();
        }
    }
    0
}

fn is_power_of_two_at(x: &Scalar, k: u32) -> bool {
    bit_length(x) == k + 1 && x.iter().map(|w| w.count_ones()).sum::<u32>() == 1
}

/// Every entry below `2^(RESPONSE_BITS + 2)` and the matrix normalised.
fn check_matrix(params: &VerifyParams, m: &[[Scalar; 2]; 2]) -> bool {
    let mid_bits = params.response_bits + HD_EXTRA_TORSION - 1;
    let mut normal = false;
    for row in m.iter() {
        for x in row.iter() {
            let bits = bit_length(x);
            if bits > mid_bits + 1 {
                return false;
            }
            if normal {
                continue;
            }
            if bits == mid_bits + 1 && !is_power_of_two_at(x, mid_bits) {
                return false;
            }
            normal = bits <= mid_bits && bits != 0;
        }
    }
    true
}

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

/// `det M = y + 2^RESPONSE_BITS mod 2^(RESPONSE_BITS + 1)` and the matrix
/// from the three transmitted entries (compressed::recover's tail).
fn recover_entry<L: Level>(
    params: &VerifyParams,
    c: &CompressedSignature<L>,
    y: &Scalar,
) -> Option<[[Scalar; 2]; 2]> {
    let r = params.response_bits;
    let t = r + 1;
    let mut two_r = [0u64; MAX_ORDER_WORDS];
    two_r[(r / 64) as usize] |= 1u64 << (r % 64);
    let det = add_mod_2exp(y, &two_r, t);
    let (m00, m01, m_var) = (c.m00, c.m01, c.m_var);
    let (m10, m11) = if m00[0] & 1 == 1 {
        let num = add_mod_2exp(&det, &mul_mod_2exp(&m01, &m_var, t), t);
        (m_var, mul_mod_2exp(&num, &inv_mod_2exp(&m00, t)?, t))
    } else {
        let num = sub_mod_2exp(&mul_mod_2exp(&m00, &m_var, t), &det, t);
        (mul_mod_2exp(&num, &inv_mod_2exp(&m01, t)?, t), m_var)
    };
    Some([[m00, m01], [m10, m11]])
}

// ---------------------------------------------------------------------------
// The shared part of both paths, step by step
// ---------------------------------------------------------------------------

/// The curve and basis of `E_pk` a prepared key holds.
fn pk_basis<L: Level>(params: &VerifyParams, pk: &PublicKey<L>) -> (EcCurve<L>, EcBasis<L>) {
    let mut curve = pk.curve.clone();
    let f_chall = params.challenge_bits + EXTRA_TORSION;
    let basis = ec_curve_to_basis_2f_from_hint(&mut curve, f_chall, pk.hint_pk).unwrap();
    (curve, basis)
}

fn challenge_curve<L: Level>(
    params: &VerifyParams,
    e_pk: &EcCurve<L>,
    b_pk: &EcBasis<L>,
    chall: &Scalar,
) -> EcCurve<L> {
    let f_chall = params.challenge_bits + EXTRA_TORSION;
    let kernel = ec_ladder3pt(&chall[..L::ORDER_WORDS], &b_pk.p, &b_pk.q, &b_pk.pmq, e_pk).unwrap();
    iso_isogeny_2chain_from(e_pk, &kernel, f_chall).unwrap()
}

fn basis_change<L: Level>(
    params: &VerifyParams,
    m: &[[Scalar; 2]; 2],
    b: &EcBasis<L>,
    e: &mut EcCurve<L>,
) -> EcBasis<L> {
    let n = (params.response_bits + HD_EXTRA_TORSION) as usize;
    let w = L::ORDER_WORDS;
    let mut out = b.clone();
    out.p = ec_biscalar_mul_verif(&m[0][0][..w], &m[1][0][..w], n, b, e).unwrap();
    out.q = ec_biscalar_mul_verif(&m[0][1][..w], &m[1][1][..w], n, b, e).unwrap();
    out
}

fn theta<L: Level>(
    params: &VerifyParams,
    e_chall: &EcCurve<L>,
    e_aux: &EcCurve<L>,
    b_chall: &EcBasis<L>,
    b_aux: &EcBasis<L>,
) -> EcCurve<L> {
    let mut e12 = ThetaCoupleCurve {
        e1: e_chall.clone(),
        e2: e_aux.clone(),
    };
    let ker = ThetaKernelCouplePoints::from_bases(b_chall, b_aux);
    let out = theta_chain_compute_and_eval(
        params.response_bits as u16,
        &mut e12,
        &ker,
        &[],
        ChainMode::Verify,
    )
    .unwrap();
    out.codomain.e1
}

/// The per-vector medians of one path's components, in table order.
#[derive(Clone, Copy, Default)]
struct Comps {
    decode: f64,
    aux_curve: f64,
    challenge: f64,
    bases: f64,
    pairings: f64,
    dlog: f64,
    recovery: f64,
    check: f64,
    change: f64,
    chain: f64,
    hash: f64,
    full: f64,
    /// the Tate alternative to the pairing stage (compressed path only)
    tate_bases: f64,
    tate_ladders: f64,
    tate_reduce: f64,
}

impl Comps {
    fn sum(&self) -> f64 {
        self.decode
            + self.aux_curve
            + self.challenge
            + self.bases
            + self.pairings
            + self.dlog
            + self.recovery
            + self.check
            + self.change
            + self.chain
            + self.hash
    }
}

/// The components shared by both paths, each timed on the state its
/// predecessors produced along that path, from the path's own fields.
fn shared<L: Level>(
    params: &VerifyParams,
    pk: &PublicKey<L>,
    e_pk: &EcCurve<L>,
    b_pk: &EcBasis<L>,
    e_aux_a: &Fp2<L>,
    chall: &Scalar,
    hint_chall: u8,
    hint_aux: u8,
    m: &[[Scalar; 2]; 2],
    msg: &[u8],
    c: &mut Comps,
) {
    let n = params.response_bits + HD_EXTRA_TORSION;
    assert!(check_matrix(params, m));
    c.check = measure(|| (), |_| check_matrix(params, m)).median_us;
    c.aux_curve = measure(|| (), |_| EcCurve::<L>::from_a(e_aux_a).unwrap()).median_us;
    let e_aux0 = EcCurve::<L>::from_a(e_aux_a).unwrap();
    c.challenge = measure(|| (), |_| challenge_curve(params, e_pk, b_pk, chall)).median_us;
    let e_chall0 = challenge_curve(params, e_pk, b_pk, chall);
    c.bases = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut ec, mut ea)| {
            (
                ec_curve_to_basis_2f_from_hint(&mut ec, n, hint_chall).unwrap(),
                ec_curve_to_basis_2f_from_hint(&mut ea, n, hint_aux).unwrap(),
            )
        },
    )
    .median_us;
    let mut e_chall = e_chall0.clone();
    let mut e_aux = e_aux0.clone();
    let b_chall = ec_curve_to_basis_2f_from_hint(&mut e_chall, n, hint_chall).unwrap();
    let b_aux = ec_curve_to_basis_2f_from_hint(&mut e_aux, n, hint_aux).unwrap();
    c.change = measure(
        || e_chall.clone(),
        |mut ec| basis_change(params, m, &b_chall, &mut ec),
    )
    .median_us;
    let mut e_chall2 = e_chall.clone();
    let b_chall2 = basis_change(params, m, &b_chall, &mut e_chall2);
    c.chain = measure(
        || (),
        |_| theta(params, &e_chall2, &e_aux, &b_chall2, &b_aux),
    )
    .median_us;
    let e_com = theta(params, &e_chall2, &e_aux, &b_chall2, &b_aux);
    assert_eq!(hash_to_challenge(params, pk, &e_com, msg), *chall);
    c.hash = measure(
        || (),
        |_| hash_to_challenge(params, pk, &e_com, msg) == *chall,
    )
    .median_us;
}

/// Both paths on one vector.
fn vector_breakdown<L: Level>(params: &VerifyParams, e: &Entry) -> (Comps, Comps) {
    let msg: &[u8] = &e.msg;
    let pk = public_key_from_bytes::<L>(params, &e.pk).unwrap();
    let sig_bytes = &e.sm[..params.signature_bytes];
    let sig = signature_from_bytes::<L>(params, sig_bytes).unwrap();
    let comp = compress(params, &sig).unwrap();
    let mut cbuf = [0u8; MAX_COMPRESSED_BYTES];
    let clen = compressed_to_bytes(params, &comp, &mut cbuf).unwrap();
    let cbytes = &cbuf[..clen];
    let key = PreparedPublicKey::new(params, &pk).unwrap();
    let (e_pk, b_pk) = pk_basis(params, &pk);
    let n = params.response_bits + HD_EXTRA_TORSION;

    // ---- the standard path, on the decoded signature's state
    let mut s = Comps::default();
    s.decode = measure(
        || (),
        |_| signature_from_bytes::<L>(params, sig_bytes).unwrap(),
    )
    .median_us;
    shared(
        params,
        &pk,
        &e_pk,
        &b_pk,
        &sig.e_aux_a,
        &sig.chall_coeff,
        sig.hint_chall,
        sig.hint_aux,
        &sig.mat_bchall_can_to_bchall,
        msg,
        &mut s,
    );
    s.full = measure(
        || (),
        |_| {
            let s = signature_from_bytes::<L>(params, sig_bytes).unwrap();
            assert!(verify_prepared(params, &key, &s, msg))
        },
    )
    .median_us;

    // ---- the compressed path, on the decoded compressed signature's state:
    // the pairings, the discrete logarithm and the recovery come first,
    // the shared components then run on the recovered matrix
    let mut c = Comps::default();
    c.decode = measure(
        || (),
        |_| compressed_from_bytes::<L>(params, cbytes).unwrap(),
    )
    .median_us;
    let cs = compressed_from_bytes::<L>(params, cbytes).unwrap();
    let mut ec = challenge_curve(params, &e_pk, &b_pk, &cs.chall_coeff);
    let mut ea = EcCurve::<L>::from_a(&cs.e_aux_a).unwrap();
    let bc = ec_curve_to_basis_2f_from_hint(&mut ec, n, cs.hint_chall).unwrap();
    let ba = ec_curve_to_basis_2f_from_hint(&mut ea, n, cs.hint_aux).unwrap();
    c.pairings = measure(
        || (ec.clone(), ea.clone()),
        |(mut e1, mut e2)| {
            (
                weil::<L>(n, &bc.p, &bc.q, &bc.pmq, &mut e1),
                weil::<L>(n, &ba.p, &ba.q, &ba.pmq, &mut e2),
            )
        },
    )
    .median_us;
    // the Tate alternative: bases at 2^m*, one ladder per curve, and the
    // two points of each basis brought down to order 2^N with the
    // difference point at that order (the verifier's bases)
    let fc = L::TWO_ADIC_EXPONENT;
    let m_star = (fc + n - 1).div_ceil(2);
    let e_chall0 = challenge_curve(params, &e_pk, &b_pk, &cs.chall_coeff);
    let e_aux0 = EcCurve::<L>::from_a(&cs.e_aux_a).unwrap();
    c.tate_bases = measure(
        || (e_chall0.clone(), e_aux0.clone()),
        |(mut e1, mut e2)| {
            (
                ec_curve_to_basis_2f_from_hint(&mut e1, m_star, cs.hint_chall).unwrap(),
                ec_curve_to_basis_2f_from_hint(&mut e2, m_star, cs.hint_aux).unwrap(),
            )
        },
    )
    .median_us;
    let (mut e3, mut e4) = (e_chall0.clone(), e_aux0.clone());
    let bc_m = ec_curve_to_basis_2f_from_hint(&mut e3, m_star, cs.hint_chall).unwrap();
    let ba_m = ec_curve_to_basis_2f_from_hint(&mut e4, m_star, cs.hint_aux).unwrap();
    c.tate_ladders = measure(
        || (e3.clone(), e4.clone()),
        |(mut e1, mut e2)| {
            (
                reduced_tate::<L>(m_star, &bc_m.p, &bc_m.q, &bc_m.pmq, &mut e1),
                reduced_tate::<L>(m_star, &ba_m.p, &ba_m.q, &ba_m.pmq, &mut e2),
            )
        },
    )
    .median_us;
    let reduce = |bm: &EcBasis<L>, e: &mut EcCurve<L>| {
        let p = ec_dbl_iter(&bm.p, (m_star - n) as usize, e);
        let q = ec_dbl_iter(&bm.pmq, (m_star - n) as usize, e);
        let r = difference_point(&p, &q, e);
        EcBasis::new(p, r, q)
    };
    c.tate_reduce = measure(
        || (e3.clone(), e4.clone()),
        |(mut e1, mut e2)| (reduce(&bc_m, &mut e1), reduce(&ba_m, &mut e2)),
    )
    .median_us;
    let (rc, ra) = (
        reduce(&bc_m, &mut e3.clone()),
        reduce(&ba_m, &mut e4.clone()),
    );
    let same = |a: &EcBasis<L>, b: &EcBasis<L>| {
        bool::from(a.p.ct_equal(&b.p) & a.q.ct_equal(&b.q) & a.pmq.ct_equal(&b.pmq))
    };
    assert!(
        same(&rc, &bc) && same(&ra, &ba),
        "the reduction gives the verifier's bases"
    );
    let (mut e1, mut e2) = (ec.clone(), ea.clone());
    let g = weil::<L>(n, &bc.p, &bc.q, &bc.pmq, &mut e1);
    let h = weil::<L>(n, &ba.p, &ba.q, &ba.pmq, &mut e2);
    c.dlog = measure(
        || [0u64; MAX_ORDER_WORDS],
        |mut y| {
            fp2_dlog_2e::<L>(&mut y[..L::ORDER_WORDS], &h.inv(), &g.inv(), n).unwrap();
            y
        },
    )
    .median_us;
    let mut y = [0u64; MAX_ORDER_WORDS];
    fp2_dlog_2e::<L>(&mut y[..L::ORDER_WORDS], &h.inv(), &g.inv(), n).unwrap();
    c.recovery = measure(|| (), |_| recover_entry(params, &cs, &y).unwrap()).median_us;
    let m_rec = recover_entry(params, &cs, &y).unwrap();
    // the recovered matrix equals the signer's below bit RESPONSE_BITS + 1
    let t = params.response_bits + 1;
    for (a, b) in m_rec
        .iter()
        .flatten()
        .zip(sig.mat_bchall_can_to_bchall.iter().flatten())
    {
        assert_eq!(reduced(a, t), reduced(b, t), "recovery");
    }
    shared(
        params,
        &pk,
        &e_pk,
        &b_pk,
        &cs.e_aux_a,
        &cs.chall_coeff,
        cs.hint_chall,
        cs.hint_aux,
        &m_rec,
        msg,
        &mut c,
    );
    c.full = measure(
        || (),
        |_| {
            let c = compressed_from_bytes::<L>(params, cbytes).unwrap();
            assert!(verify_compressed_prepared(params, &key, &c, msg))
        },
    )
    .median_us;
    (s, c)
}

/// Median, minimum and maximum over the vectors of a per-vector median.
fn over_vectors(v: &mut [f64]) -> (f64, f64, f64) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    let med = if n % 2 == 0 {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    } else {
        v[n / 2]
    };
    (med, v[0], v[n - 1])
}

fn breakdown<L: Level>(name: &str, rsp: &str) -> String {
    let params = &L::PARAMS;
    let entries = parse_rsp(rsp);
    assert_eq!(entries.len(), 100);
    let mut s_all = Vec::new();
    let mut c_all = Vec::new();
    for (i, e) in entries.iter().enumerate() {
        assert_eq!(e.count, i, "KAT entries in order");
        let (s, c) = vector_breakdown::<L>(params, e);
        s_all.push(s);
        c_all.push(c);
    }
    let mut out = String::new();
    let w = &mut out;
    writeln!(
        w,
        "==== {name}: {} KAT vectors, e_rsp = {}, N = {}, signature {} bytes, compressed {} bytes; per vector the median of {ITERS} iterations after 5 warm-ups, then over the vectors the median (min, max); kernels: {} ====",
        entries.len(),
        params.response_bits,
        params.response_bits + HD_EXTRA_TORSION,
        params.signature_bytes,
        compressed_bytes_of(params),
        sqisign_kat::arithmetic_from_env()
    )
    .unwrap();
    writeln!(
        w,
        "{:<50} {:>8} {:>8} {:>8}   {:>8} {:>8} {:>8}",
        "component (microseconds)", "S med", "S min", "S max", "C med", "C min", "C max"
    )
    .unwrap();
    let rows: [(&str, fn(&Comps) -> f64, bool); 11] = [
        ("decoding", |c| c.decode, true),
        ("E_aux from A_aux", |c| c.aux_curve, true),
        ("challenge derivation E_pk -> E_chl", |c| c.challenge, true),
        (
            "bases from the hints at 2^N, both curves",
            |c| c.bases,
            true,
        ),
        ("two Weil pairings at 2^N", |c| c.pairings, false),
        (
            "two inversions and the discrete logarithm",
            |c| c.dlog,
            false,
        ),
        ("det M and the dropped entry", |c| c.recovery, false),
        ("matrix range and normal-form check", |c| c.check, true),
        ("basis change, apply M", |c| c.change, true),
        ("(2^e_rsp, 2^e_rsp) theta chain to E_com", |c| c.chain, true),
        ("challenge hash and compare", |c| c.hash, true),
    ];
    let mut med_s = Comps::default();
    let mut med_c = Comps::default();
    let mut sum_s = 0.0;
    let mut sum_c = 0.0;
    let agg = |f: fn(&Comps) -> f64, all: &[Comps]| {
        let mut v: Vec<f64> = all.iter().map(f).collect();
        over_vectors(&mut v)
    };
    for (label, f, both) in rows {
        let (cm, cmin, cmax) = agg(f, &c_all);
        sum_c += cm;
        if both {
            let (sm, smin, smax) = agg(f, &s_all);
            sum_s += sm;
            writeln!(
                w,
                "{:<50} {:>8.1} {:>8.1} {:>8.1}   {:>8.1} {:>8.1} {:>8.1}",
                label, sm, smin, smax, cm, cmin, cmax
            )
            .unwrap();
        } else {
            writeln!(
                w,
                "{:<50} {:>8} {:>8} {:>8}   {:>8.1} {:>8.1} {:>8.1}",
                label, "-", "-", "-", cm, cmin, cmax
            )
            .unwrap();
        }
        match label {
            "decoding" => {
                med_c.decode = cm;
                med_s.decode = agg(f, &s_all).0;
            }
            "two Weil pairings at 2^N" => med_c.pairings = cm,
            "two inversions and the discrete logarithm" => med_c.dlog = cm,
            "det M and the dropped entry" => med_c.recovery = cm,
            _ => {}
        }
    }
    let tate_rows: [(&str, fn(&Comps) -> f64); 3] = [
        ("Tate: bases from the hints at 2^m*, both curves", |c| {
            c.tate_bases
        }),
        ("Tate: two ladders at 2^m*", |c| c.tate_ladders),
        ("Tate: reduction of the bases to 2^N", |c| c.tate_reduce),
    ];
    let mut tate_stage = 0.0;
    for (label, f) in tate_rows {
        let (cm, cmin, cmax) = agg(f, &c_all);
        tate_stage += cm;
        writeln!(
            w,
            "{:<50} {:>8} {:>8} {:>8}   {:>8.1} {:>8.1} {:>8.1}",
            label, "-", "-", "-", cm, cmin, cmax
        )
        .unwrap();
    }
    let (full_s, fs_min, fs_max) = agg(|c| c.full, &s_all);
    let (full_c, fc_min, fc_max) = agg(|c| c.full, &c_all);
    writeln!(
        w,
        "{:<50} {:>8.1} {:>8} {:>8}   {:>8.1} {:>8} {:>8}",
        "sum of component medians", sum_s, "", "", sum_c, "", ""
    )
    .unwrap();
    writeln!(
        w,
        "{:<50} {:>8.1} {:>8.1} {:>8.1}   {:>8.1} {:>8.1} {:>8.1}",
        "FULL verification from bytes, prepared key",
        full_s,
        fs_min,
        fs_max,
        full_c,
        fc_min,
        fc_max
    )
    .unwrap();
    let added = med_c.pairings + med_c.dlog + med_c.recovery + med_c.decode - med_s.decode;
    writeln!(
        w,
        "integrated overhead, the components the compressed path adds (pairings {:.1} + inversions and dlog {:.1} + recovery {:.1} + decoding {:.1} - {:.1}): {:.1} us, {:.1}% of the standard verification ({:.1} us); full medians differ by {:.1} us ({:.1}%)",
        med_c.pairings,
        med_c.dlog,
        med_c.recovery,
        med_c.decode,
        med_s.decode,
        added,
        100.0 * added / full_s,
        full_s,
        full_c - full_s,
        100.0 * (full_c - full_s) / full_s
    )
    .unwrap();
    let weil_stage = agg(|c| c.bases, &c_all).0 + med_c.pairings;
    let saving = weil_stage - tate_stage;
    writeln!(
        w,
        "pairing stage, medians over the vectors: Weil (bases at 2^N + two Weil pairings, the two table rows) {weil_stage:.1} us; Tate (bases at 2^m* + two ladders + reduction to 2^N, m* = {}) {tate_stage:.1} us; saving {saving:.1} us ({:.1}%); overhead with the Tate ladders {:.1} us, {:.1}% of the standard verification",
        (L::TWO_ADIC_EXPONENT + params.response_bits + 1).div_ceil(2),
        100.0 * saving / weil_stage,
        added - saving,
        100.0 * (added - saving) / full_s
    )
    .unwrap();
    let _ = med_s.sum();
    out
}

fn compressed_bytes_of(params: &VerifyParams) -> usize {
    sqisign_verify::compressed::compressed_bytes(params)
}

#[test]
fn paper_timing_level1_breakdown() {
    print!("{}", breakdown::<Level1>("level I (p324_3)", KAT_L1));
}
