//! Statistics over the official NIST round-2 KAT vectors for the compression
//! paper (Task 1). One test per level; each prints tables (a) to (g).
//!
//! Run with:
//!   CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 \
//!   cargo test -p sqisign-kat --release --features kat-compat \
//!     --test paper_stats -- --nocapture --test-threads=1

#![cfg(feature = "kat-compat")]

use std::fmt::Write as _;

use rayon::prelude::*;
use sqisign_kat::kat_parser::{self, KatEntry};
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::EcCurve;
use sqisign_verify::theta::HD_EXTRA_TORSION;
use sqisign_verify::verify::compute_challenge_curve;
use sqisign_verify::{
    CompressedSignature, ExpandedSignature, FpBackend, Level1, Level3, Level5, LevelPrecomp,
    PublicKey, Signature, Verifier,
};

const KAT_L1: &str = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
const KAT_L3: &str = include_str!("../../../reference/KAT/PQCsignKAT_529_SQIsign_lvl3.rsp");
const KAT_L5: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

/// Candidate bit positions for a "leading 1 in the expected position" reading
/// of M[0][0]; the per-vector positions are reported alongside the fixed ones.
const CANDIDATE_LABELS: [&str; 6] = [
    "E_RSP+1 (fixed, = N-1 when n_bt=0)",
    "E_RSP-1 (fixed)",
    "E_RSP-2 (fixed)",
    "N-1 (per vector, N = E_RSP+2-n_bt)",
    "N-2 (per vector)",
    "g-1 (per vector, g = pow_dim2+2 = E_RSP-r'-n_bt+2)",
];

/// Bit `i` of a little-endian limb array.
fn bit(d: &[u64], i: usize) -> bool {
    (d[i / 64] >> (i % 64)) & 1 == 1
}

/// Position of the leading one plus one (0 for zero).
fn bitlen(d: &[u64]) -> usize {
    for i in (0..d.len()).rev() {
        if d[i] != 0 {
            return i * 64 + 64 - d[i].leading_zeros() as usize;
        }
    }
    0
}

struct VecStats {
    count: usize,
    bt: u8,
    trl: u8,
    /// N = E_RSP + HD_EXTRA_TORSION - bt, the modulus width the verifier
    /// enforces on every matrix entry (check_canonical_basis_change_matrix).
    n_bits: usize,
    m00_even: bool,
    /// det_hint as the compressor wrote it into the metadata byte.
    det_hint_wire: u8,
    /// det_hint recomputed directly from bits [N-2, N) of the dropped entry.
    det_hint_direct: u8,
    m00_top_bit: bool,
    any_top_bit: bool,
    /// bit length of each entry, row-major.
    bitlen: [usize; 4],
    /// M[0][0] bit set at each candidate position in CANDIDATE_LABELS order.
    m00_cand: [bool; 6],
    /// Any entry has the bit set at each candidate position.
    any_cand: [bool; 6],
    /// Number of bt' in 0..=3 whose width N(bt') bounds every entry.
    consistent_bt: usize,
    /// Number of (bt', trl') pairs in 0..=3 x 0..=15 consistent with the matrix.
    consistent_pairs: usize,
    hint_aux_match: bool,
    hint_chall_match: bool,
    /// n_bt fits the 2-bit metadata field (n_bt <= 3).
    meta_bt_ok: bool,
    /// r' fits the 4-bit metadata field (r' <= 15).
    meta_trl_ok: bool,
    std_ok: bool,
    comp_len_ok: bool,
    comp_matrix_identical: bool,
    comp_bytes_identical: bool,
    comp_ok: bool,
    exp_len_ok: bool,
    exp_ok: bool,
}

fn analyze<L: FpBackend + LevelPrecomp>(entry: &KatEntry) -> VecStats {
    let pk = PublicKey::<L>::from_bytes(&entry.pk)
        .unwrap_or_else(|_| panic!("pk deserialization failed for count={}", entry.count));
    let sig_len = entry.sm.len() - entry.msg.len();
    let sig_bytes = &entry.sm[..sig_len];
    let sig = Signature::<L>::from_bytes(sig_bytes)
        .unwrap_or_else(|_| panic!("sig deserialization failed for count={}", entry.count));

    let bt = sig.backtracking();
    let trl = sig.two_resp_length();
    let e_rsp = L::E_RSP as usize;
    let hd = HD_EXTRA_TORSION as usize;
    let n_bits = e_rsp + hd - bt as usize;
    let pow_dim2 = e_rsp - trl as usize - bt as usize;
    let det_precision = pow_dim2 + trl as usize;

    let mat = sig.mat();
    let entries: [&[u64]; 4] = [
        mat[0][0].digits(),
        mat[0][1].digits(),
        mat[1][0].digits(),
        mat[1][1].digits(),
    ];
    let m00_even = entries[0][0] & 1 == 0;
    let bitlen_arr = [
        bitlen(entries[0]),
        bitlen(entries[1]),
        bitlen(entries[2]),
        bitlen(entries[3]),
    ];
    let m00_top_bit = bit(entries[0], n_bits - 1);
    let any_top_bit = entries.iter().any(|e| bit(e, n_bits - 1));
    // Candidate "expected bit positions" for the leading 1 (see CANDIDATE_LABELS).
    let g = pow_dim2 + hd;
    let cand_pos = [e_rsp + hd - 1, e_rsp - 1, e_rsp - 2, n_bits - 1, n_bits - 2, g - 1];
    let mut m00_cand = [false; 6];
    let mut any_cand = [false; 6];
    for (k, &p) in cand_pos.iter().enumerate() {
        m00_cand[k] = bit(entries[0], p);
        any_cand[k] = entries.iter().any(|e| bit(e, p));
    }

    // (b) det_hint, both from the compressor's wire byte and recomputed here.
    let comp = sig.compress();
    let cbytes = comp.to_bytes();
    let fp2_len = sig.e_aux_a().encode().len();
    let meta = cbytes[fp2_len];
    let det_hint_wire = (meta >> 2) & 0x03;
    // The metadata byte holds only 2 bits of n_bt and 4 bits of r'; record
    // whether this vector's values survive the packing.
    let meta_bt_ok = (meta & 0x03) == bt;
    let meta_trl_ok = ((meta >> 4) & 0x0F) == trl;
    let dropped = if m00_even { entries[2] } else { entries[3] };
    let det_hint_direct =
        (bit(dropped, det_precision) as u8) | ((bit(dropped, det_precision + 1) as u8) << 1);

    // (e) parameter inference from the matrix alone.
    let fits = |n: usize| entries.iter().all(|e| bitlen(e) <= n);
    let mut consistent_bt = 0;
    for btp in 0..=3usize {
        if fits(e_rsp + hd - btp) {
            consistent_bt += 1;
        }
    }
    let mut consistent_pairs = 0;
    for btp in 0..=3usize {
        for trlp in 0..=15usize {
            let pow = e_rsp as i32 - trlp as i32 - btp as i32;
            if pow <= 1 {
                continue;
            }
            if fits(e_rsp + hd - btp) {
                consistent_pairs += 1;
            }
        }
    }

    // (f) hint recomputation from the curves.
    let to_hint = |curve: &mut EcCurve<L>| -> u8 {
        let (_, h) = ec_curve_to_basis_2f_to_hint(
            curve,
            L::F_CHR,
            L::basis_e0_px_bytes(),
            L::basis_e0_qx_bytes(),
            L::p_cofactor_for_2f(),
            L::p_cofactor_for_2f_bitlength() as usize,
            L::torsion_even_power(),
        )
        .expect("basis to hint");
        h
    };
    let mut e_aux = EcCurve::<L>::from_a(sig.e_aux_a()).expect("E_aux from A_aux");
    let hint_aux_match = to_hint(&mut e_aux) == sig.hint_aux();
    let mut e_chall =
        compute_challenge_curve(sig.chall_coeff(), sig.backtracking(), pk.curve(), pk.hint_pk())
            .expect("challenge curve");
    let hint_chall_match = to_hint(&mut e_chall) == sig.hint_chall();

    // (g) format round trips.
    let std_ok = pk.verify(&entry.msg, &sig).is_ok();
    let expected_comp_len = fp2_len + 3 * ((e_rsp + 9) / 8) + L::LAMBDA as usize / 8 + 1;
    let comp_len_ok =
        cbytes.len() == CompressedSignature::<L>::WIRE_BYTES && cbytes.len() == expected_comp_len;
    let comp2 = CompressedSignature::<L>::from_bytes(cbytes.as_slice()).expect("compressed decode");
    let (comp_matrix_identical, comp_bytes_identical) = match comp2.decompress(&pk) {
        Ok(dec) => {
            let dm = dec.mat();
            let mat_same = (0..2).all(|i| (0..2).all(|j| dm[i][j].digits() == mat[i][j].digits()));
            (mat_same, dec.to_bytes().as_slice() == sig_bytes)
        }
        Err(_) => (false, false),
    };
    let comp_ok = pk.verify(&entry.msg, &comp2).is_ok();
    let (exp_len_ok, exp_ok) = match sig.expand(&pk) {
        Ok(exp) => {
            let eb = exp.to_bytes();
            let len_ok = eb.len() == ExpandedSignature::<L>::WIRE_BYTES;
            let ok = match ExpandedSignature::<L>::from_bytes(eb.as_slice()) {
                Ok(exp2) => pk.verify(&entry.msg, &exp2).is_ok(),
                Err(_) => false,
            };
            (len_ok, ok)
        }
        Err(_) => (false, false),
    };

    VecStats {
        count: entry.count,
        bt,
        trl,
        n_bits,
        m00_even,
        det_hint_wire,
        det_hint_direct,
        m00_top_bit,
        any_top_bit,
        bitlen: bitlen_arr,
        m00_cand,
        any_cand,
        consistent_bt,
        consistent_pairs,
        hint_aux_match,
        hint_chall_match,
        meta_bt_ok,
        meta_trl_ok,
        std_ok,
        comp_len_ok,
        comp_matrix_identical,
        comp_bytes_identical,
        comp_ok,
        exp_len_ok,
        exp_ok,
    }
}

fn report<L: FpBackend + LevelPrecomp>(name: &str, rsp: &str) -> String {
    let entries = kat_parser::parse_rsp(rsp);
    assert!(!entries.is_empty(), "no KAT entries parsed");
    let stats: Vec<VecStats> = entries.par_iter().map(analyze::<L>).collect();
    let n = stats.len();
    let e_rsp = L::E_RSP as usize;
    let hd = HD_EXTRA_TORSION as usize;

    let mut out = String::new();
    let w = &mut out;
    writeln!(w, "==== {name}: {n} KAT vectors, E_RSP={e_rsp}, HD_EXTRA_TORSION={hd}, LAMBDA={}, F_CHR={} ====", L::LAMBDA, L::F_CHR).unwrap();

    // (a)
    let even = stats.iter().filter(|s| s.m00_even).count();
    writeln!(w, "(a) M[0][0] parity: even={even} odd={} of {n}", n - even).unwrap();

    // (b)
    let mut hist = [0usize; 4];
    for s in &stats {
        assert_eq!(s.det_hint_wire, s.det_hint_direct, "count={}: wire det_hint != direct", s.count);
        hist[s.det_hint_wire as usize] += 1;
    }
    let exp = n as f64 / 4.0;
    let chi: f64 = hist.iter().map(|&o| (o as f64 - exp).powi(2) / exp).sum();
    writeln!(w, "(b) det_hint histogram {{0,1,2,3}} = {:?}; chi-square vs uniform (3 dof) = {chi:.3}", hist).unwrap();
    writeln!(w, "    (compressor wire byte and direct recomputation from the dropped entry agree on all {n})").unwrap();
    let max_bt = stats.iter().map(|s| s.bt).max().unwrap() as usize;
    let max_trl = stats.iter().map(|s| s.trl).max().unwrap() as usize;
    writeln!(w, "    contingency det_hint (rows 0..3) x n_bt (cols 0..={max_bt}):").unwrap();
    for h in 0..4 {
        let row: Vec<usize> = (0..=max_bt)
            .map(|b| stats.iter().filter(|s| s.det_hint_wire as usize == h && s.bt as usize == b).count())
            .collect();
        writeln!(w, "      hint={h}: {row:?}").unwrap();
    }
    writeln!(w, "    contingency det_hint (rows 0..3) x r' (cols 0..={max_trl}):").unwrap();
    for h in 0..4 {
        let row: Vec<usize> = (0..=max_trl)
            .map(|t| stats.iter().filter(|s| s.det_hint_wire as usize == h && s.trl as usize == t).count())
            .collect();
        writeln!(w, "      hint={h}: {row:?}").unwrap();
    }
    writeln!(w, "    contingency det_hint (rows 0..3) x M[0][0] parity (cols [even, odd]):").unwrap();
    for h in 0..4 {
        let ev = stats.iter().filter(|s| s.det_hint_wire as usize == h && s.m00_even).count();
        let od = stats.iter().filter(|s| s.det_hint_wire as usize == h && !s.m00_even).count();
        writeln!(w, "      hint={h}: [{ev}, {od}]").unwrap();
    }

    // (c)
    let bt_hist: Vec<usize> = (0..=max_bt).map(|b| stats.iter().filter(|s| s.bt as usize == b).count()).collect();
    let trl_hist: Vec<usize> = (0..=max_trl).map(|t| stats.iter().filter(|s| s.trl as usize == t).count()).collect();
    writeln!(w, "(c) n_bt histogram (0..={max_bt}) = {bt_hist:?}").unwrap();
    writeln!(w, "    r'   histogram (0..={max_trl}) = {trl_hist:?}").unwrap();
    writeln!(w, "    joint (n_bt rows 0..={max_bt}, r' cols 0..={max_trl}):").unwrap();
    for b in 0..=max_bt {
        let row: Vec<usize> = (0..=max_trl)
            .map(|t| stats.iter().filter(|s| s.bt as usize == b && s.trl as usize == t).count())
            .collect();
        writeln!(w, "      n_bt={b}: {row:?}").unwrap();
    }

    // (d)
    writeln!(w, "(d) N = E_RSP + HD_EXTRA_TORSION - n_bt (verifier bound: every entry < 2^N); matrix entry wire width = {} bits", 8 * ((e_rsp + 9) / 8)).unwrap();
    for b in 0..=max_bt {
        let cnt = stats.iter().filter(|s| s.bt as usize == b).count();
        if cnt > 0 {
            writeln!(w, "    n_bt={b}: N={} ({cnt} vectors)", e_rsp + hd - b).unwrap();
        }
    }
    let m00_top = stats.iter().filter(|s| s.m00_top_bit).count();
    let any_top = stats.iter().filter(|s| s.any_top_bit).count();
    writeln!(w, "    bit N-1 of M[0][0] set: {m00_top} of {n}").unwrap();
    writeln!(w, "    bit N-1 of at least one entry set: {any_top} of {n}").unwrap();
    for (k, label) in CANDIDATE_LABELS.iter().enumerate() {
        let m = stats.iter().filter(|s| s.m00_cand[k]).count();
        let a = stats.iter().filter(|s| s.any_cand[k]).count();
        writeln!(w, "    candidate position {label}: M[0][0] bit set in {m} of {n}; any entry bit set in {a} of {n}").unwrap();
    }
    let mut gap_m00 = [0usize; 5];
    let mut gap_max = [0usize; 5];
    for s in &stats {
        let g = s.n_bits.saturating_sub(s.bitlen[0]).min(4);
        gap_m00[g] += 1;
        let maxlen = *s.bitlen.iter().max().unwrap();
        let g2 = s.n_bits.saturating_sub(maxlen).min(4);
        gap_max[g2] += 1;
    }
    writeln!(w, "    N - bitlen(M[0][0]) histogram [0,1,2,3,>=4] = {gap_m00:?}").unwrap();
    writeln!(w, "    N - max bitlen over 4 entries histogram [0,1,2,3,>=4] = {gap_max:?}").unwrap();
    let over = stats.iter().filter(|s| s.bitlen.iter().any(|&l| l > s.n_bits)).count();
    writeln!(w, "    entries exceeding 2^N (would be rejected): {over}").unwrap();

    // (e)
    let mut pairs_hist = std::collections::BTreeMap::new();
    let mut bt_cons_hist = [0usize; 5];
    for s in &stats {
        *pairs_hist.entry(s.consistent_pairs).or_insert(0usize) += 1;
        bt_cons_hist[s.consistent_bt] += 1;
    }
    let unique = stats.iter().filter(|s| s.consistent_pairs == 1).count();
    writeln!(w, "(e) consistent (n_bt, r') pairs per vector (n_bt 0..=3, r' 0..=15, consistent iff all entries < 2^(E_RSP+2-n_bt) and E_RSP-r'-n_bt > 1):").unwrap();
    writeln!(w, "    histogram {{count: vectors}} = {pairs_hist:?}").unwrap();
    writeln!(w, "    vectors with exactly one consistent pair: {unique} of {n}").unwrap();
    writeln!(w, "    consistent n_bt values per vector, histogram [0,1,2,3,4] = {bt_cons_hist:?} (N does not depend on r', so pairs = 16 x this)").unwrap();
    let bt_unique_correct = stats.iter().filter(|s| s.consistent_bt == 1 && s.bt == 0).count();
    writeln!(w, "    vectors whose n_bt alone is uniquely determined: {} (all of them have n_bt=0: {bt_unique_correct})", bt_cons_hist[1]).unwrap();

    // (f)
    let ha = stats.iter().filter(|s| s.hint_aux_match).count();
    let hc = stats.iter().filter(|s| s.hint_chall_match).count();
    writeln!(w, "(f) recomputed hint_aux matches transmitted: {ha} of {n}; recomputed hint_chall matches: {hc} of {n}").unwrap();

    // (g)
    let c = |f: fn(&VecStats) -> bool| stats.iter().filter(|s| f(s)).count();
    let which = |f: fn(&VecStats) -> bool| -> Vec<usize> {
        stats.iter().filter(|s| !f(s)).map(|s| s.count).collect()
    };
    let fp2_len = Signature::<L>::default().e_aux_a().encode().len();
    writeln!(w, "(g) standard verifies: {} of {n}", c(|s| s.std_ok)).unwrap();
    writeln!(w, "    n_bt fits the 2-bit metadata field: {} of {n} (failing counts: {:?}); r' fits the 4-bit field: {} of {n} (failing counts: {:?})", c(|s| s.meta_bt_ok), which(|s| s.meta_bt_ok), c(|s| s.meta_trl_ok), which(|s| s.meta_trl_ok)).unwrap();
    writeln!(w, "    compressed length == WIRE_BYTES == Fp2 + 3*entry + chall + 1: {} of {n}  (WIRE_BYTES={}, Fp2={}, entry={}, chall={}, meta=1)", c(|s| s.comp_len_ok), CompressedSignature::<L>::WIRE_BYTES, fp2_len, (e_rsp + 9) / 8, L::LAMBDA as usize / 8).unwrap();
    writeln!(w, "    decompressed matrix identical to original: {} of {n} (failing counts: {:?})", c(|s| s.comp_matrix_identical), which(|s| s.comp_matrix_identical)).unwrap();
    writeln!(w, "    decompressed signature byte-identical to original standard bytes: {} of {n} (failing counts: {:?})", c(|s| s.comp_bytes_identical), which(|s| s.comp_bytes_identical)).unwrap();
    writeln!(w, "    compressed verifies: {} of {n} (failing counts: {:?})", c(|s| s.comp_ok), which(|s| s.comp_ok)).unwrap();
    writeln!(w, "    expanded length == WIRE_BYTES ({}): {} of {n}; expanded verifies after byte round trip: {} of {n}", ExpandedSignature::<L>::WIRE_BYTES, c(|s| s.exp_len_ok), c(|s| s.exp_ok)).unwrap();

    // Per-vector listing for the record.
    writeln!(w, "    per-vector: count bt trl N m00_even det_hint bitlen[m00,m01,m10,m11] consistent_bt").unwrap();
    for s in &stats {
        writeln!(w, "      {:3} {} {} {} {} {} {:?} {}", s.count, s.bt, s.trl, s.n_bits, s.m00_even as u8, s.det_hint_wire, s.bitlen, s.consistent_bt).unwrap();
    }

    for s in &stats {
        assert!(s.std_ok, "count={} standard verification failed", s.count);
        assert!(s.hint_aux_match && s.hint_chall_match, "count={} hint mismatch", s.count);
        assert!(s.exp_ok && s.exp_len_ok, "count={} expanded round trip failed", s.count);
        // The compressed round trip can only be expected to succeed when the
        // metadata byte can hold n_bt and r'.
        if s.meta_bt_ok && s.meta_trl_ok {
            assert!(s.comp_ok && s.comp_bytes_identical && s.comp_len_ok, "count={} compressed round trip failed", s.count);
        } else {
            assert!(!s.comp_ok, "count={} compressed verification unexpectedly succeeded with truncated metadata", s.count);
        }
    }
    out
}

#[test]
fn paper_stats_level1() {
    print!("{}", report::<Level1>("Level 1", KAT_L1));
}

#[test]
fn paper_stats_level3() {
    print!("{}", report::<Level3>("Level 3", KAT_L3));
}

#[test]
fn paper_stats_level5() {
    print!("{}", report::<Level5>("Level 5", KAT_L5));
}
