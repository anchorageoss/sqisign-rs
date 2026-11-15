//! Compressed-format KAT validation at all three security levels.
//!
//! For every official KAT vector (100 per level) the standard signature is
//! compressed, and the test asserts that: compression succeeds, the encoding
//! has the level's compressed length, decoding and re-encoding it is
//! byte-identical, decompressing it with the public key gives back the
//! original standard signature byte for byte, and the compressed signature
//! verifies on the entry's message and is rejected on a modified message.
//!
//! The metadata byte `[r:4 | h:1 | n_bt:3]` is also read back from the wire
//! and checked against the standard signature: `n_bt` and `r` must match its
//! fields and `h` must be bit `N-1` of the dropped matrix entry, where
//! `N = E_RSP + 2 - n_bt`.
//!
//! Per level the test prints the round-trip count, the histogram of `h`, the
//! histograms of `n_bt` and `r`, and the contingency counts of `h` against
//! `n_bt`, against `r` and against the parity of `M[0][0]`.
//!
//! Run with: cargo test -p sqisign-kat --release --features kat-compat --test kat_compressed -- --nocapture

#![cfg(feature = "kat-compat")]

use rayon::prelude::*;
use sqisign_kat::kat_parser::{self, KatEntry};
use sqisign_verify::{
    CompressedSignature, Error, FpBackend, Level1, Level3, Level5, LevelPrecomp, PublicKey,
    SecurityLevel, Signature, Verifier,
};

const KAT_LVL1: &str = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
const KAT_LVL3: &str = include_str!("../../../reference/KAT/PQCsignKAT_529_SQIsign_lvl3.rsp");
const KAT_LVL5: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

/// Expected wire sizes per level: (standard, compressed).
const L1_SIZES: (usize, usize) = (148, 129);
const L3_SIZES: (usize, usize) = (224, 196);
const L5_SIZES: (usize, usize) = (292, 257);

/// What one round-tripped vector contributes to the statistics.
#[derive(Clone, Copy)]
struct Observation {
    n_bt: u8,
    r: u8,
    h: u8,
    m00_odd: bool,
}

/// Round-trip one KAT vector through the compressed format, asserting every
/// step, and return the metadata the statistics need.
fn roundtrip_entry<L: FpBackend + LevelPrecomp>(
    entry: &KatEntry,
    sizes: (usize, usize),
    label: &str,
) -> Observation {
    let (standard_len, compressed_len) = sizes;
    let idx = entry.count;

    let pk = PublicKey::<L>::from_bytes(&entry.pk)
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: pk decode failed: {e:?}"));
    let sig_len = entry.sm.len() - entry.msg.len();
    assert_eq!(
        sig_len, standard_len,
        "{label} KAT {idx}: standard signature length"
    );
    let sig_bytes = &entry.sm[..sig_len];
    let sig = Signature::<L>::from_bytes(sig_bytes)
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: sig decode failed: {e:?}"));

    // Compression must succeed for every honest signature.
    let compressed = sig
        .compress()
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: compress failed: {e:?}"));

    // The encoding has the level's compressed length.
    let wire = compressed.to_bytes();
    assert_eq!(
        wire.len(),
        compressed_len,
        "{label} KAT {idx}: compressed wire length"
    );
    assert_eq!(
        CompressedSignature::<L>::WIRE_BYTES,
        compressed_len,
        "{label}: CompressedSignature::WIRE_BYTES"
    );

    // Decoding and re-encoding is byte-identical.
    let decoded = CompressedSignature::<L>::from_bytes(&wire)
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: compressed decode failed: {e:?}"));
    assert_eq!(
        decoded.to_bytes().as_slice(),
        wire.as_slice(),
        "{label} KAT {idx}: compressed re-encode differs"
    );

    // Decompression gives back the original standard signature byte for byte.
    let recovered = decoded
        .decompress(&pk)
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: decompress failed: {e:?}"));
    assert_eq!(
        recovered.to_bytes().as_slice(),
        sig_bytes,
        "{label} KAT {idx}: decompressed signature differs from the original (n_bt={}, r={})",
        sig.backtracking(),
        sig.two_resp_length(),
    );

    // The compressed signature verifies on the entry's message and not on a
    // modified one.
    pk.verify(&entry.msg, &decoded)
        .unwrap_or_else(|e| panic!("{label} KAT {idx}: compressed verify failed: {e:?}"));
    let mut wrong_msg = entry.msg.clone();
    wrong_msg[0] ^= 0x01;
    assert!(
        pk.verify(&wrong_msg, &decoded).is_err(),
        "{label} KAT {idx}: compressed verify accepted a modified message"
    );

    // Metadata byte [r:4 | h:1 | n_bt:3] directly after the Fp2 coefficient.
    // The public key is that same Fp2 plus one hint byte, which gives the
    // Fp2 encoded length without reaching into the type-level sizes.
    let fp2_len = entry.pk.len() - 1;
    let meta = wire[fp2_len];
    let n_bt = meta & 0x07;
    let h = (meta >> 3) & 0x01;
    let r = meta >> 4;
    assert_eq!(
        n_bt,
        sig.backtracking(),
        "{label} KAT {idx}: n_bt field of the metadata byte"
    );
    assert_eq!(
        r,
        sig.two_resp_length(),
        "{label} KAT {idx}: r field of the metadata byte"
    );

    // h is bit N-1 of the dropped entry (M11 when M00 is odd, else M10).
    let m00_odd = sig.mat()[0][0].digits()[0] & 1 != 0;
    let dropped = if m00_odd {
        &sig.mat()[1][1]
    } else {
        &sig.mat()[1][0]
    };
    let top_bit = <L as SecurityLevel>::E_RSP as usize + 2 - n_bt as usize - 1;
    let expected_h = ((dropped.digits()[top_bit / 64] >> (top_bit % 64)) & 1) as u8;
    assert_eq!(
        h, expected_h,
        "{label} KAT {idx}: h is not bit N-1 of the dropped entry"
    );

    Observation {
        n_bt,
        r,
        h,
        m00_odd,
    }
}

/// Round-trip every vector of one level and print the statistics.
fn kat_compressed_level<L: FpBackend + LevelPrecomp>(
    rsp: &str,
    sizes: (usize, usize),
    label: &str,
) -> Vec<Observation> {
    let entries = kat_parser::parse_rsp(rsp);
    assert_eq!(entries.len(), 100, "{label}: expected 100 KAT entries");
    for (i, entry) in entries.iter().enumerate() {
        assert_eq!(entry.count, i, "{label}: KAT entries out of order");
    }

    let observations: Vec<Observation> = entries
        .par_iter()
        .map(|entry| roundtrip_entry::<L>(entry, sizes, label))
        .collect();
    assert_eq!(observations.len(), 100, "{label}: round-trip count");

    let mut h_hist = [0u32; 2];
    let mut n_bt_hist = [0u32; 8];
    let mut r_hist = [0u32; 16];
    let mut h_by_n_bt = [[0u32; 2]; 8];
    let mut h_by_r = [[0u32; 2]; 16];
    let mut h_by_parity = [[0u32; 2]; 2];
    for obs in &observations {
        let h = obs.h as usize;
        h_hist[h] += 1;
        n_bt_hist[obs.n_bt as usize] += 1;
        r_hist[obs.r as usize] += 1;
        h_by_n_bt[obs.n_bt as usize][h] += 1;
        h_by_r[obs.r as usize][h] += 1;
        h_by_parity[obs.m00_odd as usize][h] += 1;
    }

    println!(
        "[{label}] compressed round-trips: {}/100 (standard {} B -> compressed {} B)",
        observations.len(),
        sizes.0,
        sizes.1
    );
    println!(
        "[{label}] h histogram: h=0: {}, h=1: {}",
        h_hist[0], h_hist[1]
    );
    let n_bt_line: Vec<String> = n_bt_hist
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(n_bt, c)| format!("n_bt={n_bt}: {c}"))
        .collect();
    println!("[{label}] n_bt histogram: {}", n_bt_line.join(", "));
    let r_line: Vec<String> = r_hist
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(r, c)| format!("r={r}: {c}"))
        .collect();
    println!("[{label}] r histogram: {}", r_line.join(", "));
    let h_n_bt_line: Vec<String> = h_by_n_bt
        .iter()
        .enumerate()
        .filter(|(_, row)| row[0] + row[1] > 0)
        .map(|(n_bt, row)| format!("n_bt={n_bt}: h=0: {}, h=1: {}", row[0], row[1]))
        .collect();
    println!("[{label}] h vs n_bt: {}", h_n_bt_line.join("; "));
    let h_r_line: Vec<String> = h_by_r
        .iter()
        .enumerate()
        .filter(|(_, row)| row[0] + row[1] > 0)
        .map(|(r, row)| format!("r={r}: h=0: {}, h=1: {}", row[0], row[1]))
        .collect();
    println!("[{label}] h vs r: {}", h_r_line.join("; "));
    println!(
        "[{label}] h vs M[0][0] parity: even: h=0: {}, h=1: {}; odd: h=0: {}, h=1: {}",
        h_by_parity[0][0], h_by_parity[0][1], h_by_parity[1][0], h_by_parity[1][1]
    );

    observations
}

#[test]
fn kat_compressed_level1() {
    kat_compressed_level::<Level1>(KAT_LVL1, L1_SIZES, "L1");
}

#[test]
fn kat_compressed_level3() {
    let observations = kat_compressed_level::<Level3>(KAT_LVL3, L3_SIZES, "L3");

    // Level III entry 97 has n_bt = 4, which the old 2-bit backtracking field
    // silently truncated to 0. It must now be representable and round-trip.
    let obs = observations[97];
    assert_eq!(obs.n_bt, 4, "L3 KAT 97: expected n_bt = 4");
    println!(
        "[L3] entry 97: n_bt={} r={} h={} M[0][0] {}: round-trip OK",
        obs.n_bt,
        obs.r,
        obs.h,
        if obs.m00_odd { "odd" } else { "even" }
    );
}

#[test]
fn kat_compressed_level5() {
    kat_compressed_level::<Level5>(KAT_LVL5, L5_SIZES, "L5");
}

/// `compress()` rejects metadata that does not fit the byte instead of
/// truncating it: n_bt needs 3 bits (0..8) and r needs 4 bits (0..16).
#[test]
fn compress_rejects_out_of_range_metadata_on_kat_signature() {
    let entries = kat_parser::parse_rsp(KAT_LVL1);
    let entry = &entries[0];
    let sig_bytes = &entry.sm[..entry.sm.len() - entry.msg.len()];
    let sig = Signature::<Level1>::from_bytes(sig_bytes).expect("KAT 0 decodes");
    assert!(sig.compress().is_ok(), "KAT 0 compresses as is");

    let mut too_much_backtracking = sig.clone();
    too_much_backtracking.set_backtracking(8);
    assert_eq!(
        too_much_backtracking.compress().err(),
        Some(Error::MalformedInput),
        "backtracking = 8 must be rejected"
    );

    let mut too_long_response = sig.clone();
    too_long_response.set_two_resp_length(16);
    assert_eq!(
        too_long_response.compress().err(),
        Some(Error::MalformedInput),
        "two_resp_length = 16 must be rejected"
    );

    // The field maxima still fit and survive the wire.
    let mut at_the_limits = sig.clone();
    at_the_limits.set_backtracking(7);
    at_the_limits.set_two_resp_length(15);
    let wire = at_the_limits
        .compress()
        .expect("n_bt = 7, r = 15 fit the metadata byte")
        .to_bytes();
    let meta = wire[entry.pk.len() - 1];
    assert_eq!(meta & 0x07, 7);
    assert_eq!(meta >> 4, 15);
}
