//! The compressed signature format (`sqisign_verify::compressed`, not in
//! the specification) measured on the NIST known-answer tests of the
//! round-3 reference (`KAT/*.rsp` of the-sqisign at tag `nist-v3`, 100
//! entries per level). For each entry the public key and the signature
//! (the first `SIGNATURE_BYTES` of `sm`; the message follows) are decoded
//! and the signature must verify; its compressed form must encode to the
//! level's length (`compressed_bytes`: 174 / 267 / 352), decode back to
//! the same bytes, verify on the entry's message and fail on a modified
//! one, and decompress to the signer's signature, byte for byte.
//!
//! Per level one line reports what the format's description predicts and
//! the paper tabulates: the sizes, how many entries take the pivot path
//! (`m00` even, so the third transmitted entry is `m11`), the histogram of
//! bit `RESPONSE_BITS` of the three transmitted entries, and that every
//! decompressed signature equals the signer's.
//!
//! All 100 entries run in release; a debug build runs the first two unless
//! `SQISIGN_KAT_ENTRIES` says otherwise. `SQISIGN_FORCE_PORTABLE=1` runs
//! the portable field kernels of the x86-64 backends.

use rayon::prelude::*;
use sqisign_verify::compressed::compressed_bytes;
use sqisign_verify::sqisign::Scalar;
use sqisign_verify::{CompressedSignature, Level, Level1, Level3, Level5, PublicKey, Signature};

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

fn entries_to_run() -> usize {
    std::env::var("SQISIGN_KAT_ENTRIES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(if cfg!(debug_assertions) {
            2
        } else {
            usize::MAX
        })
}

/// What one entry contributes to the counts.
struct Observed {
    /// `m00` is even: the compressed signature carries `m11`, not `m10`.
    pivot: bool,
    /// Bit `RESPONSE_BITS` of `m00`, `m01` and `m_var`, in bits 0 to 2.
    bits: u8,
    /// The decompressed signature equals the signer's byte for byte.
    exact: bool,
}

fn run<L: Level>(rsp: &str, name: &str, compressed_len: usize) {
    let backend = sqisign_kat::arithmetic_from_env();
    let params = &L::PARAMS;
    let siglen = params.signature_bytes;
    assert_eq!(
        compressed_bytes(params),
        compressed_len,
        "{name}: compressed length"
    );
    let entries = parse_rsp(rsp);
    assert_eq!(entries.len(), 100, "{name}: KAT entries");
    let n = entries_to_run().min(entries.len());
    let observed: Vec<Observed> = entries[..n]
        .par_iter()
        .map(|e| {
            let what = format!("{name} count {}", e.count);
            let pk = PublicKey::<L>::from_bytes(&e.pk).expect("pk decode");
            assert_eq!(e.sm.len(), siglen + e.msg.len(), "{what}: smlen");
            let sig = Signature::<L>::from_bytes(&e.sm[..siglen]).expect("sig decode");
            assert!(
                pk.verify(&e.msg, &sig).is_ok(),
                "{what}: verify ({backend})"
            );
            // compression: the level's length, and the bytes round-trip
            let c = sig.compress().expect("compress");
            let cb = c.to_bytes();
            assert_eq!(cb.len(), compressed_len, "{what}: compressed length");
            assert_eq!(
                cb.len(),
                compressed_bytes(params),
                "{what}: compressed_bytes"
            );
            let c2 = CompressedSignature::<L>::from_bytes(&cb[..]).expect("compressed decode");
            assert_eq!(
                hex::encode(&cb[..]),
                hex::encode(&c2.to_bytes()[..]),
                "{what}: compressed round trip"
            );
            // verification of the compressed form; a modified message fails
            assert!(
                pk.verify_compressed(&e.msg, &c2).is_ok(),
                "{what}: verify compressed ({backend})"
            );
            let mut bad = e.msg.clone();
            match bad.first_mut() {
                Some(b) => *b ^= 1,
                None => bad.push(0),
            }
            assert!(
                pk.verify_compressed(&bad, &c2).is_err(),
                "{what}: a modified message verified compressed"
            );
            // decompression: the signer's signature, byte for byte, and
            // the standard verifier accepts it
            let dec = c2.decompress(&pk).expect("decompress");
            assert!(
                pk.verify(&e.msg, &dec).is_ok(),
                "{what}: verify decompressed ({backend})"
            );
            let exact = dec.to_bytes()[..] == e.sm[..siglen];
            assert!(
                exact,
                "{what}: the decompressed signature is not the signer's"
            );
            let orig = &sig.inner().mat_bchall_can_to_bchall;
            let pivot = orig[0][0][0] & 1 == 0;
            assert_eq!(c2.inner().m00[0] & 1 == 0, pivot, "{what}: pivot");
            let r = params.response_bits;
            let top = |x: &Scalar| ((x[(r / 64) as usize] >> (r % 64)) & 1) as u8;
            let c = c2.inner();
            let bits = top(&c.m00) | top(&c.m01) << 1 | top(&c.m_var) << 2;
            Observed { pivot, bits, exact }
        })
        .collect();
    let verified = observed.len();
    let pivot = observed.iter().filter(|o| o.pivot).count();
    let exact = observed.iter().filter(|o| o.exact).count();
    let mut hist = [0usize; 8];
    for o in &observed {
        hist[usize::from(o.bits)] += 1;
    }
    let saved = siglen - compressed_len;
    let pct = 100.0 * saved as f64 / siglen as f64;
    println!(
        "{name}: standard {siglen} B, compressed {compressed_len} B, saved {saved} B ({pct:.1}%); \
         verified compressed {verified}/{n}; pivot path (m00 even) {pivot}; \
         decompressed signature equal to the signer's {exact}/{n}; histogram of bit RESPONSE_BITS of (m00, m01, m_var) {hist:?}"
    );
}

#[test]
fn kat_compressed_level1() {
    run::<Level1>(
        include_str!("../kat/PQCsignKAT_270_SQIsign_p324_3.rsp"),
        "p324_3",
        174,
    );
}

#[test]
fn kat_compressed_level3() {
    run::<Level3>(
        include_str!("../kat/PQCsignKAT_417_SQIsign_p500_27.rsp"),
        "p500_27",
        267,
    );
}

#[test]
fn kat_compressed_level5() {
    run::<Level5>(
        include_str!("../kat/PQCsignKAT_549_SQIsign_p664_17.rsp"),
        "p664_17",
        352,
    );
}
