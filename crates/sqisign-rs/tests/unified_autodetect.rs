//! Unified, dimension-agnostic autodetection - now that the dim-4 SQIsignHD
//! verifier lives inside `sqisign-verify` as the `hd` module, HD is just
//! another arm of the crate's existing `AnySignature` (no wrapper types):
//!
//! 1. dim-2 (standard / expanded / compressed) autodetect and verify
//!    (regression: the dim-2 path is unchanged).
//! 2. SQIsignHD signatures autodetect at 108 bytes and verify; the 64-byte HD
//!    public key autodetects through `PublicKey::from_bytes`.
//! 3. Each wire length routes to the correct arm; bad lengths reject.
//! 4. Malformed HD input and cross-dimension mismatches reject cleanly.

use crypto_bigint::U256;
use serde_json::Value;
use sqisign_rs::keygen::keypair;
use sqisign_rs::params::Level1;
use sqisign_rs::sign::sign;
use sqisign_rs::{AnySignature, CompactPublicKey, Fp2, PublicKey, Verifier};

const HD_VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/test_vectors_l1.json"
);
/// The message the Level-1 SQIsignHD oracle vectors were generated against.
const HD_MSG: [u8; 32] = [0u8; 32];

fn is_compact(sig: &AnySignature<Level1>) -> bool {
    matches!(sig, AnySignature::Compact(_))
}

// dim-2: autodetect + verify is the unchanged verifier

#[test]
fn dim2_all_formats_autodetect_and_verify() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<Level1>(&mut rng);
    let msg = b"unified autodetect";
    let sig = sign::<Level1>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    // 65-byte key autodetects as dim-2.
    let pk = PublicKey::<Level1>::from_bytes(&pk.to_bytes()).expect("pk parse");

    // Standard (148 B).
    let std_bytes = sig.to_bytes();
    let any = AnySignature::<Level1>::from_bytes(&std_bytes).expect("standard parse");
    assert!(!is_compact(&any), "148 B must route to dim-2");
    pk.verify(msg, &any).expect("standard must verify");

    // Expanded (212 B).
    let exp = sig.expand(&pk).expect("expand");
    let any = AnySignature::<Level1>::from_bytes(&exp.to_bytes()).expect("expanded parse");
    assert!(!is_compact(&any), "212 B must route to dim-2");
    pk.verify(msg, &any).expect("expanded must verify");

    // Compressed (129 B).
    let cmp = sig.compress();
    let any = AnySignature::<Level1>::from_bytes(&cmp.to_bytes()).expect("compressed parse");
    assert!(!is_compact(&any), "129 B must route to dim-2");
    pk.verify(msg, &any).expect("compressed must verify");

    // Wrong message must reject through the unified path too.
    let any = AnySignature::<Level1>::from_bytes(&std_bytes).unwrap();
    assert!(
        pk.verify(b"different message", &any).is_err(),
        "wrong message must reject"
    );
}

// HD vector plumbing (encode the oracle vectors to wire bytes)

fn le32(hexstr: &str) -> [u8; 32] {
    let s = hexstr.trim_start_matches("0x");
    let s = if s.len() % 2 == 1 {
        format!("0{s}")
    } else {
        s.to_string()
    };
    let be = hex::decode(&s).expect("valid hex");
    let mut le = [0u8; 32];
    for (i, b) in be.iter().rev().enumerate() {
        le[i] = *b;
    }
    le
}
fn parse_fp2(pair: &Value) -> Fp2<Level1> {
    let mut buf = [0u8; 64];
    buf[..32].copy_from_slice(&le32(pair[0].as_str().unwrap()));
    buf[32..].copy_from_slice(&le32(pair[1].as_str().unwrap()));
    Fp2::<Level1>::decode(&buf).expect("coordinate in [0,p)")
}
fn dec_u256(s: &str) -> U256 {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    U256::from_words(limbs)
}
fn dec_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}

/// Encode all five oracle vectors into `(sig_bytes, pk_bytes)` HD wire pairs.
fn hd_wire_vectors() -> Vec<(Vec<u8>, Vec<u8>)> {
    let doc: Value =
        serde_json::from_str(&std::fs::read_to_string(HD_VECTORS).expect("read HD vectors"))
            .expect("valid json");
    doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            let sig = &v["signature"];
            let s = sqisign_rs::hd::encode_signature(
                &parse_fp2(&sig["A_com"]),
                dec_i128(&sig["a"]),
                dec_i128(&sig["b"]),
                dec_i128(&sig["c_or_d"]),
                &dec_u256(sig["q"].as_str().unwrap()),
                sig["hint_com_P"].as_u64().unwrap() as u32,
                sig["hint_com_Q"].as_u64().unwrap() as u32,
            )
            .expect("encode HD signature");
            let p = sqisign_rs::hd::encode_public_key(
                &parse_fp2(&v["public_key"]["A_pk"]),
                v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32,
                v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32,
            )
            .expect("encode HD public key");
            (s.to_vec(), p.to_vec())
        })
        .collect()
}

#[test]
fn hd_autodetects_and_verifies() {
    let wires = hd_wire_vectors();
    assert_eq!(wires.len(), 5);
    for (i, (s, p)) in wires.iter().enumerate() {
        assert_eq!(s.len(), 108, "HD signature must be 108 bytes");
        assert_eq!(p.len(), 64, "HD public key must be 64 bytes");

        let any = AnySignature::<Level1>::from_bytes(s).expect("HD signature parse");
        assert!(
            is_compact(&any),
            "vec {i}: 108 B must route to the compact arm"
        );

        // The 64-byte compact public key is verified by a `CompactPublicKey`.
        let cpk = CompactPublicKey::<Level1>::from_bytes(p).expect("compact pk parse");
        cpk.verify(&HD_MSG, &any)
            .unwrap_or_else(|_| panic!("vec {i}: valid compact signature must verify"));

        // Tampering the message must break the recomputed Fiat-Shamir challenge.
        let mut bad = HD_MSG;
        bad[0] ^= 1;
        assert!(
            cpk.verify(&bad, &any).is_err(),
            "vec {i}: tampered message must reject"
        );

        // Cross-scheme: a dim-2 public key must NOT verify a compact signature,
        // even though `PublicKey::from_bytes` still accepts the 64-byte key.
        let dim2_pk = PublicKey::<Level1>::from_bytes(p).expect("dim-2 pk parse (legacy 64-byte)");
        assert!(
            dim2_pk.verify(&HD_MSG, &any).is_err(),
            "vec {i}: a dim-2 key must not verify a compact signature"
        );
    }
}

// routing / negative cases

#[test]
fn bad_lengths_reject_cleanly() {
    // Signature lengths that match no (dimension, format) at Level 1: rejected,
    // no panic. (108/129/148/212 are the valid signature lengths.)
    for bad_len in [0usize, 1, 64, 107, 109, 130, 149, 200, 213, 300] {
        assert!(
            AnySignature::<Level1>::from_bytes(&vec![0u8; bad_len]).is_err(),
            "len {bad_len} must be rejected as a signature"
        );
    }
    // Public-key lengths that are neither 64 (HD) nor 65 (dim-2): rejected.
    for bad_len in [0usize, 1, 63, 66, 100, 108] {
        assert!(
            PublicKey::<Level1>::from_bytes(&vec![0u8; bad_len]).is_err(),
            "len {bad_len} must be rejected as a public key"
        );
    }
}

#[test]
fn malformed_hd_rejects_cleanly() {
    // 108 bytes whose A_com (first 64) is out of range: Fp2::decode rejects
    // after the packed hint bits are masked off.
    let mut bad = [0u8; 108];
    for b in bad[..64].iter_mut() {
        *b = 0xFF;
    }
    assert!(
        AnySignature::<Level1>::from_bytes(&bad).is_err(),
        "out-of-range A_com must reject at parse"
    );

    // 64-byte pk with out-of-range A_pk: rejects.
    let bad_pk = [0xFFu8; 64];
    assert!(
        PublicKey::<Level1>::from_bytes(&bad_pk).is_err(),
        "out-of-range A_pk must reject at parse"
    );
}

#[test]
fn dimension_mismatch_rejects() {
    // A real dim-2 key + a real HD signature: cross-dimension must reject.
    let mut rng = rand::thread_rng();
    let (pk, _sk) = keypair::<Level1>(&mut rng);
    let pk_dim2 = PublicKey::<Level1>::from_bytes(&pk.to_bytes()).unwrap();

    let (hd_sig_bytes, hd_pk_bytes) = hd_wire_vectors().remove(0);
    let hd_sig = AnySignature::<Level1>::from_bytes(&hd_sig_bytes).unwrap();
    let hd_pk = PublicKey::<Level1>::from_bytes(&hd_pk_bytes).unwrap();

    // dim-2 key + HD signature: the HD verify runs against the wrong curve.
    assert!(
        pk_dim2.verify(&HD_MSG, &hd_sig).is_err(),
        "dim-2 key must reject an HD signature"
    );

    // HD key + dim-2 signature.
    let sk_msg = b"mismatch test";
    let (pk2, sk2) = keypair::<Level1>(&mut rng);
    let dim2_sig = AnySignature::<Level1>::from_bytes(
        &sign::<Level1>(&sk2, &pk2, sk_msg, &mut rng)
            .unwrap()
            .to_bytes(),
    )
    .unwrap();
    assert!(
        hd_pk.verify(sk_msg, &dim2_sig).is_err(),
        "HD key must reject a dim-2 signature"
    );
}
