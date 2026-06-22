//! The compact (108-byte-signature) public API: round-trip, tampering, the
//! unified autodetect path, and the cross-scheme verification rules.
//!
//! Compact keys (`generate_compact` → `CompactSigningKey` / `CompactPublicKey` /
//! `CompactSignature`) are a distinct scheme from the dim-2 keys: a compact key
//! verifies only compact signatures, and a dim-2 key verifies only dim-2
//! signatures. `AnySignature` autodetects the wire format, but the public key
//! must match the scheme.

use rand::rngs::StdRng;
use rand::SeedableRng;

use sqisign_rs::{generate, generate_compact, AnySignature, CompactSignature, Level1, Verifier};

/// The headline clean-API test from the task description.
#[test]
fn compact_clean_api() {
    let mut rng = StdRng::seed_from_u64(0x00C0_FFEE);
    let (pk, sk) = generate_compact(&mut rng);
    let msg = b"hello compact";
    let sig = sk.sign(msg, &mut rng).expect("signing must succeed");
    assert_eq!(sig.to_bytes().len(), 108, "compact signatures are 108 bytes");
    pk.verify(msg, &sig).expect("a fresh compact signature must verify");
}

/// One key, several messages: every fresh compact signature verifies.
#[test]
fn compact_roundtrip_multiple_messages() {
    let mut rng = StdRng::seed_from_u64(1);
    let (pk, sk) = generate_compact(&mut rng);
    let messages: [&[u8]; 5] = [
        b"",
        b"hello world",
        b"compact round-trip",
        &[0u8; 32],
        &[0xABu8; 200],
    ];
    for (i, msg) in messages.iter().enumerate() {
        let sig = sk.sign(msg, &mut rng).expect("signing must succeed");
        assert_eq!(sig.to_bytes().len(), 108, "sig {i}: wire size");
        pk.verify(msg, &sig).unwrap_or_else(|_| panic!("sig {i} must verify"));
    }
}

/// Several independent keypairs (different seeds) each round-trip.
#[test]
fn compact_roundtrip_multiple_seeds() {
    for seed in 0..3u64 {
        let mut rng = StdRng::seed_from_u64(0xA5A5_0000 + seed);
        let (pk, sk) = generate_compact(&mut rng);
        let msg = b"per-seed compact message";
        let sig = sk.sign(msg, &mut rng).expect("signing must succeed");
        pk.verify(msg, &sig)
            .unwrap_or_else(|_| panic!("seed {seed} must verify"));
    }
}

/// `AnySignature::from_bytes` autodetects the 108-byte compact format, and a
/// `CompactPublicKey` verifies it.
#[test]
fn compact_anysignature_autodetect() {
    let mut rng = StdRng::seed_from_u64(42);
    let (pk, sk) = generate_compact(&mut rng);
    let msg = b"autodetect";
    let sig = sk.sign(msg, &mut rng).expect("signing must succeed");

    let bytes = sig.to_bytes();
    assert_eq!(bytes.len(), 108);
    let any = AnySignature::<Level1>::from_bytes(&bytes).expect("must autodetect compact format");
    assert!(matches!(any, AnySignature::Compact(_)), "108 B must route to the compact arm");
    pk.verify(msg, &any).expect("compact pk must verify the compact AnySignature");
}

/// Tampering the message, the signature bytes, or the public key all reject.
#[test]
fn compact_rejects_tampering() {
    let mut rng = StdRng::seed_from_u64(7);
    let (pk, sk) = generate_compact(&mut rng);
    let msg: &[u8] = b"authentic compact message";
    let sig = sk.sign(msg, &mut rng).expect("signing must succeed");

    // Baseline accept.
    pk.verify(msg, &sig).expect("baseline accept");

    // Wrong message → recomputed Fiat-Shamir challenge differs.
    assert!(
        pk.verify(b"different message", &sig).is_err(),
        "wrong message must reject"
    );

    // A flipped byte in each field (A_com, q, a, b, c_or_d): parsing or
    // verification must fail (never accept).
    let bytes = sig.to_bytes();
    for &pos in &[5usize, 70, 85, 95, 107] {
        let mut bad = bytes;
        bad[pos] ^= 0x01;
        let rejected = match CompactSignature::<Level1>::from_bytes(&bad) {
            Ok(s) => pk.verify(msg, &s).is_err(),
            Err(_) => true, // non-canonical bytes rejected at parse time
        };
        assert!(rejected, "flipped byte at {pos} must reject");
    }

    // A different compact public key.
    let (pk2, _sk2) = generate_compact(&mut rng);
    assert!(pk2.verify(msg, &sig).is_err(), "wrong public key must reject");
}

/// Cross-scheme rules: a compact key verifies compact signatures only; a dim-2
/// key verifies dim-2 signatures only. Routing is by the public key's scheme,
/// not just the autodetected format.
#[test]
fn cross_scheme_verification() {
    let mut rng = StdRng::seed_from_u64(99);
    let msg: &[u8] = b"cross-scheme";

    // A compact keypair + signature.
    let (cpk, csk) = generate_compact(&mut rng);
    let csig = csk.sign(msg, &mut rng).expect("compact sign");
    let csig_any = AnySignature::<Level1>::from_bytes(&csig.to_bytes()).unwrap();

    // A dim-2 keypair + (standard) signature.
    let (dpk, dsk) = generate::<Level1>(&mut rng);
    let dsig = dsk.sign(msg, &mut rng).expect("dim-2 sign");
    let dsig_any = AnySignature::<Level1>::from_bytes(&dsig.to_bytes()).unwrap();

    // Native verification works for each scheme.
    cpk.verify(msg, &csig).expect("compact pk verifies compact sig");
    cpk.verify(msg, &csig_any).expect("compact pk verifies compact AnySignature");
    dpk.verify(msg, &dsig).expect("dim-2 pk verifies dim-2 sig");
    dpk.verify(msg, &dsig_any).expect("dim-2 pk verifies dim-2 AnySignature");

    // Cross-scheme verification is rejected (right format, wrong key scheme).
    assert!(
        cpk.verify(msg, &dsig_any).is_err(),
        "compact pk must NOT verify a dim-2 signature"
    );
    assert!(
        dpk.verify(msg, &csig_any).is_err(),
        "dim-2 pk must NOT verify a compact signature"
    );

    // Type-level separation: `cpk.verify(msg, &dsig)` and `dpk.verify(msg, &csig)`
    // do not even compile - `CompactPublicKey` has no `Verifier<Signature>` impl
    // and `PublicKey` has no `Verifier<CompactSignature>` impl.
}

/// Timing/throughput bench (run explicitly:
/// `cargo test -p sqisign-rs --release --test compact_api compact_measure -- --nocapture --ignored`).
#[test]
#[ignore]
fn compact_measure() {
    use std::time::Instant;
    let mut rng = StdRng::seed_from_u64(0xBEEF);
    let t = Instant::now();
    let (pk, sk) = generate_compact(&mut rng);
    eprintln!("[compact] generate_compact: {:?}", t.elapsed());
    let msg: &[u8] = b"timing";
    for i in 0..5 {
        let t = Instant::now();
        let sig = sk.sign(msg, &mut rng).unwrap();
        let dt = t.elapsed();
        let t2 = Instant::now();
        let ok = pk.verify(msg, &sig).is_ok();
        eprintln!("[compact] sign #{i}: {dt:?}   verify: {:?}  (ok={ok})", t2.elapsed());
    }
}
