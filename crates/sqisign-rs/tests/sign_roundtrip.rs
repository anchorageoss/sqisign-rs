//! End-to-end round-trip test: keygen → sign → verify.

use sqisign_rs::keygen::keypair;
use sqisign_rs::params::Level1;
use sqisign_rs::sign::sign;
use sqisign_rs::Verifier;

type L1 = Level1;

#[test]
fn sign_verify_roundtrip() {
    let mut rng = rand::thread_rng();

    let (pk, sk) = keypair::<L1>(&mut rng);

    let msg = b"SQIsign round-trip test message";
    let sig = sign::<L1>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    assert!(pk.verify(msg, &sig).is_ok(), "valid signature must verify");
}

#[test]
fn sign_verify_empty_message() {
    let mut rng = rand::thread_rng();

    let (pk, sk) = keypair::<L1>(&mut rng);

    let sig = sign::<L1>(&sk, &pk, b"", &mut rng).expect("signing must succeed");

    assert!(
        pk.verify(b"", &sig).is_ok(),
        "empty message signature must verify"
    );
}

#[test]
fn sign_verify_wrong_message_fails() {
    let mut rng = rand::thread_rng();

    let (pk, sk) = keypair::<L1>(&mut rng);
    let msg = b"correct message";
    let wrong_msg = b"wrong message";

    let sig = sign::<L1>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    assert!(
        pk.verify(wrong_msg, &sig).is_err(),
        "signature must not verify under wrong message"
    );
}

#[test]
fn sign_verify_wrong_key_fails() {
    let mut rng = rand::thread_rng();

    let (pk1, _sk1) = keypair::<L1>(&mut rng);
    let (pk2, _sk2) = keypair::<L1>(&mut rng);

    let msg = b"test message";
    let sig = sign::<L1>(&_sk1, &pk1, msg, &mut rng).expect("signing must succeed");

    assert!(
        pk2.verify(msg, &sig).is_err(),
        "signature must not verify under wrong public key"
    );
}

// ---------------------------------------------------------------------------
// signature crate trait-based API
// ---------------------------------------------------------------------------

#[test]
fn trait_randomized_signer_roundtrip() {
    use sqisign_rs::signature::{RandomizedSigner, Verifier};

    let mut rng = rand::thread_rng();
    let (_pk, signing_key) = sqisign_rs::generate::<L1>(&mut rng);

    let msg = b"trait-based signing test";
    let sig: sqisign_rs::Signature<L1> = signing_key.sign_with_rng(&mut rng, msg);

    let pk = sqisign_rs::signature::Keypair::verifying_key(&signing_key);
    pk.verify(msg, &sig)
        .expect("trait-based verify must succeed");
}

#[test]
fn trait_verifier_rejects_wrong_message() {
    use sqisign_rs::signature::{RandomizedSigner, Verifier};

    let mut rng = rand::thread_rng();
    let (_pk, signing_key) = sqisign_rs::generate::<L1>(&mut rng);

    let sig: sqisign_rs::Signature<L1> = signing_key.sign_with_rng(&mut rng, b"correct");

    let pk = signing_key.public_key();
    assert!(
        pk.verify(b"wrong", &sig).is_err(),
        "trait-based verify must reject wrong message"
    );
}

#[test]
fn trait_signature_encoding_roundtrip() {
    use sqisign_rs::signature::RandomizedSigner;

    let mut rng = rand::thread_rng();
    let (_pk, signing_key) = sqisign_rs::generate::<L1>(&mut rng);

    let sig: sqisign_rs::Signature<L1> = signing_key.sign_with_rng(&mut rng, b"encoding test");

    let bytes = sig.to_bytes();
    let bytes_slice: &[u8] = bytes.as_slice();
    let sig2 = sqisign_rs::Signature::<L1>::try_from(bytes_slice)
        .expect("SignatureEncoding round-trip must succeed");

    let enc1 = sig.to_bytes();
    let enc2 = sig2.to_bytes();
    assert_eq!(
        enc1.as_slice(),
        enc2.as_slice(),
        "encoded signatures must match"
    );
}
