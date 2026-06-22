//! Round-trip tests for serialization/deserialization of keys and signatures,
//! and end-to-end sign-then-verify through the serialization boundary.

use sqisign_rs::keygen::keypair;
use sqisign_rs::params::Level1;
use sqisign_rs::sign::sign;
use sqisign_rs::types::{PublicKey, Signature};
use sqisign_rs::Verifier;

type L = Level1;

// ---------------------------------------------------------------------------
// PublicKey round-trip
// ---------------------------------------------------------------------------

#[test]
fn pk_to_bytes_from_bytes_roundtrip() {
    let mut rng = rand::thread_rng();
    let (pk, _sk) = keypair::<L>(&mut rng);

    let enc = pk.to_bytes();
    let pk2 = PublicKey::<L>::from_bytes(&enc).expect("pk deserialization must succeed");

    assert_eq!(pk.to_bytes(), pk2.to_bytes(), "pk round-trip mismatch");
}

#[test]
fn pk_try_from_slice() {
    let mut rng = rand::thread_rng();
    let (pk, _sk) = keypair::<L>(&mut rng);

    let enc = pk.to_bytes();
    let pk2: PublicKey<L> = (&enc[..])
        .try_into()
        .expect("TryFrom<&[u8]> must succeed for valid pk bytes");

    assert_eq!(pk.to_bytes(), pk2.to_bytes());
}

#[test]
fn pk_from_bytes_rejects_short_input() {
    let short = [0u8; 4];
    assert!(
        PublicKey::<L>::from_bytes(&short).is_err(),
        "short input must be rejected"
    );
}

#[test]
fn pk_from_bytes_rejects_empty_input() {
    assert!(
        PublicKey::<L>::from_bytes(&[]).is_err(),
        "empty input must be rejected"
    );
}

// ---------------------------------------------------------------------------
// Signature round-trip
// ---------------------------------------------------------------------------

#[test]
fn sig_to_bytes_from_bytes_roundtrip() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let msg = b"signature serialization test";

    let sig = sign::<L>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    let enc = sig.to_bytes();
    let sig2 = Signature::<L>::from_bytes(&enc).expect("sig deserialization must succeed");

    assert_eq!(sig.to_bytes(), sig2.to_bytes(), "sig round-trip mismatch");
}

#[test]
fn sig_try_from_slice() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let msg = b"try_from test";

    let sig = sign::<L>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    let enc = sig.to_bytes();
    let sig2: Signature<L> = (&enc[..])
        .try_into()
        .expect("TryFrom<&[u8]> must succeed for valid sig bytes");

    assert_eq!(sig.to_bytes(), sig2.to_bytes());
}

#[test]
fn sig_from_bytes_rejects_short_input() {
    let short = [0u8; 4];
    assert!(
        Signature::<L>::from_bytes(&short).is_err(),
        "short input must be rejected"
    );
}

#[test]
fn sig_deserialized_still_verifies() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);
    let msg = b"deserialized verification test";

    let sig = sign::<L>(&sk, &pk, msg, &mut rng).expect("signing must succeed");

    let enc = sig.to_bytes();
    let sig2 = Signature::<L>::from_bytes(&enc).expect("sig deserialization must succeed");

    assert!(
        pk.verify(msg, &sig2).is_ok(),
        "deserialized signature must still verify"
    );
}

// ---------------------------------------------------------------------------
// SecretKey round-trip
// ---------------------------------------------------------------------------

#[test]
fn sk_to_bytes_from_bytes_roundtrip() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);

    let enc = sk.to_bytes().expect("sk encoding must succeed");
    let mut sk2 = sqisign_rs::sign::SecretKey::<L>::from_bytes(&enc)
        .expect("sk deserialization must succeed");
    sk2.populate_from_pk(&pk);

    let enc2 = sk2.to_bytes().expect("re-encoding must succeed");
    assert_eq!(&enc[..], &enc2[..], "sk round-trip mismatch");
}

#[test]
fn sk_from_bytes_rejects_short_input() {
    let short = [0u8; 4];
    assert!(
        sqisign_rs::sign::SecretKey::<L>::from_bytes(&short).is_err(),
        "short sk input must be rejected"
    );
}

#[test]
fn sk_from_bytes_rejects_non_invertible_matrix() {
    use sqisign_rs::quaternion::types::IbzMat2x2;
    let mut rng = rand::thread_rng();
    let (_pk, mut sk) = keypair::<L>(&mut rng);

    // A genuine key deserializes.
    let good = sk.to_bytes().expect("sk encoding must succeed");
    assert!(
        sqisign_rs::sign::SecretKey::<L>::from_bytes(&good).is_ok(),
        "a valid secret key must deserialize"
    );

    // The basis-change matrix must be invertible mod 2^TORSION_EVEN_POWER. A
    // zero (non-invertible) matrix must be rejected at deserialization rather
    // than accepted and then failing when the key is first used to sign.
    sk.mat_ba_can_to_ba0_two = IbzMat2x2::default();
    let bad = sk.to_bytes().expect("sk encoding must succeed");
    assert!(
        sqisign_rs::sign::SecretKey::<L>::from_bytes(&bad).is_err(),
        "a non-invertible basis-change matrix must be rejected"
    );
}

// ---------------------------------------------------------------------------
// Full end-to-end: keygen → serialize → deserialize → sign → verify
// ---------------------------------------------------------------------------

#[test]
fn full_roundtrip_through_serialization() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);

    // Serialize both keys
    let pk_bytes = pk.to_bytes();
    let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");

    // Deserialize
    let pk2 = PublicKey::<L>::from_bytes(&pk_bytes).expect("pk deserialization must succeed");
    let mut sk2 = sqisign_rs::sign::SecretKey::<L>::from_bytes(&sk_bytes)
        .expect("sk deserialization must succeed");
    sk2.populate_from_pk(&pk2);

    // Sign with deserialized keys
    let msg = b"full serialization round-trip test";
    let sig = sign::<L>(&sk2, &pk2, msg, &mut rng).expect("signing must succeed");

    // Serialize and deserialize the signature
    let sig_bytes = sig.to_bytes();
    let sig2 = Signature::<L>::from_bytes(&sig_bytes).expect("sig deserialization must succeed");

    // Verify with deserialized public key and signature
    assert!(
        pk2.verify(msg, &sig2).is_ok(),
        "full round-trip through serialization must verify"
    );
}

#[test]
fn sign_with_deserialized_sk_wrong_message_fails() {
    let mut rng = rand::thread_rng();
    let (pk, sk) = keypair::<L>(&mut rng);

    let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
    let mut sk2 = sqisign_rs::sign::SecretKey::<L>::from_bytes(&sk_bytes)
        .expect("sk deserialization must succeed");
    sk2.populate_from_pk(&pk);

    let sig = sign::<L>(&sk2, &pk, b"correct", &mut rng).expect("signing must succeed");

    assert!(
        pk.verify(b"wrong", &sig).is_err(),
        "wrong message must not verify"
    );
}

// ---------------------------------------------------------------------------
// SigningKey round-trip
// ---------------------------------------------------------------------------

#[test]
fn signing_key_to_bytes_from_bytes_roundtrip() {
    let mut rng = rand::thread_rng();
    let (_pk, signing_key) = sqisign_rs::generate::<L>(&mut rng);

    let enc = signing_key
        .to_bytes()
        .expect("signing key encoding must succeed");
    let signing_key2 = sqisign_rs::SigningKey::<L>::from_bytes(&enc)
        .expect("signing key deserialization must succeed");

    let enc2 = signing_key2.to_bytes().expect("re-encoding must succeed");
    assert_eq!(enc, enc2, "signing key round-trip mismatch");
}

#[test]
fn signing_key_from_bytes_rejects_short_input() {
    let short = [0u8; 4];
    assert!(
        sqisign_rs::SigningKey::<L>::from_bytes(&short).is_err(),
        "short signing key input must be rejected"
    );
}

#[test]
fn signing_key_deserialized_can_sign_and_verify() {
    let mut rng = rand::thread_rng();
    let (pk, signing_key) = sqisign_rs::generate::<L>(&mut rng);

    let enc = signing_key
        .to_bytes()
        .expect("signing key encoding must succeed");
    let signing_key2 = sqisign_rs::SigningKey::<L>::from_bytes(&enc)
        .expect("signing key deserialization must succeed");

    let msg = b"deserialized signing key test";
    let sig = signing_key2
        .sign(msg, &mut rng)
        .expect("signing with deserialized key must succeed");

    assert!(
        pk.verify(msg, &sig).is_ok(),
        "signature from deserialized key must verify against original pk"
    );
    assert!(
        signing_key2.public_key().verify(msg, &sig).is_ok(),
        "signature must verify against deserialized pk"
    );
}
