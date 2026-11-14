//! The compact (dimension-4) format at levels I, III and V: key
//! generation, signing, verification, the wire sizes, and the negatives
//! (wrong message, tampered bytes, wrong key). Feature `compact`.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::{hd_verify_bytes, HdLevel, HdReject};
use sqisign_verify::{P324_3, P500_27, P664_17};

const N1: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const N3: usize = sqisign_rs::precomp::p500_27::IBZ_NLIMBS;
const N5: usize = sqisign_rs::precomp::p664_17::IBZ_NLIMBS;

/// Public key and signature, in bytes, with the byte-aligned layout.
const SIZES: [(usize, usize); 3] = [(84, 142), (130, 214), (170, 280)];

#[test]
fn sizes() {
    assert_eq!(compact_public_key_bytes::<P324_3>(), SIZES[0].0);
    assert_eq!(compact_signature_bytes::<P324_3>(), SIZES[0].1);
    assert_eq!(compact_public_key_bytes::<P500_27>(), SIZES[1].0);
    assert_eq!(compact_signature_bytes::<P500_27>(), SIZES[1].1);
    assert_eq!(compact_public_key_bytes::<P664_17>(), SIZES[2].0);
    assert_eq!(compact_signature_bytes::<P664_17>(), SIZES[2].1);
}

fn round_trip<L: HdLevel, const N: usize>(params: &Params<N>, label: &[u8]) {
    let (pk_len, sig_len) = (
        compact_public_key_bytes::<L>(),
        compact_signature_bytes::<L>(),
    );
    let mut rng = ShakeRng::new(label, b"tst");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut rng).expect("keygen");
    let mut pkb = vec![0u8; pk_len];
    assert_eq!(compact_public_key_to_bytes(&pk, &mut pkb), Some(pk_len));
    let msg = b"compact round three";
    let mut sig = vec![0u8; sig_len];
    let m = compact_sign::<L, N>(params, &pk, &sk, msg, &mut rng, &mut sig).expect("sign");
    assert_eq!(m, sig_len);
    assert_eq!(hd_verify_bytes::<L>(&sig, &pkb, msg), Ok(()), "round trip");

    // negatives
    assert!(hd_verify_bytes::<L>(&sig, &pkb, b"another message").is_err());
    for (what, byte) in [
        ("A_com", 0usize),
        ("q", L::FP2_BYTES + 1),
        ("a", L::FP2_BYTES + L::Q_BYTES + 1),
    ] {
        let mut t = sig.clone();
        t[byte] ^= 1;
        assert!(
            hd_verify_bytes::<L>(&t, &pkb, msg).is_err(),
            "tampered {what}"
        );
    }
    assert_eq!(
        hd_verify_bytes::<L>(&sig[..sig_len - 1], &pkb, msg),
        Err(HdReject::MalformedInput)
    );
    let (pk2, _) = compact_keygen::<L, N>(params, &mut rng).expect("keygen 2");
    let mut pkb2 = vec![0u8; pk_len];
    compact_public_key_to_bytes(&pk2, &mut pkb2).unwrap();
    assert!(hd_verify_bytes::<L>(&sig, &pkb2, msg).is_err(), "wrong key");
}

#[test]
fn keygen_sign_verify_level1() {
    round_trip::<P324_3, N1>(&level1(), b"compact round trip");
}

#[test]
fn keygen_sign_verify_level3() {
    round_trip::<P500_27, N3>(&level3(), b"compact round trip III");
}

#[test]
fn keygen_sign_verify_level5() {
    round_trip::<P664_17, N5>(&level5(), b"compact round trip V");
}

#[test]
fn typed_api() {
    use signature::Verifier;
    use sqisign_rs::{generate_compact, CompactPublicKey, CompactSignature};
    let mut rng = ShakeRng::new(b"compact typed", b"tst");
    let (pk, sk) = generate_compact::<P324_3, N1>(level1(), &mut rng).expect("keygen");
    let msg = b"typed";
    let sig = sk.sign(msg, &mut rng).expect("sign");
    assert_eq!(sig.to_bytes().len(), SIZES[0].1);
    assert_eq!(pk.to_bytes().len(), SIZES[0].0);
    assert!(pk.verify(msg, &sig).is_ok());
    assert!(pk.verify(b"other", &sig).is_err());
    let pk2 = CompactPublicKey::<P324_3>::from_bytes(pk.as_bytes()).unwrap();
    let sig2 = CompactSignature::<P324_3>::from_bytes(sig.as_bytes()).unwrap();
    assert!(pk2.verify_compact(msg, &sig2).is_ok());
    assert!(CompactSignature::<P324_3>::from_bytes(&sig.as_bytes()[..SIZES[0].1 - 1]).is_err());

    let (pk, sk) = generate_compact::<P500_27, N3>(level3(), &mut rng).expect("keygen III");
    let sig = sk.sign(msg, &mut rng).expect("sign III");
    assert_eq!((pk.to_bytes().len(), sig.to_bytes().len()), SIZES[1]);
    assert!(pk.verify(msg, &sig).is_ok());
    let (pk, sk) = generate_compact::<P664_17, N5>(level5(), &mut rng).expect("keygen V");
    let sig = sk.sign(msg, &mut rng).expect("sign V");
    assert_eq!((pk.to_bytes().len(), sig.to_bytes().len()), SIZES[2]);
    assert!(pk.verify(msg, &sig).is_ok());
}
