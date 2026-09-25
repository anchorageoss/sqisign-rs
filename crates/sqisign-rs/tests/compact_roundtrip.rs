//! The compact (dimension-4) format at level I: key generation, signing,
//! verification, the wire sizes, and the negatives (wrong message,
//! tampered bytes, wrong key). Feature `compact`.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::level1;
use sqisign_verify::hd::{hd_verify_bytes, HdReject};
use sqisign_verify::P324_3;

const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;

#[test]
fn sizes() {
    assert_eq!(compact_signature_bytes::<P324_3>(), 142);
    assert_eq!(compact_public_key_bytes::<P324_3>(), 84);
}

#[test]
fn keygen_sign_verify_level1() {
    let params = level1();
    let mut rng = ShakeRng::new(b"compact round trip", b"tst");
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen");
    let mut pkb = [0u8; 84];
    let n = compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    assert_eq!(n, 84);
    let msg = b"compact round three";
    let mut sig = [0u8; 142];
    let m = compact_sign::<P324_3, N>(&params, &pk, &sk, msg, &mut rng, &mut sig).expect("sign");
    assert_eq!(m, 142);
    assert_eq!(
        hd_verify_bytes::<P324_3>(&sig, &pkb, msg),
        Ok(()),
        "round trip"
    );

    // negatives
    assert!(hd_verify_bytes::<P324_3>(&sig, &pkb, b"another message").is_err());
    let mut t = sig;
    t[0] ^= 1;
    assert!(
        hd_verify_bytes::<P324_3>(&t, &pkb, msg).is_err(),
        "tampered A_com"
    );
    let mut t = sig;
    t[82] ^= 1; // q
    assert!(
        hd_verify_bytes::<P324_3>(&t, &pkb, msg).is_err(),
        "tampered q"
    );
    let mut t = sig;
    t[82 + 22] ^= 1; // a
    assert!(
        hd_verify_bytes::<P324_3>(&t, &pkb, msg).is_err(),
        "tampered a"
    );
    assert_eq!(
        hd_verify_bytes::<P324_3>(&sig[..141], &pkb, msg),
        Err(HdReject::MalformedInput)
    );
    let (pk2, _) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen 2");
    let mut pkb2 = [0u8; 84];
    compact_public_key_to_bytes(&pk2, &mut pkb2).unwrap();
    assert!(
        hd_verify_bytes::<P324_3>(&sig, &pkb2, msg).is_err(),
        "wrong key"
    );
}

#[test]
fn typed_api_level1() {
    use signature::Verifier;
    use sqisign_rs::{generate_compact, CompactPublicKey, CompactSignature};
    let mut rng = ShakeRng::new(b"compact typed", b"tst");
    let (pk, sk) = generate_compact::<P324_3, N>(level1(), &mut rng).expect("keygen");
    let msg = b"typed";
    let sig = sk.sign(msg, &mut rng).expect("sign");
    assert_eq!(sig.to_bytes().len(), 142);
    assert_eq!(pk.to_bytes().len(), 84);
    assert!(pk.verify(msg, &sig).is_ok());
    assert!(pk.verify(b"other", &sig).is_err());
    let pk2 = CompactPublicKey::<P324_3>::from_bytes(pk.as_bytes()).unwrap();
    let sig2 = CompactSignature::<P324_3>::from_bytes(sig.as_bytes()).unwrap();
    assert!(pk2.verify_compact(msg, &sig2).is_ok());
    assert!(CompactSignature::<P324_3>::from_bytes(&sig.as_bytes()[..141]).is_err());
}
