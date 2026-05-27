//! KAT validation for SQIsign Level 5 verification.
//!
//! Deserializes the public key and signature directly from the `.rsp` file
//! (no keygen or signing involved) and verifies each signature. This tests
//! the verification path in isolation from the signing path.

use rayon::prelude::*;
use sqisign_kat::kat_parser;
use sqisign_verify::{Level5, PublicKey, Signature, Verifier};

const KAT_RSP: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

#[test]
fn kat_verify_level5() {
    let entries = kat_parser::parse_rsp(KAT_RSP);
    assert!(!entries.is_empty(), "no KAT entries parsed");

    entries.par_iter().for_each(|entry| {
        let pk = PublicKey::<Level5>::from_bytes(&entry.pk)
            .unwrap_or_else(|_| panic!("pk deserialization failed for count={}", entry.count));

        let sig_len = entry.sm.len() - entry.msg.len();
        let sig_bytes = &entry.sm[..sig_len];

        let sig = Signature::<Level5>::from_bytes(sig_bytes)
            .unwrap_or_else(|_| panic!("sig deserialization failed for count={}", entry.count));

        assert!(
            pk.verify(&entry.msg, &sig).is_ok(),
            "verification failed for count={}",
            entry.count
        );
    });
}
