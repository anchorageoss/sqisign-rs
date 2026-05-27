//! KAT validation for SQIsign Level 5 signing.
//!
//! Seeds the NIST AES-256-CTR-DRBG with each test vector's seed,
//! runs `protocols_keygen` then `protocols_sign`, and compares the
//! serialized signature against the expected bytes from the `.rsp` file.
//!
//! Also verifies the signature with `pk.verify()`.

use rayon::prelude::*;
use sqisign_kat::kat_parser;
use sqisign_kat::nist_drbg::NistDrbg;
use sqisign_rs::id2iso::sign_precomp::SigningPrecomp;
use sqisign_rs::keygen::keygen::protocols_keygen;
use sqisign_rs::sign::sign::protocols_sign;
use sqisign_rs::{Level5, Verifier};

const KAT_RSP: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02X}", b)).collect()
}

#[test]
fn kat_sign_level5() {
    let entries = kat_parser::parse_rsp(KAT_RSP);
    assert!(!entries.is_empty(), "no KAT entries parsed");

    let precomp = SigningPrecomp::<Level5>::level5();

    entries.par_iter().for_each(|entry| {
        let seed: [u8; 48] = entry
            .seed
            .as_slice()
            .try_into()
            .expect("seed must be 48 bytes");
        let mut rng = NistDrbg::new(&seed);

        // 1. Keygen
        let (pk, sk) = protocols_keygen::<Level5>(&mut rng, &precomp);

        let pk_bytes = pk.to_bytes();
        assert_eq!(
            &pk_bytes[..],
            &entry.pk[..],
            "pk mismatch for count={}",
            entry.count
        );

        // Round-trip sk through serialization (the NIST API calls
        // secret_key_to_bytes then secret_key_from_bytes between
        // keygen and signing).
        let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
        let mut sk = sqisign_rs::SecretKey::<Level5>::from_bytes(&sk_bytes)
            .expect("sk round-trip must succeed");
        sk.populate_from_pk(&pk);

        // 2. Sign
        let sig =
            protocols_sign::<Level5>(&pk, &sk, &entry.msg, &mut rng).expect("signing must succeed");
        let sig_bytes = sig.to_bytes();

        // 3. Compare signature bytes against expected sm = sig || msg
        let sig_len = sig_bytes.len();
        let expected_sig = &entry.sm[..sig_len];
        let expected_msg = &entry.sm[sig_len..];

        assert_eq!(
            expected_msg,
            &entry.msg[..],
            "sm suffix must equal msg for count={}",
            entry.count
        );

        if &sig_bytes[..] != expected_sig {
            panic!(
                "signature mismatch for count={}!\n  got:  {}\n  want: {}",
                entry.count,
                to_hex(&sig_bytes),
                to_hex(expected_sig)
            );
        }

        // 4. Verify
        assert!(
            pk.verify(&entry.msg, &sig).is_ok(),
            "signature must verify for count={}",
            entry.count
        );
    });
}
