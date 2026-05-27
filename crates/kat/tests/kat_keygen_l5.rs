//! KAT validation for SQIsign Level 5 keygen.
//!
//! Seeds the NIST AES-256-CTR-DRBG with each test vector's seed,
//! runs `protocols_keygen`, and compares the serialized pk and sk
//! against the expected bytes from the `.rsp` file.

use rayon::prelude::*;
use sqisign_kat::kat_parser;
use sqisign_kat::nist_drbg::NistDrbg;
use sqisign_rs::id2iso::sign_precomp::SigningPrecomp;
use sqisign_rs::keygen::keygen::protocols_keygen;
use sqisign_rs::Level5;

const KAT_RSP: &str = include_str!("../../../reference/KAT/PQCsignKAT_701_SQIsign_lvl5.rsp");

#[test]
fn kat_keygen_level5() {
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

        let (pk, sk) = protocols_keygen::<Level5>(&mut rng, &precomp);

        let pk_bytes = pk.to_bytes();
        assert_eq!(
            &pk_bytes[..],
            &entry.pk[..],
            "pk mismatch for count={}",
            entry.count
        );

        let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
        let pk_len = pk_bytes.len();
        let expected_sk = &entry.sk[pk_len..];
        assert_eq!(
            &sk_bytes[..expected_sk.len()],
            expected_sk,
            "sk mismatch for count={}",
            entry.count
        );
    });
}
