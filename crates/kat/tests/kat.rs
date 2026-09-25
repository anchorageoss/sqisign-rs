//! The NIST known-answer tests of the round-3 reference (`KAT/*.rsp` of
//! the-sqisign at tag `nist-v3`, 100 entries per level), reproduced byte
//! for byte by the `v3` implementation: for each entry the AES-256
//! CTR-DRBG is seeded with the entry's seed, key generation draws its root
//! seed, the encoded keys must match, signing (from the re-decoded secret
//! key, as the KAT generator does) draws its root seed, `sm = sig || msg`
//! must match, and the signature must verify (cold and under a prepared
//! key) and fail on a modified message.
//!
//! All 100 entries run in release; a debug build runs the first two unless
//! `SQISIGN_KAT_ENTRIES` says otherwise. `SQISIGN_FORCE_PORTABLE=1` runs
//! the portable field kernels of the x86-64 backends.

use rayon::prelude::*;
use sqisign_kat::nist_drbg::NistDrbg;
use sqisign_rs::sqisign::{
    keygen, level1, level3, level5, secret_key_from_bytes, secret_key_to_bytes, sign, Params,
};
use sqisign_verify::fp::FpBackend;
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::sqisign::{
    public_key_from_bytes, public_key_to_bytes, signature_from_bytes, signature_to_bytes, verify,
    verify_prepared, PreparedPublicKey, MAX_PUBLICKEY_BYTES, MAX_SIGNATURE_BYTES,
};
use sqisign_verify::{P324_3, P500_27, P664_17};

struct Entry {
    count: usize,
    seed: [u8; 48],
    msg: Vec<u8>,
    pk: Vec<u8>,
    sk: Vec<u8>,
    smlen: usize,
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
                seed: [0u8; 48],
                msg: vec![],
                pk: vec![],
                sk: vec![],
                smlen: 0,
                sm: vec![],
            });
            continue;
        }
        let e = cur.as_mut().expect("count first");
        match k {
            "seed" => e.seed = hex::decode(v).unwrap().try_into().expect("48-byte seed"),
            "msg" => e.msg = hex::decode(v).unwrap(),
            "pk" => e.pk = hex::decode(v).unwrap(),
            "sk" => e.sk = hex::decode(v).unwrap(),
            "smlen" => e.smlen = v.parse().unwrap(),
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

fn run<L: FpBackend + PrimePrecomp, const N: usize>(params: &Params<N>, rsp: &str, name: &str) {
    let backend = sqisign_kat::arithmetic_from_env();
    let entries = parse_rsp(rsp);
    assert_eq!(entries.len(), 100, "{name}: KAT entries");
    let n = entries_to_run().min(entries.len());
    entries[..n].par_iter().for_each(|e| {
        let mut drbg = NistDrbg::new(&e.seed);
        let (pk, sk) = keygen::<L, N>(params, &mut drbg).expect("keygen");
        let mut pkb = [0u8; MAX_PUBLICKEY_BYTES];
        let pklen = public_key_to_bytes(&pk, &mut pkb);
        assert_eq!(
            hex::encode(&pkb[..pklen]),
            hex::encode(&e.pk),
            "{name} count {}: pk ({backend})",
            e.count
        );
        let skb = secret_key_to_bytes(params, &sk, &pk);
        assert_eq!(
            hex::encode(&skb),
            hex::encode(&e.sk),
            "{name} count {}: sk",
            e.count
        );
        // the KAT generator signs with the decoded secret key
        let (sk2, pk2) = secret_key_from_bytes::<L, N>(params, &skb).expect("sk decode");
        let sig = sign::<L, N>(params, &pk2, &sk2, &e.msg, &mut drbg).expect("sign");
        let mut sigb = [0u8; MAX_SIGNATURE_BYTES];
        let siglen = signature_to_bytes(&params.verify, &sig, &mut sigb).expect("encode");
        let mut sm = sigb[..siglen].to_vec();
        sm.extend_from_slice(&e.msg);
        assert_eq!(sm.len(), e.smlen, "{name} count {}: smlen", e.count);
        assert_eq!(
            hex::encode(&sm),
            hex::encode(&e.sm),
            "{name} count {}: sm",
            e.count
        );
        // verification from the file's bytes, cold and prepared
        let pk3 = public_key_from_bytes::<L>(&params.verify, &e.pk).expect("pk decode");
        let sig3 = signature_from_bytes::<L>(&params.verify, &e.sm[..siglen]).expect("sig decode");
        assert!(
            verify(&params.verify, &pk3, &sig3, &e.msg),
            "{name} count {}: verify",
            e.count
        );
        let prepared = PreparedPublicKey::new(&params.verify, &pk3).expect("prepare");
        assert!(
            verify_prepared(&params.verify, &prepared, &sig3, &e.msg),
            "{name} count {}: verify, prepared key",
            e.count
        );
        let mut bad = e.msg.clone();
        match bad.first_mut() {
            Some(b) => *b ^= 1,
            None => bad.push(0),
        }
        assert!(
            !verify(&params.verify, &pk3, &sig3, &bad),
            "{name} count {}: a modified message verified",
            e.count
        );
    });
}

#[test]
fn kat_v3_level1() {
    run::<P324_3, { sqisign_rs::precomp::p324_3::IBZ_NLIMBS }>(
        &level1(),
        include_str!("../kat/PQCsignKAT_270_SQIsign_p324_3.rsp"),
        "p324_3",
    );
}

#[test]
fn kat_v3_level3() {
    run::<P500_27, { sqisign_rs::precomp::p500_27::IBZ_NLIMBS }>(
        &level3(),
        include_str!("../kat/PQCsignKAT_417_SQIsign_p500_27.rsp"),
        "p500_27",
    );
}

#[test]
fn kat_v3_level5() {
    run::<P664_17, { sqisign_rs::precomp::p664_17::IBZ_NLIMBS }>(
        &level5(),
        include_str!("../kat/PQCsignKAT_549_SQIsign_p664_17.rsp"),
        "p664_17",
    );
}
