//! Every byte of a KAT signature flipped, one position at a time, must be
//! rejected (SECENG-2438): the signature of entry 0 of each level's KAT
//! file, the byte at every position replaced by its complement, decoded
//! and verified under the file's public key and message; the test prints
//! the position count and the distribution of where the verifier stops,
//! and one line per position, so that two builds can be compared. Levels
//! III and V run in release builds only.

use sqisign_verify::fp::FpBackend;
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::sqisign::{
    level1, level3, level5, public_key_from_bytes, signature_from_bytes, verify_reason,
    PreparedPublicKey, RejectReason, VerifyParams,
};
use sqisign_verify::{P324_3, P500_27, P664_17};

struct Entry {
    msg: Vec<u8>,
    pk: Vec<u8>,
    smlen: usize,
    sm: Vec<u8>,
}

// The first entry of a KAT file (count, seed, msg, pk, sk, smlen, sm).
fn first_entry(text: &str) -> Entry {
    let mut e = Entry {
        msg: vec![],
        pk: vec![],
        smlen: 0,
        sm: vec![],
    };
    let mut seen_count = 0;
    for line in text.lines() {
        let Some((k, v)) = line.split_once('=') else {
            continue;
        };
        let (k, v) = (k.trim(), v.trim());
        match k {
            "count" => {
                seen_count += 1;
                if seen_count > 1 {
                    break;
                }
            }
            "msg" => e.msg = hex::decode(v).unwrap(),
            "pk" => e.pk = hex::decode(v).unwrap(),
            "smlen" => e.smlen = v.parse().unwrap(),
            "sm" => e.sm = hex::decode(v).unwrap(),
            _ => {}
        }
    }
    e
}

fn run<L: FpBackend + PrimePrecomp>(params: &VerifyParams, rsp: &str, name: &str) {
    let backend = sqisign_kat::arithmetic_from_env();
    let e = first_entry(rsp);
    assert_eq!(e.sm.len(), e.smlen, "{name}: smlen");
    let siglen = e.smlen - e.msg.len();
    assert_eq!(siglen, params.signature_bytes, "{name}: signature length");
    let sig_bytes = &e.sm[..siglen];
    let pk = public_key_from_bytes::<L>(params, &e.pk).expect("pk decode");
    let key = PreparedPublicKey::new(params, &pk).expect("prepare");
    let sig = signature_from_bytes::<L>(params, sig_bytes).expect("sig decode");
    assert_eq!(
        verify_reason(params, &key, &sig, &e.msg),
        Ok(()),
        "{name}: the untampered signature verifies"
    );
    let (mut decode, mut matrix, mut challenge, mut response) = (0usize, 0usize, 0usize, 0usize);
    for pos in 0..siglen {
        let mut t = sig_bytes.to_vec();
        t[pos] ^= 0xff;
        let reason = match signature_from_bytes::<L>(params, &t) {
            None => "decode",
            Some(ts) => match verify_reason(params, &key, &ts, &e.msg) {
                Ok(()) => "ACCEPTED",
                Err(RejectReason::BasisChangeMatrix) => "matrix",
                Err(RejectReason::ChallengeState) => "challenge",
                Err(RejectReason::Response) => "response",
            },
        };
        println!("tampered {name} pos {pos} {reason}");
        assert_ne!(reason, "ACCEPTED", "{name}: position {pos} accepted");
        match reason {
            "decode" => decode += 1,
            "matrix" => matrix += 1,
            "challenge" => challenge += 1,
            _ => response += 1,
        }
    }
    println!(
        "tampered {name} ({backend}) positions {siglen}: decode {decode} matrix {matrix} challenge {challenge} response {response}; all rejected"
    );
}

#[test]
fn tampered_level1() {
    run::<P324_3>(
        &level1(),
        include_str!("../kat/PQCsignKAT_270_SQIsign_p324_3.rsp"),
        "level I",
    );
}

#[test]
fn tampered_level3() {
    if cfg!(debug_assertions) {
        return;
    }
    run::<P500_27>(
        &level3(),
        include_str!("../kat/PQCsignKAT_417_SQIsign_p500_27.rsp"),
        "level III",
    );
}

#[test]
fn tampered_level5() {
    if cfg!(debug_assertions) {
        return;
    }
    run::<P664_17>(
        &level5(),
        include_str!("../kat/PQCsignKAT_549_SQIsign_p664_17.rsp"),
        "level V",
    );
}
