//! Self-generated vectors of the compact format at level I in the SQIsignHD
//! library's text format, for the Sage oracle
//! (`sqisignhd-harness/round3/verify_r3.py`): `COMPACT_VECTORS_OUT=<dir>`
//! writes `Public_keys_r3lvl1.txt`, `Signatures_r3lvl1.txt` (accepted) and
//! `Signatures_r3lvl1_bad.txt` (rejected: a tampered scalar, a wrong degree,
//! a wrong challenge, in turn). Without the variable the test only checks
//! that the vectors verify here. Deterministic seeds: the files are the
//! format's known-answer vectors.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::level1;
use sqisign_verify::hd::{challenge_limbs, hd_verify_bytes, parse_public_key, parse_signature};
use sqisign_verify::P324_3;

const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const COUNT: usize = 5;

fn hex_be(bytes_le: &[u8]) -> String {
    let mut s = String::new();
    let mut started = false;
    for b in bytes_le.iter().rev() {
        if !started && *b == 0 {
            continue;
        }
        started = true;
        s.push_str(&format!("{b:02x}"));
    }
    if s.is_empty() {
        s.push('0');
    }
    s
}

fn fp2_line(name: &str, enc: &[u8]) -> String {
    let half = enc.len() / 2;
    format!(
        "{name} = 0x{} + i*0x{}\n",
        hex_be(&enc[..half]),
        hex_be(&enc[half..])
    )
}

#[test]
fn vectors_level1() {
    let params = level1();
    let mut pks = String::new();
    let mut sigs = String::new();
    let mut bad = String::new();
    for i in 0..COUNT {
        let label = format!("compact level I vector {i}");
        let mut rng = ShakeRng::new(label.as_bytes(), b"kat");
        let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen");
        let mut pkb = [0u8; 84];
        compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
        let msg = format!("message {i}").into_bytes();
        let mut sig = [0u8; 142];
        compact_sign::<P324_3, N>(&params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
        assert_eq!(hd_verify_bytes::<P324_3>(&sig, &pkb, &msg), Ok(()));
        let s = parse_signature::<P324_3>(&sig).unwrap();
        let p = parse_public_key::<P324_3>(&pkb).unwrap();
        let chal = challenge_limbs::<P324_3>(&s, &p, &msg).unwrap();

        pks.push_str(&fp2_line("A_pk", &p.a_pk.encode()));
        pks.push_str(&format!(
            "hint_pk_P = {}\nhint_pk_Q = {}\n",
            p.hint_pk_p, p.hint_pk_q
        ));

        let scalar = |v: i128| format!("{v:x}");
        let mut qb = s.q.to_le_bytes().to_vec();
        qb.truncate(32);
        let mut cb = Vec::new();
        for l in chal.iter() {
            cb.extend_from_slice(&l.to_le_bytes());
        }
        let record = |a: i128, q_hex: &str, chal_hex: &str| {
            format!(
                "{}a = {}\nb = {}\nc_or_d = {}\nq = {}\nhint_com_P = {}\nhint_com_Q = {}\nchal = {}\n",
                fp2_line("A_com", &s.a_com.encode()),
                scalar(a),
                scalar(s.b),
                scalar(s.c_or_d),
                q_hex,
                s.hint_com_p,
                s.hint_com_q,
                chal_hex
            )
        };
        sigs.push_str(&record(s.a, &hex_be(&qb), &hex_be(&cb)));
        // negatives, one kind per vector in turn
        match i % 3 {
            0 => bad.push_str(&record(s.a ^ 1, &hex_be(&qb), &hex_be(&cb))),
            1 => {
                let mut q2 = qb.clone();
                q2[0] ^= 4; // still odd, still 3 mod 4
                bad.push_str(&record(s.a, &hex_be(&q2), &hex_be(&cb)));
            }
            _ => {
                let mut c2 = cb.clone();
                c2[0] ^= 1;
                bad.push_str(&record(s.a, &hex_be(&qb), &hex_be(&c2)));
            }
        }
    }
    if let Ok(dir) = std::env::var("COMPACT_VECTORS_OUT") {
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(format!("{dir}/Public_keys_r3lvl1.txt"), pks).unwrap();
        std::fs::write(format!("{dir}/Signatures_r3lvl1.txt"), sigs).unwrap();
        std::fs::write(format!("{dir}/Signatures_r3lvl1_bad.txt"), bad).unwrap();
    }
}
