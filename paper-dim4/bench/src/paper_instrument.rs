//! The paper's two instrumented measurements at level I (feature
//! `instrument`; never part of a timing session):
//!
//!     paper_instrument --out <dir> [--extra <n>]
//!
//! * `gram.csv`: the Gram matrix of the response lattice of every
//!   commitment, first over the 100 signatures of the loop table (the same
//!   key, messages and randomness as `paper`'s signing block, so the same
//!   signatures), then over `n` further signatures of the same stream.
//!   `gram_minima.py` computes the successive minima from it.
//! * `hints.csv`: for each of the 100 signatures, the hint pairs of the
//!   commitment curve that give a basis, and how many of the re-encodings
//!   on those bases verify when the verifier does not check the hints (the
//!   canonical pair aside). Every one of them is also given to the real
//!   verifier, which must reject it.

use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign_records,
    compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::level1;
use sqisign_verify::hd::{
    hd_verify_bytes, hd_verify_bytes_without_hint_check, parse_signature, HdReject,
};
use sqisign_verify::P324_3;

const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const SIGNS: u64 = 100;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut out = PathBuf::from(".");
    let mut extra = 0u64;
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out" => {
                i += 1;
                out = PathBuf::from(&argv[i]);
            }
            "--extra" => {
                i += 1;
                extra = argv[i].parse().expect("--extra <n>");
            }
            s => panic!("unknown argument {s}"),
        }
        i += 1;
    }
    fs::create_dir_all(&out).unwrap();
    let params = level1();
    // the stream of `paper`'s compact signing block at level I
    let mut r = ShakeRng::new(b"paper session I sign", b"sgn");
    let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut r).unwrap();
    let mut pkb = vec![0u8; compact_public_key_bytes::<P324_3>()];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let mut sig = vec![0u8; compact_signature_bytes::<P324_3>()];

    let mut gram = File::create(out.join("gram.csv")).unwrap();
    writeln!(gram, "# level I: the Gram matrix G of the response lattice of each commitment (x^T G x is twice the degree of the response x); signatures 0 to 99 are those of the loop table, the rest continue the same stream").unwrap();
    writeln!(gram, "signature,commitment,samples,tests,g00,g01,g02,g03,g10,g11,g12,g13,g20,g21,g22,g23,g30,g31,g32,g33").unwrap();
    let mut hints = File::create(out.join("hints.csv")).unwrap();
    writeln!(hints, "# level I, the 100 signatures of the loop table: hint pairs (hp, hq) in 0..20 x 0..20 of the commitment curve; with_basis = pairs whose points form a basis and give an encoding; valid_other = of those, the non-canonical pairs whose encoding verifies with the hint check off; rejected_by_verifier = of the valid_other, those the real verifier rejects as NonCanonicalHint").unwrap();
    writeln!(
        hints,
        "signature,canonical_hp,canonical_hq,with_basis,valid_other,rejected_by_verifier"
    )
    .unwrap();

    let mut buf = vec![0u8; 4096];
    for s in 0..SIGNS + extra {
        let msg = s.to_le_bytes();
        let want = s < SIGNS;
        let (_, stats, rec) =
            compact_sign_records::<P324_3, N>(&params, &pk, &sk, &msg, &mut r, &mut sig, want)
                .unwrap();
        for (c, g) in rec.gram.iter().enumerate() {
            let cells: Vec<String> =
                g.0.iter()
                    .flat_map(|row| row.iter())
                    .map(|x| x.to_hex(&mut buf).to_string())
                    .collect();
            writeln!(
                gram,
                "{s},{c},{},{},{}",
                stats.response_samples,
                stats.primality_tests,
                cells.join(",")
            )
            .unwrap();
        }
        if want {
            assert_eq!(hd_verify_bytes::<P324_3>(&sig, &pkb, &msg), Ok(()));
            let ps = parse_signature::<P324_3>(&sig).unwrap();
            let canon = (ps.hint_com_p, ps.hint_com_q);
            let mut valid_other = 0;
            let mut rejected = 0;
            let mut canon_seen = false;
            for (hp, hq, bytes) in &rec.reencodings {
                if (*hp, *hq) == canon {
                    assert_eq!(
                        bytes, &sig,
                        "the canonical pair re-encodes to the signature"
                    );
                    canon_seen = true;
                    continue;
                }
                if hd_verify_bytes_without_hint_check::<P324_3>(bytes, &pkb, &msg).is_ok() {
                    valid_other += 1;
                    if hd_verify_bytes::<P324_3>(bytes, &pkb, &msg)
                        == Err(HdReject::NonCanonicalHint)
                    {
                        rejected += 1;
                    }
                }
            }
            assert!(canon_seen);
            writeln!(
                hints,
                "{s},{},{},{},{valid_other},{rejected}",
                canon.0,
                canon.1,
                rec.reencodings.len()
            )
            .unwrap();
            eprintln!(
                "signature {s}: {} samples, {} pairs with a basis, {valid_other} other valid encodings",
                stats.response_samples,
                rec.reencodings.len()
            );
        }
    }
}
