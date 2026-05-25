use criterion::{criterion_group, criterion_main, Criterion};
use sqisign_verify::ec::basis::ec_curve_to_basis_2f_from_hint;
use sqisign_verify::hash::hash_to_challenge;
use sqisign_verify::params::{Level1, SecurityLevel};
use sqisign_verify::precomp::level1;
use sqisign_verify::types::{PublicKey, Signature};

type L1 = Level1;

const SIGNATURE_BYTES: usize = 148;

fn hex_to_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

struct KatEntry {
    pk: Vec<u8>,
    sm: Vec<u8>,
}

fn parse_first_kat_entry(content: &str) -> KatEntry {
    let mut pk = None;
    let mut sm = None;

    for line in content.lines() {
        let line = line.trim();
        if let Some(val) = line.strip_prefix("pk = ") {
            pk = Some(hex_to_bytes(val));
        } else if let Some(val) = line.strip_prefix("sm = ") {
            sm = Some(hex_to_bytes(val));
        }
        if pk.is_some() && sm.is_some() {
            break;
        }
    }

    KatEntry {
        pk: pk.unwrap(),
        sm: sm.unwrap(),
    }
}

fn bench_verify_standard(c: &mut Criterion) {
    let content = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
    let entry = parse_first_kat_entry(content);

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig_bytes = &entry.sm[..SIGNATURE_BYTES];
    let msg = &entry.sm[SIGNATURE_BYTES..];

    c.bench_function("verify_standard (Level 1)", |b| {
        b.iter(|| {
            let sig = Signature::<L1>::from_bytes(sig_bytes).unwrap();
            assert!(sig.verify(&pk, msg).is_ok());
        })
    });
}

fn bench_verify_expanded(c: &mut Criterion) {
    let content = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
    let entry = parse_first_kat_entry(content);

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig_bytes = &entry.sm[..SIGNATURE_BYTES];
    let msg = &entry.sm[SIGNATURE_BYTES..];

    let sig = Signature::<L1>::from_bytes(sig_bytes).unwrap();
    let expanded = sig.expand(&pk).expect("expand failed");

    c.bench_function("verify_expanded (Level 1)", |b| {
        b.iter(|| {
            expanded.verify(&pk, msg).expect("expanded verify failed");
        })
    });
}

fn bench_hash_to_challenge(c: &mut Criterion) {
    let content = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
    let entry = parse_first_kat_entry(content);

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let com_curve = pk.curve().clone();
    let msg = &entry.sm[SIGNATURE_BYTES..];

    c.bench_function("hash_to_challenge (Level 1)", |b| {
        b.iter(|| {
            let _chall = hash_to_challenge::<L1>(&pk, &com_curve, msg);
        })
    });
}

fn bench_verify_compressed(c: &mut Criterion) {
    let content = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
    let entry = parse_first_kat_entry(content);

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();
    let sig_bytes = &entry.sm[..SIGNATURE_BYTES];
    let msg = &entry.sm[SIGNATURE_BYTES..];

    let sig = Signature::<L1>::from_bytes(sig_bytes).unwrap();
    let compressed = sig.compress();

    c.bench_function("verify_compressed (Level 1)", |b| {
        b.iter(|| {
            compressed
                .verify(&pk, msg)
                .expect("compressed verify failed");
        })
    });
}

fn bench_basis_from_hint(c: &mut Criterion) {
    let content = include_str!("../../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
    let entry = parse_first_kat_entry(content);

    let pk = PublicKey::<L1>::from_bytes(&entry.pk).unwrap();

    c.bench_function("ec_curve_to_basis_2f_from_hint (Level 1)", |b| {
        b.iter(|| {
            let mut curve = pk.curve().clone();
            let (basis, ok) = ec_curve_to_basis_2f_from_hint(
                &mut curve,
                L1::F_CHR,
                pk.hint_pk(),
                &level1::BASIS_E0_PX_BYTES,
                &level1::BASIS_E0_QX_BYTES,
                level1::P_COFACTOR_FOR_2F,
                level1::P_COFACTOR_FOR_2F_BITLENGTH as usize,
                level1::TORSION_EVEN_POWER,
            )
            .unwrap();
            assert_eq!(ok, 1);
            std::hint::black_box(basis);
        })
    });
}

criterion_group!(
    benches,
    bench_verify_standard,
    bench_verify_expanded,
    bench_verify_compressed,
    bench_hash_to_challenge,
    bench_basis_from_hint,
);
criterion_main!(benches);
