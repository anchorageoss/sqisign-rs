//! Cross-validation tests for sqisign-verify types, serialization, and hash.
//!
//! Compares Rust output byte-for-byte against the reference implementation.
//! Reference output captured from tools/c-validate/verify_types_cval.

use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::EcCurve;
use sqisign_verify::fp::Fp2;
use sqisign_verify::hash::hash_to_challenge;
use sqisign_verify::params::{Level1, SecurityLevel};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::{PublicKey, Scalar, Signature};

type L1 = Level1;

fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

fn scalar_hex(s: &Scalar<L1>) -> String {
    let mut out = String::new();
    for &d in s.digits() {
        for b in d.to_le_bytes() {
            out.push_str(&format!("{:02x}", b));
        }
    }
    out
}

// --- Section 1: Public key serialization round-trip ---

#[test]
fn test_public_key_serialization() {
    let mut curve = EcCurve::<L1> {
        a: Fp2::from_small(6),
        c: Fp2::one(),
        ..EcCurve::default()
    };
    let (_, canonical_hint) = ec_curve_to_basis_2f_to_hint::<L1>(
        &mut curve,
        L1::F_CHR,
        L1::basis_e0_px_bytes(),
        L1::basis_e0_qx_bytes(),
        L1::p_cofactor_for_2f(),
        L1::p_cofactor_for_2f_bitlength() as usize,
        L1::torsion_even_power(),
    )
    .unwrap();

    let pk = PublicKey::<L1>::new(
        EcCurve {
            a: Fp2::from_small(6),
            c: Fp2::one(),
            ..EcCurve::default()
        },
        canonical_hint,
    );

    let enc = pk.to_bytes();
    let expected_hex = format!(
        "0600000000000000000000000000000000000000000000000000000000000000\
         0000000000000000000000000000000000000000000000000000000000000000\
         {:02x}",
        canonical_hint
    );
    assert_eq!(to_hex(&enc), expected_hex);

    // Round-trip
    let pk2 = PublicKey::<L1>::from_bytes(&enc).unwrap();
    let enc2 = pk2.to_bytes();
    assert_eq!(to_hex(&enc), to_hex(&enc2), "public key round-trip failed");
}

// --- Section 2: Signature serialization ---

#[test]
fn test_signature_serialization() {
    let mut sig = Signature::<L1>::default();
    sig.set_e_aux_a(Fp2::from_small(42));
    sig.set_backtracking(3);
    sig.set_two_resp_length(7);

    Signature::<L1>::scalar_digits_mut(&mut sig.mat_mut()[0][0])[0] = 0x0102030405060708u64;
    Signature::<L1>::scalar_digits_mut(&mut sig.mat_mut()[0][1])[0] = 0x1112131415161718u64;
    Signature::<L1>::scalar_digits_mut(&mut sig.mat_mut()[1][0])[0] = 0x2122232425262728u64;
    Signature::<L1>::scalar_digits_mut(&mut sig.mat_mut()[1][1])[0] = 0x3132333435363738u64;

    Signature::<L1>::scalar_digits_mut(sig.chall_coeff_mut())[0] = 0xAABBCCDDEEFF0011u64;
    sig.set_hint_aux(0xAA);
    sig.set_hint_chall(0xBB);

    let enc = sig.to_bytes();
    assert_eq!(
        to_hex(&enc),
        "2a00000000000000000000000000000000000000000000000000000000000000\
         00000000000000000000000000000000000000000000000000000000000000000307\
         0807060504030201000000000000000018171615141312110000000000000000\
         2827262524232221000000000000000038373635343332310000000000000000\
         1100ffeeddccbbaa0000000000000000aabb"
    );

    // Round-trip
    let sig2 = Signature::<L1>::from_bytes(&enc).unwrap();
    let enc2 = sig2.to_bytes();
    assert_eq!(to_hex(&enc), to_hex(&enc2), "signature round-trip failed");
}

// --- Section 3: hash_to_challenge ---

#[test]
fn test_hash_to_challenge() {
    let pk = PublicKey::<L1>::new(
        EcCurve {
            a: Fp2::from_small(6),
            c: Fp2::one(),
            ..EcCurve::default()
        },
        0,
    );

    let com = EcCurve::<L1> {
        a: Fp2::from_small(10),
        c: Fp2::one(),
        ..EcCurve::default()
    };

    let msg = b"SQIsign cross-validation test";

    // Verify j-invariant encoding matches
    let j1 = pk.curve().j_inv();
    let j2 = com.j_inv();
    assert_eq!(
        to_hex(&j1.encode()),
        "08630400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        to_hex(&j2.encode()),
        "adcdcfaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa010000000000000000000000000000000000000000000000000000000000000000"
    );

    let challenge = hash_to_challenge::<L1>(&pk, &com, msg);
    assert!(
        challenge.is_ok(),
        "hash_to_challenge must succeed for valid inputs"
    );
    let challenge = challenge.unwrap();
    assert_eq!(
        scalar_hex(&challenge),
        "37cddaa14a7e2bb37ca34595a4b67c0100000000000000000000000000000000"
    );
}

// --- Section 4: hash_to_challenge with zero inputs ---

#[test]
fn test_hash_to_challenge_zero_inputs() {
    let pk = PublicKey::<L1>::new(
        EcCurve {
            a: Fp2::from_small(0),
            c: Fp2::one(),
            ..EcCurve::default()
        },
        0,
    );

    let com = EcCurve::<L1> {
        a: Fp2::from_small(0),
        c: Fp2::one(),
        ..EcCurve::default()
    };

    let challenge = hash_to_challenge::<L1>(&pk, &com, b"");
    assert!(
        challenge.is_ok(),
        "hash_to_challenge must succeed for zero inputs"
    );
    let challenge = challenge.unwrap();
    assert_eq!(
        scalar_hex(&challenge),
        "c038bae2407518b728f6155d04476f0000000000000000000000000000000000"
    );
}

// --- Section 5: Signature from_bytes rejects truncated input ---

#[test]
fn test_signature_rejects_short_input() {
    let bytes = [0u8; 100]; // too short for SIGNATURE_BYTES=148
    assert!(Signature::<L1>::from_bytes(&bytes).is_err());
}

// --- Section 6: PublicKey from_bytes rejects truncated input ---

#[test]
fn test_public_key_rejects_short_input() {
    let bytes = [0u8; 10]; // too short for PUBLICKEY_BYTES=65
    assert!(PublicKey::<L1>::from_bytes(&bytes).is_err());
}
