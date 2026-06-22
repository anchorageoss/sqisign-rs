//! Phase 6: signature / public-key wire-format parsing.
//!
//! 1. Round-trip: encode each reference vector's structured fields to wire
//!    bytes, parse them back, and confirm they match (scalars up to the
//!    canonical mod-`2^r` reduction). The reference has no binary format, so the
//!    bytes are produced by our encoder from the same structured fields the
//!    Phase-5b tests use.
//! 2. End-to-end from bytes: `hd_verify_bytes_l1` accepts all 5 valid
//!    signatures and rejects tampering of the signature, public key, or message.
//! 3. Strict deserialization: trailing bytes, truncation, out-of-range field
//!    elements, and non-canonical scalars all return `Err` (never panic).

mod hd_common;
use hd_common::{load, parse_fp2, PHASE0_VECTORS, F};

use crypto_bigint::U256;
use serde_json::Value;
use sqisign_verify::hd::{
    encode_public_key, encode_signature, hd_verify_bytes_l1, hd_verify_l1, parse_public_key,
    parse_signature, HdReject, HdSignatureL1, PK_WIRE_BYTES, SIG_WIRE_BYTES,
};

const MSG: [u8; 32] = [0u8; 32];
const R_BITS: u32 = 70;

fn dec_u256(s: &str) -> U256 {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    U256::from_words(limbs)
}
fn dec_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}
fn chal_limbs_of(sig: &Value) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for ch in sig["chal"].as_str().unwrap().trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    limbs
}
fn limbs_to_le_bytes(limbs: &[u64; 4]) -> [u8; 32] {
    let mut b = [0u8; 32];
    for (i, &l) in limbs.iter().enumerate() {
        b[i * 8..i * 8 + 8].copy_from_slice(&l.to_le_bytes());
    }
    b
}

/// Structured fields for one vector (mirrors the Phase-5b test inputs).
struct Fields {
    a_pk: F,
    a_com: F,
    a: i128,
    b: i128,
    c_or_d: i128,
    q: U256,
    hint_com_p: u32,
    hint_com_q: u32,
    hint_pk_p: u32,
    hint_pk_q: u32,
    chal_limbs: [u64; 4],
}
fn fields_of(v: &Value) -> Fields {
    let sig = &v["signature"];
    Fields {
        a_pk: parse_fp2(&v["public_key"]["A_pk"]),
        a_com: parse_fp2(&sig["A_com"]),
        a: dec_i128(&sig["a"]),
        b: dec_i128(&sig["b"]),
        c_or_d: dec_i128(&sig["c_or_d"]),
        q: dec_u256(sig["q"].as_str().unwrap()),
        hint_com_p: sig["hint_com_P"].as_u64().unwrap() as u32,
        hint_com_q: sig["hint_com_Q"].as_u64().unwrap() as u32,
        hint_pk_p: v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32,
        hint_pk_q: v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32,
        chal_limbs: chal_limbs_of(sig),
    }
}
fn sig_wire(f: &Fields) -> [u8; SIG_WIRE_BYTES] {
    encode_signature(
        &f.a_com,
        f.a,
        f.b,
        f.c_or_d,
        &f.q,
        f.hint_com_p,
        f.hint_com_q,
    )
    .unwrap()
}
fn pk_wire(f: &Fields) -> [u8; PK_WIRE_BYTES] {
    encode_public_key(&f.a_pk, f.hint_pk_p, f.hint_pk_q).unwrap()
}
fn red(x: i128) -> i128 {
    x.rem_euclid(1i128 << R_BITS)
}

#[test]
fn wire_sizes() {
    assert_eq!(SIG_WIRE_BYTES, 108);
    assert_eq!(PK_WIRE_BYTES, 64);
}

#[test]
fn roundtrip_matches_structured_fields() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let f = fields_of(v);

        let ps = parse_signature(&sig_wire(&f)).expect("parse sig");
        assert!(bool::from(ps.a_com.ct_equal(&f.a_com)), "vec {vi}: A_com");
        assert_eq!(ps.q, f.q, "vec {vi}: q");
        assert_eq!(ps.a, red(f.a), "vec {vi}: a (mod 2^r)");
        assert_eq!(ps.b, red(f.b), "vec {vi}: b (mod 2^r)");
        assert_eq!(ps.c_or_d, red(f.c_or_d), "vec {vi}: c_or_d (mod 2^r)");
        assert_eq!(ps.hint_com_p, f.hint_com_p, "vec {vi}: hint_com_P");
        assert_eq!(ps.hint_com_q, f.hint_com_q, "vec {vi}: hint_com_Q");

        let pp = parse_public_key(&pk_wire(&f)).expect("parse pk");
        assert!(bool::from(pp.a_pk.ct_equal(&f.a_pk)), "vec {vi}: A_pk");
        assert_eq!(pp.hint_pk_p, f.hint_pk_p, "vec {vi}: hint_pk_P");
        assert_eq!(pp.hint_pk_q, f.hint_pk_q, "vec {vi}: hint_pk_Q");
        n += 1;
    }
    assert_eq!(n, 5);
    println!("wire round-trip: parsed fields match the structured fields for all {n} vectors");
}

#[test]
fn bytes_path_matches_structured_path() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let f = fields_of(v);

        // Structured path (recorded chal).
        let chal_bytes = limbs_to_le_bytes(&f.chal_limbs);
        let structured = hd_verify_l1(&HdSignatureL1 {
            a_pk: f.a_pk.clone(),
            a_com: f.a_com.clone(),
            hint_pk_p: f.hint_pk_p,
            hint_pk_q: f.hint_pk_q,
            hint_com_p: f.hint_com_p,
            hint_com_q: f.hint_com_q,
            message: &MSG,
            chal_limbs: &f.chal_limbs,
            claimed_chal: &chal_bytes,
            resp_a: f.a,
            resp_b: f.b,
            resp_c_or_d: f.c_or_d,
            q: f.q,
        });

        // Bytes path (recomputed chal).
        let from_bytes = hd_verify_bytes_l1(&sig_wire(&f), &pk_wire(&f), &MSG);

        assert_eq!(structured, Ok(()), "vec {vi}: structured path must accept");
        assert_eq!(from_bytes, Ok(()), "vec {vi}: bytes path must accept");
        assert_eq!(structured, from_bytes, "vec {vi}: paths must agree");
        n += 1;
    }
    assert_eq!(n, 5);
    println!("bytes path == structured path (both accept) for all {n} vectors");
}

#[test]
fn bytes_path_rejects_tampering() {
    let doc = load(PHASE0_VECTORS);
    let mut n = 0;
    for v in doc["test_vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();
        let f = fields_of(v);
        let sig = sig_wire(&f);
        let pk = pk_wire(&f);

        assert!(hd_verify_bytes_l1(&sig, &pk, &MSG).is_ok(), "vec {vi}: valid must accept");

        // Tamper a response-scalar byte (low byte of `a`): wrong response.
        let mut bad_sig = sig;
        bad_sig[F2_OFFSET_A] ^= 1;
        assert!(
            hd_verify_bytes_l1(&bad_sig, &pk, &MSG).is_err(),
            "vec {vi}: tampered signature scalar must reject"
        );

        // Tamper the commitment curve A_com (byte 0): wrong curve / challenge.
        let mut bad_sig2 = sig;
        bad_sig2[0] ^= 1;
        assert!(
            hd_verify_bytes_l1(&bad_sig2, &pk, &MSG).is_err(),
            "vec {vi}: tampered A_com must reject"
        );

        // Tamper the public key A_pk (byte 0).
        let mut bad_pk = pk;
        bad_pk[0] ^= 1;
        assert!(
            hd_verify_bytes_l1(&sig, &bad_pk, &MSG).is_err(),
            "vec {vi}: tampered public key must reject"
        );

        // Tamper a *packed hint* bit (bit 3 of A_pk's top re-byte). This leaves
        // the curve unchanged (the canonical low bits are untouched) but selects
        // a different torsion basis, so verification must still reject.
        let mut bad_pk_hint = pk;
        bad_pk_hint[31] ^= 0x08;
        assert!(
            hd_verify_bytes_l1(&sig, &bad_pk_hint, &MSG).is_err(),
            "vec {vi}: tampered packed pk hint must reject"
        );

        // Tamper the message.
        let mut bad_msg = MSG;
        bad_msg[0] ^= 1;
        assert!(
            hd_verify_bytes_l1(&sig, &pk, &bad_msg).is_err(),
            "vec {vi}: tampered message must reject"
        );
        n += 1;
    }
    assert_eq!(n, 5);
    println!("bytes path rejects tampering (signature / A_com / pk / hint / message) for all {n} vectors");
}

/// Offset of the first response scalar `a` (after A_com ‖ q).
const F2_OFFSET_A: usize = 64 + 17;

#[test]
fn strict_deserialization() {
    let doc = load(PHASE0_VECTORS);
    let f = fields_of(&doc["test_vectors"][0]);
    let sig = sig_wire(&f);
    let pk = pk_wire(&f);

    // Trailing byte: too long.
    let mut long = sig.to_vec();
    long.push(0);
    assert!(matches!(parse_signature(&long), Err(HdReject::MalformedInput)));
    assert_eq!(hd_verify_bytes_l1(&long, &pk, &MSG), Err(HdReject::MalformedInput));

    // Truncated: one byte short.
    assert!(matches!(parse_signature(&sig[..SIG_WIRE_BYTES - 1]), Err(HdReject::MalformedInput)));
    assert!(matches!(parse_public_key(&pk[..PK_WIRE_BYTES - 1]), Err(HdReject::MalformedInput)));

    // Empty.
    assert!(matches!(parse_signature(&[]), Err(HdReject::MalformedInput)));
    assert!(matches!(parse_public_key(&[]), Err(HdReject::MalformedInput)));

    // Out-of-range A_com (all 0xFF ≥ p): Fp2::decode must reject.
    let mut bad_fp2 = sig;
    for b in bad_fp2[..64].iter_mut() {
        *b = 0xFF;
    }
    assert!(matches!(parse_signature(&bad_fp2), Err(HdReject::MalformedInput)));

    // Out-of-range A_pk in the public key.
    let mut bad_pk_fp2 = pk;
    for b in bad_pk_fp2[..64].iter_mut() {
        *b = 0xFF;
    }
    assert!(matches!(parse_public_key(&bad_pk_fp2), Err(HdReject::MalformedInput)));

    // Non-canonical scalar: top two bits of the 9th byte of `a` set (value ≥ 2^r).
    let mut bad_scalar = sig;
    bad_scalar[F2_OFFSET_A + 8] |= 0x40;
    assert!(matches!(parse_signature(&bad_scalar), Err(HdReject::MalformedInput)));
    bad_scalar[F2_OFFSET_A + 8] = 0x80;
    assert!(matches!(parse_signature(&bad_scalar), Err(HdReject::MalformedInput)));

    println!("strict deserialization: trailing/truncated/empty/out-of-range/non-canonical all rejected (no panic)");
}

#[test]
fn hint_packing_roundtrips_and_preserves_curve() {
    // Use a real vector's A_com / A_pk as a canonical Fp2 value.
    let doc = load(PHASE0_VECTORS);
    let f = fields_of(&doc["test_vectors"][0]);

    // Every (hP, hQ) in the 5-bit range packs and unpacks losslessly, and the
    // recovered curve coefficient is unchanged (hint bits live above it).
    for &(hp, hq) in &[(0u32, 0u32), (19, 19), (1, 31), (31, 1), (7, 12)] {
        let sig = encode_signature(&f.a_com, f.a, f.b, f.c_or_d, &f.q, hp, hq).unwrap();
        let ps = parse_signature(&sig).expect("parse");
        assert_eq!((ps.hint_com_p, ps.hint_com_q), (hp, hq), "sig hints round-trip");
        assert!(bool::from(ps.a_com.ct_equal(&f.a_com)), "A_com preserved under packing");

        let pk = encode_public_key(&f.a_pk, hp, hq).unwrap();
        let pp = parse_public_key(&pk).expect("parse pk");
        assert_eq!((pp.hint_pk_p, pp.hint_pk_q), (hp, hq), "pk hints round-trip");
        assert!(bool::from(pp.a_pk.ct_equal(&f.a_pk)), "A_pk preserved under packing");
    }

    // A hint that does not fit in 5 bits is rejected by the encoder.
    assert!(encode_signature(&f.a_com, f.a, f.b, f.c_or_d, &f.q, 32, 0).is_none());
    assert!(encode_public_key(&f.a_pk, 0, 32).is_none());
}
