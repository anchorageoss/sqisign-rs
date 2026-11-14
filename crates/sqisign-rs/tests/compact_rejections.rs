//! The compact verifier's rejections before the chain, at levels I, III
//! and V: a degree that is not good (the goodness test repeated before
//! Cornacchia), hints that are not the canonical ones of their curves, a
//! pairing of order below `2^r`, `a` and `b` both even; and the signer's
//! restart on a fresh commitment. Feature `compact`.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{
    compact_keygen, compact_public_key_bytes, compact_public_key_to_bytes, compact_sign,
    compact_sign_stats_with_budget, compact_signature_bytes,
};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::{level1, level3, level5, Params};
use sqisign_verify::hd::{
    canonical_hints, hd_torsion_basis, hd_verify_bytes, hd_verify_parsed, is_good_degree,
    parse_public_key, parse_signature, HdLevel, HdReject, ParsedPublicKey, ParsedSignature, U6, W,
};
use sqisign_verify::{P324_3, P500_27, P664_17};

const N1: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const N3: usize = sqisign_rs::precomp::p500_27::IBZ_NLIMBS;
const N5: usize = sqisign_rs::precomp::p664_17::IBZ_NLIMBS;

fn signed<L: HdLevel, const N: usize>(
    params: &Params<N>,
    label: &[u8],
) -> (ParsedPublicKey<L>, ParsedSignature<L>, Vec<u8>) {
    let mut rng = ShakeRng::new(label, b"rej");
    let (pk, sk) = compact_keygen::<L, N>(params, &mut rng).expect("keygen");
    let mut pkb = vec![0u8; compact_public_key_bytes::<L>()];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let msg = b"rejections".to_vec();
    let mut sig = vec![0u8; compact_signature_bytes::<L>()];
    compact_sign::<L, N>(params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
    assert_eq!(hd_verify_bytes::<L>(&sig, &pkb, &msg), Ok(()));
    (
        parse_public_key::<L>(&pkb).unwrap(),
        parse_signature::<L>(&sig).unwrap(),
        msg,
    )
}

fn rejections<L: HdLevel, const N: usize>(params: &Params<N>, label: &[u8]) {
    let (pk, sig, msg) = signed::<L, N>(params, label);
    assert_eq!(hd_verify_parsed::<L>(&sig, &pk, &msg), Ok(()));
    let e = L::E_EMBED;
    assert!(is_good_degree(e, &sig.q));

    // The degree. q >= 2^e; q even (2^e - q even); q = 1 mod 4 (2^e - q = 3
    // mod 4); q = 3 mod 4 with 2^e - q composite.
    let with_q = |q: U6| {
        let mut s = sig.clone();
        s.q = q;
        hd_verify_parsed::<L>(&s, &pk, &msg)
    };
    assert_eq!(with_q(U6::pow2(e)), Err(HdReject::DegreeNotGood));
    assert_eq!(
        with_q(sig.q.wrapping_add(&U6::ONE)),
        Err(HdReject::DegreeNotGood)
    );
    assert_eq!(
        with_q(sig.q.wrapping_add(&U6::from_u64(2))),
        Err(HdReject::DegreeNotGood)
    );
    let mut composite = sig.q;
    loop {
        composite = composite.wrapping_add(&U6::from_u64(4));
        if !is_good_degree(e, &composite) {
            break;
        }
    }
    assert_eq!(with_q(composite), Err(HdReject::DegreeNotGood));
    // a good degree that is not the signature's passes the test and fails later
    let mut other = sig.q;
    loop {
        other = other.wrapping_add(&U6::from_u64(4));
        if is_good_degree(e, &other) {
            break;
        }
    }
    let verdict = with_q(other);
    assert!(verdict.is_err() && verdict != Err(HdReject::DegreeNotGood));

    // The hints: each of the four moved to another table entry.
    assert_eq!(
        canonical_hints::<L>(&pk.a_pk),
        Some((pk.hint_pk_p, pk.hint_pk_q))
    );
    assert_eq!(
        canonical_hints::<L>(&sig.a_com),
        Some((sig.hint_com_p, sig.hint_com_q))
    );
    for which in 0..4 {
        for delta in [1u32, 2, 7] {
            let (mut p, mut s) = (pk.clone(), sig.clone());
            let h = match which {
                0 => &mut p.hint_pk_p,
                1 => &mut p.hint_pk_q,
                2 => &mut s.hint_com_p,
                _ => &mut s.hint_com_q,
            };
            *h += delta;
            assert_eq!(
                hd_verify_parsed::<L>(&s, &p, &msg),
                Err(HdReject::NonCanonicalHint),
                "hint {which} + {delta}"
            );
        }
    }

    // a and b both even: the rule for the dropped entry is not defined.
    let mut s = sig.clone();
    s.a = s.a & (W::mask(L::R) ^ W::ONE);
    s.b = s.b & (W::mask(L::R) ^ W::ONE);
    assert!(hd_verify_parsed::<L>(&s, &pk, &msg).is_err());
}

#[test]
fn rejections_level1() {
    rejections::<P324_3, N1>(&level1(), b"rejections I");
}

#[test]
fn rejections_level3() {
    rejections::<P500_27, N3>(&level3(), b"rejections III");
}

#[test]
fn rejections_level5() {
    rejections::<P664_17, N5>(&level5(), b"rejections V");
}

/// A hint past the tables whose point is in `2 E`: the "basis" it gives has
/// order below `2^f`, the pairing of its `2^r`-torsion multiples has order
/// below `2^r`, and the response recovery reports it.
#[test]
fn pairing_of_small_order_is_reported() {
    use sqisign_verify::ec::basis::is_on_curve;
    use sqisign_verify::ec::EcCurve;
    use sqisign_verify::hd::{challenge_limbs, NoProfiler, ResponseScalars};
    use sqisign_verify::hd::{recover_challenge, recover_response_prof, ResponseError};
    use sqisign_verify::Fp2;
    type L = P324_3;
    let (pk, sig, msg) = signed::<L, N1>(&level1(), b"pairing order");
    let chal_limbs = challenge_limbs::<L>(&sig, &pk, &msg).unwrap();
    let chal =
        recover_challenge::<L>(&pk.a_pk, pk.hint_pk_p, pk.hint_pk_q, &chal_limbs).expect("chal");
    let scalars = ResponseScalars {
        a: sig.a,
        b: sig.b,
        c_or_d: sig.c_or_d,
        q: sig.q.low_w(),
    };
    // the honest hints: fine
    assert!(recover_response_prof::<L, _>(
        &chal,
        &sig.a_com,
        sig.hint_com_p,
        sig.hint_com_q,
        scalars,
        &mut NoProfiler
    )
    .is_ok());
    // x = h + i past the table, a square (so the point is a double) on the curve
    let curve = EcCurve::<L>::from_a(&sig.a_com).unwrap();
    let i = Fp2::<L>::i_element();
    let h = (20u32..2000)
        .find(|&h| {
            let x = Fp2::<L>::from_small(h as u64).add(&i);
            bool::from(x.is_square()) && bool::from(is_on_curve(&x, &curve))
        })
        .expect("a square x on the curve");
    assert!(hd_torsion_basis::<L>(&sig.a_com, h, sig.hint_com_q).is_some());
    let got = recover_response_prof::<L, _>(
        &chal,
        &sig.a_com,
        h,
        sig.hint_com_q,
        scalars,
        &mut NoProfiler,
    );
    assert_eq!(got.err(), Some(ResponseError::PairingOrder));
}

/// With a budget of two samples per commitment the signer restarts on
/// fresh commitments until a good degree is found, and the signatures
/// verify.
#[test]
fn restart_on_a_fresh_commitment() {
    type L = P324_3;
    let params = level1();
    let mut rng = ShakeRng::new(b"restart", b"rej");
    let (pk, sk) = compact_keygen::<L, N1>(&params, &mut rng).expect("keygen");
    let mut pkb = vec![0u8; compact_public_key_bytes::<L>()];
    compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
    let mut restarts = 0;
    for n in 0..4u8 {
        let msg = [n; 3];
        let mut sig = vec![0u8; compact_signature_bytes::<L>()];
        let (len, stats) =
            compact_sign_stats_with_budget::<L, N1>(&params, &pk, &sk, &msg, &mut rng, &mut sig, 2)
                .expect("sign with restarts");
        assert_eq!(len, sig.len());
        assert_eq!(hd_verify_bytes::<L>(&sig, &pkb, &msg), Ok(()));
        assert!(stats.response_samples <= 2 * (stats.restarts + 1));
        restarts += stats.restarts;
    }
    assert!(
        restarts > 0,
        "a budget of two samples restarted at least once"
    );
}
