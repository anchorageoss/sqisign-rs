//! De-risking the 108-byte tightening: confirm that the transmitted basis
//! hints (`hint_com_P/Q`, `hint_pk_P/Q`) are exactly the canonical hints
//! recomputed from the curve coefficient alone. If this holds for every
//! oracle vector, the hints are redundant on the wire and can be dropped.

mod hd_common;
use hd_common::{load, parse_fp2, PHASE0_VECTORS};
use sqisign_verify::hd::canonical_hints_l1;

#[test]
fn transmitted_hints_match_recomputed_canonical_hints() {
    let doc = load(PHASE0_VECTORS);
    let vectors = doc["test_vectors"].as_array().unwrap();
    assert_eq!(vectors.len(), 5);

    for (i, v) in vectors.iter().enumerate() {
        // Commitment curve hints.
        let a_com = parse_fp2(&v["signature"]["A_com"]);
        let (hp, hq) = canonical_hints_l1(&a_com).expect("A_com must be a valid curve");
        let want_p = v["signature"]["hint_com_P"].as_u64().unwrap() as u32;
        let want_q = v["signature"]["hint_com_Q"].as_u64().unwrap() as u32;
        assert_eq!(
            (hp, hq),
            (want_p, want_q),
            "vector {i}: recomputed com hints {:?} != transmitted {:?}",
            (hp, hq),
            (want_p, want_q)
        );

        // Public-key curve hints.
        let a_pk = parse_fp2(&v["public_key"]["A_pk"]);
        let (hp, hq) = canonical_hints_l1(&a_pk).expect("A_pk must be a valid curve");
        let want_p = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let want_q = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        assert_eq!(
            (hp, hq),
            (want_p, want_q),
            "vector {i}: recomputed pk hints {:?} != transmitted {:?}",
            (hp, hq),
            (want_p, want_q)
        );
    }
}
