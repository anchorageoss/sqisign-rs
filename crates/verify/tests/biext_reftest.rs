//! Biextension pairing and discrete log tests (seed=0).
//!
//! Replicates a fixed test sequence using hardcoded DRBG-derived scalars.
//! All intermediate values were captured from `tools/c-validate/biext_reftest_trace`.

use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::pairing::{ec_dlog_2_tate, ec_dlog_2_weil, reduced_tate, weil};
use sqisign_verify::ec::point::{ec_biscalar_mul, ec_dbl_iter, xadd, xdbl_a24};
use sqisign_verify::ec::{EcBasis, EcCurve};
use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level1;
use sqisign_verify::precomp::level1::*;

type L1 = Level1;

fn fp2_hex(v: &Fp2<L1>) -> String {
    v.encode().iter().map(|b| format!("{:02x}", b)).collect()
}

fn digit_hex(arr: &[u64]) -> String {
    let mut s = String::new();
    for &w in arr {
        for b in w.to_le_bytes() {
            s.push_str(&format!("{:02x}", b));
        }
    }
    s
}

fn make_e1() -> EcCurve<L1> {
    let mut e1 = EcCurve::<L1> {
        a: Fp2::from_small(6),
        c: Fp2::one(),
        ..Default::default()
    };
    e1.normalize_a24();
    e1
}

fn make_basis(e1: &mut EcCurve<L1>) -> EcBasis<L1> {
    let (basis, _) = ec_curve_to_basis_2f_to_hint(
        e1,
        TORSION_EVEN_POWER,
        &BASIS_E0_PX_BYTES,
        &BASIS_E0_QX_BYTES,
        P_COFACTOR_FOR_2F,
        P_COFACTOR_FOR_2F_BITLENGTH as usize,
        TORSION_EVEN_POWER,
    )
    .unwrap();
    basis
}

/// Biextension pairing and dlog test with seed=0 DRBG trace.
#[test]
fn biextension_reftest_seed0() {
    let e = TORSION_EVEN_POWER;
    let mut curve = make_e1();
    let even_torsion = make_basis(&mut curve);

    // ---- Verify basis point order ----
    let tmp = ec_dbl_iter(&even_torsion.p, e as usize, &mut curve);
    assert!(bool::from(tmp.is_zero()), "P must have order 2^e");
    let tmp = ec_dbl_iter(&even_torsion.q, e as usize, &mut curve);
    assert!(bool::from(tmp.is_zero()), "Q must have order 2^e");
    let tmp = ec_dbl_iter(&even_torsion.pmq, e as usize, &mut curve);
    assert!(bool::from(tmp.is_zero()), "PmQ must have order 2^e");

    // ---- Weil pairing ----
    let pq = xadd(&even_torsion.p, &even_torsion.q, &even_torsion.pmq);
    let weil_r = weil(e, &even_torsion.p, &even_torsion.q, &pq, &mut curve);
    assert_eq!(
        fp2_hex(&weil_r),
        "271be5a19d04fff6217381ddcc2d8e6cac2790b9a7fbd476228b6c9e02bbbc01efeb7b4f5b870fe3d09eb900b5cac653abb5b9429f2b1af47f330a42d542e904",
        "weil pairing must match expected value"
    );

    // ---- Tate pairing ----
    let tate_r = reduced_tate(
        e,
        &even_torsion.p,
        &even_torsion.q,
        &pq,
        &mut curve,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    );
    assert_eq!(
        fp2_hex(&tate_r),
        "2dfce0ead4191843183be24317eab42ff5e094aaf88f344ffc069bb84e72a002a81cf83d7a5fd5e215140cb765b18f058663ff39ac35a455d838fa8831455f01",
        "tate pairing must match expected value"
    );

    // ---- Weil pairing order check ----
    let mut tmp_fp2 = weil_r.clone();
    for _ in 0..e {
        tmp_fp2 = tmp_fp2.sqr();
    }
    assert!(bool::from(tmp_fp2.ct_is_one()), "weil^(2^e) must be 1");

    let mut tmp_fp2 = weil_r.clone();
    for _ in 0..e - 1 {
        tmp_fp2 = tmp_fp2.sqr();
    }
    assert!(
        !bool::from(tmp_fp2.ct_is_one()),
        "weil^(2^(e-1)) must NOT be 1"
    );

    // ---- Tate pairing order check ----
    let mut tmp_fp2 = tate_r.clone();
    for _ in 0..e {
        tmp_fp2 = tmp_fp2.sqr();
    }
    assert!(bool::from(tmp_fp2.ct_is_one()), "tate^(2^e) must be 1");

    let mut tmp_fp2 = tate_r.clone();
    for _ in 0..e - 1 {
        tmp_fp2 = tmp_fp2.sqr();
    }
    assert!(
        !bool::from(tmp_fp2.ct_is_one()),
        "tate^(2^(e-1)) must NOT be 1"
    );

    // ---- Bilinearity: e(P, Q, P-Q) = 1/e(P, Q, P+Q) ----
    let weil_r2 = weil(
        e,
        &even_torsion.p,
        &even_torsion.q,
        &even_torsion.pmq,
        &mut curve,
    );
    let weil_r2_inv = weil_r2.inv();
    assert!(
        bool::from(weil_r.ct_equal(&weil_r2_inv)),
        "bilinearity: e(P,Q,P-Q)^-1 == e(P,Q,P+Q)"
    );

    // ---- Double-bilinearity: e(2P, Q) == e(P, 2Q) == e(P,Q)^2 ----
    let a24 = curve.a24.clone();
    let pp = xdbl_a24(&even_torsion.p, &a24, false);
    let qq = xdbl_a24(&even_torsion.q, &a24, false);
    let ppq = xadd(&pq, &even_torsion.p, &even_torsion.q);
    let pqq = xadd(&pq, &even_torsion.q, &even_torsion.p);
    let w2 = weil(e, &pp, &even_torsion.q, &ppq, &mut curve);
    let w3 = weil(e, &even_torsion.p, &qq, &pqq, &mut curve);
    assert!(bool::from(w2.ct_equal(&w3)), "e(2P,Q) must equal e(P,2Q)");
    let rr1 = weil_r.sqr();
    assert!(bool::from(rr1.ct_equal(&w2)), "e(2P,Q) must equal e(P,Q)^2");

    // ---- dlog tests: replicate seed=0 DRBG scalars ----
    // Exact scalars produced by CTR-DRBG(seed=0) after the standard
    // fixup logic (mask to torsion order, set low bit).
    let scal_d1: [u64; 4] = [
        0x20948f9ae98f6190,
        0xa0275b736f247b49,
        0xa0b2a63c9d8a0719,
        0,
    ];
    let scal_d2: [u64; 4] = [
        0xf7fd3eba_ac326779,
        0xf5ddee2a_247cbf31,
        0x6ae390da_31b1a5eb,
        0,
    ];
    let scal_s1: [u64; 4] = [
        0x4d8b1e3b_4928fdd5,
        0xc6704f65_f315ee77,
        0x9d19d1bd_7354a64a,
        0,
    ];
    let scal_s2: [u64; 4] = [
        0x549544ed_19fdad5b,
        0x85d87915_e46b6388,
        0xe9c877f2_abc72ac1,
        0,
    ];

    // r1 = d1 + s1, r2 = d2 + s2
    let mut scal_r1 = [0u64; 4];
    let mut scal_r2 = [0u64; 4];
    mp_add_scalar(&mut scal_r1, &scal_d1, &scal_s1);
    mp_add_scalar(&mut scal_r2, &scal_d2, &scal_s2);

    assert_eq!(
        digit_hex(&scal_r1),
        "655fb832d6ad1f6ec0693a62d9aa976664adde10fa77cc3d0100000000000000"
    );
    assert_eq!(
        digit_hex(&scal_r2),
        "d41430c6a783924cba22e8084067b67badd078ddcc08ac540100000000000000"
    );

    let bpq = even_torsion.clone();

    let brs_p = ec_biscalar_mul(&scal_r1, &scal_r2, e as usize, &bpq, &curve).unwrap();
    let brs_q = ec_biscalar_mul(&scal_s1, &scal_s2, e as usize, &bpq, &curve).unwrap();
    let brs_pmq = ec_biscalar_mul(&scal_d1, &scal_d2, e as usize, &bpq, &curve).unwrap();
    let brs = EcBasis::new(brs_p, brs_q, brs_pmq);

    // Verify BRS point coordinates match expected values
    assert_eq!(
        fp2_hex(&brs.p.x),
        "8729fedc91ab96e73d6967b1344aa850d02c31a8bebab54c7d0de267c43c7300573c66de2b799d424bc2ce7e0600bc1d5afa01ecb6361562372f6fe1d116e603"
    );

    // ---- Weil dlog ----
    let mut r1 = [0u64; 4];
    let mut r2 = [0u64; 4];
    let mut s1 = [0u64; 4];
    let mut s2 = [0u64; 4];

    ec_dlog_2_weil(
        &mut r1, &mut r2, &mut s1, &mut s2, &bpq, &brs, &mut curve, e,
    )
    .unwrap();

    assert_eq!(
        digit_hex(&r1),
        "655fb832d6ad1f6ec0693a62d9aa976664adde10fa77cc3d0100000000000000"
    );
    assert_eq!(
        digit_hex(&r2),
        "d41430c6a783924cba22e8084067b67badd078ddcc08ac540100000000000000"
    );
    assert_eq!(
        digit_hex(&s1),
        "d5fd28493b1e8b4d77ee15f3654f70c64aa65473bdd1199d0000000000000000"
    );
    assert_eq!(
        digit_hex(&s2),
        "5badfd19ed44955488636be41579d885c12ac7abf277c8e90000000000000000"
    );

    // Verify: [r1]P + [r2]Q == R
    let check = ec_biscalar_mul(&r1, &r2, e as usize, &bpq, &curve).unwrap();
    assert!(
        bool::from(check.ct_equal(&brs.p)),
        "weil dlog: R = [r1]P + [r2]Q"
    );
    let check = ec_biscalar_mul(&s1, &s2, e as usize, &bpq, &curve).unwrap();
    assert!(
        bool::from(check.ct_equal(&brs.q)),
        "weil dlog: S = [s1]P + [s2]Q"
    );

    // ---- Tate dlog ----
    ec_dlog_2_tate(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &bpq,
        &brs,
        &mut curve,
        e,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    )
    .unwrap();

    assert_eq!(
        digit_hex(&r1),
        "655fb832d6ad1f6ec0693a62d9aa976664adde10fa77cc3d0100000000000000"
    );
    assert_eq!(
        digit_hex(&r2),
        "d41430c6a783924cba22e8084067b67badd078ddcc08ac540100000000000000"
    );
    assert_eq!(
        digit_hex(&s1),
        "d5fd28493b1e8b4d77ee15f3654f70c64aa65473bdd1199d0000000000000000"
    );
    assert_eq!(
        digit_hex(&s2),
        "5badfd19ed44955488636be41579d885c12ac7abf277c8e90000000000000000"
    );

    let check = ec_biscalar_mul(&r1, &r2, e as usize, &bpq, &curve).unwrap();
    assert!(
        bool::from(check.ct_equal(&brs.p)),
        "tate dlog: R = [r1]P + [r2]Q"
    );
    let check = ec_biscalar_mul(&s1, &s2, e as usize, &bpq, &curve).unwrap();
    assert!(
        bool::from(check.ct_equal(&brs.q)),
        "tate dlog: S = [s1]P + [s2]Q"
    );

    // ---- Tate partial dlog (e=126) ----
    let e_partial = 126u32;
    let e_diff = (TORSION_EVEN_POWER - e_partial) as usize;

    let brs_partial = EcBasis::new(
        ec_dbl_iter(&brs.p, e_diff, &mut curve),
        ec_dbl_iter(&brs.q, e_diff, &mut curve),
        ec_dbl_iter(&brs.pmq, e_diff, &mut curve),
    );

    ec_dlog_2_tate(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &bpq,
        &brs_partial,
        &mut curve,
        e_partial,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    )
    .unwrap();

    assert_eq!(
        digit_hex(&r1),
        "655fb832d6ad1f6ec0693a62d9aa972600000000000000000000000000000000"
    );
    assert_eq!(
        digit_hex(&r2),
        "d41430c6a783924cba22e8084067b63b00000000000000000000000000000000"
    );
    assert_eq!(
        digit_hex(&s1),
        "d5fd28493b1e8b4d77ee15f3654f700600000000000000000000000000000000"
    );
    assert_eq!(
        digit_hex(&s2),
        "5badfd19ed44955488636be41579d80500000000000000000000000000000000"
    );

    // Verify partial dlog round-trip
    let check = ec_biscalar_mul(&r1, &r2, e as usize, &bpq, &curve).unwrap();
    let check = ec_dbl_iter(&check, e_diff, &mut curve);
    assert!(
        bool::from(check.ct_equal(&brs_partial.p)),
        "partial tate: R"
    );
    let check = ec_biscalar_mul(&s1, &s2, e as usize, &bpq, &curve).unwrap();
    let check = ec_dbl_iter(&check, e_diff, &mut curve);
    assert!(
        bool::from(check.ct_equal(&brs_partial.q)),
        "partial tate: S"
    );
}

fn mp_add_scalar(c: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]) {
    let mut carry = 0u64;
    for i in 0..4 {
        let (s1, c1) = a[i].overflowing_add(b[i]);
        let (s2, c2) = s1.overflowing_add(carry);
        c[i] = s2;
        carry = c1 as u64 + c2 as u64;
    }
}
