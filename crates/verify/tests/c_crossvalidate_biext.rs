//! Cross-validation tests for the biextension module (pairings and dlogs).
//!
//! Expected hex values were captured from `tools/c-validate/biext_cval`.

use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::pairing::{
    clear_cofac, ec_dlog_2_tate, ec_dlog_2_weil, fp2_frob, reduced_tate, weil,
};
use sqisign_verify::ec::point::{ec_biscalar_mul, ec_dbl_iter, xadd};
use sqisign_verify::ec::{EcBasis, EcCurve};
use sqisign_verify::fp::{Fp, Fp2};
use sqisign_verify::params::Level1;
use sqisign_verify::precomp::level1::*;

type L1 = Level1;

fn fp2_from_small(re: u64, im: u64) -> Fp2<L1> {
    Fp2 {
        re: Fp::<L1>::from_small(re),
        im: Fp::<L1>::from_small(im),
    }
}

fn fp2_hex(v: &Fp2<L1>) -> String {
    let bytes = v.encode();
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

fn digit_array_hex(arr: &[u64]) -> String {
    let mut s = String::new();
    for &w in arr {
        let bytes = w.to_le_bytes();
        for b in &bytes {
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

fn make_e1_basis(e1: &mut EcCurve<L1>) -> EcBasis<L1> {
    let (basis, _hint) = ec_curve_to_basis_2f_to_hint(
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

// ---- Test 1: Weil pairing on E1 (A=6) ----
#[test]
fn test_weil_pairing_e1() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);
    let pq = xadd(&basis.p, &basis.q, &basis.pmq);
    let w = weil(TORSION_EVEN_POWER, &basis.p, &basis.q, &pq, &mut e1);
    assert_eq!(
        fp2_hex(&w),
        "271be5a19d04fff6217381ddcc2d8e6cac2790b9a7fbd476228b6c9e02bbbc01efeb7b4f5b870fe3d09eb900b5cac653abb5b9429f2b1af47f330a42d542e904"
    );
}

// ---- Test 2: Reduced Tate pairing on E1 (A=6) ----
#[test]
fn test_tate_pairing_e1() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);
    let pq = xadd(&basis.p, &basis.q, &basis.pmq);
    let t = reduced_tate(
        TORSION_EVEN_POWER,
        &basis.p,
        &basis.q,
        &pq,
        &mut e1,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    );
    assert_eq!(
        fp2_hex(&t),
        "2dfce0ead4191843183be24317eab42ff5e094aaf88f344ffc069bb84e72a002a81cf83d7a5fd5e215140cb765b18f058663ff39ac35a455d838fa8831455f01"
    );
}

// ---- Test 3: Weil pairing order check ----
#[test]
fn test_weil_pairing_order() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);
    let pq = xadd(&basis.p, &basis.q, &basis.pmq);
    let w = weil(TORSION_EVEN_POWER, &basis.p, &basis.q, &pq, &mut e1);

    // w^(2^248) should be 1
    let mut tmp = w.clone();
    for _ in 0..TORSION_EVEN_POWER {
        tmp = tmp.sqr();
    }
    assert_eq!(
        fp2_hex(&tmp),
        "01000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );

    // w^(2^247) should NOT be 1 (should be -1 = p-1)
    let mut tmp = w;
    for _ in 0..TORSION_EVEN_POWER - 1 {
        tmp = tmp.sqr();
    }
    assert_eq!(
        fp2_hex(&tmp),
        "feffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
}

// ---- Test 4: Tate pairing order check ----
#[test]
fn test_tate_pairing_order() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);
    let pq = xadd(&basis.p, &basis.q, &basis.pmq);
    let t = reduced_tate(
        TORSION_EVEN_POWER,
        &basis.p,
        &basis.q,
        &pq,
        &mut e1,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    );

    let mut tmp = t.clone();
    for _ in 0..TORSION_EVEN_POWER {
        tmp = tmp.sqr();
    }
    assert_eq!(
        fp2_hex(&tmp),
        "01000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    );

    let mut tmp = t;
    for _ in 0..TORSION_EVEN_POWER - 1 {
        tmp = tmp.sqr();
    }
    assert_eq!(
        fp2_hex(&tmp),
        "feffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000"
    );
}

// ---- Test 5: Weil bilinearity ----
#[test]
fn test_weil_bilinearity() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);
    // Use P-Q as difference (instead of P+Q)
    let w_neg = weil(TORSION_EVEN_POWER, &basis.p, &basis.q, &basis.pmq, &mut e1);
    let w_neg_inv = w_neg.inv();
    assert_eq!(
        fp2_hex(&w_neg_inv),
        "271be5a19d04fff6217381ddcc2d8e6cac2790b9a7fbd476228b6c9e02bbbc01efeb7b4f5b870fe3d09eb900b5cac653abb5b9429f2b1af47f330a42d542e904"
    );
}

// ---- Test 6: Dlog with Weil pairing ----
#[test]
fn test_weil_dlog() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);

    let known_r1: [u64; 4] = [7, 0, 0, 0];
    let known_r2: [u64; 4] = [12, 0, 0, 0];
    let known_s1: [u64; 4] = [3, 0, 0, 0];
    let known_s2: [u64; 4] = [11, 0, 0, 0];
    let diff_d1: [u64; 4] = [4, 0, 0, 0];
    let diff_d2: [u64; 4] = [1, 0, 0, 0];

    let brs_p = ec_biscalar_mul(
        &known_r1,
        &known_r2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_q = ec_biscalar_mul(
        &known_s1,
        &known_s2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_pmq =
        ec_biscalar_mul(&diff_d1, &diff_d2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();
    let brs = EcBasis::new(brs_p, brs_q, brs_pmq);

    let mut r1 = [0u64; 4];
    let mut r2 = [0u64; 4];
    let mut s1 = [0u64; 4];
    let mut s2 = [0u64; 4];

    ec_dlog_2_weil(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &basis,
        &brs,
        &mut e1,
        TORSION_EVEN_POWER,
    )
    .unwrap();

    assert_eq!(
        digit_array_hex(&r1),
        "f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&r2),
        "f4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&s1),
        "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&s2),
        "f5ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
}

// ---- Test 7: Dlog with Tate pairing ----
#[test]
fn test_tate_dlog() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);

    let known_r1: [u64; 4] = [7, 0, 0, 0];
    let known_r2: [u64; 4] = [12, 0, 0, 0];
    let known_s1: [u64; 4] = [3, 0, 0, 0];
    let known_s2: [u64; 4] = [11, 0, 0, 0];
    let diff_d1: [u64; 4] = [4, 0, 0, 0];
    let diff_d2: [u64; 4] = [1, 0, 0, 0];

    let brs_p = ec_biscalar_mul(
        &known_r1,
        &known_r2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_q = ec_biscalar_mul(
        &known_s1,
        &known_s2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_pmq =
        ec_biscalar_mul(&diff_d1, &diff_d2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();
    let brs = EcBasis::new(brs_p, brs_q, brs_pmq);

    let mut r1 = [0u64; 4];
    let mut r2 = [0u64; 4];
    let mut s1 = [0u64; 4];
    let mut s2 = [0u64; 4];

    ec_dlog_2_tate(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &basis,
        &brs,
        &mut e1,
        TORSION_EVEN_POWER,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    )
    .unwrap();

    assert_eq!(
        digit_array_hex(&r1),
        "f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&r2),
        "f4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&s1),
        "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
    assert_eq!(
        digit_array_hex(&s2),
        "f5ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff00"
    );
}

// ---- Test 8: Tate partial dlog (e=126) ----
#[test]
fn test_tate_partial_dlog() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);

    let known_r1: [u64; 4] = [7, 0, 0, 0];
    let known_r2: [u64; 4] = [12, 0, 0, 0];
    let known_s1: [u64; 4] = [3, 0, 0, 0];
    let known_s2: [u64; 4] = [11, 0, 0, 0];
    let diff_d1: [u64; 4] = [4, 0, 0, 0];
    let diff_d2: [u64; 4] = [1, 0, 0, 0];

    let brs_p = ec_biscalar_mul(
        &known_r1,
        &known_r2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_q = ec_biscalar_mul(
        &known_s1,
        &known_s2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_pmq =
        ec_biscalar_mul(&diff_d1, &diff_d2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();

    let e_diff = TORSION_EVEN_POWER - 126;
    let brs_partial = EcBasis::new(
        ec_dbl_iter(&brs_p, e_diff as usize, &mut e1),
        ec_dbl_iter(&brs_q, e_diff as usize, &mut e1),
        ec_dbl_iter(&brs_pmq, e_diff as usize, &mut e1),
    );

    let mut r1 = [0u64; 4];
    let mut r2 = [0u64; 4];
    let mut s1 = [0u64; 4];
    let mut s2 = [0u64; 4];

    ec_dlog_2_tate(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &basis,
        &brs_partial,
        &mut e1,
        126,
        TORSION_EVEN_POWER,
        P_COFACTOR_FOR_2F,
    )
    .unwrap();

    assert_eq!(
        digit_array_hex(&r1),
        "0700000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        digit_array_hex(&r2),
        "0c00000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        digit_array_hex(&s1),
        "0300000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        digit_array_hex(&s2),
        "0b00000000000000000000000000000000000000000000000000000000000000"
    );
}

// ---- Test 9: clear_cofac ----
#[test]
fn test_clear_cofac() {
    let input = fp2_from_small(3, 7);
    let result = clear_cofac::<L1>(&input, P_COFACTOR_FOR_2F);
    assert_eq!(
        fp2_hex(&result),
        "f45900000000000000000000000000000000000000000000000000000000000023d4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04"
    );
}

// ---- Test 10: fp2_frob ----
#[test]
fn test_fp2_frob() {
    let input = fp2_from_small(5, 11);
    let result = fp2_frob::<L1>(&input);
    assert_eq!(
        fp2_hex(&result),
        "0500000000000000000000000000000000000000000000000000000000000000f4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04"
    );
}

// ---- Ported quality check: Weil dlog round-trip ----
#[test]
fn test_weil_dlog_roundtrip() {
    let mut e1 = make_e1();
    let basis = make_e1_basis(&mut e1);

    let known_r1: [u64; 4] = [7, 0, 0, 0];
    let known_r2: [u64; 4] = [12, 0, 0, 0];
    let known_s1: [u64; 4] = [3, 0, 0, 0];
    let known_s2: [u64; 4] = [11, 0, 0, 0];
    let diff_d1: [u64; 4] = [4, 0, 0, 0];
    let diff_d2: [u64; 4] = [1, 0, 0, 0];

    let brs_p = ec_biscalar_mul(
        &known_r1,
        &known_r2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_q = ec_biscalar_mul(
        &known_s1,
        &known_s2,
        TORSION_EVEN_POWER as usize,
        &basis,
        &e1,
    )
    .unwrap();
    let brs_pmq =
        ec_biscalar_mul(&diff_d1, &diff_d2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();
    let brs = EcBasis::new(brs_p, brs_q, brs_pmq);

    let mut r1 = [0u64; 4];
    let mut r2 = [0u64; 4];
    let mut s1 = [0u64; 4];
    let mut s2 = [0u64; 4];

    ec_dlog_2_weil(
        &mut r1,
        &mut r2,
        &mut s1,
        &mut s2,
        &basis,
        &brs,
        &mut e1,
        TORSION_EVEN_POWER,
    )
    .unwrap();

    // Verify: [r1]P + [r2]Q should equal R
    let check_r = ec_biscalar_mul(&r1, &r2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();
    assert!(
        bool::from(check_r.ct_equal(&brs.p)),
        "R = [r1]P + [r2]Q must hold"
    );

    // Verify: [s1]P + [s2]Q should equal S
    let check_s = ec_biscalar_mul(&s1, &s2, TORSION_EVEN_POWER as usize, &basis, &e1).unwrap();
    assert!(
        bool::from(check_s.ct_equal(&brs.q)),
        "S = [s1]P + [s2]Q must hold"
    );
}
