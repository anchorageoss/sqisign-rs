//! Byte-for-byte cross-validation of theta foundation (Groups 1, 2, 3)
//! against the reference output from `tools/c-validate/theta_foundation_cval`.

use sqisign_verify::ec::basis::ec_curve_to_basis_2f_to_hint;
use sqisign_verify::ec::EcCurve;
use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level1;
use sqisign_verify::precomp::level1::*;
use sqisign_verify::theta::basis_change::{
    apply_isomorphism, apply_isomorphism_general, base_change_matrix_multiplication,
    set_base_change_matrix_from_precomp,
};
use sqisign_verify::theta::couple::{double_couple_point, double_couple_point_iter};
use sqisign_verify::theta::theta_structure::{
    double_iter, double_point, hadamard, is_product_theta_point, theta_precomputation,
    to_squared_theta,
};
use sqisign_verify::theta::{
    PrecompBasisChangeMatrix, ThetaCoupleCurve, ThetaCouplePoint, ThetaPoint, ThetaStructure,
};

type L1 = Level1;

fn fp2_hex(v: &Fp2<L1>) -> String {
    v.encode().iter().map(|b| format!("{:02x}", b)).collect()
}

fn fp2_constants() -> [Fp2<L1>; 5] {
    [
        Fp2::zero(),                  // 0
        Fp2::one(),                   // 1
        Fp2::i_element(),             // i
        Fp2::<L1>::one().neg(),       // -1
        Fp2::<L1>::i_element().neg(), // -i
    ]
}

fn precomp_matrix(idx: usize) -> PrecompBasisChangeMatrix {
    PrecompBasisChangeMatrix {
        m: SPLITTING_TRANSFORMS[idx],
    }
}

fn make_e0() -> EcCurve<L1> {
    let mut e = EcCurve::<L1> {
        a: Fp2::from_small(6),
        c: Fp2::one(),
        ..Default::default()
    };
    e.normalize_a24();
    e
}

// --- Section 2: Hadamard ---
#[test]
fn test_hadamard() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let h = hadamard(&p);
    assert_eq!(fp2_hex(&h.x), "1a000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h.y), "f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h.z), "f5ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h.t), "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// Hadamard is an involution up to scaling by 4
#[test]
fn test_hadamard_involution() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let h2 = hadamard(&hadamard(&p));
    // hadamard(hadamard(p)) = 4 * p
    assert_eq!(fp2_hex(&h2.x), "0c000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h2.y), "14000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h2.z), "1c000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&h2.t), "2c000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// --- Section 3: to_squared_theta ---
#[test]
fn test_to_squared_theta() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let sq = to_squared_theta(&p);
    assert_eq!(fp2_hex(&sq.x), "cc000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&sq.y), "a7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&sq.z), "77ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&sq.t), "38000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// --- Section 4: theta_precomputation ---
#[test]
fn test_theta_precomputation() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let mut ts = ThetaStructure::<L1> {
        null_point: p,
        ..Default::default()
    };
    theta_precomputation(&mut ts);
    assert!(ts.precomputation);

    assert_eq!(fp2_hex(&ts.cap_xyz0), "00412500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.cap_yzt0), "003a0a00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.cap_xzt0), "ff4ae8ffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.cap_xyt0), "ffa8f0ffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.xyz0), "69000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.yzt0), "81010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.xzt0), "e7000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&ts.xyt0), "a5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// --- Section 5: double_point ---
#[test]
fn test_double_point() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let mut ts = ThetaStructure::<L1> {
        null_point: p.clone(),
        ..Default::default()
    };
    let dbl = double_point(&p, &mut ts);
    assert_eq!(fp2_hex(&dbl.x), "00e02f35b90100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&dbl.y), "0020a558df0200000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&dbl.z), "00601a7c050400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&dbl.t), "00e004c3510600000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// --- Section 6: double_iter ---
#[test]
fn test_double_iter() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let mut ts = ThetaStructure::<L1> {
        null_point: p.clone(),
        ..Default::default()
    };
    let iter3 = double_iter(&p, &mut ts, 3);
    assert_eq!(fp2_hex(&iter3.x), "b0680556a714076182f298aac97f7f38b4ae0b6663e7f6d395a727cd8ac301030000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&iter3.y), "26595e8fc177b6a12e94a971a57fd408d7cd68fffad6f00b4f1742ab3c9bad010000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&iter3.z), "9c49b7c8dbda65e2da35ba38817f29d9f9ecc59892c6ea4308875c89ee7259000000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&iter3.t), "872a693b10a1c4633379dbc6387fd3793f2b80cbc1a5deb37a6691455222b1020000000000000000000000000000000000000000000000000000000000000000");
}

// double_iter consistency: double_iter(p, 3) == double(double(double(p)))
#[test]
fn test_double_iter_consistency() {
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let mut ts = ThetaStructure::<L1> {
        null_point: p.clone(),
        ..Default::default()
    };
    let iter3 = double_iter(&p, &mut ts, 3);
    let mut manual = double_point(&p, &mut ts);
    manual = double_point(&manual, &mut ts);
    manual = double_point(&manual, &mut ts);
    assert_eq!(fp2_hex(&iter3.x), fp2_hex(&manual.x));
    assert_eq!(fp2_hex(&iter3.y), fp2_hex(&manual.y));
    assert_eq!(fp2_hex(&iter3.z), fp2_hex(&manual.z));
    assert_eq!(fp2_hex(&iter3.t), fp2_hex(&manual.t));
}

// --- Section 7: is_product_theta_point ---
#[test]
fn test_is_product_theta_point() {
    // (3,5,7,11): 3*11=33 != 5*7=35 => not product
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    assert!(!bool::from(is_product_theta_point(&p)));

    // (2,3,4,6): 2*6=12=3*4 => product
    let q = ThetaPoint::<L1> {
        x: Fp2::from_small(2),
        y: Fp2::from_small(3),
        z: Fp2::from_small(4),
        t: Fp2::from_small(6),
    };
    assert!(bool::from(is_product_theta_point(&q)));
}

// --- Section 8: double_couple_point ---
#[test]
fn test_double_couple_point() {
    let mut e0 = make_e0();
    let (basis, _) = ec_curve_to_basis_2f_to_hint(
        &mut e0,
        TORSION_EVEN_POWER,
        &BASIS_E0_PX_BYTES,
        &BASIS_E0_QX_BYTES,
        P_COFACTOR_FOR_2F,
        P_COFACTOR_FOR_2F_BITLENGTH as usize,
        TORSION_EVEN_POWER,
    )
    .unwrap();
    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0,
    };
    let cp = ThetaCouplePoint {
        p1: basis.p.clone(),
        p2: basis.q.clone(),
    };
    let dbl = double_couple_point(&cp, &e12);
    assert_eq!(fp2_hex(&dbl.p1.x), "97ca5495d9ff6aef3dce6e473bcb014200133c9b0be538f099fe13aaee9d2602829e89dc776722e333de9fd60a3eff7c094f6d3964b27598877959b6d979d503");
    assert_eq!(fp2_hex(&dbl.p1.z), "a8da31cacac8d9539c7df0314cb3c0e8491c8879535dec5cbef23d0c38140d03b25be7d19bbdb0be9059285cd060a3ee723b27af93d5a8276535a59aebbac900");
    assert_eq!(fp2_hex(&dbl.p2.x), "22683b2f6645d7ec118704edb575d0f11ca0cebaf853ecffed7a8b04dc75e30267cbd184b0f75e1eec366032354c0d240b4e563bd0df18cb2a3470f151443503");
    assert_eq!(fp2_hex(&dbl.p2.z), "a7d9c24dcde340047cd1336c1e4d725161f5dd2ab0c8a118842af5bce0cd350422a32aa2069318b4d057c0c88981931aeffb4d279a01445e6457afdeaafdec02");
}

// --- Section 9: double_couple_point_iter ---
#[test]
fn test_double_couple_point_iter() {
    let mut e0 = make_e0();
    let (basis, _) = ec_curve_to_basis_2f_to_hint(
        &mut e0,
        TORSION_EVEN_POWER,
        &BASIS_E0_PX_BYTES,
        &BASIS_E0_QX_BYTES,
        P_COFACTOR_FOR_2F,
        P_COFACTOR_FOR_2F_BITLENGTH as usize,
        TORSION_EVEN_POWER,
    )
    .unwrap();
    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0,
    };
    let cp = ThetaCouplePoint {
        p1: basis.p.clone(),
        p2: basis.q.clone(),
    };
    let iter3 = double_couple_point_iter(&cp, 3, &e12);
    assert_eq!(fp2_hex(&iter3.p1.x), "bfeeb01ef0c67f8ec968d395ed833a5e645393405bdd727bacc5917a7deffd030d1b9d8e432ebd34e67c0e53ceadb9efb556d765d06ee844c6310f3238e6b303");
    assert_eq!(fp2_hex(&iter3.p1.z), "aae944d4ca6c08a48a246a4cc6725f203b6c608849d7e09b16c5d2a2ce2587011e37acbf9a8e89c5e251bbe779a0c9d42fde1c844c88122deeecdf7bab29ba00");
    assert_eq!(fp2_hex(&iter3.p2.x), "f49772a602910f9dddb36b4898a89f8d4d6b72628469ab3c55028bd7be932e0166ffbd2810bec0312f8615567239edcb4848e46f959bc116af152dec7c4bb502");
    assert_eq!(fp2_hex(&iter3.p2.z), "b3a66fe1486fcbb312b5d2d8a0b5f2c2dac137396b9079ca1f2be09fe687ca048b89059fd9702c86ad1861d66235e7c33a7ec236eafbbc2552e4ae1a71ebd801");
}

// --- Section 10: matrix multiplication ---
#[test]
fn test_base_change_matrix_multiplication() {
    let fc = fp2_constants();
    let m1 = set_base_change_matrix_from_precomp(&precomp_matrix(0), &fc);
    let m2 = set_base_change_matrix_from_precomp(&precomp_matrix(4), &fc);
    let prod = base_change_matrix_multiplication(&m1, &m2);

    let expected: [[&str; 4]; 4] = [
        [
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
            "00000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000000000000000000",
        ],
        [
            "0000000000000000000000000000000000000000000000000000000000000000fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
            "00000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000000000000000000",
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        ],
        [
            "00000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "02000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        ],
        [
            "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000",
            "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040000000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
            "00000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000000000000000000",
        ],
    ];

    for (i, row) in expected.iter().enumerate() {
        for (j, exp) in row.iter().enumerate() {
            assert_eq!(
                fp2_hex(&prod.m[i][j]),
                *exp,
                "M_prod[{}][{}] mismatch",
                i,
                j
            );
        }
    }
}

// --- Section 11: apply_isomorphism ---
#[test]
fn test_apply_isomorphism() {
    let fc = fp2_constants();
    let m1 = set_base_change_matrix_from_precomp(&precomp_matrix(0), &fc);
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let out = apply_isomorphism(&m1, &p);
    assert_eq!(fp2_hex(&out.x), "0a000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&out.y), "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040600000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&out.z), "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04");
    assert_eq!(fp2_hex(&out.t), "f5ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff041000000000000000000000000000000000000000000000000000000000000000");
}

// --- Section 12: apply_isomorphism_general (t=0) ---
#[test]
fn test_apply_isomorphism_general_t_zero() {
    let fc = fp2_constants();
    let m1 = set_base_change_matrix_from_precomp(&precomp_matrix(0), &fc);
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::zero(),
    };
    let out = apply_isomorphism_general(&m1, &p, false);
    assert_eq!(fp2_hex(&out.x), "0a000000000000000000000000000000000000000000000000000000000000000500000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&out.y), "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04faffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04");
    assert_eq!(fp2_hex(&out.z), "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040500000000000000000000000000000000000000000000000000000000000000");
    assert_eq!(fp2_hex(&out.t), "f5ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff040500000000000000000000000000000000000000000000000000000000000000");
}

// --- Additional algebraic property tests ---

// Matrix multiplication associativity: (M0*M4)*M9 == M0*(M4*M9)
#[test]
fn test_matrix_mul_associativity() {
    let fc = fp2_constants();
    let m0 = set_base_change_matrix_from_precomp(&precomp_matrix(0), &fc);
    let m4 = set_base_change_matrix_from_precomp(&precomp_matrix(4), &fc);
    let m9 = set_base_change_matrix_from_precomp(&precomp_matrix(9), &fc);

    let left = base_change_matrix_multiplication(&base_change_matrix_multiplication(&m0, &m4), &m9);
    let right =
        base_change_matrix_multiplication(&m0, &base_change_matrix_multiplication(&m4, &m9));

    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(
                fp2_hex(&left.m[i][j]),
                fp2_hex(&right.m[i][j]),
                "associativity failed at [{}][{}]",
                i,
                j
            );
        }
    }
}

// Identity matrix (SPLITTING_TRANSFORMS[9]) is the identity
#[test]
fn test_identity_matrix() {
    let fc = fp2_constants();
    let id = set_base_change_matrix_from_precomp(&precomp_matrix(9), &fc);
    let p = ThetaPoint::<L1> {
        x: Fp2::from_small(3),
        y: Fp2::from_small(5),
        z: Fp2::from_small(7),
        t: Fp2::from_small(11),
    };
    let out = apply_isomorphism(&id, &p);
    assert_eq!(fp2_hex(&out.x), fp2_hex(&p.x));
    assert_eq!(fp2_hex(&out.y), fp2_hex(&p.y));
    assert_eq!(fp2_hex(&out.z), fp2_hex(&p.z));
    assert_eq!(fp2_hex(&out.t), fp2_hex(&p.t));
}

// Precomp tables decode correctly: FP2_CONSTANTS round-trips
#[test]
fn test_fp2_constants_values() {
    let fc = fp2_constants();
    // 0
    assert!(bool::from(fc[0].ct_is_zero()));
    // 1
    assert!(bool::from(fc[1].ct_is_one()));
    // i: re=0, im=1
    assert!(bool::from(fc[2].re.ct_is_zero()));
    // -1
    let neg1 = Fp2::<L1>::one().neg();
    assert!(bool::from(fc[3].ct_equal(&neg1)));
    // -i
    let neg_i = Fp2::<L1>::i_element().neg();
    assert!(bool::from(fc[4].ct_equal(&neg_i)));
}
