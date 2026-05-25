//! Byte-for-byte cross-validation of theta isogeny (Groups 4, 5)
//! against the reference output from `tools/c-validate/theta_isog_cval`.

use sqisign_verify::ec::basis::{ec_curve_to_basis_2f_to_hint, lift_basis};
use sqisign_verify::ec::point::xdbl_a24;
use sqisign_verify::ec::{EcBasis, EcCurve};
use sqisign_verify::fp::Fp2;
use sqisign_verify::params::Level1;
use sqisign_verify::precomp::level1::*;
use sqisign_verify::theta::couple::double_couple_point;
use sqisign_verify::theta::gluing::{
    gluing_compute, gluing_eval_point, gluing_eval_point_special_case, verify_two_torsion,
};
use sqisign_verify::theta::isogeny::{
    theta_isogeny_compute, theta_isogeny_compute_2, theta_isogeny_compute_4, theta_isogeny_eval,
};
use sqisign_verify::theta::theta_structure::{double_iter, double_point, theta_precomputation};
use sqisign_verify::theta::{
    ThetaCoupleCurve, ThetaCoupleJacPoint, ThetaCouplePoint, ThetaStructure,
};

type L1 = Level1;

fn fp2_hex(v: &Fp2<L1>) -> String {
    v.encode().iter().map(|b| format!("{:02x}", b)).collect()
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

fn gen_full_basis(e0: &mut EcCurve<L1>) -> EcBasis<L1> {
    let (basis, _) = ec_curve_to_basis_2f_to_hint(
        e0,
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

/// Double down an (X:Z) basis to 2^n torsion.
fn double_down_basis(bas: &EcBasis<L1>, e0: &EcCurve<L1>, target_bits: u32) -> EcBasis<L1> {
    let dbl_count = TORSION_EVEN_POWER - target_bits;
    let mut k1 = bas.p.clone();
    let mut k2 = bas.q.clone();
    let mut k1m2 = bas.pmq.clone();
    for _ in 0..dbl_count {
        k1 = xdbl_a24(&k1, &e0.a24, e0.is_a24_computed_and_normalized);
        k2 = xdbl_a24(&k2, &e0.a24, e0.is_a24_computed_and_normalized);
        k1m2 = xdbl_a24(&k1m2, &e0.a24, e0.is_a24_computed_and_normalized);
    }
    EcBasis::new(k1, k2, k1m2)
}

// ===== Section 1: gluing_compute =====

#[test]
fn test_gluing_compute() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);

    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, ok) = lift_basis(&mut bas8_mut, &mut e_mut);
    assert!(bool::from(ok));

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing =
        gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).expect("gluing_compute should succeed");

    assert_eq!(fp2_hex(&gluing.codomain.x), "0c1ee0f5c9c69338dd1a9ca0a13f9f7664542723bedc902d1eb0f8f0ddfee50163eeb9ddcb824fbb59cff073f1546d3e57aa0e716da69bfa87499c46393de603");
    assert_eq!(fp2_hex(&gluing.codomain.y), "e58253a86f23f1305ef85493a9f61ccbf00ab7fa6e327215434a20d972bdcb005ccac17ab16f89c2c860cfbe0b435535bfb6a91661742e8c63bb8e14306ed104");
    assert_eq!(fp2_hex(&gluing.codomain.z), "e58253a86f23f1305ef85493a9f61ccbf00ab7fa6e327215434a20d972bdcb005ccac17ab16f89c2c860cfbe0b435535bfb6a91661742e8c63bb8e14306ed104");
    assert_eq!(fp2_hex(&gluing.codomain.t), "bde7c65a15804e29dfd50d86b1ad9a1f7dc146d21f8853fd67e447c1077cb10456a6c917975cc3c937f2ad0926313d2c27c344bc5442c11d3f2d81e2269fbc00");

    assert_eq!(fp2_hex(&gluing.image_k1_8.x), "df90e7d3b7a67d26fa18016e9f3c727e54a031b04495f33ec795bb595ff82c02730e9e0a69d8f5c0976a5e93a4adb7943ed1438e0db4970365892d158831e201");
    assert_eq!(fp2_hex(&gluing.image_k1_8.y), "220808291f9116d77443fa55a1afbb985c7284cf2521dc06f91e2c59bf576502e81d4d241f6bc793720e7bf50438c4cd9acddbe75823b3f0147127de946e6000");

    assert_eq!(fp2_hex(&gluing.precomputation.x), "9644a40b7c6176c3e22416fcc1634a2562828e80db91a2f05fd75fefc454ad02a36e7095aa464c1eb0db81c11bddd74b32bffd6ff01441dbc61845dce3453200");
    assert_eq!(fp2_hex(&gluing.precomputation.y), "934dc626ad51d1833f91a3067c24c1d5b924389427550f8ced32ec8bb5200d0303127c318d09637c48b790daf2088c04cc79322d0699363712c7069984670a02");
    assert_eq!(fp2_hex(&gluing.precomputation.z), "934dc626ad51d1833f91a3067c24c1d5b924389427550f8ced32ec8bb5200d0303127c318d09637c48b790daf2088c04cc79322d0699363712c7069984670a02");
    assert_eq!(fp2_hex(&gluing.precomputation.t), "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

// ===== Section 1b: basis change matrix =====

#[test]
fn test_gluing_basis_change_matrix() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let expected = [
        ["6817ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b5010000000000000000000000000000000000000000000000000000000000000000",
         "225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100",
         "225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100",
         "6617ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b501ec8e0dbe7c8224fa32ed52599e7b83fdb3825cd6e328928957e2fc7003b1fe01"],
        ["225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100",
         "6617ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b501ec8e0dbe7c8224fa32ed52599e7b83fdb3825cd6e328928957e2fc7003b1fe01",
         "6817ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b5010000000000000000000000000000000000000000000000000000000000000000",
         "225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100"],
        ["225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100",
         "6817ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b5010000000000000000000000000000000000000000000000000000000000000000",
         "6617ffc4ef22d05f46b47a878c32eeb88773d8173ab7f8340aab0e4bb7c5b501ec8e0dbe7c8224fa32ed52599e7b83fdb3825cd6e328928957e2fc7003b1fe01",
         "225daa41a50bf01fc2e6282d84bba43d2dd1f207bee752bc588eafc3e741670139f2554bd48ff870e941a0a627a7794a7477cf72204f17bad574d487ad6c0100"],
        ["99e8003b10dd2fa0b94b857873cd1147788c27e8c54807cbf554f1b4483a4a031371f241837ddb05cd12ada661847c024c7da3291cd76d76a81d038ffc4e0103",
         "dda255be5af40fe03d19d7d27b445bc2d22e0df84118ad43a771503c18be9803c60daab42b70078f16be5f59d85886b58b88308ddfb0e8452a8b2b785293fe04",
         "dda255be5af40fe03d19d7d27b445bc2d22e0df84118ad43a771503c18be9803c60daab42b70078f16be5f59d85886b58b88308ddfb0e8452a8b2b785293fe04",
         "97e8003b10dd2fa0b94b857873cd1147788c27e8c54807cbf554f1b4483a4a030000000000000000000000000000000000000000000000000000000000000000"],
    ];

    for (i, row) in expected.iter().enumerate() {
        for (j, exp) in row.iter().enumerate() {
            assert_eq!(
                fp2_hex(&gluing.basis_change.m[i][j]),
                *exp,
                "M[{}][{}] mismatch",
                i,
                j
            );
        }
    }
}

// ===== Section 2: gluing_eval_point =====

#[test]
fn test_gluing_eval_point() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let zero = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";

    let eval_k1 = gluing_eval_point(&xy_k1_8, &gluing);
    assert_eq!(fp2_hex(&eval_k1.x), zero);
    assert_eq!(fp2_hex(&eval_k1.y), zero);
    assert_eq!(fp2_hex(&eval_k1.z), zero);
    assert_eq!(fp2_hex(&eval_k1.t), zero);

    let eval_k2 = gluing_eval_point(&xy_k2_8, &gluing);
    assert_eq!(fp2_hex(&eval_k2.x), "11265b48fae47b3341b75210ef0be3707ebeb1f172ae1ff11e5d24d55a636204cd31ffec4c03350caa7e27748c21896f473c9ebca1a5d7c77480a74b14572900");
    assert_eq!(fp2_hex(&eval_k2.y), "4829c10498adf410978be5e9ee1525ea6aa6cc9b3a7076d720a0f4f4c2b38b01d4448292efa303ce29c97e6882110287eb8725bf07e171903fdda2c21affa000");
    assert_eq!(fp2_hex(&eval_k2.z), zero);
    assert_eq!(fp2_hex(&eval_k2.t), zero);
}

// ===== Section 3: gluing_eval_point_special_case =====

#[test]
fn test_gluing_eval_point_special_case() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let sc_pt = ThetaCouplePoint {
        p1: bas8.p.clone(),
        p2: bas8.q.clone(),
    };
    let sc_img =
        gluing_eval_point_special_case(&sc_pt, &gluing).expect("special case should succeed");

    assert_eq!(fp2_hex(&sc_img.x), "3ad34050d855c2564bc03fd363fc80b62e5d711a914f133e797b52b1bf5e51035a0b7fff712965f510c4321194944a7a2c5c88371f73bfb85bb244ff8cf45d03");
    assert_eq!(fp2_hex(&sc_img.y), "1b115f895449d94d91698f402a71d91f83cb98fd6791e6907136a33d26c582024acf941c315eb1e570c3bc4f9325d874a4248169f0f135a6dc4b42ab913ad801");
    assert_eq!(fp2_hex(&sc_img.z), "03b1823d2fc30fbb28ed20520f1ace7628c63f1fc12c461c960e0c3673d44b03c56c55c60f6d022a2f3db9716d499a90e31286643e8f536ca21ac0a8697fad04");
    assert_eq!(fp2_hex(&sc_img.t), "e4eea076abb626b26e9670bfd58e26e07c346702986e196f8ec95cc2d93a7d02b5306be3cea14e1a8f3c43b06cda278b5bdb7e960f0eca5923b4bd546ec52703");
}

// ===== Section 4: theta_isogeny_compute =====

#[test]
fn test_theta_isogeny_compute() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    // 32-torsion through gluing
    let bas32 = double_down_basis(&full_basis, &e0, 5);
    let mut bas32_mut = bas32.clone();
    let mut e_mut2 = e0.clone();
    let (xy_t1_32, xy_t2_32, _) = lift_basis(&mut bas32_mut, &mut e_mut2);

    let xy_t1_cp = ThetaCoupleJacPoint {
        p1: xy_t1_32.clone(),
        p2: xy_t2_32.clone(),
    };
    let xy_t2_cp = ThetaCoupleJacPoint {
        p1: xy_t2_32.clone(),
        p2: xy_t1_32.clone(),
    };

    let theta_t1 = gluing_eval_point(&xy_t1_cp, &gluing);
    let theta_t2 = gluing_eval_point(&xy_t2_cp, &gluing);

    assert_eq!(fp2_hex(&theta_t1.x), "dcffc9770eb89b3888e3da62bf570f7f11391ba3b626b7dd700c0ce1af862a02cfef898edeaf60dd9e1ccade3c317b721f0bc13e643c3e8b807df1047daeeb04");
    assert_eq!(fp2_hex(&theta_t2.x), "484f4de5ec0f357e9bb73d7e8c738da8fddbd5610056ab498f11e8583f079804a6e11f3187d19d45236f928bc288ff1cf1408a5388a282edc00b48d6dd928800");

    let mut ts = ThetaStructure {
        null_point: gluing.codomain.clone(),
        ..ThetaStructure::default()
    };
    theta_precomputation(&mut ts);

    let t1_8 = double_iter(&theta_t1, &mut ts, 2);
    let t2_8 = double_iter(&theta_t2, &mut ts, 2);

    assert_eq!(fp2_hex(&t1_8.x), "7cd3bd5d3ef9c1fda2858caa02d94fa9eee3a8fcdb1113abd8226cac2d3060002b5884fb981531d1b9919a00e17bc7b22d48fb4f530721d0e9782ac736519004");
    assert_eq!(fp2_hex(&t2_8.x), "5c7eae6305a916c3bf8d3140dc183e229bc626e0a6af3d72ecd1abedef5daf0075cc5639597f0aa0e6d5ff85fd33822c8fef8d18e49906e2c7f7e09570f50902");

    let isog = theta_isogeny_compute(&ts, &t1_8, &t2_8, false, false, false)
        .expect("theta_isogeny_compute should succeed");

    let v = "9d23708e839ed47b89cc54a9e90368984ca927dd5c2f7d1d982ac15c3aa0e5048443bf4eb379e950da83689efc691ef6055cec67b62bfc633f1d603b3d127d02";
    assert_eq!(fp2_hex(&isog.codomain.null_point.x), v);
    assert_eq!(fp2_hex(&isog.codomain.null_point.y), v);
    assert_eq!(fp2_hex(&isog.codomain.null_point.z), v);
    assert_eq!(fp2_hex(&isog.codomain.null_point.t), v);
    assert_eq!(fp2_hex(&isog.precomputation.x), v);
    assert_eq!(fp2_hex(&isog.precomputation.y), v);
    assert_eq!(fp2_hex(&isog.precomputation.z), v);
    assert_eq!(fp2_hex(&isog.precomputation.t), v);
}

// ===== Section 5: theta_isogeny_eval =====

#[test]
fn test_theta_isogeny_eval() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let bas32 = double_down_basis(&full_basis, &e0, 5);
    let mut bas32_mut = bas32.clone();
    let mut e_mut2 = e0.clone();
    let (xy_t1_32, xy_t2_32, _) = lift_basis(&mut bas32_mut, &mut e_mut2);

    let xy_t1_cp = ThetaCoupleJacPoint {
        p1: xy_t1_32.clone(),
        p2: xy_t2_32.clone(),
    };
    let xy_t2_cp = ThetaCoupleJacPoint {
        p1: xy_t2_32.clone(),
        p2: xy_t1_32.clone(),
    };

    let theta_t1 = gluing_eval_point(&xy_t1_cp, &gluing);
    let theta_t2 = gluing_eval_point(&xy_t2_cp, &gluing);

    let mut ts = ThetaStructure {
        null_point: gluing.codomain.clone(),
        ..ThetaStructure::default()
    };
    theta_precomputation(&mut ts);

    let t1_8 = double_iter(&theta_t1, &mut ts, 2);
    let t2_8 = double_iter(&theta_t2, &mut ts, 2);

    let isog = theta_isogeny_compute(&ts, &t1_8, &t2_8, false, false, false).unwrap();

    let eval_out = theta_isogeny_eval(&isog, &theta_t1);
    assert_eq!(fp2_hex(&eval_out.x), "cf93f9e9d78110e963182157e12974b2acd2e5438769a5e68cdfbf944baeb0019c637a1bfa5c14ea1314b3c9f99bd6b7f39802db9ab8c4f1c9bb947a86d24f03");
    assert_eq!(fp2_hex(&eval_out.y), "ba00170c9e99888b5e7f975de7e3493580424d5b3b615bd6ce85bde57ed0520324178a3650c0cdd71339f3d53d3821db5b4c028e9ab2bd1447f3906518643d01");
    assert_eq!(fp2_hex(&eval_out.z), "99c0afbd674a015a72f550d1cf01d6bdc5681ba1d1f48bceda450e9ed3df990402a967acd1d5cb2d27625e055f05cf2177236e4905cec7e6a4a77578df806104");
    assert_eq!(fp2_hex(&eval_out.t), "e6e553881ad40b1fec18d0d3f7d14b322a2e896e274b3c3f56e9487dbaabc5039cdb4dfc40f27c743fa8234e8fecd7c70d73b4f4b4b31238b844038c1dc11803");
}

// ===== Section 6: theta_isogeny_compute_4 =====

#[test]
fn test_theta_isogeny_compute_4() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let bas32 = double_down_basis(&full_basis, &e0, 5);
    let mut bas32_mut = bas32.clone();
    let mut e_mut2 = e0.clone();
    let (xy_t1_32, xy_t2_32, _) = lift_basis(&mut bas32_mut, &mut e_mut2);

    let xy_t1_cp = ThetaCoupleJacPoint {
        p1: xy_t1_32.clone(),
        p2: xy_t2_32.clone(),
    };
    let xy_t2_cp = ThetaCoupleJacPoint {
        p1: xy_t2_32.clone(),
        p2: xy_t1_32.clone(),
    };

    let theta_t1 = gluing_eval_point(&xy_t1_cp, &gluing);
    let theta_t2 = gluing_eval_point(&xy_t2_cp, &gluing);

    let mut ts = ThetaStructure {
        null_point: gluing.codomain.clone(),
        ..ThetaStructure::default()
    };
    theta_precomputation(&mut ts);

    let t1_8 = double_iter(&theta_t1, &mut ts, 2);
    let t2_8 = double_iter(&theta_t2, &mut ts, 2);

    let t1_4 = double_point(&t1_8, &mut ts);
    let t2_4 = double_point(&t2_8, &mut ts);

    assert_eq!(fp2_hex(&t1_4.x), "0b56b72bf800c18d6730ac38445da3034245b2186767f156cf873ba662b521034149b2d9f6cc471476c17f2d6a8120717fd7b5fc7c5b503a9127ab7ce8e10003");
    assert_eq!(fp2_hex(&t2_4.x), "a7978af6fae81245ddcc8e5a87dab69bc9c1a9f2ce86cad2506a59fe1d34610034f17fde55892fc622fb988fa9792e491bf0d40ab80211d37392f69157fb3702");

    let isog4 = theta_isogeny_compute_4(&ts, &t1_4, &t2_4, false, false);

    assert_eq!(fp2_hex(&isog4.codomain.null_point.x), "94cc584fe8d6060618d30e9b4e86580bc8092e3602eb2f9363a45572ef8b8f025a6f420d685e917daf5670bea41425e68e9db2c52aeedc72f5435d6410706303");
    assert_eq!(fp2_hex(&isog4.codomain.null_point.y), "7cce5d4bc3c347520a324237210fe711c3800c600434335e17a7f773fbc35f024a94804a5e2e39f70e88a382ed51967a29b7d2512fa341b9b0f39bd029218b03");
    assert_eq!(fp2_hex(&isog4.codomain.null_point.z), "7cce5d4bc3c347520a324237210fe711c3800c600434335e17a7f773fbc35f024a94804a5e2e39f70e88a382ed51967a29b7d2512fa341b9b0f39bd029218b03");
    assert_eq!(fp2_hex(&isog4.codomain.null_point.t), "488256cc139bb79143cb95abad0b7ed797d595522d4cfedb2d993668e80d5e0231d38200ea42979240e52fe6969a1a50733aa40527adb8962bbe21933d1ed303");

    assert_eq!(fp2_hex(&isog4.precomputation.x), "44e6e4626ab4ab82af18c68d302d4c19ee09e32d31e2f5d0fc4c7ae673515b04c5445b1bf1ef7c6751b48d56f6a2df4d5cdcc47fa81c75d603596e5060ff8600");
    assert_eq!(fp2_hex(&isog4.precomputation.y), "671a04be4ffdefbfacecf9466bed5a192293c8d116374b0ab284bab77da22f0080e66b348f27884713b9954290f99877e290f890709fa7028eb6bccdff5ceb02");
    assert_eq!(fp2_hex(&isog4.precomputation.z), "671a04be4ffdefbfacecf9466bed5a192293c8d116374b0ab284bab77da22f0080e66b348f27884713b9954290f99877e290f890709fa7028eb6bccdff5ceb02");
    assert_eq!(fp2_hex(&isog4.precomputation.t), "65ccaf70168b8686e995364f16d5475c4f2f89e40ba57e3dbd37e440bba2ba0455d584c3210611c13906b8d0d4878f60aa1920b61d4229d1a42469bc37916701");

    let eval4_out = theta_isogeny_eval(&isog4, &theta_t1);
    assert_eq!(fp2_hex(&eval4_out.x), "a40b9a6ffac088e3068fb1c72f906daed9052953bb684fb4ddbad80180854303a83de34f261723171ee1d20eee42d30e57c461985fb255a9ae16f002851a9700");
    assert_eq!(fp2_hex(&eval4_out.y), "a4d3684cb066acf0b9e901c714174a32055fcb144fc13861991c9d538cb34d00fb4966760492c067e15b46acdfdf16756f63a99e1dd0249f17e6114e64419704");
    assert_eq!(fp2_hex(&eval4_out.z), "bf65cff1cf799b4b1c7ecd1be954dd3fef2738fc824295aaed08451efbc4a100149f8e787938b905bd5a69c9914aa20da7896d3872cfdc03973a6100b74f0e02");
    assert_eq!(fp2_hex(&eval4_out.t), "287ddf01736e8021324ef5c5229bc83877dad9487f63acabf85c960c22f67801b06449dbe037b4d44077972995cf8fab51c1ca8f845f0716d277bf63df744403");
}

// ===== Section 7b: theta_isogeny_compute_2 =====

#[test]
fn test_theta_isogeny_compute_2() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);
    let mut bas8_mut = bas8.clone();
    let mut e_mut = e0.clone();
    let (xy_k1, xy_k2, _) = lift_basis(&mut bas8_mut, &mut e_mut);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };
    let xy_k1_8 = ThetaCoupleJacPoint {
        p1: xy_k1.clone(),
        p2: xy_k2.clone(),
    };
    let xy_k2_8 = ThetaCoupleJacPoint {
        p1: xy_k2.clone(),
        p2: xy_k1.clone(),
    };

    let gluing = gluing_compute(&e12, &xy_k1_8, &xy_k2_8, true).unwrap();

    let bas32 = double_down_basis(&full_basis, &e0, 5);
    let mut bas32_mut = bas32.clone();
    let mut e_mut2 = e0.clone();
    let (xy_t1_32, xy_t2_32, _) = lift_basis(&mut bas32_mut, &mut e_mut2);

    let xy_t1_cp = ThetaCoupleJacPoint {
        p1: xy_t1_32.clone(),
        p2: xy_t2_32.clone(),
    };
    let xy_t2_cp = ThetaCoupleJacPoint {
        p1: xy_t2_32.clone(),
        p2: xy_t1_32.clone(),
    };

    let theta_t1 = gluing_eval_point(&xy_t1_cp, &gluing);
    let theta_t2 = gluing_eval_point(&xy_t2_cp, &gluing);

    let mut ts = ThetaStructure {
        null_point: gluing.codomain.clone(),
        ..ThetaStructure::default()
    };
    theta_precomputation(&mut ts);

    let t1_8 = double_iter(&theta_t1, &mut ts, 2);
    let t2_8 = double_iter(&theta_t2, &mut ts, 2);
    let t1_2 = double_iter(&t1_8, &mut ts, 2);
    let t2_2 = double_iter(&t2_8, &mut ts, 2);

    assert_eq!(fp2_hex(&t1_2.x), "11fbf8b22232be12e012b9c88c9dd3e5f7a53f2e9b8f98bf5370cc25e2f5030012d5c70c570fbeebaa37690605dc3bf6aa9b4f8ee6d2ddeddd3c0bb691ffe300");
    assert_eq!(fp2_hex(&t2_2.x), "315528b91aac5ca2f5c979c5840f303860acee777b4a2da217d6e4dde3e4ed049de1bf40a7531ae688fbe91b89ccc14b69a19d5aba9553e4d48d520d2aa80304");

    let isog2 = theta_isogeny_compute_2(&ts, &t1_2, &t2_2, false, false);

    assert_eq!(fp2_hex(&isog2.codomain.null_point.x), "3ea1a96b4193ec2f07c024b56df71122af13dfb730b72f864a1b10da1245b70014b737acaee0711d9073492b20dbf46e2314e645f37869ed89c93c71f3c05f02");
    assert_eq!(fp2_hex(&isog2.codomain.null_point.y), "3063690bbca4bd41cb021b23045894fd278f15fd96ecca445a1e9d293b1e44021a70636381b377bc3feb3d0898f73ad6f7010210e46a81230b2216c6088be101");
    assert_eq!(fp2_hex(&isog2.codomain.null_point.z), "3063690bbca4bd41cb021b23045894fd278f15fd96ecca445a1e9d293b1e44021a70636381b377bc3feb3d0898f73ad6f7010210e46a81230b2216c6088be101");
    assert_eq!(fp2_hex(&isog2.codomain.null_point.t), "4abbf19b9081c8ce2a8ff75f1aecd3f8366ad329dd09a0c5c4c5c3524afcd500ee83eb29c9506e5ffa4296690f63eef21c1cbe0a84fb2e617c122650692b7104");

    assert_eq!(fp2_hex(&isog2.precomputation.x), "51cf900ca46b85055c0bf6c3d59fba9534f7b9fddff135b14f81cab5c5120704796896a34451b96da2fcaaea39870437728f8287a69d92c9bc73c3550da3d900");
    assert_eq!(fp2_hex(&isog2.precomputation.y), "baf57d3ecb58ef61d19aca7a22665254692f1f0dee265297f592adc969da6d042adb7244b13350e6e1eb3df44da1e7ea719ee32d53bb54d15f21a85f3a484504");
    assert_eq!(fp2_hex(&isog2.precomputation.z), "baf57d3ecb58ef61d19aca7a22665254692f1f0dee265297f592adc969da6d042adb7244b13350e6e1eb3df44da1e7ea719ee32d53bb54d15f21a85f3a484504");
    assert_eq!(fp2_hex(&isog2.precomputation.t), "cf2a4168e7a112a5744c70a2803335a4c14fd80fea45020ec9c8190aa91ff403cbdb0468d2ae231043685e8dbe2b3bf5db57ce82d155edc40a7166528b8a8204");

    let eval2_out = theta_isogeny_eval(&isog2, &theta_t1);
    assert_eq!(fp2_hex(&eval2_out.x), "bfd00fc863ddba2c2b2f2763ce829810e2f7cd9d5cfddc77bcc9af93990243047f0c9f99586ab6d70bd683e863f3487d3a96370e282c338890f27dce1079c504");
    assert_eq!(fp2_hex(&eval2_out.y), "7b0cb58a565874a7ff18314e3c69bde9fda4535d7761b2a1bd03acba257049049cb6ca2b650346b34efe54e5e1bf44d0e313d1d70e4dd0625025a86d3b7b2702");
    assert_eq!(fp2_hex(&eval2_out.z), "da44a4510e1d38ebfab2812f3ebadb889815e9280e6cf1a5e1e631072f3a41016364822108f5e07355207eb9e9def40394ceb62b565dca77ed73dd94d9eb7203");
    assert_eq!(fp2_hex(&eval2_out.t), "fbdb401105bcf1d5927054e4f51ed5cf87d6f9145a7b6a802e5db5a8b7c2aa019d266e800c86cb82d8bd162c5619e27cad4e306f7cee923daad659a6ea637301");
}

// ===== Section 8: verify_two_torsion =====

#[test]
fn test_verify_two_torsion() {
    let mut e0 = make_e0();
    let full_basis = gen_full_basis(&mut e0);
    let bas8 = double_down_basis(&full_basis, &e0, 3);

    let e12 = ThetaCoupleCurve {
        e1: e0.clone(),
        e2: e0.clone(),
    };

    let k1_cp = ThetaCouplePoint {
        p1: bas8.p.clone(),
        p2: bas8.q.clone(),
    };
    let k2_cp = ThetaCouplePoint {
        p1: bas8.q.clone(),
        p2: bas8.p.clone(),
    };

    let k1_4 = double_couple_point(&k1_cp, &e12);
    let k2_4 = double_couple_point(&k2_cp, &e12);
    let k1_2 = double_couple_point(&k1_4, &e12);
    let k2_2 = double_couple_point(&k2_4, &e12);

    assert!(verify_two_torsion(&k1_2, &k2_2, &e12));
    assert!(!verify_two_torsion(&k1_2, &k1_2, &e12));
}
