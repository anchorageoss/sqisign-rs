//! Phase 5b.4: the Kani-embedding integer matrices and the sum-of-two-squares,
//! validated against the additive oracle `kani_matrices_l1.json` (produced by
//! `sqisignhd-harness/extract_kani_matrices.py`) for all 5 vectors.
//!
//! Exact integer equality (these are integer/matrix quantities, not projective):
//!
//! * the norm equation `N = 2^136 - q = a1² + a2²` (`a1` odd, `a2` even);
//! * `matrix_F`, `matrix_F_dual` (8×8 over Z/2^70);
//! * the kernel blocks `C, D` of `complete_kernel_matrix_F1` / `_F2_dual`
//!   (the columns `kernel_basis` consumes);
//! * `gluing_base_change_matrix_dim2_F1` / `_F2` (4×4 over Z/4).
//!
//! The symplectic completion (`complete_symplectic_dim4`) routes through PARI's
//! `matsolvemod` in the reference (see NOTES); here we validate it satisfies
//! the symplectic conditions on the real kernel blocks for all 5 vectors.

mod hd_common;
use hd_common::{load, PHASE0_VECTORS};

use crypto_bigint::{Integer, U256};
use serde_json::Value;
use sqisign_verify::hd::{
    complete_symplectic_dim4, gluing_dim2_f1, gluing_dim2_f2, is_symplectic_dim4, kernel_matrix_f1,
    kernel_matrix_f2_dual, matrix_f, matrix_f_dual, norm_equation_2f_minus_q,
};

const KANI_VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/kani_matrices_l1.json"
);
const MASK70: u128 = (1u128 << 70) - 1;

/// Decimal string -> U256 (values < 2^256).
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

/// Low 128 bits of a U256.
fn low_u128(x: &U256) -> u128 {
    let w = x.to_words();
    (w[0] as u128) | ((w[1] as u128) << 64)
}

fn u256_to_u128(x: &U256) -> u128 {
    let w = x.to_words();
    assert!(w[2] == 0 && w[3] == 0, "value exceeds 128 bits");
    (w[0] as u128) | ((w[1] as u128) << 64)
}

/// Parse an `r×c` matrix of decimal strings into `u128`.
fn parse_mat(node: &Value) -> Vec<Vec<u128>> {
    node.as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|e| e.as_str().unwrap().parse::<u128>().unwrap())
                .collect()
        })
        .collect()
}

fn eq_mat8(got: &[[u128; 8]; 8], want: &[Vec<u128>], label: &str) {
    for i in 0..8 {
        for j in 0..8 {
            assert_eq!(got[i][j], want[i][j], "{label}: entry ({i},{j})");
        }
    }
}

fn eq_mat4(got: &[[u128; 4]; 4], want: &[Vec<u128>], label: &str) {
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(got[i][j], want[i][j], "{label}: entry ({i},{j})");
        }
    }
}

fn eq_mat4_u8(got: &[[u8; 4]; 4], want: &[Vec<u128>], label: &str) {
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(got[i][j] as u128, want[i][j], "{label}: entry ({i},{j})");
        }
    }
}

#[test]
fn kani_matrices_match_oracle() {
    let kani = load(KANI_VECTORS);
    let _main = load(PHASE0_VECTORS); // ensures the harness vectors are present
    let mut n = 0;
    for v in kani["vectors"].as_array().unwrap() {
        let vi = v["index"].as_u64().unwrap();

        // post-swap (a1 odd, a2 even), q, and their reductions.
        let a1 = dec_u256(v["a1"].as_str().unwrap());
        let a2 = dec_u256(v["a2"].as_str().unwrap());
        let q = dec_u256(v["q"].as_str().unwrap());

        // (1) Sum-of-two-squares: N = 2^136 - q = a1² + a2².
        let (ra1, ra2) = norm_equation_2f_minus_q(136, &q).expect("Cornacchia");
        assert_eq!(ra1, a1, "vec {vi}: a1 mismatch");
        assert_eq!(ra2, a2, "vec {vi}: a2 mismatch");
        let n_val = U256::ONE.shl(136).wrapping_sub(&q);
        let chk = ra1.wrapping_mul(&ra1).wrapping_add(&ra2.wrapping_mul(&ra2));
        assert_eq!(chk, n_val, "vec {vi}: a1²+a2² != N");
        assert!(
            bool::from(ra1.is_odd()) && !bool::from(ra2.is_odd()),
            "vec {vi}: parity"
        );

        // Reduced inputs for the matrices.
        let a1m = u256_to_u128(&a1);
        let a2m = u256_to_u128(&a2);
        let qm70 = low_u128(&q) & MASK70;

        // (2) matrix_F, matrix_F_dual (8×8 over Z/2^70).
        eq_mat8(
            &matrix_f(a1m, a2m, qm70),
            &parse_mat(&v["matrix_F"]),
            &format!("vec {vi} matrix_F"),
        );
        eq_mat8(
            &matrix_f_dual(a1m, a2m, qm70),
            &parse_mat(&v["matrix_F_dual"]),
            &format!("vec {vi} matrix_F_dual"),
        );

        // (3) kernel blocks C, D (the kernel_basis columns).
        let (c1, d1) = kernel_matrix_f1(a1m, a2m, qm70);
        eq_mat4(&c1, &parse_mat(&v["ckm_F1_C"]), &format!("vec {vi} F1_C"));
        eq_mat4(&d1, &parse_mat(&v["ckm_F1_D"]), &format!("vec {vi} F1_D"));
        let (c2, d2) = kernel_matrix_f2_dual(a1m, a2m, qm70);
        eq_mat4(
            &c2,
            &parse_mat(&v["ckm_F2dual_C"]),
            &format!("vec {vi} F2dual_C"),
        );
        eq_mat4(
            &d2,
            &parse_mat(&v["ckm_F2dual_D"]),
            &format!("vec {vi} F2dual_D"),
        );

        // (4) gluing dim-2 matrices (4×4 over Z/4).
        eq_mat4_u8(
            &gluing_dim2_f1(a1m, a2m, qm70),
            &parse_mat(&v["gluing_dim2_F1"]),
            &format!("vec {vi} gluing_F1"),
        );
        eq_mat4_u8(
            &gluing_dim2_f2(a1m, a2m, qm70),
            &parse_mat(&v["gluing_dim2_F2"]),
            &format!("vec {vi} gluing_F2"),
        );

        // (5) The symplectic completion is a valid completion of the real
        // kernel blocks (property check; PARI byte-match deferred - see NOTES).
        let m1 = complete_symplectic_dim4(&c1, &d1, MASK70).expect("F1 completion");
        assert!(
            is_symplectic_dim4(&m1, MASK70),
            "vec {vi}: F1 completion not symplectic"
        );
        let m2 = complete_symplectic_dim4(&c2, &d2, MASK70).expect("F2_dual completion");
        assert!(
            is_symplectic_dim4(&m2, MASK70),
            "vec {vi}: F2_dual completion not symplectic"
        );

        n += 1;
    }
    assert_eq!(n, 5);
    println!("Kani matrices + sum-of-two-squares match the oracle for all {n} vectors (completion validated symplectic)");
}
