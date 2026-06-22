//! Phase 5b.6 (front half): the dimension-2 `(2,2)`-isogeny chain validated
//! against the oracle's recorded `dim2_codomain_null` (`product_theta_l1.json`)
//! for all 5 vectors, both half-chains (F1 and F2_dual).
//!
//! This is the inner engine of the gluing chain. The test self-derives
//! everything the chain consumes - the canonical dim-1 bases (`make_canonical`),
//! their product null (`Theta12.zero()`), the symplectic→theta base change
//! `N_dim2 = base_change_theta_dim2(M0·M1, e4)`, and the kernel `B_K_dim2` - from
//! the response basis (stages 2-3), then builds [`IsogenyChainDim2`] and compares
//! its codomain theta-null to the oracle (projectively).
//!
//! The Kani integers `a1, a2, q, m` are read from `kani_matrices_l1.json` (the
//! validated 5b.4 quantities); the full self-derivation from the norm equation
//! is exercised in `hd_verify`. The base change uses `e4 = e₄(T1,T2)⁻¹` (the
//! biextension `weil` is PARI's inverse, 5b.2), matching the oracle's PARI `e4`.
//!
//! For `m = 1` (vectors 1, 4) the chain is the gluing alone; `m ∈ {2,3,6}`
//! (vectors 0, 2, 3) exercise the plain `ThetaIsogenyDim2` steps and so confirm
//! the image map (the codomain at step `k` is computed from the 8-torsion pushed
//! through steps `0..k`).

mod hd_common;
use hd_common::{load, parse_fp2, F};

use serde_json::Value;
use sqisign_verify::ec::jacobian::{jac_add, jac_dbl};
use sqisign_verify::ec::pairing::weil;
use sqisign_verify::ec::{EcCurve, JacPoint};
use sqisign_verify::hd::{
    base_change_theta_dim2, gluing_dim2_f1, gluing_dim2_f2, jac_to_affine, make_canonical,
    product_theta_dim2, recover_challenge_l1, recover_response_l1, IsogenyChainDim2,
    ResponseScalars, ThetaStructureDim1, TuplePoint,
};
use sqisign_verify::Level1;

const PRODUCT_THETA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/product_theta_l1.json"
);
const KANI: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../sqisignhd-harness/kani_matrices_l1.json"
);

/// f - 2 = 70 - 2 doublings take an order-2^70 point to the 4-torsion.
const TO_4TORSION: u32 = 68;

// decimal parsing

fn decimal_to_le_limbs(s: &str) -> [u64; 4] {
    let mut limbs = [0u64; 4];
    for ch in s.trim().bytes() {
        let mut carry = (ch - b'0') as u128;
        for l in limbs.iter_mut() {
            let prod = (*l as u128) * 10 + carry;
            *l = prod as u64;
            carry = prod >> 64;
        }
    }
    limbs
}
fn dec_i128(v: &Value) -> i128 {
    v.as_str().unwrap().parse::<i128>().unwrap()
}
/// Parse a decimal that fits in `u128` (a1, a2 are ≈ 71 bits).
fn dec_u128(v: &Value) -> u128 {
    v.as_str().unwrap().trim().parse::<u128>().unwrap()
}
/// Parse a decimal mod `2^128` (q is ≈ 132 bits; only `q mod 4` is used).
fn dec_u128_wrap(v: &Value) -> u128 {
    let mut acc = 0u128;
    for ch in v.as_str().unwrap().trim().bytes() {
        acc = acc.wrapping_mul(10).wrapping_add((ch - b'0') as u128);
    }
    acc
}

// small EC helpers

fn jac_dbl_n(p: &JacPoint<Level1>, n: u32, curve: &EcCurve<Level1>) -> JacPoint<Level1> {
    let mut acc = p.clone();
    for _ in 0..n {
        acc = jac_dbl(&acc, curve);
    }
    acc
}

/// `[k]·P` for a small non-negative `k` (≤ 511 here), MSB-first double-and-add.
/// `jac_add`/`jac_dbl` are robust to identity / `P=Q` / `P=-Q`.
fn jac_mul(p: &JacPoint<Level1>, k: u64, curve: &EcCurve<Level1>) -> JacPoint<Level1> {
    if k == 0 {
        return JacPoint::identity();
    }
    let top = 63 - k.leading_zeros();
    let mut acc = p.clone();
    for i in (0..top).rev() {
        acc = jac_dbl(&acc, curve);
        if (k >> i) & 1 == 1 {
            acc = jac_add(&acc, p, curve);
        }
    }
    acc
}

#[inline]
fn jac_sub(
    a: &JacPoint<Level1>,
    b: &JacPoint<Level1>,
    curve: &EcCurve<Level1>,
) -> JacPoint<Level1> {
    jac_add(a, &b.neg(), curve)
}

/// `e_4(U, V)` (biextension `weil`, = PARI's inverse) for full points `U, V`.
fn weil4(u: &JacPoint<Level1>, v: &JacPoint<Level1>, curve: &mut EcCurve<Level1>) -> F {
    let uv = jac_add(u, &v.neg(), curve);
    weil(2, &u.to_xz(), &v.to_xz(), &uv.to_xz(), curve)
}

// matrices over Z/4

/// `M0` (`M_product_dim2`) from the two canonical change matrices `MT` (E_com)
/// and `MU` (E_chal): the interleaved block layout of the reference.
fn m0_from_canon(mt: &[[u8; 2]; 2], mu: &[[u8; 2]; 2]) -> [[u8; 4]; 4] {
    [
        [mt[0][0], 0, mt[1][0], 0],
        [0, mu[0][0], 0, mu[1][0]],
        [mt[0][1], 0, mt[1][1], 0],
        [0, mu[0][1], 0, mu[1][1]],
    ]
}

/// `(A·B) mod 4`, returned as `i64` for [`base_change_theta_dim2`].
fn mat4_mul_mod4(a: &[[u8; 4]; 4], b: &[[u8; 4]; 4]) -> [[i64; 4]; 4] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut s = 0u32;
            for k in 0..4 {
                s += (a[i][k] as u32) * (b[k][j] as u32);
            }
            (s % 4) as i64
        })
    })
}

// projective comparison

fn proj_eq4(a: &[F; 4], b: &[F; 4]) -> bool {
    for i in 0..4 {
        for j in (i + 1)..4 {
            if !bool::from(a[i].mul(&b[j]).ct_equal(&b[i].mul(&a[j]))) {
                return false;
            }
        }
    }
    true
}

/// The dim-1 theta null `(x+1, x-1)` of a canonical 4-torsion point (affine).
fn dim1_null(u1: &JacPoint<Level1>) -> [F; 2] {
    let (x, _) = jac_to_affine(u1);
    ThetaStructureDim1::<Level1>::from_torsion(&x, &F::one())
        .null()
        .clone()
}

#[test]
fn dim2_chain_matches_oracle() {
    let doc = load(hd_common::PHASE0_VECTORS);
    let pt = load(PRODUCT_THETA);
    let kani = load(KANI);
    let mut n = 0;

    for ((v, vp), vk) in doc["test_vectors"]
        .as_array()
        .unwrap()
        .iter()
        .zip(pt["vectors"].as_array().unwrap().iter())
        .zip(kani["vectors"].as_array().unwrap().iter())
    {
        let vi = v["index"].as_u64().unwrap();
        assert_eq!(v["index"], vp["index"]);
        assert_eq!(v["index"], vk["index"]);
        let s1 = &v["stage1_recover_basis"];
        let sig = &v["signature"];

        // stages 2-3: self-derive the response basis
        let a_pk = parse_fp2(&s1["E_pk_A"]);
        let hp = v["public_key"]["hint_pk_P"].as_u64().unwrap() as u32;
        let hq = v["public_key"]["hint_pk_Q"].as_u64().unwrap() as u32;
        let chal_limbs = decimal_to_le_limbs(sig["chal"].as_str().unwrap());
        let chal = recover_challenge_l1(&a_pk, hp, hq, &chal_limbs).unwrap();

        let a_com = parse_fp2(&s1["E_com_A"]);
        let hcp = sig["hint_com_P"].as_u64().unwrap() as u32;
        let hcq = sig["hint_com_Q"].as_u64().unwrap() as u32;
        let q_wrap = dec_u128_wrap(&sig["q"]);
        let rsp = recover_response_l1(
            &chal,
            &a_com,
            hcp,
            hcq,
            ResponseScalars {
                a: dec_i128(&sig["a"]),
                b: dec_i128(&sig["b"]),
                c_or_d: dec_i128(&sig["c_or_d"]),
                q: q_wrap,
            },
        )
        .unwrap();

        // Kani integers (validated 5b.4 quantities)
        let a1 = dec_u128(&vk["a1"]);
        let a2 = dec_u128(&vk["a2"]);
        let m = vk["m"].as_u64().unwrap() as usize;
        let q4 = dec_u128_wrap(&vk["q"]); // only q mod 4 is used downstream
        let lamb = q4 & 3; // q is odd ⇒ inverse_mod(q,4) = q mod 4

        // the two curves; E1 = E_com, E2 = E_chal
        let mut e_com = EcCurve::<Level1>::from_a(&a_com).unwrap();
        e_com.normalize_a24();
        let mut e_chal = chal.e_chal.clone();
        e_chal.normalize_a24();

        // canonical bases on E_com and E_chal
        // (P1_4, Q1_4) = 2^68·(R_com, S_com); (R2_4, lamb·S2_4) = 2^68·(images).
        let p1_4 = jac_dbl_n(&rsp.r_com, TO_4TORSION, &e_com);
        let q1_4 = jac_dbl_n(&rsp.s_com, TO_4TORSION, &e_com);
        let (t1, t2, mt) = make_canonical(&p1_4, &q1_4, &mut e_com).unwrap();

        let r2_4 = jac_dbl_n(&rsp.phi_rsp_r_com, TO_4TORSION, &e_chal);
        let mut s2_4 = jac_dbl_n(&rsp.phi_rsp_s_com, TO_4TORSION, &e_chal);
        if lamb == 3 {
            s2_4 = jac_add(&jac_dbl(&s2_4, &e_chal), &s2_4, &e_chal); // 3·S2_4
        }
        let (u1, _u2, mu) = make_canonical(&r2_4, &s2_4, &mut e_chal).unwrap();

        // e4 = e₄(T1, T2) in PARI convention = biextension weil inverse.
        let e4 = weil4(&t1, &t2, &mut e_com).inv();

        // product theta null (Theta12.zero()), self-derived
        let null_com = dim1_null(&t1);
        let null_chal = dim1_null(&u1);
        let zero12 = product_theta_dim2(&null_com, &null_chal);
        // sanity: equals the oracle's recorded product null (both half-chains).
        for hc in vp["half_chains"].as_array().unwrap() {
            let want: [F; 4] = core::array::from_fn(|i| parse_fp2(&hc["theta12_null"][i]));
            assert!(proj_eq4(&zero12, &want), "vec {vi}: theta12_null mismatch");
        }

        let m0 = m0_from_canon(&mt, &mu);

        // points of order 2^(m+3): 2^(67-m)·(response basis)
        let k = 67 - m as u32;
        let p1_m = jac_dbl_n(&rsp.r_com, k, &e_com);
        let q1_m = jac_dbl_n(&rsp.s_com, k, &e_com);
        let r2_m = jac_dbl_n(&rsp.phi_rsp_r_com, k, &e_chal);
        let s2_m = jac_dbl_n(&rsp.phi_rsp_s_com, k, &e_chal);

        let two_mp2 = 1u128 << (m + 2);
        let a1r = (a1 % two_mp2) as u64;
        let a2r = (a2 % two_mp2) as u64;
        let (s1c, s2c) = (2 * a1r, 2 * a2r);

        for (hc_idx, hc) in vp["half_chains"].as_array().unwrap().iter().enumerate() {
            let dual = hc["dual"].as_bool().unwrap();
            assert_eq!(dual, hc_idx == 1, "vec {vi}: F2_dual expected at index 1");

            // M1 = gluing matrix; M10 = M0·M1 mod 4; N_dim2 base-changes theta.
            let m1 = if !dual {
                gluing_dim2_f1(a1, a2, q4)
            } else {
                gluing_dim2_f2(a1, a2, q4)
            };
            let m10 = mat4_mul_mod4(&m0, &m1);
            let n_dim2 = base_change_theta_dim2(&m10, &e4);

            // B_K_dim2 (kernel of the 2^m-isogeny chain in dimension 2).
            let s1p = jac_mul(&p1_m, s1c, &e_com);
            let s2p = jac_mul(&p1_m, s2c, &e_com);
            let s1q = jac_mul(&q1_m, s1c, &e_com);
            let s2q = jac_mul(&q1_m, s2c, &e_com);
            let two_r2 = jac_mul(&r2_m, 2, &e_chal);
            let two_s2 = jac_mul(&s2_m, 2, &e_chal);

            let (tp0, tp1): (TuplePoint<Level1>, TuplePoint<Level1>) = if !dual {
                (
                    TuplePoint::new(jac_sub(&s1p, &s2q, &e_com), two_r2.clone()),
                    TuplePoint::new(jac_add(&s1q, &s2p, &e_com), two_s2.clone()),
                )
            } else {
                (
                    TuplePoint::new(jac_add(&s1p, &s2q, &e_com), two_r2.neg()),
                    TuplePoint::new(jac_sub(&s1q, &s2p, &e_com), two_s2.neg()),
                )
            };

            let chain = IsogenyChainDim2::new(&tp0, &tp1, &zero12, &n_dim2, m, &e_com, &e_chal)
                .expect("oracle vector: dim-2 gluing chain is well-formed");

            let want: [F; 4] = core::array::from_fn(|i| parse_fp2(&hc["dim2_codomain_null"][i]));
            let tag = if dual { "F2_dual" } else { "F1" };
            assert!(
                proj_eq4(chain.codomain_null(), &want),
                "vec {vi} {tag} (m={m}): dim2 chain codomain null mismatch"
            );
        }

        n += 1;
    }
    assert_eq!(n, 5);
    println!(
        "dim-2 (2,2)-isogeny chain reproduced the oracle's dim2_codomain_null \
         for all {n} vectors, both half-chains (m ∈ {{1,2,3,6}})"
    );
}
