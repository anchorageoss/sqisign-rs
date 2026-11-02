//! The two base changes of the dimension-4 gluing in their fast forms
//! against the dense products (SECENG-2437): the sparse dimension-4 base
//! change (`BaseChangeDim4`) and the folded dimension-2 base change
//! (`FoldedBaseChangeDim2`) give the dense product's field elements on the
//! matrices of the five level-I vectors, on random matrices of each detected
//! shape, for both fourth roots of unity, and the dense fall-back is taken
//! on a matrix of another shape; the inverses up to a scalar undo the base
//! changes on the vectors' matrices and are absent where the rows are not
//! orthogonal. Feature `compact`.
#![cfg(feature = "compact")]

use sqisign_rs::compact::{compact_keygen, compact_public_key_to_bytes, compact_sign};
use sqisign_rs::mp::ShakeRng;
use sqisign_rs::sqisign::level1;
use sqisign_verify::hd::dim2::{montgomery_to_theta_matrix_dim2, FoldedBaseChangeDim2};
use sqisign_verify::hd::params::HdLevel;
use sqisign_verify::hd::{
    apply_base_change_theta_dim4, apply_mat4, base_change_theta_dim2, base_change_theta_dim4,
    challenge_limbs, gluing_bc_dim4_f1, gluing_bc_dim4_f2, gluing_dim2_f1, gluing_dim2_f2,
    hd_kernel_setup, parse_public_key, parse_signature, starting_two_symplectic_matrices,
    BaseChangeDim4, HdSignature,
};
use sqisign_verify::{Fp, Fp2, P324_3};

type F = Fp2<P324_3>;
const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
const PK_BYTES: usize = 84;
const SIG_BYTES: usize = 142;

// A pseudo-random field element from a seed (squarings and additions).
fn elt(seed: u64) -> F {
    let mut x = F::from_small(seed | 1);
    let s = F::from_small(seed ^ 0x5555);
    for _ in 0..40 {
        x = x.sqr().add(&s);
    }
    x
}

fn same(a: &F, b: &F) -> bool {
    a.encode() == b.encode()
}

fn roots() -> [F; 2] {
    let i = F {
        re: Fp::<P324_3>::zero(),
        im: Fp::<P324_3>::one(),
    };
    [i.clone(), i.neg()]
}

fn pow(e4: &F, k: usize) -> F {
    let mut z = F::one();
    for _ in 0..k {
        z = z.mul(e4);
    }
    z
}

// A small deterministic generator for shapes and exponents.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> usize {
        (self.next() % n) as usize
    }
}

// `a` and `b` are the same projective point: a[j] b[0] = a[0] b[j] for all j.
fn proportional<const K: usize>(a: &[F; K], b: &[F; K]) -> bool {
    (0..K).all(|j| same(&a[j].mul(&b[0]), &a[0].mul(&b[j])))
}

fn check_dim4(n: &[[F; 16]; 16], e4: &F, sparse_expected: bool, seed: u64, label: &str) {
    let bc = BaseChangeDim4::new(n.clone(), e4);
    assert_eq!(bc.is_sparse(), sparse_expected, "{label}: shape detection");
    for t in 0..4u64 {
        let p: [F; 16] = core::array::from_fn(|j| elt(seed + 100 * t + j as u64));
        let fast = bc.apply(&p);
        let dense = apply_base_change_theta_dim4(n, &p);
        for j in 0..16 {
            assert!(
                same(&fast[j], &dense[j]),
                "{label}: point {t} coordinate {j}"
            );
        }
        // the inverse up to a scalar, where the form has one, undoes the base change
        if let Some(back) = bc.apply_inverse_scaled(&fast) {
            assert!(proportional(&back, &p), "{label}: inverse of point {t}");
        }
    }
}

// The base change of a theta matrix has the inverse form.
fn check_dim4_inverse_present(n: &[[F; 16]; 16], e4: &F, label: &str) {
    let bc = BaseChangeDim4::new(n.clone(), e4);
    let p: [F; 16] = core::array::from_fn(|j| elt(90_000 + j as u64));
    assert!(
        bc.apply_inverse_scaled(&p).is_some(),
        "{label}: no inverse form"
    );
}

fn check_dim2_inverse_present(zero12: &[F; 4], n: &[[F; 4]; 4], label: &str) {
    let f = FoldedBaseChangeDim2::new(zero12, n).expect(label);
    let v: [F; 4] = core::array::from_fn(|j| elt(91_000 + j as u64));
    assert!(
        f.apply_inverse_scaled(&v).is_some(),
        "{label}: no inverse form"
    );
}

fn check_dim2(zero12: &[F; 4], n: &[[F; 4]; 4], folded_expected: bool, seed: u64, label: &str) {
    let folded = FoldedBaseChangeDim2::new(zero12, n);
    assert_eq!(
        folded.is_some(),
        folded_expected,
        "{label}: shape detection"
    );
    let (m2, _) = montgomery_to_theta_matrix_dim2(zero12, n);
    if let Some(f) = folded {
        for t in 0..4u64 {
            let v: [F; 4] = core::array::from_fn(|j| elt(seed + 100 * t + j as u64));
            let fast = f.apply(&v);
            let dense = apply_mat4(&m2, &v);
            for j in 0..4 {
                assert!(
                    same(&fast[j], &dense[j]),
                    "{label}: point {t} coordinate {j}"
                );
            }
            if let Some(back) = f.apply_inverse_scaled(&fast) {
                assert!(proportional(&back, &v), "{label}: inverse of point {t}");
            }
        }
    }
}

// The dimension-4 shape: four nonzero entries per row, each 4 e4^k.
fn random_dim4(g: &mut Lcg, e4: &F) -> [[F; 16]; 16] {
    let four = F::one().mul_small(4);
    core::array::from_fn(|_| {
        let mut row: [F; 16] = core::array::from_fn(|_| F::zero());
        let mut cols = [0usize; 4];
        let mut filled = 0;
        while filled < 4 {
            let c = g.below(16);
            if !cols[..filled].contains(&c) {
                cols[filled] = c;
                filled += 1;
            }
        }
        for &c in cols.iter() {
            row[c] = four.mul(&pow(e4, g.below(4)));
        }
        row
    })
}

// The dimension-2 shapes: four entries per row, each e4^k (`two_per_row`
// false), or two entries per row, each 2 e4^k (`two_per_row` true).
fn random_dim2(g: &mut Lcg, e4: &F, two_per_row: bool) -> [[F; 4]; 4] {
    core::array::from_fn(|_| {
        let mut row: [F; 4] = core::array::from_fn(|_| F::zero());
        if two_per_row {
            let c0 = g.below(4);
            let c1 = (c0 + 1 + g.below(3)) % 4;
            let two = F::one().mul_small(2);
            row[c0] = two.mul(&pow(e4, g.below(4)));
            row[c1] = two.mul(&pow(e4, g.below(4)));
        } else {
            for r in row.iter_mut() {
                *r = pow(e4, g.below(4));
            }
        }
        row
    })
}

#[test]
fn random_matrices_of_each_shape_and_the_fall_back() {
    let mut g = Lcg(2437);
    for (r, e4) in roots().iter().enumerate() {
        for t in 0..8u64 {
            let n4 = random_dim4(&mut g, e4);
            check_dim4(
                &n4,
                e4,
                true,
                10_000 + 1000 * r as u64 + 50 * t,
                &format!("dim4 root {r} matrix {t}"),
            );
            let zero12: [F; 4] =
                core::array::from_fn(|j| elt(20_000 + 1000 * r as u64 + 50 * t + j as u64));
            let n2a = random_dim2(&mut g, e4, false);
            check_dim2(
                &zero12,
                &n2a,
                true,
                30_000 + 1000 * r as u64 + 50 * t,
                &format!("dim2 four-per-row root {r} matrix {t}"),
            );
            let n2b = random_dim2(&mut g, e4, true);
            check_dim2(
                &zero12,
                &n2b,
                true,
                40_000 + 1000 * r as u64 + 50 * t,
                &format!("dim2 two-per-row root {r} matrix {t}"),
            );
        }
        // a matrix of the shape whose rows are not orthogonal has no inverse form
        let four = F::one().mul_small(4);
        let mut n4u: [[F; 16]; 16] = core::array::from_fn(|_| core::array::from_fn(|_| F::zero()));
        for (i, row) in n4u.iter_mut().enumerate() {
            for c in 0..4 {
                row[(i + c) % 16] = four.clone();
            }
        }
        let bcu = BaseChangeDim4::new(n4u.clone(), e4);
        assert!(
            bcu.is_sparse(),
            "dim4 root {r}: the shifted-ones matrix has the shape"
        );
        let p: [F; 16] = core::array::from_fn(|j| elt(80_000 + j as u64));
        assert!(
            bcu.apply_inverse_scaled(&p).is_none(),
            "dim4 root {r}: no inverse form expected"
        );
        let n2u: [[F; 4]; 4] = core::array::from_fn(|_| core::array::from_fn(|_| F::one()));
        let zero12u: [F; 4] = core::array::from_fn(|j| elt(81_000 + j as u64));
        let fu = FoldedBaseChangeDim2::new(&zero12u, &n2u).expect("the all-ones matrix folds");
        let v: [F; 4] = core::array::from_fn(|j| elt(82_000 + j as u64));
        assert!(
            fu.apply_inverse_scaled(&v).is_none(),
            "dim2 root {r}: no inverse form expected"
        );
        // another shape: one entry 3 (not 1, 2 or 4 times a root) in each matrix
        let mut n4 = random_dim4(&mut g, e4);
        n4[5][7] = F::one().mul_small(3);
        check_dim4(
            &n4,
            e4,
            false,
            50_000 + 1000 * r as u64,
            &format!("dim4 fall-back root {r}"),
        );
        let mut n4b = random_dim4(&mut g, e4);
        n4b[0][1] = elt(77);
        check_dim4(
            &n4b,
            e4,
            false,
            51_000 + 1000 * r as u64,
            &format!("dim4 fall-back random entry root {r}"),
        );
        let zero12: [F; 4] = core::array::from_fn(|j| elt(60_000 + j as u64));
        let mut n2 = random_dim2(&mut g, e4, false);
        n2[2][2] = F::one().mul_small(3);
        check_dim2(
            &zero12,
            &n2,
            false,
            61_000,
            &format!("dim2 fall-back root {r}"),
        );
        let mut n2b = random_dim2(&mut g, e4, false);
        n2b[1][3] = elt(78);
        check_dim2(
            &zero12,
            &n2b,
            false,
            62_000,
            &format!("dim2 fall-back random entry root {r}"),
        );
    }
}

// The matrices of one signature: both halves of the dimension-4 and the
// dimension-2 gluing base changes, as the verifier builds them.
fn check_vector(tag: &str, pkb: &[u8; PK_BYTES], sig: &[u8; SIG_BYTES], msg: &[u8], seed: u64) {
    let s = parse_signature::<P324_3>(sig).unwrap();
    let p = parse_public_key::<P324_3>(pkb).unwrap();
    let chal = challenge_limbs::<P324_3>(&s, &p, msg).unwrap();
    let hdsig = HdSignature {
        a_pk: p.a_pk.clone(),
        a_com: s.a_com.clone(),
        hint_pk_p: p.hint_pk_p,
        hint_pk_q: p.hint_pk_q,
        hint_com_p: s.hint_com_p,
        hint_com_q: s.hint_com_q,
        chal_limbs: &chal,
        resp_a: s.a,
        resp_b: s.b,
        resp_c_or_d: s.c_or_d,
        q: s.q,
    };
    let k = hd_kernel_setup::<P324_3>(&hdsig).unwrap();
    let q4 = s.q.low_u128();
    let mask = (1u128 << <P324_3 as HdLevel>::R) - 1;
    let (m1, m2) = starting_two_symplectic_matrices(k.a1, k.a2, q4, mask).unwrap();
    let mg4 = [
        gluing_bc_dim4_f1(k.a1, k.a2, q4, k.m, &m1),
        gluing_bc_dim4_f2(k.a1, k.a2, q4, k.m, &m2),
    ];
    for (h, m) in mg4.iter().enumerate() {
        let n4 = base_change_theta_dim4::<P324_3>(m, &k.e4);
        check_dim4(
            &n4,
            &k.e4,
            true,
            seed + 10 * h as u64,
            &format!("{tag} dim4 half {}", h + 1),
        );
        check_dim4_inverse_present(&n4, &k.e4, &format!("{tag} dim4 half {}", h + 1));
    }
    let mg2 = [
        gluing_dim2_f1(k.a1, k.a2, q4),
        gluing_dim2_f2(k.a1, k.a2, q4),
    ];
    for (h, m) in mg2.iter().enumerate() {
        let m10: [[i64; 4]; 4] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let acc: u32 = k.m0[i]
                    .iter()
                    .zip(m.iter())
                    .map(|(a, row)| (*a as u32) * (row[j] as u32))
                    .sum();
                (acc % 4) as i64
            })
        });
        let n2 = base_change_theta_dim2::<P324_3>(&m10, &k.e4);
        check_dim2(
            &k.zero12,
            &n2,
            true,
            seed + 20 + 10 * h as u64,
            &format!("{tag} dim2 half {}", h + 1),
        );
        check_dim2_inverse_present(&k.zero12, &n2, &format!("{tag} dim2 half {}", h + 1));
    }
}

#[test]
fn the_vectors_matrices() {
    let params = level1();
    for i in 0..5u64 {
        let label = format!("compact level I vector {i}");
        let mut rng = ShakeRng::new(label.as_bytes(), b"kat");
        let (pk, sk) = compact_keygen::<P324_3, N>(&params, &mut rng).expect("keygen");
        let mut pkb = [0u8; PK_BYTES];
        compact_public_key_to_bytes(&pk, &mut pkb).unwrap();
        let msg = format!("message {i}").into_bytes();
        let mut sig = [0u8; SIG_BYTES];
        compact_sign::<P324_3, N>(&params, &pk, &sk, &msg, &mut rng, &mut sig).expect("sign");
        check_vector(&format!("vector {i}"), &pkb, &sig, &msg, 70_000 + 100 * i);
    }
}
