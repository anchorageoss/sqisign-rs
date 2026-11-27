//! Items 1 and 2 of paper-shortiso: Corollary A of "Shortest isogenies
//! between supersingular elliptic curves with complex multiplication" at
//! the three SQIsign primes, on the crate's own field, curve, basis of
//! `E0[2^f]`, stored endomorphism action, quaternion layer and
//! ideal-to-isogeny translation.
//!
//! Item 1 (curve side). The kernel is `G = theta'(E0[2m])` with
//! `theta' = iota + pi - t m`, `t = D mod 2`, `m = 2^e`. On `E0[m]` this is
//! `(iota + pi)/2 - t m/2`, and `(iota + pi)/2` is the third stored
//! generator of `O0`: the matrix used is `ACTION_GEN3`, nothing else. The
//! isogeny `E0 -> E0/G` is then walked three ways:
//!   (a) 2-isogenies one at a time (x-only formulas written here on the
//!       crate's `Fp2`), which gives every j-invariant of the chain;
//!   (b) the crate's 4-isogeny (`iso_xisog_4`, `iso_xeval_4`) one step at a
//!       time, the kernel generator pushed through each step;
//!   (c) the crate's strategy-based chain `iso_isogeny_2chain_from`, fed
//!       with a point two doublings above the generator as the crate
//!       itself feeds it (it walks `length - 2` steps).
//! The generic module of item 3 is compared with it: same kernel points,
//! same chain of j-invariants.
//!
//! Item 2 (quaternion side). The left ideal `I = O0 theta + O0 m`,
//! `theta = (t m + i + j)/2`, its norm, primitivity, right order, the
//! element `theta / m` of discriminant `-D` in it, and the crate's
//! ideal-to-isogeny on an equivalent ideal of odd norm.

mod shortiso_common;

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use sha3::{Digest, Sha3_256};
use shortiso_common::curve::Curve;
use shortiso_common::field::{Field, F2};
use shortiso_common::{
    check_record, class_polynomial_checked, cor_b, corollary_a_kernel, j_class_number_one, pow2,
    walk,
};
use sqisign_rs::id2iso::{
    e0_basis_even, e0_curve, ideal_to_isogeny_qlapoty, kernel_dlogs_to_ideal_even,
    matrix_application_even_basis, E0Actions, Id2IsoStats,
};
use sqisign_rs::mp::{Ibz, ShakeRng};
use sqisign_rs::quat::integers::generate_random_prime;
use sqisign_rs::quat::{Mat2x2, QuatAlg, QuatAlgElem, QuatIdeal, QuatLattice, Vec2, Vec4};
use sqisign_verify::ec::isogeny::{iso_isogeny_2chain_from, iso_xeval_4, iso_xisog_4};
use sqisign_verify::ec::point::{
    ec_biscalar_mul, ec_dbl_iter_basis, ec_mul, test_point_order_twof, xdbl_a24,
};
use sqisign_verify::ec::{EcCurve, EcPoint, MAX_ORDER_WORDS};
use sqisign_verify::fp::{Fp, Fp2, FpBackend};
use sqisign_verify::params::{Prime, P324_3, P500_27, P664_17};
use sqisign_verify::precomp::PrimePrecomp;
use std::collections::HashSet;

// ---------------------------------------------------------------- glue

fn fp2_to_f2<L: FpBackend>(x: &Fp2<L>) -> F2 {
    let b = x.encode();
    let n = b.len() / 2;
    F2 {
        re: BigUint::from_bytes_le(&b[..n]),
        im: BigUint::from_bytes_le(&b[n..]),
    }
}

fn fp_from_big<L: FpBackend>(v: &BigUint) -> Fp<L> {
    let n = Fp2::<L>::zero().encode().len() / 2;
    let mut bytes = v.to_bytes_le();
    bytes.resize(n, 0);
    Fp::<L>::decode(&bytes).expect("reduced residue")
}

fn fp2_from_int<L: FpBackend>(v: &BigInt, p: &BigUint) -> Fp2<L> {
    let pm = BigInt::from_biguint(Sign::Plus, p.clone());
    Fp2 {
        re: fp_from_big(&v.mod_floor(&pm).to_biguint().unwrap()),
        im: Fp::zero(),
    }
}

fn eq2<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

fn affine_x<L: FpBackend>(p: &EcPoint<L>) -> Fp2<L> {
    p.x.mul(&p.z.inv())
}

fn digits<const N: usize>(x: &Ibz<N>) -> [u64; MAX_ORDER_WORDS] {
    let mut d = [0u64; MAX_ORDER_WORDS];
    if !x.is_zero() {
        x.to_digits(&mut d);
    }
    d
}

fn ibz_to_big<const N: usize>(x: &Ibz<N>) -> BigInt {
    if x.is_zero() {
        return BigInt::zero();
    }
    let neg = x.cmp_i32(0) == std::cmp::Ordering::Less;
    let a = x.abs();
    let mut d = [0u64; 64];
    a.to_digits(&mut d);
    let mut bytes = Vec::new();
    for w in d {
        bytes.extend_from_slice(&w.to_le_bytes());
    }
    let v = BigInt::from_bytes_le(Sign::Plus, &bytes);
    if neg {
        -v
    } else {
        v
    }
}

fn big_to_ibz<const N: usize>(x: &BigInt) -> Ibz<N> {
    let v = Ibz::<N>::from_le_bytes(&x.magnitude().to_bytes_le());
    if x.is_negative() {
        v.neg()
    } else {
        v
    }
}

// ------------------------------------------------- item 1: curve side

/// `(A + 2C : 4C)` of the 2-isogeny with kernel `<T>`, `T = (x : z)` of
/// order 2, `x != 0`.
fn two_isogeny_codomain<L: FpBackend>(t: &EcPoint<L>) -> EcPoint<L> {
    let x2 = t.x.sqr();
    let z2 = t.z.sqr();
    EcPoint::new(z2.sub(&x2), z2)
}

/// The image of `Q` under that 2-isogeny.
fn two_isogeny_eval<L: FpBackend>(q: &EcPoint<L>, t: &EcPoint<L>) -> EcPoint<L> {
    let t0 = t.x.add(&t.z).mul(&q.x.sub(&q.z));
    let t1 = t.x.sub(&t.z).mul(&q.x.add(&q.z));
    EcPoint::new(q.x.mul(&t0.add(&t1)), q.z.mul(&t0.sub(&t1)))
}

fn dbl_n<L: FpBackend>(p: &EcPoint<L>, a24: &EcPoint<L>, n: u32) -> EcPoint<L> {
    let mut r = p.clone();
    for _ in 0..n {
        r = xdbl_a24(&r, a24, false);
    }
    r
}

fn j_of_a24<L: FpBackend>(a24: &EcPoint<L>) -> Fp2<L> {
    EcCurve::from_a24(a24).j_inv()
}

struct CurveSide<L: Prime, const N: usize> {
    /// The kernel generator (the image of `P'` or `Q'`).
    kernel: EcPoint<L>,
    /// Its coordinates `(a, b)` on the basis of `E0[m]`: `K = [a] P' + [b] Q'`.
    coords: Vec2<N>,
    /// The images of the basis `(P', Q')` of `E0[m]`.
    images: (EcPoint<L>, EcPoint<L>),
    /// `j_0, ..., j_e`.
    js: Vec<Fp2<L>>,
    which: &'static str,
}

fn curve_side<L: FpBackend + PrimePrecomp, const N: usize>(
    name: &str,
    act: &E0Actions<N>,
    d: u64,
    e: u32,
) -> CurveSide<L, N> {
    let f = L::TWO_ADIC_EXPONENT;
    assert_eq!(cor_b(L::COFACTOR, f), (d, e), "{name}: Corollary B");
    let t = d % 2;
    let m = Ibz::<N>::one().mul_2exp(e);
    let w = L::ORDER_WORDS;

    // theta' / 2 = (iota + pi)/2 - t m/2 on E0[m]: ACTION_GEN3 - t (m/2) Id
    let mut mat: Mat2x2<N> = act.action_gen3;
    if t == 1 {
        let half = Ibz::<N>::one().mul_2exp(e - 1);
        mat.0[0][0] = mat.0[0][0].sub(&half);
        mat.0[1][1] = mat.0[1][1].sub(&half);
    }
    for row in mat.0.iter_mut() {
        for x in row.iter_mut() {
            *x = x.modulo(&m);
        }
    }
    // integers first: the image of (Z/m)^2 is cyclic of order m exactly
    // when the determinant vanishes mod m and an entry is odd
    let det = mat.0[0][0]
        .mul(&mat.0[1][1])
        .sub(&mat.0[0][1].mul(&mat.0[1][0]))
        .modulo(&m);
    assert!(
        det.is_zero(),
        "{name}: det of the action on E0[m] is 0 mod m"
    );
    let col_odd = |c: usize| mat.0[0][c].is_odd() || mat.0[1][c].is_odd();
    assert!(
        col_odd(0) || col_odd(1),
        "{name}: the action on E0[m] has an odd entry"
    );

    // the images of the basis of E0[m] = [2^(f-e)] E0[2^f]
    let mut e0 = e0_curve::<L>();
    let b0 = e0_basis_even::<L>();
    let mut bas = ec_dbl_iter_basis(&b0, (f - e) as usize, &mut e0);
    assert!(
        test_point_order_twof(&bas.p, &e0, e as usize)
            && test_point_order_twof(&bas.q, &e0, e as usize),
        "{name}: basis of E0[m]"
    );
    let mut mm = mat;
    matrix_application_even_basis(&mut bas, &e0, &mut mm, e, false)
        .expect("matrix application on E0[m]");
    let (k1, k2) = (bas.p.clone(), bas.q.clone());

    // on the curve: one image has order exactly m, the other is the
    // multiple of it the matrix predicts
    let (gen_col, which) = if col_odd(0) {
        (0, "the image of P")
    } else {
        (1, "the image of Q")
    };
    let (kernel, other) = if gen_col == 0 {
        (k1.clone(), k2.clone())
    } else {
        (k2.clone(), k1.clone())
    };
    assert!(
        test_point_order_twof(&kernel, &e0, e as usize),
        "{name}: the generator of G has order exactly m"
    );
    let row = if mat.0[0][gen_col].is_odd() { 0 } else { 1 };
    let inv = mat.0[row][gen_col].invmod(&m).expect("odd entry");
    let lambda = mat.0[row][1 - gen_col].mul(&inv).modulo(&m);
    for r in 0..2 {
        assert_eq!(
            mat.0[r][gen_col].mul(&lambda).modulo(&m),
            mat.0[r][1 - gen_col],
            "{name}: the second column is lambda times the first"
        );
    }
    let lam_k = ec_mul(&kernel, &digits(&lambda)[..w], e as usize, &mut e0);
    if lambda.is_zero() {
        assert!(bool::from(other.is_zero()), "{name}: second image is zero");
    } else {
        assert!(
            bool::from(lam_k.ct_equal(&other)),
            "{name}: the second image is [lambda] the generator: G is cyclic of order m"
        );
    }

    // (a) 2-isogenies one at a time
    let a24_e0 = EcPoint::<L>::new(Fp2::from_small(2), Fp2::from_small(4));
    let mut a24 = a24_e0.clone();
    let mut k = kernel.clone();
    let mut js = vec![j_of_a24(&a24)];
    assert!(eq2(&js[0], &Fp2::from_small(1728)));
    for step in 0..e {
        let t2 = dbl_n(&k, &a24, e - 1 - step);
        assert!(
            !bool::from(t2.is_zero())
                && bool::from(xdbl_a24(&t2, &a24, false).is_zero())
                && !bool::from(t2.x.ct_is_zero()),
            "{name}: step {step}: a kernel point of order 2 other than (0, 0)"
        );
        let next = two_isogeny_codomain(&t2);
        k = two_isogeny_eval(&k, &t2);
        a24 = next;
        js.push(j_of_a24(&a24));
    }
    assert!(bool::from(k.is_zero()), "{name}: kernel exhausted");

    // (b) the crate's 4-isogeny, one step at a time (an odd e starts with
    // one 2-isogeny)
    let mut a24 = a24_e0.clone();
    let mut k = kernel.clone();
    let mut done = 0u32;
    if e % 2 == 1 {
        let t2 = dbl_n(&k, &a24, e - 1);
        let next = two_isogeny_codomain(&t2);
        k = two_isogeny_eval(&k, &t2);
        a24 = next;
        done = 1;
    }
    while done < e {
        let t4 = dbl_n(&k, &a24, e - done - 2);
        let (kps, next) = iso_xisog_4(&t4);
        k = iso_xeval_4(&k, &kps);
        a24 = next;
        done += 2;
        assert!(
            eq2(&j_of_a24(&a24), &js[done as usize]),
            "{name}: 4-isogeny walk agrees with the 2-isogeny walk at step {done}"
        );
    }
    assert!(bool::from(k.is_zero()));

    // (c) the crate's strategy-based chain. As the crate calls it
    // (`EXTRA_TORSION`), `iso_isogeny_2chain_from` with `length` takes a
    // point of order `2^length` and walks `length - 2` steps: the last
    // point of order 4 only normalises the model. So it is given a point
    // K~ with [4] K~ = K and length e + 2 (at an odd e: one 2-isogeny
    // first, then length e + 1).
    let mut e0c = e0_curve::<L>();
    let above = ec_dbl_iter_basis(&b0, (f - e - 2) as usize, &mut e0c);
    let k_above = ec_biscalar_mul(
        &digits(&mat.0[0][gen_col])[..w],
        &digits(&mat.0[1][gen_col])[..w],
        e as usize + 2,
        &above,
        &e0c,
    )
    .expect("a point above the generator");
    assert!(
        bool::from(dbl_n(&k_above, &a24_e0, 2).ct_equal(&kernel)),
        "{name}: [4] K~ = K"
    );
    let cod = if e % 2 == 0 {
        iso_isogeny_2chain_from(&e0_curve::<L>(), &k_above, e + 2)
    } else {
        let t2 = dbl_n(&k_above, &a24_e0, e + 1);
        let a24_1 = two_isogeny_codomain(&t2);
        let k1_above = two_isogeny_eval(&k_above, &t2);
        iso_isogeny_2chain_from(&EcCurve::from_a24(&a24_1), &k1_above, e + 1)
    }
    .expect("strategy-based chain");
    assert!(
        eq2(&cod.j_inv(), &js[e as usize]),
        "{name}: the strategy-based chain ends on the same j"
    );

    CurveSide {
        kernel,
        coords: Vec2([mat.0[0][gen_col], mat.0[1][gen_col]]),
        images: (k1, k2),
        js,
        which,
    }
}

/// The assertions on the chain of j-invariants.
fn chain_assertions<L: FpBackend>(name: &str, d: u64, e: u32, js: &[Fp2<L>]) -> String {
    let p = BigUint::from_bytes_le(L::prime_le_bytes());
    assert_eq!(js.len(), e as usize + 1);
    for (i, j) in js.iter().enumerate() {
        assert!(
            bool::from(j.im.ct_is_zero()),
            "{name}: j at step {i} has a non-zero imaginary part"
        );
    }
    let set: HashSet<Vec<u8>> = js.iter().map(|j| j.encode().to_vec()).collect();
    assert_eq!(set.len(), js.len(), "{name}: a j-invariant repeats");
    let j_end = &js[e as usize];
    match d {
        3 => {
            assert!(eq2(j_end, &Fp2::zero()), "{name}: final j = 0");
            assert!(
                eq2(&js[e as usize - 1], &Fp2::from_small(54000)),
                "{name}: j = 54000 at step {}",
                e - 1
            );
            format!("j = 0 after {e} steps; j = 54000 at step {}", e - 1)
        }
        68 => {
            let h = class_polynomial_checked(68);
            let mut acc = Fp2::<L>::zero();
            for c in h.iter().rev() {
                acc = acc.mul(j_end).add(&fp2_from_int::<L>(c, &p));
            }
            assert!(
                bool::from(acc.ct_is_zero()),
                "{name}: final j is a root of H_-68"
            );
            format!("root of H_-68 after {e} steps")
        }
        _ => {
            let j = j_class_number_one(d).expect("class number one");
            assert!(
                eq2(j_end, &fp2_from_int::<L>(&j, &p)),
                "{name}: final j = {j}"
            );
            format!("j = {j} after {e} steps")
        }
    }
}

/// The generic module against the crate: the same kernel points from
/// `iota` and `pi` evaluated on points, the same chain of j-invariants.
fn generic_cross_check<L: FpBackend + PrimePrecomp, const N: usize>(
    name: &str,
    d: u64,
    e: u32,
    side: &CurveSide<L, N>,
) {
    let f = L::TWO_ADIC_EXPONENT;
    let p = BigUint::from_bytes_le(L::prime_le_bytes());
    let fld = Field::new(&p);
    let cur = Curve::e0(&fld);
    let b0 = e0_basis_even::<L>();
    let px = fp2_to_f2(&affine_x(&b0.p));
    let qx = fp2_to_f2(&affine_x(&b0.q));
    let pmqx = fp2_to_f2(&affine_x(&b0.pmq));
    let pp = cur.lift_x(&fld, &px).expect("P0 on E0");
    let mut qq = cur.lift_x(&fld, &qx).expect("Q0 on E0");
    let x_of = |pt: &shortiso_common::curve::Pt| cur.to_affine(&fld, pt).unwrap().0;
    if x_of(&cur.sub(&fld, &pp, &qq)) != pmqx {
        qq = cur.neg(&fld, &qq);
    }
    assert_eq!(x_of(&cur.sub(&fld, &pp, &qq)), pmqx, "{name}: x(P0 - Q0)");
    // a basis of E0[2m] above the crate's basis of E0[m]
    let rp = cur.dbl_iter(&fld, &pp, (f - e - 1) as u64);
    let rq = cur.dbl_iter(&fld, &qq, (f - e - 1) as u64);
    let m = pow2(e);
    let ker = corollary_a_kernel(&fld, d, &m, &[(2, e)], &(rp, rq));
    assert_eq!(
        x_of(&ker.images.0),
        fp2_to_f2(&affine_x(&side.images.0)),
        "{name}: theta'(R_P) from points equals the stored matrix on P'"
    );
    assert_eq!(
        x_of(&ker.images.1),
        fp2_to_f2(&affine_x(&side.images.1)),
        "{name}: theta'(R_Q) from points equals the stored matrix on Q'"
    );
    let (js, _) = walk(&fld, &ker.gen, &vec![2; e as usize]);
    let crate_js: Vec<F2> = side.js.iter().map(fp2_to_f2).collect();
    assert_eq!(js, crate_js, "{name}: the two chains of j-invariants");
}

// --------------------------------------------- item 2: quaternion side

type Row = [BigInt; 4];

/// The lattice generated by integer vectors over a common denominator, in
/// upper-triangular Hermite normal form (columns `v0..v3`, `v_k` zero
/// below row `k`), as `(denominator, columns)` with the content removed.
fn hnf(gens: &[(BigInt, Row)]) -> (BigInt, [Row; 4]) {
    let denom = gens
        .iter()
        .fold(BigInt::one(), |acc, (d, _)| acc.lcm(&d.abs()));
    let mut vecs: Vec<Row> = gens
        .iter()
        .map(|(d, v)| {
            let s = &denom / d;
            [&v[0] * &s, &v[1] * &s, &v[2] * &s, &v[3] * &s]
        })
        .collect();
    let mut basis: [Row; 4] = Default::default();
    for row in (0..4).rev() {
        // Euclid on the entries of this row
        loop {
            vecs.retain(|v| v.iter().any(|x| !x.is_zero()));
            let nz: Vec<usize> = (0..vecs.len())
                .filter(|&i| !vecs[i][row].is_zero())
                .collect();
            assert!(!nz.is_empty(), "full rank");
            if nz.len() == 1 {
                let mut piv = vecs.swap_remove(nz[0]);
                if piv[row].is_negative() {
                    for x in piv.iter_mut() {
                        *x = -&*x;
                    }
                }
                basis[row] = piv;
                break;
            }
            let min = *nz.iter().min_by_key(|&&i| vecs[i][row].abs()).unwrap();
            let piv = vecs[min].clone();
            for &i in &nz {
                if i != min {
                    let q = vecs[i][row].div_floor(&piv[row]);
                    for c in 0..4 {
                        let s = &piv[c] * &q;
                        vecs[i][c] -= s;
                    }
                }
            }
        }
        assert!(vecs.iter().all(|v| v[row].is_zero()));
    }
    // reduce the entries above each pivot
    for col in 1..4 {
        for row in (0..col).rev() {
            let q = basis[col][row].div_floor(&basis[row][row]);
            if !q.is_zero() {
                let piv = basis[row].clone();
                for c in 0..4 {
                    basis[col][c] -= &piv[c] * &q;
                }
            }
        }
    }
    let mut g = denom.clone();
    for v in &basis {
        for x in v {
            g = g.gcd(x);
        }
    }
    for v in basis.iter_mut() {
        for x in v.iter_mut() {
            *x = &*x / &g;
        }
    }
    (denom / g, basis)
}

fn elem_to_big<const N: usize>(x: &QuatAlgElem<N>) -> (BigInt, Row) {
    (
        ibz_to_big(&x.denom),
        [
            ibz_to_big(&x.coord.0[0]),
            ibz_to_big(&x.coord.0[1]),
            ibz_to_big(&x.coord.0[2]),
            ibz_to_big(&x.coord.0[3]),
        ],
    )
}

fn lattice_from_big<const N: usize>(denom: &BigInt, cols: &[Row; 4]) -> QuatLattice<N> {
    let mut lat = QuatLattice::<N>::init();
    lat.denom = big_to_ibz(denom);
    for (c, col) in cols.iter().enumerate() {
        for r in 0..4 {
            lat.basis.0[r][c] = big_to_ibz(&col[r]);
        }
    }
    lat
}

fn column<const N: usize>(lat: &QuatLattice<N>, c: usize) -> QuatAlgElem<N> {
    QuatAlgElem {
        denom: lat.denom,
        coord: Vec4([
            lat.basis.0[0][c],
            lat.basis.0[1][c],
            lat.basis.0[2][c],
            lat.basis.0[3][c],
        ]),
    }
}

/// `det(basis) / denom^4` as a reduced fraction (the bases are triangular).
fn covolume<const N: usize>(lat: &QuatLattice<N>) -> (BigInt, BigInt) {
    let num = (0..4).fold(BigInt::one(), |acc, i| acc * ibz_to_big(&lat.basis.0[i][i]));
    let den = ibz_to_big(&lat.denom).pow(4);
    let g = num.gcd(&den);
    (num / &g, den / g)
}

fn scaled<const N: usize>(x: &QuatAlgElem<N>, denom_factor: &Ibz<N>) -> QuatAlgElem<N> {
    let mut y = *x;
    y.denom = y.denom.mul(denom_factor);
    y
}

struct QuatSide<const N: usize> {
    ideal: QuatIdeal<N>,
    right_order: QuatLattice<N>,
    record: String,
}

fn quaternion_side<const N: usize>(name: &str, alg: &QuatAlg<N>, d: u64, e: u32) -> QuatSide<N> {
    let t = (d % 2) as i32;
    let m = Ibz::<N>::one().mul_2exp(e);
    let m_big = BigInt::one() << e;
    let p_big = ibz_to_big(&alg.p);
    assert_eq!(
        BigInt::from(d) * &m_big * &m_big - 1,
        p_big,
        "p = D m^2 - 1"
    );

    // theta = (t m + i + j) / 2
    let mut theta = QuatAlgElem::<N>::set(2, 0, 1, 1, 0);
    if t == 1 {
        theta.coord.0[0] = m;
    }
    assert!(alg.o0.contains(&theta).0, "{name}: theta in O0");
    let (tn, td) = theta.norm(alg);
    assert!(td.is_one());
    assert_eq!(
        ibz_to_big(&tn),
        &m_big * &m_big * BigInt::from(d + t as u64) / 4,
        "{name}: n(theta) = m^2 (D + t) / 4"
    );

    // I = O0 theta + O0 m by the crate, and independently from generators
    let ideal = QuatIdeal::create_o0_pow_two(&theta, &m, alg)
        .expect("the ideal has an inert Hermite form: it is not inside O0 (1 + i)");
    assert_eq!(ideal.norm, m, "{name}: the crate's ideal has norm m");
    let lat = ideal.to_lattice();
    let o0_basis: Vec<QuatAlgElem<N>> = (0..4).map(|c| column(&alg.o0, c)).collect();
    let m_elem = QuatAlgElem::<N>::scalar(&m, &Ibz::one());
    let mut gens = Vec::new();
    for b in &o0_basis {
        for g in [&theta, &m_elem] {
            let prod = b.mul(g, alg);
            assert!(
                lat.contains(&prod).0,
                "{name}: generator in the crate's ideal"
            );
            gens.push(elem_to_big(&prod));
        }
    }
    let (gd, gcols) = hnf(&gens);
    let gen_lat = lattice_from_big::<N>(&gd, &gcols);
    assert_eq!(
        covolume(&gen_lat),
        covolume(&lat),
        "{name}: O0 theta + O0 m is the crate's ideal (same index, one inside the other)"
    );
    assert!(gen_lat.equals(&lat), "{name}: the two lattices are equal");

    // a left O0-ideal of index m^2 in O0, that is of norm m
    let cols: Vec<QuatAlgElem<N>> = (0..4).map(|c| column(&lat, c)).collect();
    for b in &o0_basis {
        for c in &cols {
            assert!(lat.contains(&b.mul(c, alg)).0, "{name}: O0 I inside I");
        }
    }
    let (o_num, o_den) = covolume(&alg.o0);
    let (i_num, i_den) = covolume(&lat);
    assert_eq!(
        &i_num * &o_den,
        &o_num * &i_den * &m_big * &m_big,
        "{name}: [O0 : I] = m^2"
    );
    // the norm as the gcd of the norms of I: n(c_k) and n(c_k + c_l)
    let mut g = BigInt::zero();
    for a in 0..4 {
        for b in a..4 {
            let x = if a == b {
                cols[a]
            } else {
                cols[a].add(&cols[b])
            };
            let (n, dn) = x.norm(alg);
            assert!(dn.is_one());
            g = g.gcd(&ibz_to_big(&n));
        }
    }
    assert_eq!(g, m_big, "{name}: gcd of the norms of I is m");
    // cyclic: not inside l O0 for any prime l; the norm is a power of two,
    // so l = 2 is the only case
    let mut two_o0 = alg.o0;
    two_o0.denom = Ibz::one();
    assert!(
        !lat.inclusion(&two_o0),
        "{name}: I is not contained in 2 O0"
    );

    // right order: the lattice of conj(c_k) c_l / m
    let mut rgens = Vec::new();
    for a in &cols {
        for b in &cols {
            rgens.push(elem_to_big(&scaled(&a.conj().mul(b, alg), &m)));
        }
    }
    let (rd, rcols) = hnf(&rgens);
    let ro = lattice_from_big::<N>(&rd, &rcols);
    let rbasis: Vec<QuatAlgElem<N>> = (0..4).map(|c| column(&ro, c)).collect();
    assert!(
        ro.contains(&QuatAlgElem::set(1, 1, 0, 0, 0)).0,
        "{name}: 1 in O_R"
    );
    for a in &rbasis {
        for b in &rbasis {
            assert!(ro.contains(&a.mul(b, alg)).0, "{name}: O_R is a ring");
        }
        for c in &cols {
            assert!(lat.contains(&c.mul(a, alg)).0, "{name}: I O_R inside I");
        }
    }
    assert_eq!(
        covolume(&ro),
        covolume(&alg.o0),
        "{name}: O_R has the covolume of a maximal order"
    );

    // omega = theta / m in O_R: trace t, norm (D + t) / 4, discriminant -D
    let omega = scaled(&theta, &m);
    assert!(ro.contains(&omega).0, "{name}: theta / m in O_R");
    let (on, od) = omega.norm(alg);
    let trace_num = ibz_to_big(&omega.coord.0[0]) * 2;
    let trace_den = ibz_to_big(&omega.denom);
    assert_eq!(trace_num, BigInt::from(t) * &trace_den, "{name}: trace t");
    assert!(od.is_one() && ibz_to_big(&on) == BigInt::from((d + t as u64) / 4));
    let disc = BigInt::from(t * t) - 4 * ibz_to_big(&on);
    assert_eq!(disc, -BigInt::from(d), "{name}: discriminant -D");
    // primitive: no order of Q(omega) strictly containing Z[omega] is in
    // O_R. Such an order contains (t' + sqrt(-D)/l)/2 for a prime l with
    // -D/l^2 a discriminant, t' = D/l^2 mod 2, and sqrt(-D) = (i + j)/m.
    let mut over = Vec::new();
    for l in 2..=d {
        if d % (l * l) != 0 || !(2..l).all(|q| l % q != 0) {
            continue;
        }
        let d2 = d / (l * l);
        if d2 % 4 != 0 && d2 % 4 != 3 {
            continue;
        }
        let ml = m.mul(&Ibz::set(l as i64, 32));
        let mut w = QuatAlgElem::<N>::set(2, 0, 1, 1, 0);
        w.denom = ml.mul(&Ibz::two());
        if d2 % 2 == 1 {
            w.coord.0[0] = ml;
        }
        assert!(
            !ro.contains(&w).0,
            "{name}: the order of discriminant -{d2} is not in O_R"
        );
        over.push(format!("-{d2} (index {l})"));
    }

    let two_adic: u64 = (&p_big + 1u32).trailing_zeros().unwrap();
    let mut record = format!(
        "# Right order O_R(I) of I = O0*theta + O0*m, theta = (t*m + i + j)/2, at p = {}*2^{} - 1.\n\
         # It is the endomorphism ring of the endpoint E0/G of Corollary A (D = {d}, m = 2^{e}, t = {t}).\n\
         # Basis: four elements (c0 + c1*i + c2*j + c3*k)/denom of B = (-1, -p), i^2 = -1, j^2 = -p, k = i*j;\n\
         # i is iota: (x, y) -> (-x, i*y), j is the p-power Frobenius of E0: y^2 = x^3 + x.\n\
         # omega = theta/m = ({t}*m + i + j)/(2*m) lies in it: trace {t}, norm {}, discriminant -{d}; Z[omega] is\n\
         # optimally embedded (orders checked and found absent: {}).\n\
         # Produced by crates/sqisign-rs/tests/shortiso_sqisign.rs (Hermite normal form over (1, i, j, k)).\n",
        (&p_big + 1u32) >> two_adic,
        two_adic,
        (d + t as u64) / 4,
        if over.is_empty() { "none exists".to_string() } else { over.join(", ") },
    );
    record.push_str(&format!("denom = {rd}\n"));
    for (k, col) in rcols.iter().enumerate() {
        record.push_str(&format!(
            "b{k} = [{}, {}, {}, {}]\n",
            col[0], col[1], col[2], col[3]
        ));
    }
    QuatSide {
        ideal,
        right_order: ro,
        record,
    }
}

/// What the crate's ideal-to-isogeny layer says about `I`.
struct Id2Iso {
    /// Attempts of `ideal_to_isogeny_qlapoty` on equivalent ideals of odd
    /// prime norm, and how many returned a curve.
    direct_attempts: u32,
    direct_successes: u32,
    /// Whether every success gave the j-invariant of the chain.
    direct_j_ok: bool,
    /// Attempt at which the comparison through a second ideal succeeded.
    pushed_attempt: u32,
}

/// The crate's ideal-to-isogeny layer on `I`, three ways.
///
/// (i) `kernel_dlogs_to_ideal_even` (KernelToIdeal) on the coordinates of
/// the generator of `G` must return `I` itself: `E0[I] = G`, so `E0/I` is
/// the endpoint of item 1.
///
/// (ii) `ideal_to_isogeny_qlapoty` on equivalent ideals `J` of odd prime
/// norm (the equivalence is verified). Qlapoty starts from the smallest
/// ideal of the class and needs its norm odd; here the smallest ideal is
/// `I` itself, of norm `m = 2^e` (asserted below with the crate's
/// `shortest_equivalent`: no isogeny from `E0` to the endpoint is shorter
/// than `m`), so the translation is expected to return `None`. Whatever it
/// returns is recorded, and a returned curve must have the j of item 1.
///
/// (iii) `ideal_to_isogeny_qlapoty` on `I ∩ L` for a random ideal `L` of
/// odd prime norm, against the curve side: `E0 / (I ∩ L)` must be
/// `E_L / phi_L(G)`, the kernel `G` pushed through `phi_L: E0 -> E_L` by
/// the basis images the translation of `L` returns, then walked by the
/// crate's strategy-based chain.
fn ideal_to_isogeny_checks<L: FpBackend + PrimePrecomp, const N: usize>(
    name: &str,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
    ideal: &QuatIdeal<N>,
    side: &CurveSide<L, N>,
    e: u32,
) -> Id2Iso {
    let f = L::TWO_ADIC_EXPONENT;
    let w = L::ORDER_WORDS;
    let m = Ibz::<N>::one().mul_2exp(e);
    let j_end = &side.js[e as usize];
    let lat = ideal.to_lattice();

    // (i) kernel to ideal
    let (from_kernel, split) = kernel_dlogs_to_ideal_even(&side.coords, e, alg, act)
        .expect("the generator of G has full order");
    let (sn, sd) = split.norm(alg);
    assert!(
        sn.is_one() && sd.is_one(),
        "{name}: KernelToIdeal: no split part"
    );
    assert!(
        from_kernel.norm == ideal.norm && from_kernel.x == ideal.x && from_kernel.y == ideal.y,
        "{name}: KernelToIdeal(G) = I"
    );

    // (ii) the direct translation
    let mut out = Id2Iso {
        direct_attempts: 0,
        direct_successes: 0,
        direct_j_ok: true,
        pushed_attempt: 0,
    };
    for attempt in 0..8u32 {
        let mut rng = ShakeRng::new(format!("shortiso/{name}/{attempt}").as_bytes(), b"i2i");
        let (gen, equiv) = ideal
            .small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut rng)
            .expect("an equivalent ideal of odd prime norm");
        // equiv = I conj(gen) / n(I): gen in I, n(gen) = n(I) n(equiv), and
        // every c conj(gen) / n(I) lies in equiv
        assert!(equiv.norm.is_odd(), "{name}: equivalent ideal of odd norm");
        assert!(lat.contains(&gen).0, "{name}: gen in I");
        let (gn, gd) = gen.norm(alg);
        assert!(gd.is_one());
        assert_eq!(
            gn,
            ideal.norm.mul(&equiv.norm),
            "{name}: n(gen) = n(I) n(J)"
        );
        let elat = equiv.to_lattice();
        for c in 0..4 {
            let x = scaled(&column(&lat, c).mul(&gen.conj(), alg), &ideal.norm);
            assert!(elat.contains(&x).0, "{name}: I conj(gen) / n(I) inside J");
        }
        // the smallest ideal of the class has norm m
        let (_, _, smallest) = equiv.shortest_equivalent(alg, &mut rng);
        assert_eq!(
            smallest.norm, m,
            "{name}: the smallest equivalent ideal has norm m = 2^{e}"
        );
        out.direct_attempts += 1;
        let mut stats = Id2IsoStats::default();
        if let Some(res) = ideal_to_isogeny_qlapoty::<L, N>(&equiv, alg, act, &mut rng, &mut stats)
        {
            out.direct_successes += 1;
            out.direct_j_ok &= eq2(&res.codomain.j_inv(), j_end);
        }
    }
    assert!(
        out.direct_j_ok,
        "{name}: ideal-to-isogeny returned a curve other than the endpoint"
    );

    // (iii) through a second ideal
    let it = alg.primality_num_iter;
    for attempt in 0..48u32 {
        let mut rng = ShakeRng::new(format!("shortiso/{name}/L/{attempt}").as_bytes(), b"i2i");
        let n = generate_random_prime::<N>(true, alg.p.bitsize() as u32 - 20, it, &mut rng)
            .expect("a prime");
        let l_ideal = QuatIdeal::random_given_prime_norm(&n, alg, &mut rng).expect("an ideal");
        let mut stats = Id2IsoStats::default();
        let Some(phi_l) =
            ideal_to_isogeny_qlapoty::<L, N>(&l_ideal, alg, act, &mut rng, &mut stats)
        else {
            continue;
        };
        // I ∩ L, checked: inside both, norm the product
        let both = QuatIdeal::intersect_o0(&l_ideal, ideal, None);
        let blat = both.to_lattice();
        assert!(
            blat.inclusion(&lat) && blat.inclusion(&l_ideal.to_lattice()),
            "{name}: I ∩ L inside I and L"
        );
        assert_eq!(both.norm, ideal.norm.mul(&l_ideal.norm), "{name}: n(I ∩ L)");
        let Some(quot) = ideal_to_isogeny_qlapoty::<L, N>(&both, alg, act, &mut rng, &mut stats)
        else {
            continue;
        };
        // phi_L(G) on E_L, from the images of the basis of E0[2^f]; a point
        // two doublings above its generator for the strategy-based chain
        let mut el = phi_l.codomain.clone();
        let above = ec_dbl_iter_basis(&phi_l.basis, (f - e - 2) as usize, &mut el);
        let Some(k_above) = ec_biscalar_mul(
            &digits(&side.coords.0[0])[..w],
            &digits(&side.coords.0[1])[..w],
            e as usize + 2,
            &above,
            &el,
        ) else {
            continue;
        };
        let cod = if e % 2 == 0 {
            iso_isogeny_2chain_from(&el, &k_above, e + 2)
        } else {
            el.normalize_a24();
            let a24 = el.a24.clone();
            let t2 = dbl_n(&k_above, &a24, e + 1);
            if bool::from(t2.x.ct_is_zero()) {
                continue; // first kernel point is (0, 0) on this model
            }
            let a24_1 = two_isogeny_codomain(&t2);
            let k1_above = two_isogeny_eval(&k_above, &t2);
            iso_isogeny_2chain_from(&EcCurve::from_a24(&a24_1), &k1_above, e + 1)
        };
        let Some(cod) = cod else {
            continue; // kernel above (0, 0) on this model: another L
        };
        assert!(
            eq2(&cod.j_inv(), &quot.codomain.j_inv()),
            "{name}: E0 / (I ∩ L) is E_L / phi_L(G)"
        );
        assert!(
            !eq2(&cod.j_inv(), &phi_l.codomain.j_inv()),
            "{name}: the pushed chain moved"
        );
        // negative control: another kernel of the same order (first
        // coordinate shifted by 4) does not give that curve
        if e % 2 == 0 {
            let a_bad = side.coords.0[0].add(&Ibz::set(4, 4));
            let wrong = ec_biscalar_mul(
                &digits(&a_bad)[..w],
                &digits(&side.coords.0[1])[..w],
                e as usize + 2,
                &above,
                &el,
            )
            .and_then(|k| iso_isogeny_2chain_from(&el, &k, e + 2));
            if let Some(wrong) = wrong {
                assert!(
                    !eq2(&wrong.j_inv(), &quot.codomain.j_inv()),
                    "{name}: a different kernel gives a different curve"
                );
            }
        }
        out.pushed_attempt = attempt;
        return out;
    }
    panic!("{name}: no ideal L gave a comparison in 48 attempts");
}

fn digest<L: FpBackend>(js: &[Fp2<L>]) -> String {
    let mut h = Sha3_256::new();
    for j in js {
        h.update(j.encode());
    }
    h.finalize()[..8]
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn run_level<L: FpBackend + PrimePrecomp, const N: usize>(
    name: &str,
    file: &str,
    alg: &QuatAlg<N>,
    act: &E0Actions<N>,
    d: u64,
    e: u32,
) {
    let start = std::time::Instant::now();
    // item 1
    let side = curve_side::<L, N>(name, act, d, e);
    let endpoint = chain_assertions(name, d, e, &side.js);
    println!(
        "ROW|item1|{name}|{}*2^{}-1|D={d}|m=2^{e}|kernel: ACTION_GEN3{} on E0[m], generator {}, cyclic of order m|{endpoint}|all j in F_p, none repeats|2-steps = 4-steps = strategy chain|chain digest {}|PASS",
        L::COFACTOR,
        L::TWO_ADIC_EXPONENT,
        if d % 2 == 1 { " - (m/2)" } else { "" },
        side.which,
        digest(&side.js),
    );
    generic_cross_check(name, d, e, &side);
    println!("ROW|item3-cross|{name}|generic module: same kernel points as the stored matrix, same {} j-invariants|PASS", side.js.len());
    let _ = &side.kernel;

    // item 2
    let quat = quaternion_side(name, alg, d, e);
    let j_end = &side.js[e as usize];
    let i2i = ideal_to_isogeny_checks::<L, N>(name, alg, act, &quat.ideal, &side, e);
    let p = BigUint::from_bytes_le(L::prime_le_bytes());
    let fld = Field::new(&p);
    let record = format!(
        "{}j_endpoint = {}\n",
        quat.record,
        fld.hex(&fp2_to_f2(j_end))
    );
    check_record(&format!("endpoints/{file}.txt"), &record);
    let _ = &quat.right_order;
    println!(
        "ROW|item2|{name}|I = O0 theta + O0 m: norm m, not in 2 O0|right order maximal, contains theta/m of discriminant -{d}, primitive|KernelToIdeal(G) = I|smallest equivalent ideal has norm m|ideal_to_isogeny_qlapoty on I: {} of {} attempts returned a curve|on I ∩ L: E0/(I ∩ L) = E_L/phi_L(G) (attempt {})|endpoints/{file}.txt|PASS|{:.1}s",
        i2i.direct_successes,
        i2i.direct_attempts,
        i2i.pushed_attempt,
        start.elapsed().as_secs_f64()
    );
}

#[test]
fn level_1() {
    run_level::<P324_3, { sqisign_rs::precomp::p324_3::IBZ_NLIMBS }>(
        "level I",
        "level1",
        &sqisign_rs::quat::params::level1(),
        &sqisign_rs::id2iso::params::level1(),
        3,
        162,
    );
}

#[test]
fn level_3() {
    run_level::<P500_27, { sqisign_rs::precomp::p500_27::IBZ_NLIMBS }>(
        "level III",
        "level3",
        &sqisign_rs::quat::params::level3(),
        &sqisign_rs::id2iso::params::level3(),
        27,
        250,
    );
}

#[test]
fn level_5() {
    run_level::<P664_17, { sqisign_rs::precomp::p664_17::IBZ_NLIMBS }>(
        "level V",
        "level5",
        &sqisign_rs::quat::params::level5(),
        &sqisign_rs::id2iso::params::level5(),
        68,
        331,
    );
}
