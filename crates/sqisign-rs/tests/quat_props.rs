//! Algebraic properties of the quaternion layer that hold independently of
//! the reference: membership, norms, the norm equations of Qlapoty.

mod quat_common;

use sqisign_rs::mp::{Ibz, ShakeRng};
use sqisign_rs::quat::integers::generate_random_prime;
use sqisign_rs::quat::lll::applications::reduce_basis;
use sqisign_rs::quat::qlapoty::qlapoty;
use sqisign_rs::quat::{QuatAlg, QuatAlgElem, QuatIdeal, Vec4};

fn rng(label: &str) -> ShakeRng {
    ShakeRng::new(label.as_bytes(), b"prp")
}

fn check<const N: usize>(alg: &QuatAlg<N>, label: &str) {
    let mut r = rng(label);
    // O0 contains 1, i, (i+j)/2, (1+k)/2 and their products
    let gens = [
        QuatAlgElem::<N>::set(1, 1, 0, 0, 0),
        QuatAlgElem::set(1, 0, 1, 0, 0),
        QuatAlgElem::set(2, 0, 1, 1, 0),
        QuatAlgElem::set(2, 1, 0, 0, 1),
    ];
    for a in &gens {
        for b in &gens {
            let prod = a.mul(b, alg);
            assert!(
                alg.o0.contains(&prod).0,
                "O0 is closed under multiplication"
            );
        }
    }
    assert!(
        !alg.o0.contains(&QuatAlgElem::set(2, 1, 0, 0, 0)).0,
        "1/2 is not in O0"
    );
    // norm is multiplicative
    let a = QuatAlgElem::<N>::set(3, 5, -7, 2, 1);
    let b = QuatAlgElem::<N>::set(1, -4, 3, 1, -6);
    let (na, da) = a.norm(alg);
    let (nb, db) = b.norm(alg);
    let (nab, dab) = a.mul(&b, alg).norm(alg);
    assert_eq!(na.mul(&nb).mul(&dab), nab.mul(&da).mul(&db));
    // conj(a) a = n(a)
    let (n, d) = a.norm(alg);
    let aa = a.conj().mul(&a, alg);
    assert!(aa.equals(&QuatAlgElem::scalar(&n, &d)));

    for it in 0..2 {
        let bits = alg.p_bits as u32 + 30 + 50 * it;
        let norm = generate_random_prime::<N>(false, bits, alg.primality_num_iter, &mut r).unwrap();
        let ideal = QuatIdeal::random_given_prime_norm(&norm, alg, &mut r).unwrap();
        let lat = ideal.to_lattice();
        assert!(lat.basis.is_hnf());
        let gen = ideal.odd_inert_gen();
        assert!(lat.contains(&gen).0, "the generator lies in the ideal");
        let (gn, gd) = gen.norm(alg);
        assert!(
            gd.is_one() && gn.modulo(&ideal.norm).is_zero(),
            "generator norm divisible by the ideal norm"
        );
        assert!(lat.inclusion(&alg.o0), "ideal inside O0");
        // reduced basis spans the same lattice and its vectors have norm about sqrt(p) N
        let red = reduce_basis(&ideal, alg);
        assert!(red.inclusion(&lat), "reduced vectors lie in the ideal");
        let (det_red, _) = red.basis.inv_with_det_as_denom();
        let (det_lat, _) = lat.basis.inv_with_det_as_denom();
        assert!(
            det_red.abs() == det_lat.abs() && red.denom == lat.denom,
            "reduction preserves the lattice"
        );
        for c in 0..4 {
            let v = QuatAlgElem {
                denom: red.denom,
                coord: Vec4([
                    red.basis.0[0][c],
                    red.basis.0[1][c],
                    red.basis.0[2][c],
                    red.basis.0[3][c],
                ]),
            };
            let (n, d) = v.norm(alg);
            assert!(d.is_one());
            let (q, rem) = n.div(&ideal.norm);
            assert!(rem.is_zero(), "norms in the ideal are multiples of N");
            assert!(q.bitsize() <= alg.p_bits + 4, "reduced vectors are short");
        }
        // shortest equivalent ideal has small norm
        let (g, equiv, red_ideal) = ideal.shortest_equivalent(alg, &mut r);
        assert!(
            red_ideal.norm.bitsize() <= alg.p_bits / 2 + 8,
            "shortest equivalent norm about sqrt(p)"
        );
        let (en, ed) = equiv.norm(alg);
        assert!(
            ed.is_one() && en == ideal.norm.mul(&red_ideal.norm),
            "equiv has norm N * n(red)"
        );
        let (gn, gd) = g.norm(alg);
        assert!(gd.is_one() && gn.modulo(&red_ideal.norm).is_zero());
        assert!(red_ideal.to_lattice().contains(&g).0);
        // Qlapoty: beta1 in the ideal, norm d1 N, theta of norm d1 (2^e - d1)
        let (beta1, d1, theta) = qlapoty(&ideal, alg, &mut r).expect("qlapoty");
        assert!(lat.contains(&beta1).0, "beta1 in the ideal");
        let (bn, bd) = beta1.norm(alg);
        assert!(bd.is_one() && bn == d1.mul(&ideal.norm));
        assert!(d1.is_odd());
        let two_e = Ibz::<N>::one().mul_2exp(alg.qlapoty_used_power_of_two);
        let d2 = two_e.sub(&d1);
        let (tn, td) = theta.norm(alg);
        assert!(td.is_one() && tn == d1.mul(&d2), "theta has norm d1 d2");
        assert!(alg.o0.contains(&theta).0, "theta is an endomorphism");
        // small equivalent coprime ideal of prime norm
        let (g2, eq) = ideal
            .small_equivalent_coprime(Some(&Ibz::zero()), alg, &mut r)
            .unwrap();
        assert!(eq.norm.probab_prime(alg.primality_num_iter, &mut r));
        let (gn2, gd2) = g2.norm(alg);
        assert!(gd2.is_one() && gn2 == ideal.norm.mul(&eq.norm));
    }
}

#[test]
fn level1() {
    check(&sqisign_rs::quat::params::level1(), "level1");
}

#[test]
fn level3() {
    check(&sqisign_rs::quat::params::level3(), "level3");
}

#[test]
fn level5() {
    check(&sqisign_rs::quat::params::level5(), "level5");
}
