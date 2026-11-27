//! Test-only support for the checks of "Shortest isogenies between
//! supersingular elliptic curves with complex multiplication"
//! (paper-shortiso/): arithmetic over an arbitrary prime, the kernel of
//! Corollary A, isogeny chains by Velu's formulas and class polynomials.
//! Nothing here is used by the library.

#![allow(dead_code)]

pub mod classpoly;
pub mod curve;
pub mod field;
pub mod quat;

use curve::{dlog_smooth, frobenius, has_exact_order, iota, Curve, Isogeny, Pt};
use field::{Field, F2};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use std::collections::HashSet;
use std::path::PathBuf;

/// `paper-shortiso/` at the workspace root.
pub fn paper_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../paper-shortiso")
}

/// Compare `content` with the checked-in record at `paper-shortiso/<rel>`;
/// with `SHORTISO_WRITE=1` the record is (re)written instead.
pub fn check_record(rel: &str, content: &str) {
    let path = paper_dir().join(rel);
    if std::env::var("SHORTISO_WRITE").is_ok_and(|v| v == "1") {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, content).unwrap();
        return;
    }
    let have = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("{} missing; run with SHORTISO_WRITE=1", path.display()));
    assert_eq!(have, content, "{} differs from the computed record", rel);
}

/// Corollary B: `(D, e)` with `m = 2^e` for `p = c 2^f - 1`, `c` odd.
pub fn cor_b(c: u64, f: u32) -> (u64, u32) {
    assert!(c % 2 == 1);
    if f % 2 == 0 {
        if c % 4 == 3 {
            (c, f / 2)
        } else {
            (4 * c, f / 2 - 1)
        }
    } else {
        (8 * c, (f - 3) / 2)
    }
}

/// The j-invariants of the curves with complex multiplication by the
/// orders of class number one that the checks meet.
pub fn j_class_number_one(d: u64) -> Option<BigInt> {
    let j: i128 = match d {
        3 => 0,
        4 => 1728,
        7 => -3375,
        8 => 8000,
        11 => -32768,
        12 => 54000,
        16 => 287496,
        19 => -884736,
        27 => -12288000,
        28 => 16581375,
        43 => -884736000,
        67 => -147197952000,
        163 => -262537412640768000,
        _ => return None,
    };
    Some(BigInt::from(j))
}

/// The kernel of Corollary A on `E: y^2 = x^3 + x`.
pub struct Kernel {
    /// A generator of `G`, of exact order `m`.
    pub gen: Pt,
    /// `theta'(P)` and `theta'(Q)` for the basis `(P, Q)` of `E[2m]` given,
    /// `theta' = iota + pi - t m`, `t = D mod 2`.
    pub images: (Pt, Pt),
    /// Which combination `theta'(P) + c theta'(Q)` the generator is
    /// (`c = u32::MAX` for `theta'(Q)` alone).
    pub combo: u32,
}

/// `G = (iota + pi - t m)(E[2m])`, `t = D mod 2`, from a basis `(P, Q)` of
/// `E[2m]`. Asserts that `G` is cyclic of order exactly `m`: a generator
/// of exact order `m` is exhibited and both images are shown to lie in the
/// group it generates, by a verified discrete logarithm.
pub fn corollary_a_kernel(
    f: &Field,
    d: u64,
    m: &BigUint,
    fac_m: &[(u32, u32)],
    basis: &(Pt, Pt),
) -> Kernel {
    let e = Curve::e0(f);
    let t = d % 2;
    let theta = |r: &Pt| {
        let s = e.add(f, &iota(f, r), &frobenius(f, r));
        if t == 1 {
            e.sub(f, &s, &e.mul(f, m, r))
        } else {
            s
        }
    };
    let k1 = theta(&basis.0);
    let k2 = theta(&basis.1);
    assert!(e.on_curve(f, &k1) && e.on_curve(f, &k2));
    let mut found = None;
    if has_exact_order(&e, f, &k1, m, fac_m) {
        found = Some((k1.clone(), 0));
    } else if has_exact_order(&e, f, &k2, m, fac_m) {
        found = Some((k2.clone(), u32::MAX));
    } else {
        let mut cand = k1.clone();
        for c in 1..64 {
            cand = e.add(f, &cand, &k2);
            if has_exact_order(&e, f, &cand, m, fac_m) {
                found = Some((cand.clone(), c));
                break;
            }
        }
    }
    let (gen, combo) = found.expect("G has an element of exact order m");
    assert!(
        dlog_smooth(&e, f, &gen, &k1, m, fac_m).is_some()
            && dlog_smooth(&e, f, &gen, &k2, m, fac_m).is_some(),
        "G is cyclic: both images lie in the group of the generator"
    );
    Kernel {
        gen,
        images: (k1, k2),
        combo,
    }
}

/// Walk the isogeny with cyclic kernel `<k>` on `E: y^2 = x^3 + x` as a
/// chain of prime-degree steps (`steps`, in the order taken). Returns the
/// j-invariants `j_0 = 1728, j_1, ..., j_n` and the last curve. The kernel
/// generator is pushed through every step and must die exactly at the end.
pub fn walk(f: &Field, k: &Pt, steps: &[u32]) -> (Vec<F2>, Curve) {
    let mut e = Curve::e0(f);
    let mut k = k.clone();
    let mut rem = steps
        .iter()
        .fold(BigUint::one(), |acc, &l| acc * BigUint::from(l));
    let mut js = vec![e.j(f)];
    for &l in steps {
        rem /= l;
        let kl = e.mul(f, &rem, &k);
        let phi = Isogeny::new(&e, f, &kl, l);
        k = phi.eval(f, &k);
        e = phi.codomain;
        assert!(e.on_curve(f, &k), "pushed point on the codomain");
        js.push(e.j(f));
    }
    assert!(rem.is_one() && e.is_inf(f, &k), "kernel exhausted");
    (js, e)
}

/// The prime steps of a factorisation, smallest primes first.
pub fn steps_ascending(fac: &[(u32, u32)]) -> Vec<u32> {
    fac.iter()
        .flat_map(|&(l, v)| std::iter::repeat(l).take(v as usize))
        .collect()
}

/// What a chain satisfied, beyond the assertions that abort.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainFacts {
    pub all_in_fp: bool,
    pub no_repeat: bool,
}

pub fn chain_facts(f: &Field, js: &[F2]) -> ChainFacts {
    let set: HashSet<&F2> = js.iter().collect();
    ChainFacts {
        all_in_fp: js.iter().all(|j| f.in_fp(j)),
        no_repeat: set.len() == js.len(),
    }
}

/// How the endpoint was identified.
pub fn endpoint_check(f: &Field, d: u64, j_end: &F2) -> (bool, String) {
    if let Some(j) = j_class_number_one(d) {
        (f.int(&j) == *j_end, format!("j = {j}"))
    } else {
        let h = class_polynomial_checked(d);
        (
            f.is_zero(&f.eval_poly(&h, j_end)),
            format!("root of H_-{d} (degree {})", h.len() - 1),
        )
    }
}

/// `H_{-D}`, recomputed and compared with the checked-in record
/// `paper-shortiso/classpoly/H_-D.txt`.
pub fn class_polynomial_checked(d: u64) -> Vec<BigInt> {
    let h = classpoly::class_polynomial(d);
    check_record(&format!("classpoly/H_-{d}.txt"), &classpoly::to_text(d, &h));
    h
}

pub fn pow2(e: u32) -> BigUint {
    BigUint::one() << e
}

pub fn is_zero_big(x: &BigUint) -> bool {
    x.is_zero()
}
