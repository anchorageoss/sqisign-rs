//! Item 7 of paper-shortiso: the minimality half of Corollary A, proved by
//! computation at given primes.
//!
//! What is computed. Let `(delta; X, c, e)` be an admissible quadruple,
//! `theta = (t delta + X i + c j + e k)/2` in `O0 = End(E_1728)` and
//! `I = O0 theta + O0 delta`, the left ideal whose kernel is the subgroup
//! `G` of the paper (for `p = D m^2 - 1` the quadruple is `(m; 1, 1, 0)`).
//! For a left ideal `J` of `O0` the map `beta -> nrd(beta) / nrd(J)` on `J`
//! is a positive definite integral quaternary form, and its values are
//! exactly the degrees of the isogenies between `E_1728` and the curve of
//! `J` (Deuring correspondence: `Hom(E_J, E_1728)` is `J`, the degree is
//! `nrd(beta) / nrd(J)`). So the least degree of an isogeny from `E_1728`
//! to that curve is the minimum of the form, and the minimum of a form is
//! a finite computation: reduce, bound the coordinates of any vector not
//! longer than a basis vector, evaluate every vector of the box. All of
//! it is exact integer arithmetic (`shortiso_common::quat`).
//!
//! Which ideals. `gamma = (X i + c j + e k)/delta` has square `-D` and
//! lies in the right order of `I`. For every order `O'` of discriminant
//! `-D/g^2` containing `O_D`, and every reduced primitive form `(a, b, c)`
//! of that discriminant, `J = I (Z a + Z (-b + gamma/g)/2)`. The test
//! checks, exactly, that each `J` is a left `O0`-module, that `O'`
//! multiplies it into itself on the right, and that no larger order does.
//! By the theory of oriented curves (Onuki 2021, Proposition 3.3 and
//! Theorem 3.4: the class group of `O'` acts freely and transitively on
//! the reductions of the curves with complex multiplication by `O'`, and
//! every curve with `O'` optimally embedded is one of these or a
//! Frobenius conjugate), every curve admitting `O_D` is the curve of one
//! of these `J` or its Frobenius conjugate, and conjugation does not
//! change the least degree from `E_1728`, which is defined over `F_p`. The
//! number of classes is checked against the class numbers.
//!
//! So: the minimum of each form is a proof (given Deuring's
//! correspondence) of the least degree to that one curve, and the minimum
//! over all the forms is a proof (given also the cited statement on the
//! class group action) of the least degree to any curve admitting `O_D`.
//!
//! The crate's own routine is compared at the three SQIsign primes: its
//! `reduce_basis` (the first vector of which `shortest_equivalent`
//! returns) is a constant-time Lagrange-Gauss reduction over `Z[i]` with a
//! fixed number of rounds and `sqrt(p)` replaced by its integer part. It
//! reduces; it does not certify. Its first vector is checked here to have
//! the minimal norm found by the enumeration.

mod shortiso_common;

use num_bigint::{BigInt, Sign};
use num_integer::{Integer, Roots};
use num_traits::{One, Zero};
use shortiso_common::classpoly::{class_number, reduced_forms};
use shortiso_common::cor_b;
use shortiso_common::quat::{exact_minimum, norm_form, o0_basis, Elem, Lat, Minimum};
use sqisign_rs::mp::Ibz;
use sqisign_rs::quat::lll::applications::reduce_basis;
use sqisign_rs::quat::{QuatAlg, QuatAlgElem, QuatIdeal};

fn big(x: i64) -> BigInt {
    BigInt::from(x)
}

/// The orders containing `O_D`: `(g, D / g^2)`.
fn overorders(d: u64) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut g = 1u64;
    while g * g <= d {
        if d % (g * g) == 0 {
            let dd = d / (g * g);
            if dd % 4 == 0 || dd % 4 == 3 {
                out.push((g, dd));
            }
        }
        g += 1;
    }
    out
}

struct Class {
    /// the order is `O_{D/g^2}`
    g: u64,
    disc: u64,
    form: (i64, i64, i64),
    minimum: Minimum,
}

struct Quad {
    delta: BigInt,
    x: BigInt,
    c: BigInt,
    e: BigInt,
}

/// `I = O0 theta + O0 delta`, the conductor index `g0` of the order
/// `O_{D/g0^2}` that `I` carries (the largest order of `Q(gamma)` in its
/// right order), and the classes of the curves admitting that order built
/// from `I`, each with the exact minimum of its norm form. When `g0 = 1`
/// these are the classes of all the curves admitting `O_D`.
fn classes(p: &BigInt, d: u64, q: &Quad) -> (Lat, u64, Vec<Class>) {
    let o0 = o0_basis();
    let o0_lat = Lat::from_gens(&o0);
    let t = big((d % 2) as i64);
    let delta = &q.delta;
    assert_eq!(
        BigInt::from(d) * delta * delta,
        &q.x * &q.x + p * (&q.c * &q.c + &q.e * &q.e),
        "D delta^2 = X^2 + p (c^2 + e^2)"
    );
    let theta = Elem::new(&big(2), [&(&t * delta), &q.x, &q.c, &q.e]);
    assert!(o0_lat.contains(&theta), "theta in O0 (the two congruences)");
    let (tn, td) = theta.norm(p);
    assert_eq!(
        tn,
        td * (BigInt::from(d) + &t) / big(4) * delta * delta,
        "nrd(theta) = n delta^2"
    );

    let mut gens = Vec::new();
    for b in &o0 {
        gens.push(b.mul(&theta, p));
        gens.push(b.mul(&Elem::scalar(delta), p));
    }
    let ideal = Lat::from_gens(&gens);
    assert!(ideal.stable_left(&o0, p), "I is a left O0-module");
    let (_, kappa) = norm_form(&ideal, p);
    assert_eq!(kappa, &ideal.den * &ideal.den * delta, "nrd(I) = delta");

    // gamma = (X i + c j + e k) / delta, gamma^2 = -D
    let gamma_sq = Elem::new(delta, [&BigInt::zero(), &q.x, &q.c, &q.e]);
    let sq = gamma_sq.mul(&gamma_sq, p);
    assert!(
        sq.c[0] == -BigInt::from(d) * &sq.den && sq.c[1..].iter().all(|v| v.is_zero()),
        "gamma^2 = -D"
    );

    // generator of the order of discriminant -D/g^2: (t' + gamma / g) / 2
    let omega = |gg: u64, dd: u64| -> Elem {
        let s = BigInt::from(dd % 2) * delta * BigInt::from(gg);
        Elem::new(&(big(2) * delta * BigInt::from(gg)), [&s, &q.x, &q.c, &q.e])
    };
    // The order that I carries: the largest order of Q(gamma) in the right
    // order of I. It contains O_D (theta / delta is in the right order);
    // for p = D m^2 - 1 it is O_D itself, for other primes it can be
    // larger (Proposition 5 of the paper).
    let all_orders = overorders(d);
    let stable: Vec<u64> = all_orders
        .iter()
        .filter(|&&(g, disc)| ideal.stable_right(&[omega(g, disc)], p))
        .map(|&(g, _)| g)
        .collect();
    let g0 = *stable.iter().max().expect("O_D is in the right order of I");
    for &(g, _) in &all_orders {
        assert_eq!(
            stable.contains(&g),
            g0 % g == 0,
            "the orders in the right order of I form a chain below one order"
        );
    }
    // the orders containing the one I carries
    let orders: Vec<(u64, u64)> = all_orders
        .iter()
        .copied()
        .filter(|&(g, _)| g % g0 == 0)
        .collect();
    let mut out = Vec::new();
    for &(g, disc) in &orders {
        let gb = BigInt::from(g);
        let den = big(2) * delta * &gb;
        let om = omega(g, disc);
        let forms = reduced_forms(disc);
        assert_eq!(forms.len(), class_number(disc));
        for &(a, b, c) in &forms {
            assert_eq!(b * b - 4 * a * c, -(disc as i64));
            // the proper O'-ideal Z a + Z beta, beta = (-b + gamma / g) / 2
            let beta = Elem::new(&den, [&(big(-b) * delta * &gb), &q.x, &q.c, &q.e]);
            let rel = beta
                .mul(&beta, p)
                .add(&Elem::scalar(&big(b)).mul(&beta, p))
                .add(&Elem::scalar(&big(a * c)));
            assert!(rel.is_zero(), "beta^2 + b beta + a c = 0");
            let j = ideal.times(&[Elem::scalar(&big(a)), beta], p);
            if g == g0 && a == 1 {
                assert_eq!(j, ideal, "I O' = I for the order O' that I carries");
            }
            assert!(j.stable_left(&o0, p), "J is a left O0-module");
            assert!(
                j.stable_right(std::slice::from_ref(&om), p),
                "J O' inside J"
            );
            // no larger order: the minimal overorders of O' have prime index
            for &(g2, disc2) in &orders {
                if g2 > g && g2 % g == 0 && is_prime_small(g2 / g) {
                    assert!(
                        !j.stable_right(&[omega(g2, disc2)], p),
                        "the order of discriminant -{disc2} does not stabilise J: O' is optimally embedded"
                    );
                }
            }
            let (bmat, _) = norm_form(&j, p);
            let minimum = exact_minimum(&bmat);
            out.push(Class {
                g,
                disc,
                form: (a, b, c),
                minimum,
            });
        }
    }
    let expected: usize = orders.iter().map(|&(_, dd)| class_number(dd)).sum();
    assert_eq!(out.len(), expected);
    (ideal, g0, out)
}

fn is_prime_small(n: u64) -> bool {
    n >= 2 && (2..n).take_while(|k| k * k <= n).all(|k| n % k != 0)
}

/// The least value over the classes of the orders containing `O_{D/g^2}`.
fn least_over(classes: &[Class], g: u64) -> BigInt {
    classes
        .iter()
        .filter(|c| c.g % g == 0)
        .map(|c| c.minimum.min.clone())
        .min()
        .unwrap()
}

/// `value / unit` to six decimals, by integer arithmetic (for the log
/// only; no decision uses it).
fn ratio(value: &BigInt, unit: &BigInt) -> String {
    let scaled = (value * big(1_000_000)).div_floor(unit);
    let (ip, fp) = scaled.div_rem(&big(1_000_000));
    format!("{ip}.{fp:06}")
}

// ---------------------------------------------------------- the crate

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

/// The crate's ideal `O0 theta + O0 m` as a `Lat`, and
/// `nrd(v1) / nrd(I)` for the first vector `v1` of the crate's reduced
/// basis, which is the norm `shortest_equivalent` returns.
fn crate_side<const N: usize>(alg: &QuatAlg<N>, d: u64, e: u32) -> (Lat, BigInt) {
    let m = Ibz::<N>::one().mul_2exp(e);
    let mut theta = QuatAlgElem::<N>::set(2, 0, 1, 1, 0);
    if d % 2 == 1 {
        theta.coord.0[0] = m;
    }
    let ideal = QuatIdeal::create_o0_pow_two(&theta, &m, alg).expect("inert Hermite form");
    assert_eq!(ideal.norm, m);
    let to_elem = |lat: &sqisign_rs::quat::QuatLattice<N>, c: usize| Elem {
        den: ibz_to_big(&lat.denom),
        c: std::array::from_fn(|r| ibz_to_big(&lat.basis.0[r][c])),
    };
    let lat = ideal.to_lattice();
    let cols: Vec<Elem> = (0..4).map(|c| to_elem(&lat, c)).collect();
    let red = reduce_basis(&ideal, alg);
    let v1 = to_elem(&red, 0);
    let p = ibz_to_big(&alg.p);
    let (n, den) = v1.norm(&p);
    let (q, r) = n.div_rem(&(den * (BigInt::one() << e)));
    assert!(r.is_zero());
    // the reduced basis spans the same lattice
    let red_lat = Lat::from_gens(&(0..4).map(|c| to_elem(&red, c)).collect::<Vec<_>>());
    let ideal_lat = Lat::from_gens(&cols);
    assert_eq!(
        red_lat, ideal_lat,
        "the crate's reduced basis spans the ideal"
    );
    (ideal_lat, q)
}

// ------------------------------------------------- the nine primes

struct Summary {
    classes: usize,
    evaluated: u64,
}

fn crypto_case(group: &str, c: u64, f: u32, crate_ideal: Option<(Lat, BigInt)>) -> Summary {
    let start = std::time::Instant::now();
    let (d, e) = cor_b(c, f);
    let m = BigInt::one() << e;
    let p = BigInt::from(c) * (BigInt::one() << f) - 1;
    assert_eq!(BigInt::from(d) * &m * &m - 1, p);
    let label = format!("{c}*2^{f}-1");
    // p is inert in Q(sqrt(-D)): (-D)^((p-1)/2) = -1 mod p
    let minus_d: BigInt = &p - BigInt::from(d);
    assert_eq!(
        minus_d.modpow(&((&p - BigInt::one()) / big(2)), &p),
        &p - 1,
        "{label}: p is inert in Q(sqrt(-D))"
    );
    let quad = Quad {
        delta: m.clone(),
        x: BigInt::one(),
        c: BigInt::one(),
        e: BigInt::zero(),
    };
    let (ideal, g0, cls) = classes(&p, d, &quad);
    assert_eq!(
        g0, 1,
        "{label}: E_1728/G carries O_D itself and no larger order"
    );

    let mut crate_note = String::from("crate: not applicable (no parameter set)");
    if let Some((crate_lat, crate_norm)) = crate_ideal {
        assert_eq!(
            crate_lat, ideal,
            "{label}: the crate's ideal is O0 theta + O0 m"
        );
        crate_note = format!(
            "crate: same ideal; first vector of reduce_basis has norm {} m",
            ratio(&crate_norm, &m)
        );
        assert_eq!(
            crate_norm, cls[0].minimum.min,
            "{label}: the crate's reduced first vector is a shortest vector"
        );
    }

    // the endpoint: the class of I itself
    let endpoint = &cls[0];
    assert!(endpoint.g == 1 && endpoint.form.0 == 1);
    assert_eq!(
        endpoint.minimum.min, m,
        "{label}: the least degree from E_1728 to E_1728/G is m"
    );
    // every curve admitting O_D
    for cl in &cls {
        assert!(
            cl.minimum.min >= m,
            "{label}: no curve admitting O_D is closer than m"
        );
    }
    // every curve admitting the order of discriminant -D/g^2: Corollary A
    // for D/g^2, where p = (D/g^2) (g m)^2 - 1
    for (g, disc) in overorders(d) {
        assert_eq!(
            least_over(&cls, g),
            BigInt::from(g) * &m,
            "{label}: least degree to a curve admitting O_{disc} is {g} m"
        );
    }
    let at_m = cls.iter().filter(|c| c.minimum.min == m).count();
    let evaluated: u64 = cls.iter().map(|c| c.minimum.evaluated).sum();
    for cl in &cls {
        println!(
            "CLASS|{label}|order of discriminant -{}|form ({}, {}, {})|least degree = {} m = {}|minimal vectors {}|box half sides {:?}|vectors evaluated {}",
            cl.disc,
            cl.form.0,
            cl.form.1,
            cl.form.2,
            ratio(&cl.minimum.min, &m),
            cl.minimum.min,
            cl.minimum.count,
            cl.minimum.half_sides,
            cl.minimum.evaluated,
        );
    }
    let orders: Vec<String> = overorders(d)
        .iter()
        .map(|&(_, dd)| format!("h(-{dd}) = {}", class_number(dd)))
        .collect();
    println!(
        "ROW|item7|{group}|{label}|D={d}|m=2^{e}|classes of curves admitting O_D: {} ({})|least degree to E_1728/G = m exactly|least degree over all classes = m, attained by {at_m} class(es)|{crate_note}|vectors evaluated {evaluated}|PASS|{:.1}s",
        cls.len(),
        orders.join(", "),
        start.elapsed().as_secs_f64()
    );
    Summary {
        classes: cls.len(),
        evaluated,
    }
}

const SCOPE: &str = "\
SCOPE|What the crate's routine proves: no minimum. QuatIdeal::shortest_equivalent returns the first vector of reduce_basis, a constant-time Lagrange-Gauss reduction over Z[i] with a fixed number of rounds, run with sqrt(p) replaced by its integer part; the result is not checked. It reduces and does not certify a minimum.
SCOPE|What this test proves, per class: the exact minimum of the norm form of a left ideal J of O0, by integer enumeration of a box that contains every vector not longer than a basis vector. By the Deuring correspondence that minimum is the least degree of an isogeny from E_1728 to the curve of J.
SCOPE|For the class of I = O0 theta + O0 m the minimum is m: no isogeny from E_1728 to the endpoint E_1728/G has degree below m. (That G is the kernel of I is item 2.)
SCOPE|For D = 3 the endpoint j = 0 is the only curve admitting O_3, so this is all of the minimality of Corollary A at 3*2^324-1.
SCOPE|For the other primes the curves admitting O_D are enumerated by acting on I with the class group of every order containing O_D (one more curve, j = 0, at 27*2^500-1; the other three roots of H_-68 at 17*2^664-1), and the minimum is computed for each: every one is at least m, and only the endpoint attains m. This is all of the minimality of Corollary A at these primes, given that the class group of each order acts transitively on the curves carrying it up to Frobenius conjugation (Onuki 2021, Proposition 3.3 and Theorem 3.4). That each J is a left O0-ideal carrying exactly the stated order, and the number of classes, are checked here exactly.
SCOPE|Observed and not claimed in the paper: at each of the nine primes the least degree to the curve of the reduced form (a, b, c) of discriminant -D/g^2 is exactly g a m.";

#[test]
fn sqisign_primes_exact_minimum() {
    println!("{SCOPE}");
    let l1 = crate_side(&sqisign_rs::quat::params::level1(), 3, 162);
    let l3 = crate_side(&sqisign_rs::quat::params::level3(), 27, 250);
    let l5 = crate_side(&sqisign_rs::quat::params::level5(), 68, 331);
    let a = crypto_case("sqisign", 3, 324, Some(l1));
    let b = crypto_case("sqisign", 27, 500, Some(l3));
    let c = crypto_case("sqisign", 17, 664, Some(l5));
    assert_eq!((a.classes, b.classes, c.classes), (1, 2, 4));
    println!(
        "SUMMARY|sqisign|3 primes|{} classes|{} vectors evaluated|least degree to any curve admitting O_D is m at each",
        a.classes + b.classes + c.classes,
        a.evaluated + b.evaluated + c.evaluated
    );
}

#[test]
fn mike_primes_exact_minimum() {
    let mut classes = 0;
    let mut evaluated = 0;
    for (c, f) in [
        (117, 374),
        (77, 566),
        (41, 758),
        (633, 308),
        (593, 474),
        (317, 628),
    ] {
        let s = crypto_case("mike", c, f, None);
        classes += s.classes;
        evaluated += s.evaluated;
    }
    assert_eq!(classes, 10 + 8 + 8 + 20 + 24 + 10);
    println!(
        "SUMMARY|mike|6 primes|{classes} classes|{evaluated} vectors evaluated|least degree to any curve admitting O_D is m at each"
    );
}

// ---------------------------------------------- small primes, any p

/// The least admissible quadruple for `D`, by search.
fn least_quadruple(p: u64, d: u64) -> Quad {
    let t = d % 2;
    for delta in 1u64.. {
        let top = d * delta * delta;
        let mut c = 0u64;
        while p * c * c <= top {
            let mut e = 0u64;
            while p * (c * c + e * e) <= top {
                let x2 = top - p * (c * c + e * e);
                let x = x2.sqrt();
                if x * x == x2 && (x + c) % 2 == 0 && (e + t * delta) % 2 == 0 {
                    return Quad {
                        delta: BigInt::from(delta),
                        x: BigInt::from(x),
                        c: BigInt::from(c),
                        e: BigInt::from(e),
                    };
                }
                e += 1;
            }
            c += 1;
        }
    }
    unreachable!()
}

fn is_inert(p: u64, d: u64) -> bool {
    // the Legendre symbol of -D at the odd prime p
    let a = (p - d % p) % p;
    if a == 0 {
        return false;
    }
    let mut acc = 1u128;
    let mut base = a as u128;
    let mut k = (p - 1) / 2;
    while k > 0 {
        if k & 1 == 1 {
            acc = acc * base % p as u128;
        }
        base = base * base % p as u128;
        k >>= 1;
    }
    acc == (p - 1) as u128
}

/// The same machinery at small primes of any shape, against an independent
/// search. For every non-square `D <= 100` and every prime `p = 3 mod 4`,
/// `3 < p < 1000`, inert in `Q(sqrt(-D))`: the least admissible quadruple
/// is found by search, `I` is built from it, and the least value of the
/// norm form of `I` must be its `delta`. Then, for every order `O'`
/// containing the order that `I` carries, the least value over the classes
/// of the orders containing `O'` must be the least `delta` of an
/// admissible quadruple for the discriminant of `O'`, found by its own
/// search (Proposition 3 of the paper for that discriminant). Here
/// `E_1728/G` may carry an order larger than `O_D`; the curves on which
/// `O_D` is then optimally embedded are not of the form above and are not
/// enumerated.
#[test]
fn small_primes_any_shape() {
    let start = std::time::Instant::now();
    let primes: Vec<u64> = (7u64..1000)
        .filter(|&n| n % 4 == 3 && is_prime_small(n))
        .collect();
    let mut pairs = 0u64;
    let mut total_classes = 0u64;
    let mut evaluated = 0u64;
    let mut max_h = 0usize;
    let mut discs = 0u64;
    let mut larger = 0u64;
    let mut compared = 0u64;
    for d in 3u64..=100 {
        if !(d % 4 == 0 || d % 4 == 3) || d.sqrt().pow(2) == d {
            continue;
        }
        discs += 1;
        for &p in &primes {
            if !is_inert(p, d) {
                continue;
            }
            let pb = BigInt::from(p);
            let quad = least_quadruple(p, d);
            let (_, g0, cls) = classes(&pb, d, &quad);
            assert_eq!(
                cls[0].minimum.min, quad.delta,
                "D = {d}, p = {p}: the class of I"
            );
            if g0 > 1 {
                larger += 1;
            }
            for (g, disc) in overorders(d) {
                if g % g0 != 0 {
                    continue;
                }
                let expected = least_quadruple(p, disc).delta;
                assert_eq!(
                    least_over(&cls, g),
                    expected,
                    "D = {d}, p = {p}: least degree to a curve admitting O_{disc}"
                );
                compared += 1;
            }
            pairs += 1;
            total_classes += cls.len() as u64;
            evaluated += cls.iter().map(|c| c.minimum.evaluated).sum::<u64>();
            max_h = max_h.max(cls.iter().filter(|c| c.g == g0).count());
        }
    }
    println!(
        "ROW|item7|small primes|{discs} discriminants -D, D <= 100 not a square|primes p = 3 mod 4, 3 < p < 1000, inert|{pairs} pairs (D, p), {larger} of them with E_1728/G carrying an order larger than O_D|{total_classes} classes, class numbers up to {max_h}|least value over the classes = least admissible delta found by search: {compared} comparisons, one for each order containing the one carried|vectors evaluated {evaluated}|PASS|{:.1}s",
        start.elapsed().as_secs_f64()
    );
}
