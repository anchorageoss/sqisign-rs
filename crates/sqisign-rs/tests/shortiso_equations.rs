//! Item 4 of paper-shortiso: the norm equations behind "Shortest isogenies
//! between supersingular elliptic curves with complex multiplication", in
//! machine integers only.
//!
//! An isogeny of degree `delta` from `E: y^2 = x^3 + x` over `F_p`,
//! `p = 3 mod 4`, to a curve admitting the order of discriminant `-D` gives
//! the element `delta (t + sqrt(-D)) / 2 = (t delta + X i + c j + e k) / 2`
//! of `End(E) = <1, i, (i + j)/2, (1 + k)/2>`, `t = D mod 2`, that is an
//! admissible quadruple `(delta; X, c, e)`:
//!
//!   `D delta^2 = X^2 + p (c^2 + e^2)`, `X = c mod 2`, `e = t delta mod 2`.
//!
//! For `D` in `{3, 8, 11, 12, 27}` and every prime `p = 3 mod 4` with
//! `3 < p < 3 * 10^6` inert in `Q(sqrt(-D))`, the least `delta` with
//! `D delta^2 - p` a perfect square is computed (`c^2 + e^2 = 1`) and every
//! `delta` below it is searched exhaustively for an admissible quadruple
//! with any `c^2 + e^2`: none may exist. For `D` in `{7, 28}` a second
//! value of `c^2 + e^2` competes and the minimum is taken over the two
//! conditions stated at `two_conditions`.
//!
//! The paper assumes `p > 3` throughout, so `p = 3` is left out. It is
//! inert in `Q(sqrt(-7))` only among these fields (it divides 3, 12 and 27
//! and splits for 8 and 11), so it would count for `D = 7` and `D = 28`.

use std::collections::BTreeMap;

const LIMIT: u64 = 3_000_000;

fn primes_below(n: u64) -> Vec<u64> {
    let n = n as usize;
    let mut sieve = vec![true; n];
    sieve[0] = false;
    sieve[1] = false;
    let mut i = 2;
    while i * i < n {
        if sieve[i] {
            let mut k = i * i;
            while k < n {
                sieve[k] = false;
                k += i;
            }
        }
        i += 1;
    }
    (0..n).filter(|&i| sieve[i]).map(|i| i as u64).collect()
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    r
}

fn square_root(n: u64) -> Option<u64> {
    let r = isqrt(n);
    (r * r == n).then_some(r)
}

/// The Legendre symbol `(a / p)` for an odd prime `p`, by Euler.
fn legendre(a: i64, p: u64) -> i32 {
    let a = a.rem_euclid(p as i64) as u64;
    if a == 0 {
        return 0;
    }
    let mut r = 1u128;
    let mut b = a as u128;
    let mut e = (p - 1) / 2;
    let p = p as u128;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    if r == 1 {
        1
    } else {
        -1
    }
}

/// Whether `(delta; X, c, e)` satisfies the parity conditions for `D`.
fn parities(d: u64, delta: u64, x: u64, c: u64, e: u64) -> bool {
    x % 2 == c % 2 && e % 2 == (d % 2 * delta) % 2
}

/// Whether the quadruple comes from a larger order: for a prime `l` with
/// `-D / l^2` a discriminant, `(delta; X/l, c/l, e/l)` is an admissible
/// quadruple for `D / l^2`. Then `delta (t' + sqrt(-D)/l) / 2` lies in
/// `End(E)` and the order of discriminant `-D` is not optimally embedded.
fn from_larger_order(d: u64, delta: u64, x: u64, c: u64, e: u64) -> bool {
    (2..=d).any(|l| {
        d % (l * l) == 0
            && matches!((d / (l * l)) % 4, 0 | 3)
            && x % l == 0
            && c % l == 0
            && e % l == 0
            && parities(d / (l * l), delta, x / l, c / l, e / l)
    })
}

/// The values `n = c^2 + e^2` for which an admissible quadruple
/// `(delta; X, c, e)` exists (with `optimal_only`: one that does not come
/// from a larger order), and the number of `(c, e)` pairs tried.
fn admissible_at(d: u64, p: u64, delta: u64, optimal_only: bool, tried: &mut u64) -> Vec<u64> {
    let lhs = d * delta * delta;
    let mut found = Vec::new();
    let mut c = 0u64;
    while p * c * c <= lhs {
        let mut e = 0u64;
        while p * (c * c + e * e) <= lhs {
            *tried += 1;
            let n = c * c + e * e;
            if let Some(x) = square_root(lhs - p * n) {
                if parities(d, delta, x, c, e)
                    && !(optimal_only && from_larger_order(d, delta, x, c, e))
                {
                    found.push(n);
                }
            }
            e += 1;
        }
        c += 1;
    }
    found.sort();
    found.dedup();
    found
}

/// The least `delta` with `D delta^2 - n p` a perfect square `X^2` for
/// which `(delta; X, c, e)` is admissible with the given `(c, e)`,
/// `c^2 + e^2 = n`; `None` if the parities of `(c, e)` are never
/// admissible for this `D` (`e` odd with `D` even).
fn least_delta(d: u64, p: u64, c: u64, e: u64, optimal_only: bool) -> Option<u64> {
    if d % 2 == 0 && e % 2 == 1 {
        return None;
    }
    let n = c * c + e * e;
    let mut delta = isqrt(n * p / d).max(1);
    loop {
        assert!(
            delta < 1_000_000,
            "no solution for D = {d}, p = {p}, n = {n}"
        );
        if d * delta * delta >= n * p {
            if let Some(x) = square_root(d * delta * delta - n * p) {
                if parities(d, delta, x, c, e)
                    && !(optimal_only && from_larger_order(d, delta, x, c, e))
                {
                    return Some(delta);
                }
            }
        }
        delta += 1;
    }
}

struct Tally {
    primes: u64,
    max_delta: u64,
    max_delta_p: u64,
    deltas_searched: u64,
    pairs_tried: u64,
    below: u64,
    /// The primes with an admissible quadruple below the claimed minimum.
    primes_below: std::collections::BTreeSet<u64>,
    /// `p = D m^2 - 1` among the primes: the least `delta` must be `m`.
    of_form_dm2: u64,
    /// `n = c^2 + e^2` at the least `delta`, with multiplicity over primes.
    n_at_min: BTreeMap<Vec<u64>, u64>,
}

/// `conditions`: the `(c, e)` whose least `delta` are computed; the claim
/// is that the least degree is the minimum over them and that nothing
/// admissible lies below.
fn run(d: u64, conditions: &[(u64, u64)], optimal_only: bool) -> Tally {
    let mut tally = Tally {
        primes: 0,
        max_delta: 0,
        max_delta_p: 0,
        deltas_searched: 0,
        pairs_tried: 0,
        below: 0,
        primes_below: Default::default(),
        of_form_dm2: 0,
        n_at_min: BTreeMap::new(),
    };
    for p in primes_below(LIMIT) {
        if p <= 3 || p % 4 != 3 || d % p == 0 || legendre(-(d as i64), p) != -1 {
            continue;
        }
        tally.primes += 1;
        let claimed = conditions
            .iter()
            .filter_map(|&(c, e)| least_delta(d, p, c, e, optimal_only))
            .min()
            .unwrap();
        // exhaustive search below the claimed minimum
        let mut delta = isqrt(p / d).max(1);
        while delta < claimed {
            tally.deltas_searched += 1;
            let hits = admissible_at(d, p, delta, optimal_only, &mut tally.pairs_tried);
            if !hits.is_empty() {
                tally.below += 1;
                tally.primes_below.insert(p);
                if tally.below <= 2 {
                    println!(
                        "BELOW|D={d}|p={p}|admissible quadruple at delta={delta} < {claimed}, c^2+e^2 in {hits:?} (a failure unless this is the D = 28 control)"
                    );
                }
            }
            delta += 1;
        }
        let at_min = admissible_at(d, p, claimed, optimal_only, &mut tally.pairs_tried);
        assert!(
            !at_min.is_empty(),
            "D = {d}, p = {p}: admissible at the minimum"
        );
        *tally.n_at_min.entry(at_min).or_insert(0) += 1;
        if claimed > tally.max_delta {
            tally.max_delta = claimed;
            tally.max_delta_p = p;
        }
        // Corollary A: p = D m^2 - 1 gives least delta m
        let m = isqrt((p + 1) / d);
        if d * m * m == p + 1 {
            tally.of_form_dm2 += 1;
            if !tally.primes_below.contains(&p) {
                assert_eq!(claimed, m, "D = {d}, p = {p} = D m^2 - 1: least delta is m");
            }
        }
    }
    tally
}

fn report(d: u64, what: &str, t: &Tally) {
    println!(
        "ROW|item4|D={d}|{what}|primes p = 3 mod 4, 3 < p < {LIMIT}, inert in Q(sqrt(-{d})): {}|largest least delta {} (p = {})|{} values of delta below the minimum searched, {} pairs (c, e) tried|admissible quadruples below the minimum: {} (at {} primes)|c^2+e^2 at the minimum, by number of primes: {:?}|of which p = D m^2 - 1: {} (least delta = m)|{}",
        t.primes,
        t.max_delta,
        t.max_delta_p,
        t.deltas_searched,
        t.pairs_tried,
        t.below,
        t.primes_below.len(),
        t.n_at_min,
        t.of_form_dm2,
        if t.below == 0 { "PASS" } else { "FAIL" },
    );
}

/// `D` in `{3, 8, 11, 12, 27}`: the least degree is the least `delta` with
/// `D delta^2 - p` a perfect square.
#[test]
fn one_condition() {
    let tallies: Vec<(u64, Tally)> = std::thread::scope(|s| {
        let handles: Vec<_> = [3u64, 8, 11, 12, 27]
            .into_iter()
            .map(|d| s.spawn(move || (d, run(d, &[(1, 0), (0, 1)], false))))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for (d, t) in &tallies {
        report(*d, "least delta with D delta^2 - p a square", t);
    }
    for (d, t) in &tallies {
        assert_eq!(
            t.below, 0,
            "D = {d}: an admissible quadruple below the minimum"
        );
        assert!(
            t.n_at_min.keys().all(|k| k.contains(&1)),
            "D = {d}: c^2 + e^2 = 1 at the minimum"
        );
    }
}

/// The two-condition statements, as reconstructed here (the paper's text
/// was not at hand; see paper-shortiso/REPORT.md).
///
/// `D = 7`: 2 splits in `Q(sqrt(-7))` and `c = e = 1` is admissible, so the
/// least degree is the least `delta` with `7 delta^2 - p` or
/// `7 delta^2 - 2p` a perfect square.
///
/// `D = 28`: a quadruple for `D = 7` doubles to one for `D = 28`
/// (`c^2 + e^2` multiplied by 4), so `c^2 + e^2 = 4` competes: the least
/// degree is the least `delta` with `28 delta^2 - p` or `7 delta^2 - p` a
/// perfect square. This is checked over all admissible quadruples, and
/// again over those only that embed the order of discriminant `-28`
/// optimally (not twice a quadruple for `D = 7`).
///
/// Control, expected to fail and not asserted to hold: the pair
/// `28 delta^2 - p`, `28 delta^2 - 8p` (the doubled second condition of
/// `D = 7`) is not a correct statement for `D = 28`.
#[test]
fn two_conditions() {
    let (t7, t28, t28_opt, t28_ctl) = std::thread::scope(|s| {
        let a = s.spawn(|| run(7, &[(1, 0), (0, 1), (1, 1)], false));
        let b = s.spawn(|| run(28, &[(1, 0), (2, 0), (0, 2)], false));
        let c = s.spawn(|| run(28, &[(1, 0), (2, 0), (0, 2)], true));
        let d = s.spawn(|| run(28, &[(1, 0), (2, 2)], false));
        (
            a.join().unwrap(),
            b.join().unwrap(),
            c.join().unwrap(),
            d.join().unwrap(),
        )
    });
    report(
        7,
        "least delta with 7 delta^2 - p or 7 delta^2 - 2p a square",
        &t7,
    );
    report(
        28,
        "least delta with 28 delta^2 - p or 7 delta^2 - p a square",
        &t28,
    );
    report(
        28,
        "the same, optimal embeddings of the order of discriminant -28 only",
        &t28_opt,
    );
    println!(
        "INFO|item4|D=28|control, expected to fail: with 28 delta^2 - p and 28 delta^2 - 8p as the two conditions, {} admissible quadruples lie below the minimum, at {} of {} primes; smallest case p = {}",
        t28_ctl.below,
        t28_ctl.primes_below.len(),
        t28_ctl.primes,
        t28_ctl.primes_below.iter().next().unwrap(),
    );
    assert_eq!(t7.below, 0, "D = 7");
    assert_eq!(t28.below, 0, "D = 28");
    assert_eq!(t28_opt.below, 0, "D = 28, optimal embeddings");
    assert!(t28_ctl.below > 0, "the control discriminates");
}
