//! Elliptic curve arithmetic tests.
//!
//! Tests: xDBL+xADD consistency, xDBLADD, xDBL variants, zero identities,
//! and Jacobian arithmetic (ADD, DBL, negation, commutativity, associativity,
//! Weierstrass round-trip).

use sqisign_verify::ec::jacobian::{jac_add, jac_dbl, jac_dbl_ws};
use sqisign_verify::ec::point::*;
use sqisign_verify::ec::{EcCurve, EcPoint, JacPoint};
use sqisign_verify::fp::{Fp, Fp2, FpBackend};
use sqisign_verify::params::Level1;

// ---------------------------------------------------------------------------
// Simple deterministic PRNG for test reproducibility
// ---------------------------------------------------------------------------

struct TestRng {
    state: [u64; 4],
}

impl TestRng {
    fn new(seed: u64) -> Self {
        Self {
            state: [
                seed ^ 0x9E3779B97F4A7C15,
                seed.wrapping_mul(0x6C62272E07BB0142),
                seed.wrapping_mul(0xBF58476D1CE4E5B9),
                seed.wrapping_mul(0x94D049BB133111EB),
            ],
        }
    }

    fn next_u64(&mut self) -> u64 {
        // xoshiro256**
        let result = self.state[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }

    fn fill_bytes(&mut self, buf: &mut [u8]) {
        let mut i = 0;
        while i < buf.len() {
            let val = self.next_u64();
            let bytes = val.to_le_bytes();
            let remaining = buf.len() - i;
            let to_copy = remaining.min(8);
            buf[i..i + to_copy].copy_from_slice(&bytes[..to_copy]);
            i += to_copy;
        }
    }
}

// ---------------------------------------------------------------------------
// Test helper: random field elements and points
// ---------------------------------------------------------------------------

fn fp_random<L: FpBackend>(rng: &mut TestRng) -> Fp<L> {
    let mut bytes = [0u8; 64]; // oversized for reduction
    rng.fill_bytes(&mut bytes);
    Fp::<L>::decode_reduce(&bytes)
}

fn fp2_random<L: FpBackend>(rng: &mut TestRng) -> Fp2<L> {
    Fp2 {
        re: fp_random::<L>(rng),
        im: fp_random::<L>(rng),
    }
}

/// Check if a projective point is on the curve `y^2 = x^3 + (A/C)x^2 + x`.
/// Tests whether `xz(C^2x^2 + zACx + z^2C^2)` is a square.
fn projective_is_on_curve<L: FpBackend>(p: &EcPoint<L>, curve: &EcCurve<L>) -> bool {
    let t0 = curve.c.mul(&p.x);
    let t1 = t0.mul(&p.z);
    let t1 = t1.mul(&curve.a);
    let t2 = curve.c.mul(&p.z);
    let t0 = t0.sqr();
    let t2 = t2.sqr();
    let t0 = t0.add(&t1);
    let t0 = t0.add(&t2);
    let t0 = t0.mul(&p.x);
    let t0 = t0.mul(&p.z);
    bool::from(t0.is_square()) || bool::from(t0.ct_is_zero())
}

/// Generate a random normalized point `(x : 1)` on the curve.
fn ec_random_normalized<L: FpBackend>(rng: &mut TestRng, curve: &EcCurve<L>) -> EcPoint<L> {
    loop {
        let p = EcPoint {
            x: fp2_random::<L>(rng),
            z: Fp2::one(),
        };
        if projective_is_on_curve(&p, curve) {
            return p;
        }
    }
}

/// Generate a random projective point on the curve (randomized Z).
fn ec_random<L: FpBackend>(rng: &mut TestRng, curve: &EcCurve<L>) -> EcPoint<L> {
    let mut p = ec_random_normalized::<L>(rng, curve);
    let z = fp2_random::<L>(rng);
    p.x = p.x.mul(&z);
    p.z = z;
    p
}

/// Compute the projective difference point P-Q.
/// Based on Proposition 3 of https://eprint.iacr.org/2017/518.pdf
fn projective_difference_point<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> EcPoint<L> {
    let t0 = p.x.mul(&q.x);
    let t1 = p.z.mul(&q.z);
    let bxx = t0.sub(&t1);
    let mut bxx = bxx.sqr();
    bxx = bxx.mul(&curve.c);

    let mut bxz = t0.add(&t1);
    let t0 = p.x.mul(&q.z);
    let t1 = p.z.mul(&q.x);
    let bzz_sum = t0.add(&t1);
    bxz = bxz.mul(&bzz_sum);
    let bzz_diff = t0.sub(&t1);
    let mut bzz = bzz_diff.sqr();
    bzz = bzz.mul(&curve.c);
    bxz = bxz.mul(&curve.c);
    let t0_t1 = t0.mul(&t1);
    let mut extra = t0_t1.mul(&curve.a);
    extra = extra.add(&extra);
    bxz = bxz.add(&extra);

    // Normalization factor: C * conj(C)^2 * conj(P.z)^2 * conj(Q.z)^2
    let c_conj = curve.c.conjugate();
    let c_conj_sq = c_conj.sqr();
    let mut norm = c_conj_sq.mul(&curve.c);
    let pz_conj = p.z.conjugate();
    let pz_conj_sq = pz_conj.sqr();
    norm = norm.mul(&pz_conj_sq);
    let qz_conj = q.z.conjugate();
    let qz_conj_sq = qz_conj.sqr();
    norm = norm.mul(&qz_conj_sq);

    bxx = bxx.mul(&norm);
    bxz = bxz.mul(&norm);
    bzz = bzz.mul(&norm);

    // Solve quadratic
    let disc = bxz.sqr();
    let prod = bxx.mul(&bzz);
    let disc = disc.sub(&prod);
    let disc = disc.sqrt();

    EcPoint {
        x: bxz.add(&disc),
        z: bzz,
    }
}

// ---------------------------------------------------------------------------
// Test 1: xDBL + xADD consistency
// ---------------------------------------------------------------------------

#[test]
fn test_xdbl_xadd() {
    let mut rng = TestRng::new(42);
    let iterations = 100;

    // Test with A=0, C=1
    let mut curve = EcCurve::<Level1>::default();
    curve.normalize_a24();
    run_xdbl_xadd(&mut rng, &curve, iterations);

    // Test with randomized C
    let mut curve2 = curve.clone();
    curve2.c = fp2_random::<Level1>(&mut rng);
    curve2.a = curve2.a.mul(&curve2.c);
    curve2.is_a24_computed_and_normalized = false;
    curve2.normalize_a24();
    run_xdbl_xadd(&mut rng, &curve2, iterations);
}

fn run_xdbl_xadd(rng: &mut TestRng, curve: &EcCurve<Level1>, n: usize) {
    for _ in 0..n {
        let p = ec_random::<Level1>(rng, curve);
        let q = ec_random::<Level1>(rng, curve);
        let pq = projective_difference_point(&p, &q, curve);

        // 2(P + Q) == 2P + 2Q
        let r1 = xadd(&p, &q, &pq);
        let r1 = ec_dbl(&r1, curve);
        let p2 = ec_dbl(&p, curve);
        let q2 = ec_dbl(&q, curve);
        let pq2 = ec_dbl(&pq, curve);
        let r2 = xadd(&p2, &q2, &pq2);
        assert!(bool::from(r1.ct_equal(&r2)), "Failed 2(P + Q) = 2P + 2Q");

        // (P+Q) + (P-Q) == 2P
        let pq_sum = xadd(&p2, &q2, &pq2);
        let q4 = ec_dbl(&q2, curve);
        let r1 = xadd(&pq_sum, &pq2, &q4);
        let p4 = ec_dbl(&p2, curve);
        assert!(bool::from(r1.ct_equal(&p4)), "Failed (P+Q) + (P-Q) = 2P");
    }
}

// ---------------------------------------------------------------------------
// Test 2: xDBLADD
// ---------------------------------------------------------------------------

#[test]
fn test_xdbladd() {
    let mut rng = TestRng::new(123);
    let iterations = 100;

    let mut curve = EcCurve::<Level1>::default();
    curve.normalize_a24();
    run_xdbladd(&mut rng, &curve, iterations);

    let mut curve2 = curve.clone();
    curve2.c = fp2_random::<Level1>(&mut rng);
    curve2.a = curve2.a.mul(&curve2.c);
    curve2.is_a24_computed_and_normalized = false;
    curve2.normalize_a24();
    run_xdbladd(&mut rng, &curve2, iterations);
}

fn run_xdbladd(rng: &mut TestRng, curve: &EcCurve<Level1>, n: usize) {
    let a24 = curve.ac_to_a24();

    for _ in 0..n {
        let p = ec_random::<Level1>(rng, curve);
        let q = ec_random::<Level1>(rng, curve);
        let pq = projective_difference_point(&p, &q, curve);

        let (r1, r2) = xdbladd(&p, &q, &pq, &a24, false);
        let expected_add = xadd(&p, &q, &pq);
        assert!(
            bool::from(r2.ct_equal(&expected_add)),
            "Failed addition in xDBLADD"
        );
        let expected_dbl = ec_dbl(&p, curve);
        assert!(
            bool::from(r1.ct_equal(&expected_dbl)),
            "Failed doubling in xDBLADD"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 3: xDBL variants
// ---------------------------------------------------------------------------

#[test]
fn test_xdbl_variants() {
    let mut rng = TestRng::new(456);
    let iterations = 100;

    let mut curve = EcCurve::<Level1>::default();
    curve.normalize_a24();
    run_xdbl_variants(&mut rng, &mut curve, iterations);

    let mut curve2 = curve.clone();
    curve2.c = fp2_random::<Level1>(&mut rng);
    curve2.a = curve2.a.mul(&curve2.c);
    curve2.is_a24_computed_and_normalized = false;
    curve2.normalize_a24();
    run_xdbl_variants(&mut rng, &mut curve2, iterations);
}

fn run_xdbl_variants(rng: &mut TestRng, curve: &mut EcCurve<Level1>, n: usize) {
    let a24 = curve.ac_to_a24();
    let mut a24norm = a24.clone();
    a24norm.normalize();

    // Randomize projective representation for the non-normalized case
    let z = fp2_random::<Level1>(rng);
    let a24_rand = EcPoint {
        x: a24.x.mul(&z),
        z: a24.z.mul(&z),
    };

    let ac = EcPoint {
        x: curve.a.clone(),
        z: curve.c.clone(),
    };

    for _ in 0..n {
        let p = ec_random::<Level1>(rng, curve);
        let r1 = xdbl(&p, &ac);
        let r2 = xdbl_a24(&p, &a24_rand, false);
        let r3 = xdbl_a24(&p, &a24norm, true);
        let r4 = xdbl_e0(&p);

        assert!(
            bool::from(r1.ct_equal(&r2)),
            "xDBL and xDBL_A24 don't match"
        );
        assert!(
            bool::from(r1.ct_equal(&r3)),
            "xDBL and xDBL_A24 normalized don't match"
        );
        assert!(bool::from(r1.ct_equal(&r4)), "xDBL and xDBL_E0 don't match");
    }
}

// ---------------------------------------------------------------------------
// Test 4: Zero identities
// ---------------------------------------------------------------------------

#[test]
fn test_zero_identities() {
    let mut rng = TestRng::new(789);
    let iterations = 100;

    let mut curve = EcCurve::<Level1>::default();
    curve.normalize_a24();
    run_zero_identities(&mut rng, &mut curve, iterations);

    let mut curve2 = curve.clone();
    curve2.c = fp2_random::<Level1>(&mut rng);
    curve2.a = curve2.a.mul(&curve2.c);
    curve2.is_a24_computed_and_normalized = false;
    curve2.normalize_a24();
    run_zero_identities(&mut rng, &mut curve2, iterations);
}

fn run_zero_identities(rng: &mut TestRng, curve: &mut EcCurve<Level1>, n: usize) {
    let ec_zero = EcPoint::<Level1>::identity();
    assert!(bool::from(ec_zero.is_zero()));

    for _ in 0..n {
        let p = ec_random::<Level1>(rng, curve);

        // 0 + 0 = 0
        let r = xadd(&ec_zero, &ec_zero, &ec_zero);
        assert!(bool::from(r.is_zero()), "Failed 0 + 0 = 0");

        // P - P = 0
        let p2 = ec_dbl(&p, curve);
        let r = xadd(&p, &p, &p2);
        assert!(bool::from(r.is_zero()), "Failed P - P = 0");

        // 2*0 = 0
        let r = ec_dbl(&ec_zero, curve);
        assert!(bool::from(r.is_zero()), "Failed 2*0 = 0");

        // P + 0 = P
        let r = xadd(&p, &ec_zero, &p);
        assert!(bool::from(r.ct_equal(&p)), "Failed P + 0 = P (first)");
        let r = xadd(&ec_zero, &p, &p);
        assert!(bool::from(r.ct_equal(&p)), "Failed P + 0 = P (second)");

        // xDBLADD with zero
        let a24 = curve.ac_to_a24();
        let (_r_dbl, r_add) = xdbladd(&p, &ec_zero, &p, &a24, false);
        assert!(
            bool::from(r_add.ct_equal(&p)),
            "Failed P + 0 = P in xDBLADD"
        );

        let (r_dbl2, r_add) = xdbladd(&ec_zero, &p, &p, &a24, false);
        assert!(
            bool::from(r_add.ct_equal(&p)),
            "Failed 0 + P = P in xDBLADD"
        );
        assert!(bool::from(r_dbl2.is_zero()), "Failed 2*0 = 0 in xDBLADD");
    }
}

// ---------------------------------------------------------------------------
// Test 5: Jacobian arithmetic
// ---------------------------------------------------------------------------

#[test]
fn test_jacobian() {
    let mut rng = TestRng::new(1337);
    let iterations = 100;

    let mut curve = EcCurve::<Level1>::default();
    curve.normalize_a24();
    run_jacobian(&mut rng, &curve, iterations);

    let mut curve2 = curve.clone();
    curve2.c = fp2_random::<Level1>(&mut rng);
    curve2.a = curve2.a.mul(&curve2.c);
    curve2.is_a24_computed_and_normalized = false;
    curve2.normalize_a24();
    run_jacobian(&mut rng, &curve2, iterations);
}

fn run_jacobian(rng: &mut TestRng, curve: &EcCurve<Level1>, n: usize) {
    let jac_zero = JacPoint::<Level1>::identity();

    for _ in 0..n {
        let p_xz = ec_random_normalized::<Level1>(rng, curve);
        let q_xz = ec_random_normalized::<Level1>(rng, curve);

        // Lift to Jacobian
        let (py, _) = curve.recover_y(&p_xz.x);
        let s = JacPoint::new(p_xz.x, py, Fp2::one());
        let (qy, _) = curve.recover_y(&q_xz.x);
        let t = JacPoint::new(q_xz.x, qy, Fp2::one());

        // 0 + 0 = 0
        let r = jac_add(&jac_zero, &jac_zero, curve);
        assert!(bool::from(r.ct_equal(&jac_zero)), "Failed 0 + 0 = 0 in jac");

        // 2*0 = 0
        let r = jac_dbl(&jac_zero, curve);
        assert!(bool::from(r.ct_equal(&jac_zero)), "Failed 2*0 = 0 in jac");

        // P + (-P) = 0
        let neg_s = s.neg();
        let r = jac_add(&s, &neg_s, curve);
        assert!(bool::from(r.ct_equal(&jac_zero)), "Failed P - P = 0 in jac");

        // P + 0 = P
        let r = jac_add(&s, &jac_zero, curve);
        assert!(bool::from(r.ct_equal(&s)), "Failed P + 0 = P in jac");
        let r = jac_add(&jac_zero, &s, curve);
        assert!(bool::from(r.ct_equal(&s)), "Failed 0 + P = P in jac");
        let r = jac_add(&s, &jac_zero, curve);
        assert!(
            bool::from(r.ct_equal(&s)),
            "Failed 0 + P = P in jac (third)"
        );

        // P + P = 2P
        let r_dbl = jac_dbl(&s, curve);
        let r_add = jac_add(&s, &s, curve);
        assert!(
            bool::from(r_dbl.ct_equal(&r_add)),
            "Failed P + P = 2*P in jac"
        );

        // Commutativity: T + S = S + T (overwrites T with S + T)
        let r1 = jac_add(&t, &s, curve);
        let t = jac_add(&s, &t, curve);
        assert!(bool::from(r1.ct_equal(&t)), "Failed P + Q = Q + P in jac");

        // Second commutativity with the mutated T
        let r1 = jac_add(&t, &s, curve);
        let r2 = jac_add(&s, &t, curve);
        assert!(
            bool::from(r1.ct_equal(&r2)),
            "Failed P + Q = Q + P in jac (second)"
        );

        // Associativity: (S + T) + R = R + T + S
        // where R = 2*(T+S) to ensure R differs from (T+S)
        let r = jac_add(&s, &t, curve);
        let r = jac_dbl(&r, curve);
        let lhs = {
            let st = jac_add(&s, &t, curve);
            jac_add(&st, &r, curve)
        };
        let rhs = {
            let rt = jac_add(&r, &t, curve);
            jac_add(&rt, &s, curve)
        };
        assert!(
            bool::from(lhs.ct_equal(&rhs)),
            "Failed (P + Q) + R = P + (Q + R) in jac"
        );

        // Weierstrass round-trip for identity
        let (ws_zero, _, ao3) = jac_zero.to_ws(curve);
        let rt_zero = JacPoint::from_ws(&ws_zero, &ao3, curve);
        assert!(
            bool::from(rt_zero.ct_equal(&jac_zero)),
            "Failed Weierstrass round-trip for zero"
        );

        // Weierstrass round-trip for a point
        let (ws_s, _, ao3) = s.to_ws(curve);
        let rt_s = JacPoint::from_ws(&ws_s, &ao3, curve);
        assert!(
            bool::from(rt_s.ct_equal(&s)),
            "Failed Weierstrass round-trip for S"
        );

        let s = jac_dbl(&s, curve);
        let (ws_s2, _, ao3) = s.to_ws(curve);
        let rt_s2 = JacPoint::from_ws(&ws_s2, &ao3, curve);
        assert!(
            bool::from(rt_s2.ct_equal(&s)),
            "Failed Weierstrass round-trip for 2S"
        );

        // Weierstrass doubling for identity
        let (ws_zero, t_val, ao3) = jac_zero.to_ws(curve);
        let (ws_dbl, _) = jac_dbl_ws(&ws_zero, &t_val);
        let rt = JacPoint::from_ws(&ws_dbl, &ao3, curve);
        assert!(
            bool::from(rt.ct_equal(&jac_zero)),
            "Failed 2*0 = 0 in Weierstrass"
        );

        // Weierstrass doubling matches Montgomery doubling (s = 2*S_original here)
        let (ws_s, t_val, ao3) = s.to_ws(curve);
        let (ws_dbl_s, _) = jac_dbl_ws(&ws_s, &t_val);
        let rt_dbl = JacPoint::from_ws(&ws_dbl_s, &ao3, curve);
        let mont_dbl = jac_dbl(&s, curve);
        assert!(
            bool::from(rt_dbl.ct_equal(&mont_dbl)),
            "Failed doubling in Weierstrass"
        );
    }
}
