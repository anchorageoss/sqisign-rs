//! Deterministic torsion bases: the difference point, the entangled-basis
//! search with hints (spec Algorithms 4.94 and 4.95, `TorsionBasisToHint`
//! and `TorsionBasisFromHint`), the precomputed basis of `E0`, and a basis
//! of the full rational group `E(Fp2) = E[p + 1]` for PRISM.

use super::point::{ec_dbl_iter, ec_mul, xadd, xdbl_a24, xdbl_e0};
#[cfg(feature = "compact")]
use super::JacPoint;
use super::{EcBasis, EcCurve, EcPoint, MAX_ORDER_WORDS};
use crate::fp::{Fp, Fp2, FpBackend};
use crate::precomp::PrimePrecomp;
use subtle::Choice;

/// The two candidates for `x(P +- Q)` from `x(P)` and `x(Q)`, as
/// `(B_XZ + d : B_ZZ)` and `(B_XZ - d : B_ZZ)` where `d^2 = B_XZ^2 - B_XX
/// B_ZZ` (Proposition 3 of ePrint 2017/518). The biquadratics are scaled by
/// `C conj(C)^2 conj(X_P Z_Q - Z_P X_Q)^2` first so that the deterministic
/// square root does not depend on the projective representatives. The
/// reference's basis code uses the second candidate.
#[inline]
pub fn difference_points<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> (EcPoint<L>, EcPoint<L>) {
    let t0 = p.x.mul(&q.x);
    let t1 = p.z.mul(&q.z);
    let bxx = t0.sub(&t1).sqr().mul(&curve.c);
    let bxz = t0.add(&t1);
    let t0 = p.x.mul(&q.z);
    let t1 = p.z.mul(&q.x);
    let bxz = bxz.mul(&t0.add(&t1));
    let bzz = t0.sub(&t1).sqr().mul(&curve.c);
    let bxz = bxz.mul(&curve.c);
    let two_a_term = t0.mul(&t1).mul(&curve.a);
    let bxz = bxz.add(&two_a_term.add(&two_a_term));

    // scale by C conj(C) conj(B_ZZ) = C conj(C)^2 conj(X_P Z_Q - Z_P X_Q)^2
    let scale = curve.c.conjugate().mul(&curve.c).mul(&bzz.conjugate());
    let bxx = bxx.mul(&scale);
    let bxz = bxz.mul(&scale);
    let bzz = bzz.mul(&scale);

    let disc = bxz.sqr().sub(&bxx.mul(&bzz)).sqrt();
    (
        EcPoint::new(bxz.add(&disc), bzz.clone()),
        EcPoint::new(bxz.sub(&disc), bzz),
    )
}

/// The reference's deterministic `x(P - Q)`: the second candidate of
/// [`difference_points`].
#[inline]
pub fn difference_point<L: FpBackend>(
    p: &EcPoint<L>,
    q: &EcPoint<L>,
    curve: &EcCurve<L>,
) -> EcPoint<L> {
    difference_points(p, q, curve).1
}

/// Whether `x` is the x-coordinate of a point of `E: y^2 = x^3 + A x^2 + x`
/// (`C = 1` assumed).
#[inline]
pub fn is_on_curve<L: FpBackend>(x: &Fp2<L>, curve: &EcCurve<L>) -> Choice {
    let t0 = x.add(&curve.a).mul(x).add_one().mul(x);
    t0.is_square()
}

/// Multiply away the odd cofactor and the surplus two-power so that a point
/// of order `k 2^n`, `n` maximal, becomes a point of order `2^f`.
#[inline]
fn clear_cofactor_for_maximal_even_order<L: FpBackend + PrimePrecomp>(
    p: &EcPoint<L>,
    curve: &mut EcCurve<L>,
    f: u32,
) -> EcPoint<L> {
    let cof = [L::COFACTOR];
    let mut res = ec_mul(p, &cof, L::COFACTOR_BITLENGTH, curve);
    for _ in 0..(L::TWO_ADIC_EXPONENT - f) {
        res = xdbl_a24(&res, &curve.a24, curve.is_a24_computed_and_normalized);
    }
    res
}

/// Candidates the entangled-basis searches try before failing closed
/// (`find_nqr_factor`: `b` values; `find_na_x_coord`: multiples of `A`).
/// The reference's counters wrap (16 and 8 bits) and its loops never end;
/// on an honest curve each candidate succeeds with probability about 1/4
/// and 1/2, so the bounded searches fail with probability `(3/4)^1024 <
/// 2^-425` and `2^-256`, while on a crafted curve they return `None`
/// instead of hanging, after about 2 000 quadratic-character tests (tens
/// of milliseconds at level I): for `A` a square with `Re(A^2) = 2`, `z^2 -
/// A^2 i b` lies in `F_p` for every `b` and is always a square in `F_p^2`,
/// so the reference's search has no solution at all (`docs/API_BOUNDS.md`).
pub const BASIS_SEARCH_NQR_CANDIDATES: u32 = 1024;
/// See [`BASIS_SEARCH_NQR_CANDIDATES`].
pub const BASIS_SEARCH_QR_CANDIDATES: u32 = 256;

/// The search for `x(P) = -A / (1 + i b)` with `1 + b^2` a non-square in
/// `Fp`, used when `A` is a square. Returns the projective point `(-A : 1 +
/// i b)` and the hint `b` (0 when `b` does not fit in 7 bits). The counter
/// is 16 bits wide like the reference's; unlike the reference's, the search
/// gives up after [`BASIS_SEARCH_NQR_CANDIDATES`] values of `b`.
#[inline]
fn find_nqr_factor<L: FpBackend>(curve: &EcCurve<L>, start: u8) -> Option<(EcPoint<L>, u8)> {
    let mut n: u16 = start as u16;
    let a_sq = curve.a.sqr();
    let mut tried = 0u32;
    loop {
        // b = n - 1 with 1 + b^2 a non-square in Fp
        loop {
            let tmp = Fp::<L>::from_small((n as u64) * (n as u64) + 1);
            let qr = bool::from(tmp.is_square());
            n = n.wrapping_add(1);
            tried += 1;
            if !qr {
                break;
            }
            if tried >= BASIS_SEARCH_NQR_CANDIDATES {
                return None;
            }
        }
        let b = Fp::<L>::from_small((n.wrapping_sub(1)) as u64);
        let z = Fp2 {
            re: Fp::one(),
            im: b.clone(),
        };
        let t0 = Fp2 {
            re: Fp::zero(),
            im: b,
        };
        // -A/z is on the curve iff A^2 (z - 1) - z^2 is a non-square
        let t = a_sq.mul(&t0).sub(&z.sqr());
        if !bool::from(t.is_square()) {
            let hint = if n <= 128 { (n - 1) as u8 } else { 0 };
            return Some((EcPoint::new(curve.a.neg(), z), hint));
        }
        if tried >= BASIS_SEARCH_NQR_CANDIDATES {
            return None;
        }
    }
}

/// `P = (-A : 1 + i b)` and `Q = (-(A Z_P + X_P) : Z_P)` when `A` is a square.
#[inline]
fn pq_from_nqr<L: FpBackend>(
    curve: &EcCurve<L>,
    start: u8,
) -> Option<(EcPoint<L>, EcPoint<L>, u8)> {
    let (p, hint) = find_nqr_factor(curve, start)?;
    let q = EcPoint::new(curve.a.mul(&p.z).add(&p.x).neg(), p.z.clone());
    Some((p, q, hint))
}

/// The search for `x(P) = n A` on the curve, used when `A` is a non-square.
/// The counter is 8 bits wide like the reference's; unlike the reference's,
/// the search gives up after [`BASIS_SEARCH_QR_CANDIDATES`] multiples.
#[inline]
fn find_na_x_coord<L: FpBackend>(curve: &EcCurve<L>, start: u8) -> Option<(Fp2<L>, u8)> {
    let mut n: u8 = start;
    let mut x = if n == 1 {
        curve.a.clone()
    } else {
        curve.a.mul_small(n as u32)
    };
    let mut tried = 0u32;
    while !bool::from(is_on_curve(&x, curve)) {
        x = x.add(&curve.a);
        n = n.wrapping_add(1);
        tried += 1;
        if tried >= BASIS_SEARCH_QR_CANDIDATES {
            return None;
        }
    }
    let hint = if n < 128 { n } else { 0 };
    Some((x, hint))
}

/// `P = (n A : 1)` and `Q = (-(A + x_P) : 1)` when `A` is a non-square.
#[inline]
fn pq_from_qr<L: FpBackend>(curve: &EcCurve<L>, start: u8) -> Option<(EcPoint<L>, EcPoint<L>, u8)> {
    let (x, hint) = find_na_x_coord(curve, start)?;
    let q = EcPoint::from_x(curve.a.add(&x).neg());
    Some((EcPoint::from_x(x), q, hint))
}

/// The precomputed basis of `E0[2^f]`, scaled down from `E0[2^TWO_ADIC_EXPONENT]`.
#[inline]
fn ec_basis_e0_2f<L: FpBackend + PrimePrecomp>(curve: &EcCurve<L>, f: u32) -> EcBasis<L> {
    let mut p =
        EcPoint::from_x(Fp2::<L>::decode(L::e0_basis_px()).expect("invariant: precomputed x(P0)"));
    let mut q =
        EcPoint::from_x(Fp2::<L>::decode(L::e0_basis_qx()).expect("invariant: precomputed x(Q0)"));
    for _ in 0..(L::TWO_ADIC_EXPONENT - f) {
        p = xdbl_e0(&p);
        q = xdbl_e0(&q);
    }
    let pmq = difference_point(&p, &q, curve);
    EcBasis::new(p, q, pmq)
}

/// `TorsionBasisToHint`: the deterministic basis `(P, Q)` of `E[2^f]` with
/// `[2^(f-1)] Q = (0, 0)`, and the one-byte hint that lets
/// [`ec_curve_to_basis_2f_from_hint`] rebuild it quickly. The curve is
/// normalised in place. `None` when the bounded entangled-basis search
/// finds no point ([`BASIS_SEARCH_NQR_CANDIDATES`]): never on an honest
/// curve, always on a crafted one such as `Re(A^2) = 2`.
///
/// When `e > 0` (and `e < f`), the returned basis is adjusted so that its
/// `[2^(f-e)]` multiple is exactly what `ec_curve_to_basis_2f_from_hint(curve, e,
/// hint)` returns, whose own difference-point choice would otherwise be
/// independent.
#[inline]
pub fn ec_curve_to_basis_2f_to_hint<L: FpBackend + PrimePrecomp>(
    curve: &mut EcCurve<L>,
    f: u32,
    e: u32,
) -> Option<(EcBasis<L>, u8)> {
    curve.normalize_curve_and_a24();
    if bool::from(curve.a.ct_is_zero()) {
        return Some((ec_basis_e0_2f(curve, f), 0));
    }
    let hint_a = bool::from(curve.a.is_square());
    let (p, q, hint) = if !hint_a {
        pq_from_qr(curve, 1)?
    } else {
        pq_from_nqr(curve, 1)?
    };
    let p = clear_cofactor_for_maximal_even_order(&p, curve, f);
    let q = clear_cofactor_for_maximal_even_order(&q, curve, f);

    // The basis is (P, R, Q) with R the deterministic difference: Q above
    // (0, 0) is guaranteed for R = P - Q of the reference's construction.
    let (plus, minus) = difference_points(&p, &q, curve);
    let mut basis_q = minus;
    if e > 0 {
        debug_assert!(e < f);
        let t1 = ec_dbl_iter(&basis_q, (f - e) as usize, curve);
        let t2 = ec_dbl_iter(&plus, (f - e) as usize, curve);
        // (X1 Z2 - X2 Z1) (Z1 Z2) conj(Z1 Z2)^2, compared with the square
        // root of its square: picks the candidate the e-torsion recomputation
        // will pick.
        let t0 = t1.x.mul(&t2.z).sub(&t2.x.mul(&t1.z));
        let zz = t1.z.mul(&t2.z);
        let t0 = t0.mul(&zz.mul(&zz.conjugate().sqr()));
        let root = t0.sqr().sqrt();
        basis_q = super::point::select_point(&basis_q, &plus, t0.ct_equal(&root));
    }
    let basis = EcBasis::new(p, basis_q, q);
    debug_assert!(hint < 128);
    Some((basis, (hint << 1) | hint_a as u8))
}

/// `TorsionBasisFromHint`: rebuild the basis of `E[2^f]` from its hint. The
/// curve is normalised in place. Always succeeds for a hint produced by
/// [`ec_curve_to_basis_2f_to_hint`] on the same curve; for other inputs the
/// result is some pair of points whose order the caller must check. `None`
/// only when the hint asks for the search past 128 and the bounded search
/// fails (a crafted curve, see [`ec_curve_to_basis_2f_to_hint`]).
#[inline]
pub fn ec_curve_to_basis_2f_from_hint<L: FpBackend + PrimePrecomp>(
    curve: &mut EcCurve<L>,
    f: u32,
    hint: u8,
) -> Option<EcBasis<L>> {
    curve.normalize_curve_and_a24();
    if bool::from(curve.a.ct_is_zero()) {
        return Some(ec_basis_e0_2f(curve, f));
    }
    let hint_a = hint & 1 == 1;
    let hint_p = hint >> 1;
    let (p, q) = if hint_p == 0 {
        // the signer did not find a small index; redo the search past 128
        let (p, q, _) = if !hint_a {
            pq_from_qr(curve, 128)?
        } else {
            pq_from_nqr(curve, 128)?
        };
        (p, q)
    } else if !hint_a {
        let x = curve.a.mul_small(hint_p as u32);
        let q = EcPoint::from_x(curve.a.add(&x).neg());
        (EcPoint::from_x(x), q)
    } else {
        let z = Fp2 {
            re: Fp::one(),
            im: Fp::<L>::from_small(hint_p as u64),
        };
        let p = EcPoint::new(curve.a.neg(), z.clone());
        let q = EcPoint::new(curve.a.mul(&z).add(&p.x).neg(), z);
        (p, q)
    };
    let p = clear_cofactor_for_maximal_even_order(&p, curve, f);
    let q = clear_cofactor_for_maximal_even_order(&q, curve, f);
    let basis_q = difference_point(&p, &q, curve);
    Some(EcBasis::new(p, basis_q, q))
}

/// A deterministic basis `(R, T)` of the full rational group `E(Fp2) = E[p + 1]`,
/// with `x(R - T)`, for the ECLIPSE verifier.
///
/// Candidates `x = 1 + m i`, `m = 1, 2, ...`, are tested for lying on the
/// curve; a pair generates iff for every prime `l | p + 1` the images in
/// `E[l]` are independent, which on the Kummer line means the `[(p+1)/l]`
/// multiples are non-zero and not multiples of each other (a handful of
/// comparisons since the odd part of `p + 1` is tiny). Variable time: the
/// curve is public. `None` only if the search exceeds `2^16` candidates.
#[inline]
pub fn ec_basis_p_plus_one<L: FpBackend + PrimePrecomp>(
    curve: &mut EcCurve<L>,
) -> Option<EcBasis<L>> {
    curve.normalize_curve_and_a24();
    let f = L::TWO_ADIC_EXPONENT as usize;
    let cof = L::COFACTOR;
    // odd prime factors of the cofactor
    let mut primes = [0u64; 8];
    let mut nprimes = 0;
    let mut c = cof;
    let mut d = 3u64;
    while d * d <= c {
        if c % d == 0 {
            primes[nprimes] = d;
            nprimes += 1;
            while c % d == 0 {
                c /= d;
            }
        }
        d += 2;
    }
    if c > 1 {
        primes[nprimes] = c;
        nprimes += 1;
    }
    // (p + 1) as limbs, and the bit length
    let mut n = [0u64; MAX_ORDER_WORDS];
    {
        let mut bytes = [0u8; 8 * MAX_ORDER_WORDS];
        let p = L::prime_le_bytes();
        bytes[..p.len()].copy_from_slice(p);
        for (i, w) in n.iter_mut().enumerate() {
            *w = u64::from_le_bytes(bytes[8 * i..8 * i + 8].try_into().expect("8 bytes"));
        }
        let mut carry = 1u64;
        for w in n.iter_mut() {
            let (s, o) = w.overflowing_add(carry);
            *w = s;
            carry = o as u64;
        }
    }
    let nbits = f + (64 - cof.leading_zeros() as usize);

    // [N/l] R for the relevant l; a point generates iff none is zero
    let cofactor_multiple = |pt: &EcPoint<L>, curve: &mut EcCurve<L>, l: u64| -> EcPoint<L> {
        // N / l: for l = 2 shift; for odd l divide the cofactor and keep 2^f
        if l == 2 {
            let mut m = [0u64; MAX_ORDER_WORDS];
            m[0] = cof;
            let r = ec_mul(pt, &m, 64 - cof.leading_zeros() as usize, curve);
            ec_dbl_iter(&r, f - 1, curve)
        } else {
            let mut m = [0u64; MAX_ORDER_WORDS];
            m[0] = cof / l;
            let r = ec_mul(pt, &m, 64 - (cof / l).leading_zeros() as usize, curve);
            ec_dbl_iter(&r, f, curve)
        }
    };
    let full_order = |pt: &EcPoint<L>, curve: &mut EcCurve<L>| -> bool {
        if bool::from(cofactor_multiple(pt, curve, 2).is_zero()) {
            return false;
        }
        for &l in &primes[..nprimes] {
            if bool::from(cofactor_multiple(pt, curve, l).is_zero()) {
                return false;
            }
        }
        // and [N] pt = 0, i.e. the point really is on E (not its twist)
        bool::from(ec_mul(pt, &n, nbits, curve).is_zero())
    };

    let mut m = 1u64;
    let mut first: Option<EcPoint<L>> = None;
    while m < (1 << 16) {
        let x = Fp2 {
            re: Fp::<L>::one(),
            im: Fp::<L>::from_small(m),
        };
        m += 1;
        if !bool::from(is_on_curve(&x, curve)) {
            continue;
        }
        let pt = EcPoint::from_x(x);
        if !full_order(&pt, curve) {
            continue;
        }
        match &first {
            None => first = Some(pt),
            Some(r) => {
                // independence at each prime
                let r2 = cofactor_multiple(r, curve, 2);
                let t2 = cofactor_multiple(&pt, curve, 2);
                if bool::from(r2.ct_equal(&t2)) {
                    continue;
                }
                let mut independent = true;
                for &l in &primes[..nprimes] {
                    let rl = cofactor_multiple(r, curve, l);
                    let tl = cofactor_multiple(&pt, curve, l);
                    // tl must not be +-[k] rl for 1 <= k <= (l-1)/2
                    let mut acc = rl.clone();
                    let mut prev = EcPoint::<L>::identity();
                    for _ in 0..(l - 1) / 2 {
                        if bool::from(acc.ct_equal(&tl)) {
                            independent = false;
                            break;
                        }
                        let next = xadd(&acc, &rl, &prev);
                        prev = acc;
                        acc = next;
                    }
                    if !independent {
                        break;
                    }
                }
                if !independent {
                    continue;
                }
                let pmq = difference_point(r, &pt, curve);
                return Some(EcBasis::new(r.clone(), pt, pmq));
            }
        }
    }
    None
}

// ---- full points for the dimension-4 layer (feature `compact`) ----

#[cfg(feature = "compact")]
/// Recover the y-coordinate of a point on the Montgomery curve
/// `y^2 = x^3 + Ax^2 + x` with `C = 1`.
///
/// Takes the x-coordinate `px` and the curve. Computes `y = sqrt(x^3 + Ax^2 + x)`
/// and returns `(y, is_on_curve)` indicating whether `px` is on the curve.
#[inline]
pub fn ec_recover_y<L: FpBackend>(px: &Fp2<L>, curve: &EcCurve<L>) -> (Fp2<L>, Choice) {
    let t0 = px.sqr();
    let y = t0.mul(&curve.a); // Ax^2
    let y = y.add(px); // Ax^2 + x
    let t0 = t0.mul(px);
    let mut y = y.add(&t0); // x^3 + Ax^2 + x
    let valid = y.sqrt_verify();
    (y, valid)
}

/// Lift an x-only basis `{P, Q, P-Q}` to Jacobian coordinates, assuming
/// `P.z = 1` and `E.C = 1` (normalized curve and point).
///
/// Returns `Choice(1)` if P's x-coordinate is on the curve.
#[cfg(feature = "compact")]
#[inline]
pub fn lift_basis_normalized<L: FpBackend>(
    basis: &mut EcBasis<L>,
    curve: &EcCurve<L>,
) -> (JacPoint<L>, JacPoint<L>, Choice) {
    let (py, ret) = ec_recover_y(&basis.p.x, curve);

    let jp = JacPoint::new(basis.p.x.clone(), py.clone(), Fp2::one());

    // Okeya-Sakurai y-recovery for Q
    let v1 = jp.x.mul(&basis.q.z);
    let v2 = basis.q.x.add(&v1);
    let v3 = basis.q.x.sub(&v1);
    let v3 = v3.sqr();
    let v3 = v3.mul(&basis.pmq.x);

    let v1_2a = curve.a.add(&curve.a);
    let v1_2a = v1_2a.mul(&basis.q.z);
    let v2 = v2.add(&v1_2a);

    let v4 = jp.x.mul(&basis.q.x);
    let v4 = v4.add(&basis.q.z);
    let v2 = v2.mul(&v4);

    let v1_sub = v1_2a.mul(&basis.q.z);
    let v2 = v2.sub(&v1_sub);
    let v2 = v2.mul(&basis.pmq.z);

    let qy_num = v3.sub(&v2);

    let v1_denom = py.add(&py);
    let v1_denom = v1_denom.mul(&basis.q.z);
    let v1_denom = v1_denom.mul(&basis.pmq.z);

    let qx = basis.q.x.mul(&v1_denom);
    let qz = basis.q.z.mul(&v1_denom);

    // Transform to Jacobian coordinates
    let qz_sq = qz.sqr();
    let qy_jac = qy_num.mul(&qz_sq);
    let qx_jac = qx.mul(&qz);

    let jq = JacPoint::new(qx_jac, qy_jac, qz);

    (jp, jq, ret)
}

/// Lift an x-only basis to Jacobian coordinates (general case).
///
/// Normalizes the curve and P.z first, then delegates to
/// `lift_basis_normalized`.
#[cfg(feature = "compact")]
#[inline]
pub fn lift_basis<L: FpBackend>(
    basis: &mut EcBasis<L>,
    curve: &mut EcCurve<L>,
) -> (JacPoint<L>, JacPoint<L>, Choice) {
    // Batch-invert P.z and C
    let mut inverses = [basis.p.z.clone(), curve.c.clone()];
    let mut t1 = [Fp2::<L>::zero(), Fp2::zero()];
    let mut t2 = [Fp2::<L>::zero(), Fp2::zero()];
    Fp2::batched_inv(&mut inverses, &mut t1, &mut t2);

    basis.p.x = basis.p.x.mul(&inverses[0]);
    basis.p.z = Fp2::one();
    curve.a = curve.a.mul(&inverses[1]);
    curve.c = Fp2::one();

    lift_basis_normalized(basis, curve)
}
