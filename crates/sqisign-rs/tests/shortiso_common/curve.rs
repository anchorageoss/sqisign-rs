//! Short Weierstrass curves `y^2 = x^3 + a x + b` over the generic
//! `F_{p^2}`: Jacobian points with full addition, the endomorphisms `iota`
//! and `pi` of `E: y^2 = x^3 + x` evaluated directly on points, isogenies
//! of degree 2 and of odd prime degree by Velu's formulas, the j-invariant,
//! torsion bases and discrete logarithms in smooth cyclic groups.

use super::field::{Field, F2};
use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Curve {
    pub a: F2,
    pub b: F2,
}

/// Jacobian point `(X : Y : Z)`, affine `(X / Z^2, Y / Z^3)`; `Z = 0` is
/// the point at infinity.
#[derive(Clone, Debug)]
pub struct Pt {
    pub x: F2,
    pub y: F2,
    pub z: F2,
}

impl Curve {
    /// `E: y^2 = x^3 + x`.
    pub fn e0(f: &Field) -> Self {
        Self {
            a: f.one(),
            b: f.zero(),
        }
    }

    /// `1728 * 4 a^3 / (4 a^3 + 27 b^2)`.
    pub fn j(&self, f: &Field) -> F2 {
        let a3 = f.muls(&f.mul(&f.sqr(&self.a), &self.a), 4);
        let den = f.add(&a3, &f.muls(&f.sqr(&self.b), 27));
        assert!(!f.is_zero(&den), "singular curve");
        f.mul(&f.muls(&a3, 1728), &f.inv(&den))
    }

    pub fn inf(&self, f: &Field) -> Pt {
        Pt {
            x: f.one(),
            y: f.one(),
            z: f.zero(),
        }
    }

    pub fn affine(&self, f: &Field, x: &F2, y: &F2) -> Pt {
        let p = Pt {
            x: x.clone(),
            y: y.clone(),
            z: f.one(),
        };
        assert!(self.on_curve(f, &p), "point not on the curve");
        p
    }

    pub fn is_inf(&self, f: &Field, p: &Pt) -> bool {
        f.is_zero(&p.z)
    }

    pub fn on_curve(&self, f: &Field, p: &Pt) -> bool {
        if self.is_inf(f, p) {
            return true;
        }
        // Y^2 = X^3 + a X Z^4 + b Z^6
        let z2 = f.sqr(&p.z);
        let z4 = f.sqr(&z2);
        let z6 = f.mul(&z4, &z2);
        let rhs = f.add(
            &f.add(
                &f.mul(&f.sqr(&p.x), &p.x),
                &f.mul(&f.mul(&self.a, &p.x), &z4),
            ),
            &f.mul(&self.b, &z6),
        );
        f.sqr(&p.y) == rhs
    }

    /// The affine coordinates, or `None` at infinity.
    pub fn to_affine(&self, f: &Field, p: &Pt) -> Option<(F2, F2)> {
        if self.is_inf(f, p) {
            return None;
        }
        let zi = f.inv(&p.z);
        let zi2 = f.sqr(&zi);
        Some((f.mul(&p.x, &zi2), f.mul(&p.y, &f.mul(&zi2, &zi))))
    }

    pub fn eq(&self, f: &Field, p: &Pt, q: &Pt) -> bool {
        match (self.is_inf(f, p), self.is_inf(f, q)) {
            (true, true) => true,
            (true, false) | (false, true) => false,
            _ => {
                let pz2 = f.sqr(&p.z);
                let qz2 = f.sqr(&q.z);
                f.mul(&p.x, &qz2) == f.mul(&q.x, &pz2)
                    && f.mul(&p.y, &f.mul(&qz2, &q.z)) == f.mul(&q.y, &f.mul(&pz2, &p.z))
            }
        }
    }

    pub fn neg(&self, f: &Field, p: &Pt) -> Pt {
        Pt {
            x: p.x.clone(),
            y: f.neg(&p.y),
            z: p.z.clone(),
        }
    }

    pub fn dbl(&self, f: &Field, p: &Pt) -> Pt {
        if self.is_inf(f, p) || f.is_zero(&p.y) {
            return self.inf(f);
        }
        let y2 = f.sqr(&p.y);
        let s = f.muls(&f.mul(&p.x, &y2), 4);
        let z2 = f.sqr(&p.z);
        let m = f.add(&f.muls(&f.sqr(&p.x), 3), &f.mul(&self.a, &f.sqr(&z2)));
        let x3 = f.sub(&f.sqr(&m), &f.muls(&s, 2));
        let y3 = f.sub(&f.mul(&m, &f.sub(&s, &x3)), &f.muls(&f.sqr(&y2), 8));
        let z3 = f.muls(&f.mul(&p.y, &p.z), 2);
        Pt {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Full addition: infinity, equal and opposite points are handled.
    pub fn add(&self, f: &Field, p: &Pt, q: &Pt) -> Pt {
        if self.is_inf(f, p) {
            return q.clone();
        }
        if self.is_inf(f, q) {
            return p.clone();
        }
        let pz2 = f.sqr(&p.z);
        let qz2 = f.sqr(&q.z);
        let u1 = f.mul(&p.x, &qz2);
        let u2 = f.mul(&q.x, &pz2);
        let s1 = f.mul(&p.y, &f.mul(&qz2, &q.z));
        let s2 = f.mul(&q.y, &f.mul(&pz2, &p.z));
        if u1 == u2 {
            return if s1 == s2 {
                self.dbl(f, p)
            } else {
                self.inf(f)
            };
        }
        let h = f.sub(&u2, &u1);
        let r = f.sub(&s2, &s1);
        let h2 = f.sqr(&h);
        let h3 = f.mul(&h2, &h);
        let v = f.mul(&u1, &h2);
        let x3 = f.sub(&f.sub(&f.sqr(&r), &h3), &f.muls(&v, 2));
        let y3 = f.sub(&f.mul(&r, &f.sub(&v, &x3)), &f.mul(&s1, &h3));
        let z3 = f.mul(&f.mul(&p.z, &q.z), &h);
        Pt {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn sub(&self, f: &Field, p: &Pt, q: &Pt) -> Pt {
        self.add(f, p, &self.neg(f, q))
    }

    pub fn mul(&self, f: &Field, k: &BigUint, p: &Pt) -> Pt {
        let mut r = self.inf(f);
        for i in (0..k.bits()).rev() {
            r = self.dbl(f, &r);
            if k.bit(i) {
                r = self.add(f, &r, p);
            }
        }
        r
    }

    pub fn dbl_iter(&self, f: &Field, p: &Pt, n: u64) -> Pt {
        let mut r = p.clone();
        for _ in 0..n {
            r = self.dbl(f, &r);
        }
        r
    }

    /// A point with the given `x`, if `x^3 + a x + b` is a square.
    pub fn lift_x(&self, f: &Field, x: &F2) -> Option<Pt> {
        let rhs = f.add(&f.add(&f.mul(&f.sqr(x), x), &f.mul(&self.a, x)), &self.b);
        let y = f.sqrt(&rhs)?;
        Some(self.affine(f, x, &y))
    }
}

/// `iota (x, y) = (-x, i y)` on `E: y^2 = x^3 + x`.
pub fn iota(f: &Field, p: &Pt) -> Pt {
    Pt {
        x: f.neg(&p.x),
        y: f.mul_i(&p.y),
        z: p.z.clone(),
    }
}

/// The `p`-power Frobenius `pi (x, y) = (x^p, y^p)`.
pub fn frobenius(f: &Field, p: &Pt) -> Pt {
    Pt {
        x: f.frob(&p.x),
        y: f.frob(&p.y),
        z: f.frob(&p.z),
    }
}

/// A prime factorisation `[(l, v)]` of a smooth integer by trial division
/// up to `bound`; `None` if a cofactor remains.
pub fn factor_smooth(n: &BigUint, bound: u32) -> Option<Vec<(u32, u32)>> {
    let mut n = n.clone();
    let mut out = Vec::new();
    let mut l = 2u32;
    while l <= bound && !n.is_one() {
        let mut v = 0;
        while (&n % l).is_zero() {
            n /= l;
            v += 1;
        }
        if v > 0 {
            out.push((l, v));
        }
        l += 1;
    }
    n.is_one().then_some(out)
}

/// Whether `P` has order exactly `n`, for `n` with the given factorisation.
pub fn has_exact_order(e: &Curve, f: &Field, p: &Pt, n: &BigUint, fac: &[(u32, u32)]) -> bool {
    if !e.is_inf(f, &e.mul(f, n, p)) {
        return false;
    }
    fac.iter()
        .all(|&(l, _)| !e.is_inf(f, &e.mul(f, &(n / l), p)))
}

/// Whether `P, Q` of order dividing `n` form a basis of `E[n]`: for every
/// prime `l | n` the points `[n/l] P`, `[n/l] Q` are independent in `E[l]`
/// (`l` is small: the `l` multiples are enumerated).
pub fn is_basis(e: &Curve, f: &Field, p: &Pt, q: &Pt, n: &BigUint, fac: &[(u32, u32)]) -> bool {
    if !e.is_inf(f, &e.mul(f, n, p)) || !e.is_inf(f, &e.mul(f, n, q)) {
        return false;
    }
    for &(l, _) in fac {
        let pl = e.mul(f, &(n / l), p);
        let ql = e.mul(f, &(n / l), q);
        if e.is_inf(f, &pl) || e.is_inf(f, &ql) {
            return false;
        }
        let mut acc = e.inf(f);
        for _ in 0..l {
            if e.eq(f, &acc, &ql) {
                return false;
            }
            acc = e.add(f, &acc, &pl);
        }
    }
    true
}

/// A basis of `E[n]` on `E: y^2 = x^3 + x`, `n | p + 1`, found
/// deterministically: points with `x = k + 2i`, `k = 1, 2, ...`, cleared of
/// the cofactor `(p + 1) / n` (the group `E(F_{p^2})` is `(Z / (p+1))^2`).
/// (With `x = k + i` the factor `x - i` is always a square and every such
/// point lies above the same point of order 2: no basis.)
/// `skip` discards that many valid first points, to vary the basis.
pub fn torsion_basis(
    e: &Curve,
    f: &Field,
    n: &BigUint,
    fac: &[(u32, u32)],
    skip: usize,
) -> (Pt, Pt) {
    let p1 = &f.p + 1u32;
    assert!((&p1 % n).is_zero(), "n | p + 1");
    let cof = &p1 / n;
    let mut k = 0u64;
    let mut next = |want: &dyn Fn(&Pt) -> bool| loop {
        k += 1;
        assert!(k < 5000, "no torsion point found");
        let x = f.from_parts(&BigUint::from(k), &BigUint::from(2u32));
        let Some(r) = e.lift_x(f, &x) else { continue };
        let t = e.mul(f, &cof, &r);
        if want(&t) {
            return t;
        }
    };
    let mut skipped = 0;
    let p = loop {
        let t = next(&|t| has_exact_order(e, f, t, n, fac));
        if skipped == skip {
            break t;
        }
        skipped += 1;
    };
    let q = next(&|t| is_basis(e, f, &p, t, n, fac));
    (p, q)
}

/// The discrete logarithm of `h` in the cyclic group generated by `g` of
/// exact order `n` (smooth, factorisation `fac`), by Pohlig and Hellman, or
/// `None` if `h` is not in `<g>`. The answer is verified.
pub fn dlog_smooth(
    e: &Curve,
    f: &Field,
    g: &Pt,
    h: &Pt,
    n: &BigUint,
    fac: &[(u32, u32)],
) -> Option<BigUint> {
    if !e.is_inf(f, &e.mul(f, n, h)) {
        return None;
    }
    let mut x = BigUint::zero();
    let mut modulus = BigUint::one();
    for &(l, v) in fac {
        let lv = BigUint::from(l).pow(v);
        let cof = n / &lv;
        let gl = e.mul(f, &cof, g);
        let mut rem = e.mul(f, &cof, h);
        // powers[k] = [l^k] gl
        let mut powers = vec![gl.clone()];
        for k in 1..v as usize {
            let prev = powers[k - 1].clone();
            powers.push(e.mul(f, &BigUint::from(l), &prev));
        }
        let top = &powers[v as usize - 1]; // order l
        let mut xl = BigUint::zero();
        for k in 0..v {
            // digit d with [l^(v-1-k)] rem = [d] top
            let t = e.mul(f, &BigUint::from(l).pow(v - 1 - k), &rem);
            let mut acc = e.inf(f);
            let mut digit = None;
            for d in 0..l {
                if e.eq(f, &acc, &t) {
                    digit = Some(d);
                    break;
                }
                acc = e.add(f, &acc, top);
            }
            let d = digit?;
            if d != 0 {
                let step = e.mul(f, &BigUint::from(d), &powers[k as usize]);
                rem = e.sub(f, &rem, &step);
                xl += BigUint::from(d) * BigUint::from(l).pow(k);
            }
        }
        if !e.is_inf(f, &rem) {
            return None;
        }
        // CRT: x = x mod modulus, x = xl mod lv
        let inv = modulus.modpow(&(&lv - &lv / l - 1u32), &lv); // modulus^(phi(l^v) - 1)
        let diff = ((&lv + &xl) - (&x % &lv)) % &lv;
        x += &modulus * ((diff * inv) % &lv);
        modulus *= &lv;
    }
    e.eq(f, &e.mul(f, &x, g), h).then_some(x)
}

/// The isogeny with kernel `<K>`, `K` of prime order `l` (affine multiples
/// of `K` are computed here), by Velu's formulas. Returns the codomain and
/// the data to evaluate it.
pub struct Isogeny {
    pub l: u32,
    /// `(x_Q, t_Q, u_Q)` over `Q = K, 2K, ..., ((l-1)/2) K`; for `l = 2`
    /// the single entry `(x_K, t, 0)`.
    terms: Vec<(F2, F2, F2)>,
    pub codomain: Curve,
}

impl Isogeny {
    pub fn new(e: &Curve, f: &Field, k: &Pt, l: u32) -> Self {
        assert!(!e.is_inf(f, k), "kernel generator at infinity");
        assert!(
            e.is_inf(f, &e.mul(f, &BigUint::from(l), k)),
            "kernel generator of order {l}"
        );
        let mut terms = Vec::new();
        let mut t_sum = f.zero();
        let mut w_sum = f.zero();
        if l == 2 {
            let (x0, y0) = e.to_affine(f, k).unwrap();
            assert!(f.is_zero(&y0));
            let t = f.add(&f.muls(&f.sqr(&x0), 3), &e.a);
            t_sum = t.clone();
            w_sum = f.mul(&x0, &t);
            terms.push((x0, t, f.zero()));
        } else {
            assert!(l % 2 == 1);
            let mut q = k.clone();
            for _ in 0..(l - 1) / 2 {
                let (xq, yq) = e.to_affine(f, &q).unwrap();
                let gx = f.add(&f.muls(&f.sqr(&xq), 3), &e.a);
                let t = f.muls(&gx, 2);
                let u = f.muls(&f.sqr(&yq), 4);
                t_sum = f.add(&t_sum, &t);
                w_sum = f.add(&w_sum, &f.add(&u, &f.mul(&xq, &t)));
                terms.push((xq, t, u));
                q = e.add(f, &q, k);
            }
        }
        let codomain = Curve {
            a: f.sub(&e.a, &f.muls(&t_sum, 5)),
            b: f.sub(&e.b, &f.muls(&w_sum, 7)),
        };
        Self { l, terms, codomain }
    }

    /// The image of a point of the domain.
    pub fn eval(&self, f: &Field, p: &Pt) -> Pt {
        let cod = &self.codomain;
        if f.is_zero(&p.z) {
            return cod.inf(f);
        }
        if self.l == 2 {
            // x' = x + t / (x - x0), y' = y (1 - t / (x - x0)^2), in
            // Jacobian coordinates with U = X - x0 Z^2 and Z' = Z U.
            let (x0, t, _) = &self.terms[0];
            let z2 = f.sqr(&p.z);
            let z4 = f.sqr(&z2);
            let u = f.sub(&p.x, &f.mul(x0, &z2));
            let tz4 = f.mul(t, &z4);
            let x = f.mul(&f.add(&f.mul(&p.x, &u), &tz4), &u);
            let y = f.mul(&f.mul(&p.y, &f.sub(&f.sqr(&u), &tz4)), &u);
            let z = f.mul(&p.z, &u);
            return Pt { x, y, z };
        }
        // affine: x' = x + sum t/(x - xQ) + u/(x - xQ)^2,
        //         y' = y (1 - sum t/(x - xQ)^2 + 2u/(x - xQ)^3)
        let zi = f.inv(&p.z);
        let zi2 = f.sqr(&zi);
        let x = f.mul(&p.x, &zi2);
        let y = f.mul(&p.y, &f.mul(&zi2, &zi));
        let mut xs = x.clone();
        let mut ys = f.one();
        for (xq, t, u) in &self.terms {
            let d = f.sub(&x, xq);
            if f.is_zero(&d) {
                return cod.inf(f);
            }
            let di = f.inv(&d);
            let di2 = f.sqr(&di);
            let di3 = f.mul(&di2, &di);
            xs = f.add(&xs, &f.add(&f.mul(t, &di), &f.mul(u, &di2)));
            ys = f.sub(&ys, &f.add(&f.mul(t, &di2), &f.muls(&f.mul(u, &di3), 2)));
        }
        Pt {
            x: xs,
            y: f.mul(&y, &ys),
            z: f.one(),
        }
    }
}
