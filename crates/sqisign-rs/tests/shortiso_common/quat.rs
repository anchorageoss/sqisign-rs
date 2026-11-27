//! The quaternion algebra `B = (-1, -p)` over `Q` in exact integer
//! arithmetic on `num-bigint`: elements, lattices in Hermite normal form,
//! the norm form of a left ideal of `O0`, and the exact minimum of a
//! positive definite quaternary form. No floating point anywhere; nothing
//! here uses the crate.
//!
//! Conventions: `i^2 = -1`, `j^2 = -p`, `k = i j`; an element is
//! `(c0 + c1 i + c2 j + c3 k) / den`. `O0 = <1, i, (i + j)/2, (1 + k)/2>`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub type Row = [BigInt; 4];

#[derive(Clone, Debug, PartialEq)]
pub struct Elem {
    pub den: BigInt,
    pub c: Row,
}

fn big(x: i64) -> BigInt {
    BigInt::from(x)
}

impl Elem {
    pub fn new(den: &BigInt, c: [&BigInt; 4]) -> Self {
        Self {
            den: den.clone(),
            c: [c[0].clone(), c[1].clone(), c[2].clone(), c[3].clone()],
        }
    }

    pub fn small(den: i64, c: [i64; 4]) -> Self {
        Self {
            den: big(den),
            c: [big(c[0]), big(c[1]), big(c[2]), big(c[3])],
        }
    }

    pub fn scalar(x: &BigInt) -> Self {
        Self {
            den: BigInt::one(),
            c: [x.clone(), BigInt::zero(), BigInt::zero(), BigInt::zero()],
        }
    }

    pub fn mul(&self, o: &Self, p: &BigInt) -> Self {
        let (a, b) = (&self.c, &o.c);
        let c0 = &a[0] * &b[0] - &a[1] * &b[1] - p * (&a[2] * &b[2] + &a[3] * &b[3]);
        let c1 = &a[0] * &b[1] + &a[1] * &b[0] + p * (&a[2] * &b[3] - &a[3] * &b[2]);
        let c2 = &a[0] * &b[2] + &a[2] * &b[0] + &a[3] * &b[1] - &a[1] * &b[3];
        let c3 = &a[0] * &b[3] + &a[3] * &b[0] + &a[1] * &b[2] - &a[2] * &b[1];
        Self {
            den: &self.den * &o.den,
            c: [c0, c1, c2, c3],
        }
    }

    pub fn add(&self, o: &Self) -> Self {
        let mut c: Row = Default::default();
        for (r, slot) in c.iter_mut().enumerate() {
            *slot = &self.c[r] * &o.den + &o.c[r] * &self.den;
        }
        Self {
            den: &self.den * &o.den,
            c,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.c.iter().all(|x| x.is_zero())
    }

    /// `(num, den)` of the reduced norm, not reduced to lowest terms.
    pub fn norm(&self, p: &BigInt) -> (BigInt, BigInt) {
        let c = &self.c;
        (
            &c[0] * &c[0] + &c[1] * &c[1] + p * (&c[2] * &c[2] + &c[3] * &c[3]),
            &self.den * &self.den,
        )
    }
}

/// The basis `1, i, (i + j)/2, (1 + k)/2` of `O0`.
pub fn o0_basis() -> [Elem; 4] {
    [
        Elem::small(1, [1, 0, 0, 0]),
        Elem::small(1, [0, 1, 0, 0]),
        Elem::small(2, [0, 1, 1, 0]),
        Elem::small(2, [1, 0, 0, 1]),
    ]
}

/// A full lattice in `B`: `rows[r] / den`, in Hermite normal form (row `r`
/// is supported on the coordinates `0..=r`, the entries left of each pivot
/// reduced, the content divided out), so that two lattices are equal
/// exactly when their `Lat` are.
#[derive(Clone, Debug, PartialEq)]
pub struct Lat {
    pub den: BigInt,
    pub rows: [Row; 4],
}

/// Hermite normal form of the lattice spanned by `gens` (pairs
/// denominator, coordinates).
pub fn hnf(gens: &[(BigInt, Row)]) -> (BigInt, [Row; 4]) {
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

impl Lat {
    pub fn from_gens(gens: &[Elem]) -> Self {
        let pairs: Vec<(BigInt, Row)> = gens.iter().map(|g| (g.den.clone(), g.c.clone())).collect();
        let (den, rows) = hnf(&pairs);
        Self { den, rows }
    }

    pub fn basis(&self) -> [Elem; 4] {
        std::array::from_fn(|r| Elem {
            den: self.den.clone(),
            c: self.rows[r].clone(),
        })
    }

    pub fn contains(&self, x: &Elem) -> bool {
        // x = v / den with v integral, then back substitution
        let mut v: Row = Default::default();
        for r in 0..4 {
            let (q, rem) = (&x.c[r] * &self.den).div_rem(&x.den);
            if !rem.is_zero() {
                return false;
            }
            v[r] = q;
        }
        for r in (0..4).rev() {
            let (q, rem) = v[r].div_rem(&self.rows[r][r]);
            if !rem.is_zero() {
                return false;
            }
            for c in 0..=r {
                let s = &self.rows[r][c] * &q;
                v[c] -= s;
            }
        }
        v.iter().all(|x| x.is_zero())
    }

    /// The lattice spanned by the products `a b`, `a` in `self`, `b` in
    /// `gens`.
    pub fn times(&self, gens: &[Elem], p: &BigInt) -> Self {
        let mut out = Vec::new();
        for a in self.basis() {
            for b in gens {
                out.push(a.mul(b, p));
            }
        }
        Self::from_gens(&out)
    }

    /// `gens * self` inside `self`.
    pub fn stable_left(&self, gens: &[Elem], p: &BigInt) -> bool {
        self.basis()
            .iter()
            .all(|b| gens.iter().all(|g| self.contains(&g.mul(b, p))))
    }

    /// `self * gens` inside `self`.
    pub fn stable_right(&self, gens: &[Elem], p: &BigInt) -> bool {
        self.basis()
            .iter()
            .all(|b| gens.iter().all(|g| self.contains(&b.mul(g, p))))
    }
}

pub type Mat = [[BigInt; 4]; 4];

pub fn det(m: &[Vec<BigInt>]) -> BigInt {
    let n = m.len();
    if n == 1 {
        return m[0][0].clone();
    }
    let mut acc = BigInt::zero();
    for c in 0..n {
        if m[0][c].is_zero() {
            continue;
        }
        let minor: Vec<Vec<BigInt>> = (1..n)
            .map(|r| {
                (0..n)
                    .filter(|&k| k != c)
                    .map(|k| m[r][k].clone())
                    .collect()
            })
            .collect();
        let term = &m[0][c] * det(&minor);
        if c % 2 == 0 {
            acc += term;
        } else {
            acc -= term;
        }
    }
    acc
}

fn to_vecs(m: &Mat) -> Vec<Vec<BigInt>> {
    m.iter().map(|r| r.to_vec()).collect()
}

pub fn det4(m: &Mat) -> BigInt {
    det(&to_vecs(m))
}

/// The determinant of `m` with row and column `i` deleted.
pub fn principal_minor(m: &Mat, i: usize) -> BigInt {
    let v: Vec<Vec<BigInt>> = (0..4)
        .filter(|&r| r != i)
        .map(|r| {
            (0..4)
                .filter(|&c| c != i)
                .map(|c| m[r][c].clone())
                .collect()
        })
        .collect();
    det(&v)
}

/// The normalised norm form of a left ideal `J` of `O0` (integral or
/// fractional): the matrix `B` with `x^T B x = 2 nrd(sum x_r b_r) / nrd(J)`
/// on the basis `b_r` of `J`, together with `kappa = den^2 nrd(J)`.
///
/// `B` is integral with even diagonal and determinant `p^2`: the values of
/// `x^T B x / 2` are the integers `nrd(beta) / nrd(J)`, `beta` in `J`, which
/// are the norms of the left ideals of `O0` equivalent to `J`, that is the
/// degrees of the isogenies between `E0` and the curve of `J`.
///
/// `nrd(J)` is not taken on trust: `kappa` is the exact fourth root of
/// `det(A) / p^2`, `A` the integral Gram matrix of `trd(x conj(y)) den^2`,
/// which is `den^2 nrd(J)` because `[O0 : J] = nrd(J)^2` and the form
/// `trd(x conj(y))` has determinant `p^2` on `O0`. Every division is
/// checked to be exact.
pub fn norm_form(j: &Lat, p: &BigInt) -> (Mat, BigInt) {
    let mut a: Mat = Default::default();
    for r in 0..4 {
        for s in 0..4 {
            let (x, y) = (&j.rows[r], &j.rows[s]);
            a[r][s] =
                big(2) * (&x[0] * &y[0] + &x[1] * &y[1] + p * (&x[2] * &y[2] + &x[3] * &y[3]));
        }
    }
    let d = det4(&a);
    let p2 = p * p;
    let (q, rem) = d.div_rem(&p2);
    assert!(rem.is_zero() && q.is_positive(), "det(A) = kappa^4 p^2");
    let kappa = q.nth_root(4);
    assert_eq!(kappa.pow(4), q, "det(A) / p^2 is an exact fourth power");
    let mut b: Mat = Default::default();
    for r in 0..4 {
        for s in 0..4 {
            let (v, rem) = a[r][s].div_rem(&kappa);
            assert!(rem.is_zero(), "the norm form is integral");
            b[r][s] = v;
        }
        assert!(b[r][r].is_even(), "even diagonal");
    }
    assert_eq!(det4(&b), p2, "the norm form has determinant p^2");
    (b, kappa)
}

// ------------------------------------------------ exact rationals (LLL)

#[derive(Clone, Debug)]
struct Fr {
    n: BigInt,
    d: BigInt,
}

impl Fr {
    fn new(n: BigInt, d: BigInt) -> Self {
        assert!(!d.is_zero());
        let g = n.gcd(&d);
        let (mut n, mut d) = (n / &g, d / &g);
        if d.is_negative() {
            n = -n;
            d = -d;
        }
        Self { n, d }
    }
    fn int(n: &BigInt) -> Self {
        Self {
            n: n.clone(),
            d: BigInt::one(),
        }
    }
    fn sub(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.d - &o.n * &self.d, &self.d * &o.d)
    }
    fn mul(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.n, &self.d * &o.d)
    }
    fn div(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.d, &self.d * &o.n)
    }
    /// A nearest integer (halves go up; either choice gives a valid size
    /// reduction).
    fn round(&self) -> BigInt {
        (big(2) * &self.n + &self.d).div_floor(&(big(2) * &self.d))
    }
    fn lt(&self, o: &Self) -> bool {
        &self.n * &o.d < &o.n * &self.d
    }
}

/// Gram-Schmidt data of the Gram matrix `g`: `(mu, bstar)`.
fn gso(g: &Mat) -> (Vec<Vec<Fr>>, Vec<Fr>) {
    let zero = Fr::int(&BigInt::zero());
    let mut mu = vec![vec![zero.clone(); 4]; 4];
    let mut bs: Vec<Fr> = Vec::new();
    for i in 0..4 {
        for j in 0..i {
            let mut s = Fr::int(&g[i][j]);
            for k in 0..j {
                s = s.sub(&mu[j][k].mul(&mu[i][k]).mul(&bs[k]));
            }
            mu[i][j] = s.div(&bs[j]);
        }
        let mut s = Fr::int(&g[i][i]);
        for k in 0..i {
            s = s.sub(&mu[i][k].mul(&mu[i][k]).mul(&bs[k]));
        }
        bs.push(s);
    }
    (mu, bs)
}

/// LLL (parameter 3/4) on a positive definite Gram matrix, in exact
/// rational arithmetic. Returns `(G', U)` with `G' = U^T G U`; `U` is
/// built from elementary operations and its determinant and the identity
/// are checked again by the caller, so the correctness of what follows
/// does not depend on this function reducing well.
pub fn lll_gram(g: &Mat) -> (Mat, Mat) {
    let mut g = g.clone();
    let mut u: Mat = Default::default();
    for (r, row) in u.iter_mut().enumerate() {
        for (c, x) in row.iter_mut().enumerate() {
            *x = big((r == c) as i64);
        }
    }
    let three_quarters = Fr::new(big(3), big(4));
    let mut k = 1usize;
    while k < 4 {
        for j in (0..k).rev() {
            let (mu, _) = gso(&g);
            let r = mu[k][j].round();
            if !r.is_zero() {
                // b_k -= r b_j
                for c in 0..4 {
                    let s = &g[j][c] * &r;
                    g[k][c] -= s;
                }
                for row in g.iter_mut() {
                    let s = &row[j] * &r;
                    row[k] -= s;
                }
                for row in u.iter_mut() {
                    let s = &row[j] * &r;
                    row[k] -= s;
                }
            }
        }
        let (mu, bs) = gso(&g);
        let lhs = &bs[k];
        let rhs = three_quarters
            .sub(&mu[k][k - 1].mul(&mu[k][k - 1]))
            .mul(&bs[k - 1]);
        if lhs.lt(&rhs) {
            g.swap(k, k - 1);
            for row in g.iter_mut() {
                row.swap(k, k - 1);
            }
            for row in u.iter_mut() {
                row.swap(k, k - 1);
            }
            k = (k - 1).max(1);
        } else {
            k += 1;
        }
    }
    (g, u)
}

/// What the enumeration found for one form.
#[derive(Clone, Debug)]
pub struct Minimum {
    /// the least value of `x^T B x / 2` over nonzero `x`
    pub min: BigInt,
    /// the number of `x` attaining it (`x` and `-x` both counted)
    pub count: u64,
    /// the half sides of the box that was searched
    pub half_sides: [i64; 4],
    /// the number of nonzero `x` evaluated (one of each pair `x`, `-x`)
    pub evaluated: u64,
}

/// The exact minimum of the positive definite even form `B` (Gram matrix,
/// `q(x) = x^T B x / 2`) over nonzero integer vectors.
///
/// 1. `B' = U^T B U` with `U` from LLL; `det U = +-1` and the identity are
///    recomputed here, so `B'` is the same form in another basis.
/// 2. `M0 = min_r q(e_r) = min_r B'_rr / 2` is a value of the form, so the
///    minimum is at most `M0`.
/// 3. Any `x` with `q(x) <= M0` has `x_r^2 <= 2 M0 (B'^-1)_rr` (Cauchy and
///    Schwarz: `x_r = <B'^-1/2 e_r, B'^1/2 x>`). With `(B'^-1)_rr = C_rr /
///    det B'`, `C_rr` the principal minor, this is `|x_r| <= h_r`, `h_r =
///    floor(sqrt(2 M0 C_rr / det B'))`, computed with integer division and
///    integer square root, each checked.
/// 4. Every nonzero `x` of the box `|x_r| <= h_r` is evaluated in integers.
///    The least value found is the minimum of the form.
pub fn exact_minimum(b: &Mat) -> Minimum {
    let (g, u) = lll_gram(b);
    // 1. same form
    let du = det4(&u);
    assert!(du.is_one() || (-&du).is_one(), "U is unimodular");
    for r in 0..4 {
        for s in 0..4 {
            let mut acc = BigInt::zero();
            for a in 0..4 {
                for c in 0..4 {
                    acc += &u[a][r] * &b[a][c] * &u[c][s];
                }
            }
            assert_eq!(acc, g[r][s], "B' = U^T B U");
        }
    }
    let dg = det4(&g);
    assert!(dg.is_positive());
    assert_eq!(dg, det4(b));
    // positive definite: the leading principal minors are positive
    for n in 1..=4 {
        let lead: Vec<Vec<BigInt>> = (0..n).map(|r| g[r][..n].to_vec()).collect();
        assert!(det(&lead).is_positive(), "positive definite");
    }
    // 2.
    let m0 = (0..4).map(|r| &g[r][r] / big(2)).min().unwrap();
    let bound = big(2) * &m0;
    // 3.
    let mut half_sides = [0i64; 4];
    for (r, side) in half_sides.iter_mut().enumerate() {
        let c = principal_minor(&g, r);
        assert!(c.is_positive());
        let v = (&bound * &c).div_floor(&dg);
        let h = v.sqrt();
        assert!(
            &h * &h <= v && v < (&h + 1) * (&h + 1),
            "integer square root"
        );
        *side = i64::try_from(h).expect("small box");
    }
    let size: u128 = half_sides.iter().map(|&h| 2 * h as u128 + 1).product();
    assert!(size <= 20_000_000, "the box is small: {half_sides:?}");
    // 4.
    let q2 = |x: &[i64; 4]| -> BigInt {
        let mut acc = BigInt::zero();
        for r in 0..4 {
            if x[r] == 0 {
                continue;
            }
            for s in 0..4 {
                if x[s] != 0 {
                    acc += &g[r][s] * big(x[r] * x[s]);
                }
            }
        }
        acc
    };
    let mut best: Option<BigInt> = None;
    let mut count = 0u64;
    let mut evaluated = 0u64;
    let h = half_sides;
    for x0 in -h[0]..=h[0] {
        for x1 in -h[1]..=h[1] {
            for x2 in -h[2]..=h[2] {
                for x3 in -h[3]..=h[3] {
                    let x = [x0, x1, x2, x3];
                    // one of each pair x, -x: the first nonzero entry positive
                    match x.iter().find(|&&v| v != 0) {
                        None => continue,
                        Some(&v) if v < 0 => continue,
                        _ => {}
                    }
                    evaluated += 1;
                    let v = q2(&x);
                    assert!(v.is_positive() && v.is_even());
                    let v = v / big(2);
                    match &best {
                        Some(b) if &v > b => {}
                        Some(b) if &v == b => count += 2,
                        _ => {
                            best = Some(v);
                            count = 2;
                        }
                    }
                }
            }
        }
    }
    let min = best.expect("the box contains the basis vectors");
    assert!(min <= m0);
    Minimum {
        min,
        count,
        half_sides,
        evaluated,
    }
}
