//! Kernels of `(2^n, 2^n)`-isogenies from `E0 x E0`, built from
//! endomorphisms of `E0` through Kani's lemma, for testing the theta chains
//! of `sqisign_verify::theta`. This is the construction the SQIsign reference
//! uses in `dim2id2iso` with the output of Qlapoty.
//!
//! For an odd `u` and `theta` in `End(E0)` of norm `u (2^n - u)`, the group
//! `{([u] P, theta(P)) : P in E0[2^n]}` is a maximal isotropic subgroup.
//! Writing `theta = psi o phi` with `deg phi = u`, it is the graph
//! `{(phi^(Q), psi(Q)) : Q in E'[2^n]}` of two isogenies from `E'` whose
//! degrees sum to `2^n`, so by Kani's lemma the codomain is a product of
//! elliptic curves. The chain wants generators of order `2^(n + 2)`, so
//! the construction is applied to a basis of `E0[2^(n + 2)]`.
//!
//! `theta` must be taken in the full maximal order `O0 = Z<1, i, (i + j)/2,
//! (1 + k)/2>`: every element of the suborder `Z<1, i, j, k>` with odd norm
//! acts on `E0[2]` as an automorphism of `E0` (`j` acts like `i` there), so
//! the first `(2, 2)`-step would land on a product again and the theta
//! formulas break down. With `theta = a + b i + c (i + j)/2 + d (1 + k)/2`,
//! `4 N(theta) = (2a + d)^2 + (2b + c)^2 + p (c^2 + d^2)`; the search takes
//! `d` odd, `c` even, and needs `4 u (2^n - u) > p`, so lengths above about
//! `log2(p) / 2` are reachable, which covers what signature verification uses.

#![allow(dead_code)]

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Zero};
use sqisign_verify::ec::point::{ec_biscalar_mul, ec_mul};
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::precomp::PrimePrecomp;
use sqisign_verify::theta::{ThetaCouplePoint, ThetaKernelCouplePoints};

/// A `2x2` matrix over `Z / 2^f` as big integers, `[[a(P0), a(Q0)], [b(P0), b(Q0)]]`.
pub type Mat = [[BigUint; 2]; 2];

pub fn mat_from_bytes<const N: usize>(m: &[[[u8; N]; 2]; 2]) -> Mat {
    [
        [
            BigUint::from_bytes_le(&m[0][0]),
            BigUint::from_bytes_le(&m[0][1]),
        ],
        [
            BigUint::from_bytes_le(&m[1][0]),
            BigUint::from_bytes_le(&m[1][1]),
        ],
    ]
}

pub fn big_to_limbs(v: &BigUint, words: usize) -> Vec<u64> {
    let mut out = v.to_u64_digits();
    assert!(out.len() <= words, "scalar too large");
    out.resize(words, 0);
    out
}

pub fn pow2(k: u32) -> BigUint {
    BigUint::one() << k
}

/// The prime as a big integer.
pub fn prime<L: FpBackend>() -> BigUint {
    BigUint::from_bytes_le(L::prime_le_bytes())
}

/// The curve `E0: y^2 = x^3 + x` with its constant normalised.
pub fn e0<L: FpBackend>() -> EcCurve<L> {
    let mut e = EcCurve::<L>::default();
    e.normalize_a24();
    e
}

/// The precomputed basis `(P0, Q0, P0 - Q0)` of `E0[2^f]`.
pub fn e0_basis<L: FpBackend + PrimePrecomp>() -> EcBasis<L> {
    let dec = |b: &[u8]| EcPoint::from_x(Fp2::<L>::decode(b).expect("precomp x-coordinate"));
    EcBasis::new(
        dec(L::e0_basis_px()),
        dec(L::e0_basis_qx()),
        dec(L::e0_basis_pmqx()),
    )
}

/// The action matrices of the `O0` generators `i`, `(i + j)/2`, `(1 + k)/2`.
pub struct Gens {
    pub g2: Mat,
    pub g3: Mat,
    pub g4: Mat,
}

/// An endomorphism `a + b i + c (i + j)/2 + d (1 + k)/2` of `E0`.
#[derive(Clone, Debug)]
pub struct Endo {
    pub a: BigInt,
    pub b: BigInt,
    pub c: BigInt,
    pub d: BigInt,
}

fn reduce(v: &BigInt, modulus: &BigUint) -> BigUint {
    let m = BigInt::from(modulus.clone());
    v.mod_floor(&m).to_biguint().expect("non-negative")
}

impl Endo {
    /// Its action on `(P0, Q0)` modulo `2^f`.
    pub fn matrix(&self, g: &Gens, f: u32) -> Mat {
        let modulus = pow2(f);
        let (b, c, d) = (
            reduce(&self.b, &modulus),
            reduce(&self.c, &modulus),
            reduce(&self.d, &modulus),
        );
        let a = reduce(&self.a, &modulus);
        let mut out: Mat = Default::default();
        for (r, row) in out.iter_mut().enumerate() {
            for (col, entry) in row.iter_mut().enumerate() {
                let mut acc = &b * &g.g2[r][col] + &c * &g.g3[r][col] + &d * &g.g4[r][col];
                if r == col {
                    acc += &a;
                }
                *entry = acc % &modulus;
            }
        }
        out
    }

    /// `((2a + d)^2 + (2b + c)^2 + p (c^2 + d^2)) / 4`.
    pub fn norm(&self, p: &BigUint) -> BigUint {
        let x = &self.a * 2 + &self.d;
        let y = &self.b * 2 + &self.c;
        let n: BigInt =
            &x * &x + &y * &y + BigInt::from(p.clone()) * (&self.c * &self.c + &self.d * &self.d);
        let (q, r) = n.div_rem(&BigInt::from(4));
        assert!(r.is_zero(), "norm not integral");
        q.to_biguint().expect("positive norm")
    }
}

/// Miller-Rabin with fixed bases; enough for test parameter search.
pub fn is_probable_prime(m: &BigUint) -> bool {
    if m < &BigUint::from(4u32) {
        return *m == BigUint::from(2u32) || *m == BigUint::from(3u32);
    }
    if m.is_even() {
        return false;
    }
    let one = BigUint::one();
    let m1 = m - &one;
    let s = m1.trailing_zeros().unwrap_or(0);
    let d = &m1 >> s;
    'outer: for base in [
        2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    ] {
        let b = BigUint::from(base);
        if &b % m == BigUint::zero() {
            continue;
        }
        let mut x = b.modpow(&d, m);
        if x == one || x == m1 {
            continue;
        }
        for _ in 1..s {
            x = x.modpow(&BigUint::from(2u32), m);
            if x == m1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// `(x, y)` with `x^2 + y^2 = m` for a prime `m = 1 mod 4` (Cornacchia).
pub fn two_squares_prime(m: &BigUint) -> Option<(BigUint, BigUint)> {
    let one = BigUint::one();
    let m1 = m - &one;
    let half = &m1 >> 1u32;
    let quarter = &m1 >> 2u32;
    let mut c = BigUint::from(2u32);
    let r = loop {
        if c.modpow(&half, m) == m1 {
            break c.modpow(&quarter, m);
        }
        c += &one;
    };
    // Euclid until the remainder drops below sqrt(m)
    let mut a = m.clone();
    let mut b = r;
    while &b * &b > *m {
        let t = &a % &b;
        a = b;
        b = t;
    }
    let x = b;
    let y2 = m - &x * &x;
    let y = y2.sqrt();
    if &y * &y == y2 {
        Some((x, y))
    } else {
        None
    }
}

/// An odd `u` and `theta` in `O0` of norm `u (2^n - u)` for a chain of length
/// `n`, whose action on `E0[2]` is not that of an automorphism. The search
/// starts from an odd `u` near `2^(n - 1)` offset by `seed`.
pub fn find_endomorphism(n: u32, p: &BigUint, g: &Gens, f: u32, seed: u32) -> (BigUint, Endo) {
    assert!(n >= 2);
    let two_n = pow2(n);
    let mut u = pow2(n - 1) + BigUint::from(2 * seed + 1);
    let two = BigUint::from(2u32);
    let mod2 = |m: &Mat| -> [[u8; 2]; 2] {
        let f = |x: &BigUint| (x % &two).to_u64_digits().first().copied().unwrap_or(0) as u8;
        [[f(&m[0][0]), f(&m[0][1])], [f(&m[1][0]), f(&m[1][1])]]
    };
    let id2 = [[1u8, 0], [0, 1]];
    let i2 = mod2(&g.g2);
    loop {
        assert!(u < two_n, "no endomorphism found");
        let t = &u * (&two_n - &u);
        let four_t = &t << 2u32;
        assert!(
            four_t > *p,
            "length {n} too short: 4 u (2^n - u) must exceed p"
        );
        // c even, d odd: then 4T - p (c^2 + d^2) = 1 mod 4
        for (c, d) in [
            (0u32, 1u32),
            (2, 1),
            (0, 3),
            (2, 3),
            (4, 1),
            (4, 3),
            (0, 5),
            (2, 5),
        ] {
            let pcd = p * BigUint::from(c * c + d * d);
            if pcd >= four_t {
                continue;
            }
            let m = &four_t - &pcd;
            debug_assert_eq!(&m % 4u32, BigUint::one());
            if !is_probable_prime(&m) {
                continue;
            }
            let Some((x, y)) = two_squares_prime(&m) else {
                continue;
            };
            // one of x, y is odd; 2a + d = odd one, 2b + c = even one
            let (x, y) = if x.is_odd() { (x, y) } else { (y, x) };
            assert!(x.is_odd() && y.is_even());
            let e = Endo {
                a: (BigInt::from(x) - BigInt::from(d)) / 2,
                b: (BigInt::from(y) - BigInt::from(c)) / 2,
                c: BigInt::from(c),
                d: BigInt::from(d),
            };
            assert_eq!(e.norm(p), t);
            let m2 = mod2(&e.matrix(g, f));
            if m2 == id2 || m2 == i2 {
                continue;
            }
            return (u, e);
        }
        u += 2u32;
    }
}

/// Everything a test needs about a Kani kernel on `E0 x E0`.
pub struct KaniKernel<L: FpBackend> {
    pub n: u32,
    pub u: BigUint,
    pub endo: Endo,
    /// Kernel generators of order `2^(n + 2)`.
    pub ker: ThetaKernelCouplePoints<L>,
    /// The basis `(P', Q', P' - Q')` of `E0[2^(n + 2)]` used.
    pub basis_n2: EcBasis<L>,
}

/// Build the kernel `([u] P', theta(P')), ([u] Q', theta(Q'))` for the
/// basis `(P', Q') = [2^(f - n - 2)] (P0, Q0)`.
pub fn kani_kernel<L: FpBackend + PrimePrecomp>(
    n: u32,
    u: &BigUint,
    endo: &Endo,
    f: u32,
    words: usize,
    g: &Gens,
) -> KaniKernel<L> {
    assert!(n + 2 <= f);
    let mut e = e0::<L>();
    let b0 = e0_basis::<L>();
    let modulus = pow2(f);
    let s = pow2(f - n - 2);
    let m = endo.matrix(g, f);
    let fbits = f as usize;

    let e_ro = e.clone();
    let mul = |k: &BigUint, pt: &EcPoint<L>, e: &mut EcCurve<L>| {
        let k = k % &modulus;
        ec_mul(pt, &big_to_limbs(&k, words), fbits, e)
    };
    let bimul = |a: &BigUint, b: &BigUint| {
        let a = a % &modulus;
        let b = b % &modulus;
        ec_biscalar_mul(
            &big_to_limbs(&a, words),
            &big_to_limbs(&b, words),
            fbits,
            &b0,
            &e_ro,
        )
        .expect("biscalar")
    };

    let p_n2 = mul(&s, &b0.p, &mut e);
    let q_n2 = mul(&s, &b0.q, &mut e);
    let pmq_n2 = mul(&s, &b0.pmq, &mut e);

    let us = u * &s;
    let t1 = ThetaCouplePoint::new(
        mul(&us, &b0.p, &mut e),
        bimul(&(&s * &m[0][0]), &(&s * &m[1][0])),
    );
    let t2 = ThetaCouplePoint::new(
        mul(&us, &b0.q, &mut e),
        bimul(&(&s * &m[0][1]), &(&s * &m[1][1])),
    );
    KaniKernel {
        n,
        u: u.clone(),
        endo: endo.clone(),
        ker: ThetaKernelCouplePoints { t1, t2 },
        basis_n2: EcBasis::new(p_n2, q_n2, pmq_n2),
    }
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub fn fp2_hex<L: FpBackend>(a: &Fp2<L>) -> String {
    hex(&a.encode())
}

pub fn fp2_from_hex<L: FpBackend>(s: &str) -> Fp2<L> {
    let bytes: Vec<u8> = (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap())
        .collect();
    Fp2::<L>::decode(&bytes).expect("fp2 hex")
}

pub fn point_hex<L: FpBackend>(p: &EcPoint<L>) -> String {
    format!("{} {}", fp2_hex(&p.x), fp2_hex(&p.z))
}

pub fn couple_hex<L: FpBackend>(p: &ThetaCouplePoint<L>) -> String {
    format!("{} {}", point_hex(&p.p1), point_hex(&p.p2))
}

pub fn points_equal<L: FpBackend>(a: &EcPoint<L>, b: &EcPoint<L>) -> bool {
    bool::from(a.ct_equal(b))
}

pub fn curves_equal<L: FpBackend>(a: &EcCurve<L>, b: &EcCurve<L>) -> bool {
    bool::from(a.a.mul(&b.c).ct_equal(&b.a.mul(&a.c)))
}
