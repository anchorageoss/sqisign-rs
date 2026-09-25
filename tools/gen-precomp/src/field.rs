//! Build-time field arithmetic: `Fp2` through the library's own backends, and
//! the quadratic extension `Fp4 = Fp2(u)`, `u^2 = nu`, needed to halve points
//! of `E0[2^f]` (their halves live in `E0(Fp4)`).
//!
//! Everything here is variable-time and allocates freely; it only ever runs
//! on public constants inside the generator.

use num_bigint::BigUint;
use sqisign_verify::fp::{Fp, Fp2, FpBackend};
use sqisign_verify::params::Prime;

/// Field operations over some element type, so the curve code can run over
/// both `Fp2` and `Fp4` without duplication.
pub trait FieldOps {
    type E: Clone;
    fn one(&self) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn neg(&self, a: &Self::E) -> Self::E;
    fn inv(&self, a: &Self::E) -> Self::E;
    fn is_zero(&self, a: &Self::E) -> bool;
    fn eq(&self, a: &Self::E, b: &Self::E) -> bool;
    fn sqr(&self, a: &Self::E) -> Self::E {
        self.mul(a, a)
    }
    fn double(&self, a: &Self::E) -> Self::E {
        self.add(a, a)
    }
}

// ---------------------------------------------------------------------
// Fp2 helpers
// ---------------------------------------------------------------------

pub fn fp_to_biguint<L: FpBackend>(x: &Fp<L>) -> BigUint {
    BigUint::from_bytes_le(&x.encode())
}

pub fn prime<L: Prime>() -> BigUint {
    BigUint::from_bytes_le(L::prime_le_bytes())
}

/// `re + im * i` from small integers.
pub fn fp2_small<L: FpBackend>(re: u64, im: u64) -> Fp2<L> {
    Fp2 {
        re: Fp::<L>::from_small(re),
        im: Fp::<L>::from_small(im),
    }
}

pub fn fp2_eq<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    bool::from(a.ct_equal(b))
}

/// The spec's total order on `Fp2`: compare `(re, im)` lexicographically
/// using the minimal non-negative representatives.
pub fn fp2_lex_less<L: FpBackend>(a: &Fp2<L>, b: &Fp2<L>) -> bool {
    let (ar, ai) = (fp_to_biguint(&a.re), fp_to_biguint(&a.im));
    let (br, bi) = (fp_to_biguint(&b.re), fp_to_biguint(&b.im));
    (ar, ai) < (br, bi)
}

/// Variable-time `x^e` in `Fp2` for a big exponent.
pub fn fp2_pow<L: FpBackend>(x: &Fp2<L>, e: &BigUint) -> Fp2<L> {
    x.pow_vartime(&e.to_u64_digits())
}

/// Stateless `FieldOps` over `Fp2`.
pub struct Fp2Ops<L: FpBackend>(core::marker::PhantomData<L>);

impl<L: FpBackend> Fp2Ops<L> {
    pub fn new() -> Self {
        Self(core::marker::PhantomData)
    }
}

impl<L: FpBackend> FieldOps for Fp2Ops<L> {
    type E = Fp2<L>;
    fn one(&self) -> Fp2<L> {
        Fp2::one()
    }
    fn add(&self, a: &Fp2<L>, b: &Fp2<L>) -> Fp2<L> {
        a.add(b)
    }
    fn sub(&self, a: &Fp2<L>, b: &Fp2<L>) -> Fp2<L> {
        a.sub(b)
    }
    fn mul(&self, a: &Fp2<L>, b: &Fp2<L>) -> Fp2<L> {
        a.mul(b)
    }
    fn sqr(&self, a: &Fp2<L>) -> Fp2<L> {
        a.sqr()
    }
    fn neg(&self, a: &Fp2<L>) -> Fp2<L> {
        a.neg()
    }
    fn inv(&self, a: &Fp2<L>) -> Fp2<L> {
        a.inv()
    }
    fn is_zero(&self, a: &Fp2<L>) -> bool {
        bool::from(a.ct_is_zero())
    }
    fn eq(&self, a: &Fp2<L>, b: &Fp2<L>) -> bool {
        fp2_eq(a, b)
    }
}

// ---------------------------------------------------------------------
// Fp4 = Fp2(u), u^2 = nu
// ---------------------------------------------------------------------

/// An element `c0 + c1 u` of `Fp4`.
#[derive(Clone)]
pub struct Fp4<L: FpBackend> {
    pub c0: Fp2<L>,
    pub c1: Fp2<L>,
}

/// The extension's fixed data: the non-square `nu` and `nu^((p-1)/2)`, which
/// is `u^(p-1)` and drives the `p`-power Frobenius.
pub struct Fp4Ctx<L: FpBackend> {
    pub nu: Fp2<L>,
    nu_frob: Fp2<L>,
}

impl<L: FpBackend> Fp4Ctx<L> {
    /// Deterministically pick the first non-square `nu = a + i` with
    /// `a = 1, 2, ...`.
    pub fn new() -> Self {
        let p = prime::<L>();
        let half = (&p - BigUint::from(1u32)) >> 1;
        for a in 1u64..1000 {
            let nu = fp2_small::<L>(a, 1);
            if !bool::from(nu.is_square()) {
                let nu_frob = fp2_pow(&nu, &half);
                // nu^((p-1)/2) squared is nu^(p-1), which must be nu^(p-1) = conj(nu)/nu
                let check = nu_frob.sqr().mul(&nu);
                assert!(fp2_eq(&check, &nu.conjugate()), "nu^p = conj(nu)");
                return Self { nu, nu_frob };
            }
        }
        panic!("no non-square found");
    }

    pub fn lift(&self, a: &Fp2<L>) -> Fp4<L> {
        Fp4 {
            c0: a.clone(),
            c1: Fp2::zero(),
        }
    }

    /// The element as an `Fp2` element, if its `u` part vanishes.
    pub fn descend(&self, a: &Fp4<L>) -> Option<Fp2<L>> {
        if bool::from(a.c1.ct_is_zero()) {
            Some(a.c0.clone())
        } else {
            None
        }
    }

    /// `p`-power Frobenius: `(c0 + c1 u)^p = conj(c0) + conj(c1) nu^((p-1)/2) u`.
    pub fn frobenius(&self, a: &Fp4<L>) -> Fp4<L> {
        Fp4 {
            c0: a.c0.conjugate(),
            c1: a.c1.conjugate().mul(&self.nu_frob),
        }
    }

    /// Square root, or `None` if `a` is not a square. The standard descent:
    /// with `n = c0^2 - nu c1^2` (the norm to `Fp2`) and `s = sqrt(n)`, the
    /// root is `x + (c1 / 2x) u` where `x^2 = (c0 + s)/2` or `(c0 - s)/2`.
    pub fn sqrt(&self, a: &Fp4<L>) -> Option<Fp4<L>> {
        let two_inv = Fp2::<L>::one().half();
        let result = if bool::from(a.c1.ct_is_zero()) {
            if bool::from(a.c0.is_square()) {
                Fp4 {
                    c0: a.c0.sqrt(),
                    c1: Fp2::zero(),
                }
            } else {
                let t = a.c0.mul(&self.nu.inv());
                if !bool::from(t.is_square()) {
                    return None;
                }
                Fp4 {
                    c0: Fp2::zero(),
                    c1: t.sqrt(),
                }
            }
        } else {
            let n = a.c0.sqr().sub(&self.nu.mul(&a.c1.sqr()));
            if !bool::from(n.is_square()) {
                return None;
            }
            let s = n.sqrt();
            let mut t = a.c0.add(&s).mul(&two_inv);
            if !bool::from(t.is_square()) {
                t = a.c0.sub(&s).mul(&two_inv);
                if !bool::from(t.is_square()) {
                    return None;
                }
            }
            let x = t.sqrt();
            let y = a.c1.mul(&x.add(&x).inv());
            Fp4 { c0: x, c1: y }
        };
        if self.eq(&self.sqr(&result), a) {
            Some(result)
        } else {
            None
        }
    }
}

impl<L: FpBackend> FieldOps for Fp4Ctx<L> {
    type E = Fp4<L>;
    fn one(&self) -> Fp4<L> {
        self.lift(&Fp2::one())
    }
    fn add(&self, a: &Fp4<L>, b: &Fp4<L>) -> Fp4<L> {
        Fp4 {
            c0: a.c0.add(&b.c0),
            c1: a.c1.add(&b.c1),
        }
    }
    fn sub(&self, a: &Fp4<L>, b: &Fp4<L>) -> Fp4<L> {
        Fp4 {
            c0: a.c0.sub(&b.c0),
            c1: a.c1.sub(&b.c1),
        }
    }
    fn mul(&self, a: &Fp4<L>, b: &Fp4<L>) -> Fp4<L> {
        let ac = a.c0.mul(&b.c0);
        let bd = a.c1.mul(&b.c1);
        let ad_bc = a.c0.add(&a.c1).mul(&b.c0.add(&b.c1)).sub(&ac).sub(&bd);
        Fp4 {
            c0: ac.add(&self.nu.mul(&bd)),
            c1: ad_bc,
        }
    }
    fn neg(&self, a: &Fp4<L>) -> Fp4<L> {
        Fp4 {
            c0: a.c0.neg(),
            c1: a.c1.neg(),
        }
    }
    fn inv(&self, a: &Fp4<L>) -> Fp4<L> {
        // (c0 + c1 u)^-1 = (c0 - c1 u) / (c0^2 - nu c1^2)
        let n = a.c0.sqr().sub(&self.nu.mul(&a.c1.sqr())).inv();
        Fp4 {
            c0: a.c0.mul(&n),
            c1: a.c1.neg().mul(&n),
        }
    }
    fn is_zero(&self, a: &Fp4<L>) -> bool {
        bool::from(a.c0.ct_is_zero() & a.c1.ct_is_zero())
    }
    fn eq(&self, a: &Fp4<L>, b: &Fp4<L>) -> bool {
        fp2_eq(&a.c0, &b.c0) && fp2_eq(&a.c1, &b.c1)
    }
}

/// Canonical `re || im` hex of an `Fp2` element.
pub fn fp2_hex<L: FpBackend>(a: &Fp2<L>) -> String {
    a.encode().iter().map(|b| format!("{b:02x}")).collect()
}
