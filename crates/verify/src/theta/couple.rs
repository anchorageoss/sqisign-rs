use crate::ec::jacobian::{jac_add, jac_dbl, jac_dbl_ws};
use crate::ec::point::{ec_dbl, test_point_order_twof};
use crate::ec::EcBasis;
use crate::fp::FpBackend;

use super::{ThetaCoupleCurve, ThetaCoupleJacPoint, ThetaCouplePoint, ThetaKernelCouplePoints};

/// Double a couple point on E1 x E2.
///
/// `(P1, P2) -> ([2]P1, [2]P2)`
#[inline]
pub fn double_couple_point<L: FpBackend>(
    p: &ThetaCouplePoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> ThetaCouplePoint<L> {
    ThetaCouplePoint {
        p1: ec_dbl(&p.p1, &e12.e1),
        p2: ec_dbl(&p.p2, &e12.e2),
    }
}

/// Iterated doubling of a couple point: \[2ⁿ\](P1, P2).
#[inline]
pub fn double_couple_point_iter<L: FpBackend>(
    p: &ThetaCouplePoint<L>,
    n: u32,
    e12: &ThetaCoupleCurve<L>,
) -> ThetaCouplePoint<L> {
    if n == 0 {
        p.clone()
    } else {
        let mut out = double_couple_point(p, e12);
        for _ in 1..n {
            out = double_couple_point(&out, e12);
        }
        out
    }
}

/// Add two Jacobian couple points on E1 x E2.
///
/// (P1, P2) + (Q1, Q2) = (P1+Q1, P2+Q2)
#[inline]
pub fn add_couple_jac_points<L: FpBackend>(
    t1: &ThetaCoupleJacPoint<L>,
    t2: &ThetaCoupleJacPoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> ThetaCoupleJacPoint<L> {
    ThetaCoupleJacPoint {
        p1: jac_add(&t1.p1, &t2.p1, &e12.e1),
        p2: jac_add(&t1.p2, &t2.p2, &e12.e2),
    }
}

/// Double a Jacobian couple point on E1 x E2.
#[inline]
pub fn double_couple_jac_point<L: FpBackend>(
    p: &ThetaCoupleJacPoint<L>,
    e12: &ThetaCoupleCurve<L>,
) -> ThetaCoupleJacPoint<L> {
    ThetaCoupleJacPoint {
        p1: jac_dbl(&p.p1, &e12.e1),
        p2: jac_dbl(&p.p2, &e12.e2),
    }
}

/// Iterated Jacobian doubling with Weierstrass optimization for n >= 2.
#[inline]
pub fn double_couple_jac_point_iter<L: FpBackend>(
    p: &ThetaCoupleJacPoint<L>,
    n: u32,
    e12: &ThetaCoupleCurve<L>,
) -> ThetaCoupleJacPoint<L> {
    if n == 0 {
        p.clone()
    } else if n == 1 {
        double_couple_jac_point(p, e12)
    } else {
        let (mut q1, mut t1, a1) = p.p1.to_ws(&e12.e1);
        let (mut q2, mut t2, a2) = p.p2.to_ws(&e12.e2);

        for _ in 0..n {
            let (new_q1, new_t1) = jac_dbl_ws(&q1, &t1);
            q1 = new_q1;
            t1 = new_t1;
            let (new_q2, new_t2) = jac_dbl_ws(&q2, &t2);
            q2 = new_q2;
            t2 = new_t2;
        }

        ThetaCoupleJacPoint {
            p1: crate::ec::JacPoint::from_ws(&q1, &a1, &e12.e1),
            p2: crate::ec::JacPoint::from_ws(&q2, &a2, &e12.e2),
        }
    }
}

/// Convert Jacobian couple points to (X : Z) couple points.
#[inline]
pub fn couple_jac_to_xz<L: FpBackend>(p: &ThetaCoupleJacPoint<L>) -> ThetaCouplePoint<L> {
    ThetaCouplePoint {
        p1: p.p1.to_xz(),
        p2: p.p2.to_xz(),
    }
}

/// Pack two EC bases B1, B2 into a kernel triple (T1, T2, T1-T2).
///
/// T1 = (B1.P, B2.P), T2 = (B1.Q, B2.Q), T1-T2 = (B1.PmQ, B2.PmQ)
#[inline]
pub fn copy_bases_to_kernel<L: FpBackend>(
    b1: &EcBasis<L>,
    b2: &EcBasis<L>,
) -> ThetaKernelCouplePoints<L> {
    ThetaKernelCouplePoints {
        t1: ThetaCouplePoint {
            p1: b1.p.clone(),
            p2: b2.p.clone(),
        },
        t2: ThetaCouplePoint {
            p1: b1.q.clone(),
            p2: b2.q.clone(),
        },
        t1m2: ThetaCouplePoint {
            p1: b1.pmq.clone(),
            p2: b2.pmq.clone(),
        },
    }
}

/// Test if both points in a couple point have order exactly 2ᵗ.
#[inline]
pub fn test_couple_point_order_twof<L: FpBackend>(
    p: &ThetaCouplePoint<L>,
    e: &ThetaCoupleCurve<L>,
    t: usize,
) -> subtle::Choice {
    let check_p1 = test_point_order_twof(&p.p1, &e.e1, t);
    let check_p2 = test_point_order_twof(&p.p2, &e.e2, t);
    check_p1 & check_p2
}
