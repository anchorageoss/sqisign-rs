//! Algebraic properties of the theta-coordinate primitives. The chains
//! themselves are tested in `crates/prism` where the `E0` endomorphism
//! matrices needed to build kernels live.

mod common;

use common::{eq2, DetRng, ITER};
use sqisign_verify::fp::{Fp2, FpBackend};
use sqisign_verify::theta::structure::{
    hadamard, invert_point, is_product_theta_point, pointwise_product, pointwise_square,
    to_squared_theta,
};
use sqisign_verify::theta::ThetaPoint;

fn random_point<L: FpBackend>(rng: &mut DetRng) -> ThetaPoint<L> {
    ThetaPoint::new(
        rng.random_fp2(),
        rng.random_fp2(),
        rng.random_fp2(),
        rng.random_fp2(),
    )
}

fn proj_eq<L: FpBackend>(a: &ThetaPoint<L>, b: &ThetaPoint<L>) -> bool {
    eq2(&a.x.mul(&b.y), &a.y.mul(&b.x))
        && eq2(&a.x.mul(&b.z), &a.z.mul(&b.x))
        && eq2(&a.x.mul(&b.t), &a.t.mul(&b.x))
}

fn props<L: FpBackend>(label: &[u8]) {
    let mut rng = DetRng::new(label);
    let four = Fp2::<L>::from_small(4);
    for _ in 0..ITER {
        let p = random_point::<L>(&mut rng);
        // H(H(P)) = 4 P
        let hh = hadamard(&hadamard(&p));
        assert!(eq2(&hh.x, &p.x.mul(&four)) && eq2(&hh.t, &p.t.mul(&four)));
        // P * inv(P) is constant coordinate-wise
        let prod = pointwise_product(&p, &invert_point(&p));
        assert!(eq2(&prod.x, &prod.y) && eq2(&prod.x, &prod.z) && eq2(&prod.x, &prod.t));
        // inv(inv(P)) = P projectively
        assert!(proj_eq(&invert_point(&invert_point(&p)), &p));
        // squared theta = H(P^2)
        assert!(proj_eq(
            &to_squared_theta(&p),
            &hadamard(&pointwise_square(&p))
        ));
        // tensor products of two Kummer points are product theta points
        let (a, b, c, d) = (
            rng.random_fp2::<L>(),
            rng.random_fp2::<L>(),
            rng.random_fp2::<L>(),
            rng.random_fp2::<L>(),
        );
        let tensor = ThetaPoint::new(a.mul(&c), a.mul(&d), b.mul(&c), b.mul(&d));
        assert!(bool::from(is_product_theta_point(&tensor)));
        // a random point is not
        assert!(!bool::from(is_product_theta_point(&p)));
    }
}

#[test]
fn p324_3() {
    props::<sqisign_verify::params::P324_3>(b"theta-props/p324_3");
}

#[test]
fn p500_27() {
    props::<sqisign_verify::params::P500_27>(b"theta-props/p500_27");
}

#[test]
fn p664_17() {
    props::<sqisign_verify::params::P664_17>(b"theta-props/p664_17");
}
