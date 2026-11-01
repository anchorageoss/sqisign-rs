//! Field helpers of the dimension-4 layer: inversions, over the crate's
//! constant-time `Fp2` (the round-2 layer carried a variable-time level-1
//! path here; the round-3 kernels do not need one).

use crate::fp::{Fp2, FpBackend};

/// Invert every element of `x` in place with two scratch slices of the same
/// length (Montgomery's trick, one inversion).
#[inline]
pub fn batched_inv<L: FpBackend>(x: &mut [Fp2<L>], t1: &mut [Fp2<L>], t2: &mut [Fp2<L>]) {
    Fp2::batched_inv(x, t1, t2)
}

/// `1 / x` (`0` for `x = 0`).
#[inline]
pub fn inv<L: FpBackend>(x: &Fp2<L>) -> Fp2<L> {
    x.inv()
}

/// `out[i] = prod_{j != i} x[j]`: `1 / x` up to the common factor `prod_j x[j]`, with no
/// inversion (`3n - 1` products by prefix and suffix products).
#[inline]
pub fn products_of_others<L: FpBackend, const N: usize>(x: &[Fp2<L>; N]) -> [Fp2<L>; N] {
    let mut pre: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::one());
    for i in 1..N {
        pre[i] = pre[i - 1].mul(&x[i - 1]);
    }
    let mut out: [Fp2<L>; N] = core::array::from_fn(|_| Fp2::zero());
    let mut suf = Fp2::<L>::one();
    for i in (0..N).rev() {
        out[i] = pre[i].mul(&suf);
        suf.mul_assign(&x[i]);
    }
    out
}
