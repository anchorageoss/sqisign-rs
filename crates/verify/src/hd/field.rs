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
