//! Counters of the `Fp2` operations (the cargo feature `op-counters`,
//! SECENG-2434): one relaxed atomic per kind, read by the profiler at its
//! phase boundaries. Without the feature every call here is empty and
//! compiles away; the shipped builds and the benchmark sessions do not
//! have it.

#[cfg(feature = "op-counters")]
use core::sync::atomic::{AtomicU64, Ordering};

/// Multiplications (`mul`, `mul_assign`).
pub const OP_MUL: usize = 0;
/// Squarings (`sqr`, `sqr_assign`).
pub const OP_SQR: usize = 1;
/// Additions: `add`, `sub`, `neg`, `add_one`, `conjugate`, `mul_by_i` and
/// the in-place forms, each as one.
pub const OP_ADD: usize = 2;
/// Inversions (`inv`; a batched inversion counts its one inversion and its
/// multiplications).
pub const OP_INV: usize = 3;
/// Square roots (`sqrt`).
pub const OP_SQRT: usize = 4;
/// Legendre symbols (`is_square`).
pub const OP_LEGENDRE: usize = 5;
/// The small scalings `mul_small` and `half`.
pub const OP_SMALL: usize = 6;
/// The number of kinds.
pub const OPS: usize = 7;
/// The names of the kinds, in order.
pub const OP_NAMES: [&str; OPS] = ["mul", "sqr", "add", "inv", "sqrt", "legendre", "small"];
/// Whether the counters are compiled in.
pub const ENABLED: bool = cfg!(feature = "op-counters");

#[cfg(feature = "op-counters")]
static COUNTS: [AtomicU64; OPS] = [const { AtomicU64::new(0) }; OPS];

/// Count one operation of kind `k`.
#[cfg(feature = "op-counters")]
#[inline(always)]
pub fn count(k: usize) {
    COUNTS[k].fetch_add(1, Ordering::Relaxed);
}

/// Count one operation of kind `k` (nothing: the feature is off).
#[cfg(not(feature = "op-counters"))]
#[inline(always)]
pub fn count(_k: usize) {}

/// The counts so far, in the order of [`OP_NAMES`].
#[cfg(feature = "op-counters")]
pub fn snapshot() -> [u64; OPS] {
    core::array::from_fn(|i| COUNTS[i].load(Ordering::Relaxed))
}

/// The counts so far (all zero: the feature is off).
#[cfg(not(feature = "op-counters"))]
pub fn snapshot() -> [u64; OPS] {
    [0; OPS]
}
