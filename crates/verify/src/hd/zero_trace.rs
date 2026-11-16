//! Counters of all-zero vectors and points met during verification.
//!
//! A theta point, a theta null point and a Kummer point are projective: the
//! vector with every coordinate zero is none of them. This module counts,
//! site by site, how often the verifier produces or compares such a vector.
//! With the cargo feature `instrument` each site increments two atomics (the
//! visits and the all-zero visits); without it the calls are empty and the
//! verifier is unchanged. Measurement only.

use crate::fp::{Fp2, FpBackend};

/// The sites, in the order of the counters.
pub const SITES: [&str; 28] = [
    "dim2: Jacobian point (0 : y : 0) entering the gluing",
    "dim2: base-changed product point",
    "dim2: gluing codomain null",
    "dim2: gluing image",
    "dim2: doubled point",
    "dim2: plain codomain null",
    "dim2: plain image",
    "dim2: dual plain image",
    "dim2: dual gluing, product coordinates",
    "dim2: dual gluing, not a product point",
    "dim4: point entering the gluing",
    "dim4: gluing codomain null",
    "dim4: gluing image",
    "dim4: doubled point",
    "dim4: plain codomain null",
    "dim4: plain image",
    "dim4: dual plain image",
    "dim4: splitting, before the base change",
    "dim4: splitting, product coordinates",
    "dim4: splitting, a dim-2 factor",
    "compare: projective_eq, left",
    "compare: projective_eq, right",
    "compare: x_eq, left with z = 0",
    "compare: x_eq, left (0:0)",
    "compare: x_eq, right with z = 0",
    "compare: x_eq, right (0:0)",
    "image: fourth component not the identity",
    "image: fourth component (0:0)",
];

pub const D2_KUMMER_00: usize = 0;
pub const D2_BASE_CHANGE: usize = 1;
pub const D2_GLUING_NULL: usize = 2;
pub const D2_GLUING_IMAGE: usize = 3;
pub const D2_DOUBLE: usize = 4;
pub const D2_PLAIN_NULL: usize = 5;
pub const D2_PLAIN_IMAGE: usize = 6;
pub const D2_DUAL_PLAIN: usize = 7;
pub const D2_DUAL_PRODUCT: usize = 8;
pub const D2_DUAL_NOT_PRODUCT: usize = 9;
pub const D4_SPLIT_EVAL: usize = 10;
pub const D4_GLUING_NULL: usize = 11;
pub const D4_GLUING_IMAGE: usize = 12;
pub const D4_DOUBLE: usize = 13;
pub const D4_PLAIN_NULL: usize = 14;
pub const D4_PLAIN_IMAGE: usize = 15;
pub const D4_DUAL_PLAIN: usize = 16;
pub const D4_SPLIT_PRE: usize = 17;
pub const D4_SPLIT_PRODUCT: usize = 18;
pub const D4_SPLIT_FACTOR: usize = 19;
pub const EQ_LEFT: usize = 20;
pub const EQ_RIGHT: usize = 21;
pub const XEQ_LEFT_Z0: usize = 22;
pub const XEQ_LEFT_00: usize = 23;
pub const XEQ_RIGHT_Z0: usize = 24;
pub const XEQ_RIGHT_00: usize = 25;
pub const IMAGE_C3_NOT_ID: usize = 26;
pub const IMAGE_C3_00: usize = 27;

#[cfg(feature = "instrument")]
mod counters {
    use core::sync::atomic::{AtomicU64, Ordering};

    const N: usize = super::SITES.len();
    #[allow(clippy::declare_interior_mutable_const)]
    const Z: AtomicU64 = AtomicU64::new(0);
    static VISITS: [AtomicU64; N] = [Z; N];
    static HITS: [AtomicU64; N] = [Z; N];

    #[inline]
    pub fn record(site: usize, hit: bool) {
        VISITS[site].fetch_add(1, Ordering::Relaxed);
        if hit {
            HITS[site].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// `(visits, hits)` of every site since the last [`reset`].
    pub fn snapshot() -> [(u64, u64); N] {
        core::array::from_fn(|i| {
            (
                VISITS[i].load(Ordering::Relaxed),
                HITS[i].load(Ordering::Relaxed),
            )
        })
    }

    /// Set every counter to zero.
    pub fn reset() {
        for i in 0..N {
            VISITS[i].store(0, Ordering::Relaxed);
            HITS[i].store(0, Ordering::Relaxed);
        }
    }
}

#[cfg(feature = "instrument")]
pub use counters::{reset, snapshot};

/// Record a visit of `site`, and a hit when `hit` holds.
#[inline(always)]
pub fn flag(site: usize, hit: bool) {
    #[cfg(feature = "instrument")]
    counters::record(site, hit);
    #[cfg(not(feature = "instrument"))]
    let _ = (site, hit);
}

/// Record a visit of `site`, and a hit when every entry of `v` is zero.
#[inline(always)]
pub fn vector<L: FpBackend>(site: usize, v: &[Fp2<L>]) {
    #[cfg(feature = "instrument")]
    counters::record(site, v.iter().all(|c| bool::from(c.ct_is_zero())));
    #[cfg(not(feature = "instrument"))]
    let _ = (site, v);
}
