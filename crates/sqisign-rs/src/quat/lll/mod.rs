//! Constant-time lattice reduction (spec Section 4.3.3): Lagrange-Gauss
//! over `Z` and `Z[i]` with Lehmer-style windows, the dimension-4 LLL with
//! a Minkowski finish on a dual LDL profile, and their applications to
//! `O0`-ideals. Inputs of honest size are reduced with failure probability
//! below `2^-lambda`; the loop counts are fixed functions of the public
//! bounds in [`super::LllConfig`].

pub mod applications;
pub mod dim4;
pub mod gaussian;
pub mod lehmer;
pub mod lg2;

/// Dimension-2 window parameters (the reference's `DIM2_*`).
pub(crate) const DIM2_THRESHOLD: i32 = 30;
pub(crate) const DIM2_W: i32 = 2 * DIM2_THRESHOLD + 3;

pub(crate) const fn dim2_outer_its(in_bits: i32) -> i32 {
    in_bits / (DIM2_THRESHOLD - 1) + 1
}

pub(crate) const fn dim2_val_bits(in_bits: i32) -> i32 {
    in_bits + 1
}

pub(crate) const fn dim2_work_bits(in_bits: i32) -> i32 {
    dim2_val_bits(in_bits) + 32
}

pub(crate) const fn dim2_shift_max(in_bits: i32) -> i32 {
    let v = dim2_val_bits(in_bits) - DIM2_THRESHOLD;
    if v > 0 {
        v
    } else {
        0
    }
}

/// Dimension-2 over `Z[i]` parameters (the reference's `DIM2I_*`).
pub(crate) const DIM2I_THRESHOLD: i32 = 27;
pub(crate) const DIM2I_W: i32 = 2 * DIM2I_THRESHOLD + 5;
pub(crate) const DIM2I_UP_SHIFT: i32 = 3;
pub(crate) const DIM2I_INNER_ITS: i32 = DIM2I_THRESHOLD + 1;

pub(crate) fn dim2i_req_drop(in_bits: i32, sqrt_p_bits: i32) -> i32 {
    ((in_bits - sqrt_p_bits).max(0) + 1) / 2 + 1
}

pub(crate) fn dim2i_outer_its(in_bits: i32, sqrt_p_bits: i32) -> i32 {
    (dim2i_req_drop(in_bits, sqrt_p_bits) + DIM2I_THRESHOLD - 2) / (DIM2I_THRESHOLD - 1) + 4
}

pub(crate) fn dim2i_shift_max(in_bits: i32) -> i32 {
    (dim2_val_bits(in_bits) - DIM2I_THRESHOLD).max(0)
}

pub(crate) fn dim2i_herm_bits(in_bits: i32, p_bits: i32) -> i32 {
    2 * dim2_val_bits(in_bits) + p_bits + 3
}
