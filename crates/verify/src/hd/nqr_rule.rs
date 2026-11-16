//! The non-residue tables of the basis-from-hint convention at levels III
//! and V, by a deterministic rule of this crate (the SQIsignHD library's
//! tables are drawn from a seeded generator and exist for its own primes
//! only).
//!
//! The rule. Write `x_j = j + i` in `F_p^2 = F_p[i]/(i^2 + 1)` for the
//! integers `j = 1, 2, 3, ...` in increasing order.
//!
//! * `NQR_RULE[k]` is the `k`-th `j` (from `k = 0`) such that `x_j` is not a
//!   square in `F_p^2`.
//! * `Z_NQR_RULE[k]` is the `k`-th `j` such that `x_j` is a square and
//!   `x_j − 1` is not.
//!
//! Twenty entries each, as at level I. The entry of a hint is `x_j` for the
//! table's `j`; a hint past the table is `hint + i`, as in the library. The
//! constants below are the output of `compute_rule`, written by running the
//! test `rule_tables_are_the_rule` against empty tables and copying the four
//! arrays it prints; the same test recomputes them on every run.

/// Level III, `p = 27 · 2^500 − 1`: the `j` with `j + i` a non-residue.
pub(crate) const NQR_RULE_P500_27: [u16; 20] = [
    4, 5, 8, 9, 10, 12, 13, 14, 16, 17, 18, 23, 24, 27, 28, 29, 31, 32, 33, 35,
];
/// Level III: the `j` with `j + i` a square and `j + i − 1` a non-residue.
pub(crate) const Z_NQR_RULE_P500_27: [u16; 20] = [
    6, 11, 15, 19, 25, 30, 34, 36, 38, 41, 46, 52, 55, 59, 63, 68, 71, 73, 76, 81,
];
/// Level V, `p = 17 · 2^664 − 1`: the `j` with `j + i` a non-residue.
pub(crate) const NQR_RULE_P664_17: [u16; 20] = [
    5, 8, 10, 18, 20, 21, 22, 23, 24, 25, 26, 27, 30, 31, 34, 37, 42, 45, 46, 47,
];
/// Level V: the `j` with `j + i` a square and `j + i − 1` a non-residue.
pub(crate) const Z_NQR_RULE_P664_17: [u16; 20] = [
    6, 9, 11, 19, 28, 32, 35, 38, 43, 48, 54, 63, 68, 70, 72, 74, 77, 81, 85, 89,
];

/// The first twenty `j ≥ 1` of each rule, recomputed.
#[cfg(test)]
pub(crate) fn compute_rule<L: crate::fp::FpBackend>() -> ([u16; 20], [u16; 20]) {
    use crate::fp::Fp2;
    let i = Fp2::<L>::i_element();
    let one = Fp2::<L>::one();
    let (mut nqr, mut z) = ([0u16; 20], [0u16; 20]);
    let (mut n, mut m) = (0usize, 0usize);
    let mut j = 1u16;
    while n < 20 || m < 20 {
        let x = Fp2::<L>::from_small(j as u64).add(&i);
        if !bool::from(x.is_square()) {
            if n < 20 {
                nqr[n] = j;
                n += 1;
            }
        } else if !bool::from(x.sub(&one).is_square()) && m < 20 {
            z[m] = j;
            m += 1;
        }
        j += 1;
    }
    (nqr, z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{P500_27, P664_17};

    #[test]
    fn rule_tables_are_the_rule() {
        let (n3, z3) = compute_rule::<P500_27>();
        let (n5, z5) = compute_rule::<P664_17>();
        assert_eq!(
            (n3, z3, n5, z5),
            (
                NQR_RULE_P500_27,
                Z_NQR_RULE_P500_27,
                NQR_RULE_P664_17,
                Z_NQR_RULE_P664_17
            )
        );
    }
}
