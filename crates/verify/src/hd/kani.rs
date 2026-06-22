//! Phase 5b.4 - the Kani-embedding integer linear algebra and the
//! sum-of-two-squares norm equation.
//!
//! This is the conceptual heart of "HD": Kani's lemma embeds the response
//! isogeny into a dimension-4 `(2,2,2,2)`-isogeny, and the embedding is realised
//! by a family of **integer matrices mod `2^k`** (`basis_change/
//! kani_base_change.py`) together with the decomposition `N = 2^f - q =
//! a₁² + a₂²`. None of this has a dimension-2 analogue. The matrices feed the
//! gluing/strategy chain in 5b.6; getting their conventions right here is what
//! lets the derived per-step kernels match the oracle later.
//!
//! # Integer widths (chosen, and why)
//!
//! * **Norm equation.** `N = 2^f - q` is ~136 bits at Level 1 (`f = 136`,
//!   `q < 2^136`), and Cornacchia squares values up to `⌊√N⌋ ≈ 2^68`
//!   (`b² < 2^136`). A fixed 256-bit [`U256`] holds `N` and every square with
//!   room to spare; the modular exponentiation for `√(-1) mod N` uses
//!   `MontyForm` (internal wide products, no overflow). No `alloc`.
//! * **Symplectic matrices over `Z/2^f`** (`f = r = 70` inside `KaniEndoHalf`).
//!   Entries are kept as `u128` in `[0, 2^70)`. Products of two reduced entries
//!   reach ~140 bits and *wrap* `u128` (`> 2^128`), but `2^70 | 2^128`, so
//!   `wrapping_mul` followed by a mask to 70 bits is exactly the product mod
//!   `2^70`. Sums are handled the same way.
//! * **Gluing matrices over `Z/4`.** Entries fit a `u8`; arithmetic is done in
//!   `i64` and reduced into `[0,4)`.
//!
//! # What is oracle-exact here vs. what routes through PARI
//!
//! The **closed-form** matrices - [`matrix_f`], [`matrix_f_dual`], the kernel
//! blocks [`kernel_matrix_f1`]/[`kernel_matrix_f2_dual`] (the columns
//! `kernel_basis` consumes), and the dim-2 gluing matrices
//! [`gluing_dim2_f1`]/[`gluing_dim2_f2`] - are reproduced byte-for-byte against
//! the oracle. The **symplectic completion** ([`complete_symplectic_dim4`]),
//! used to assemble `M1`/`M2` (`starting_two_symplectic_matrices`), routes in
//! sage through `solve_right` over `Z/2^k`, which dispatches to **PARI's
//! `matsolvemod`**: a particular solution chosen by PARI's Hermite normal form.
//! The completion is genuinely non-unique (its freedom is a symmetric matrix:
//! `(A,B) ↦ (A,B)+(C,D)·S`), so matching PARI's *specific* representative
//! byte-for-byte requires reproducing `matsolvemod`'s HNF reduction. We
//! therefore validate the completion by the **symplectic property** (it is a
//! valid completion of the real kernel blocks for all 5 vectors).

use crypto_bigint::modular::{MontyForm, MontyParams};
use crypto_bigint::{Integer, NonZero, Odd, U256};

// Sum of two squares (Cornacchia), N = 2^f - q = a₁² + a₂², N prime ≡ 1 (mod 4).

/// Integer square root `⌊√n⌋` (Newton), for `U256`.
fn isqrt_u256(n: &U256) -> U256 {
    if *n == U256::ZERO {
        return U256::ZERO;
    }
    let nbits = 256 - n.leading_zeros();
    let mut x = U256::ONE.shl(nbits.div_ceil(2));
    loop {
        let (q, _) = n.div_rem(&NonZero::new(x).expect("Newton iterate is positive"));
        let x_new = x.wrapping_add(&q).shr(1);
        if x_new >= x {
            return x;
        }
        x = x_new;
    }
}

/// Cornacchia's algorithm: given `n` prime with `n ≡ 1 (mod 4)`, return
/// `(a1, a2)` with `a1² + a2² = n`, canonicalised so that **`a1` is odd and
/// `a2` is even** (the post-swap convention `KaniEndoHalf` uses). Returns
/// `None` if `n` is not a sum of two squares (e.g. not such a prime) or is even.
pub fn sum_of_two_squares(n: &U256) -> Option<(U256, U256)> {
    // n must be odd to be an odd prime ≡ 1 (mod 4).
    let n_odd: Odd<U256> = Option::from(Odd::new(*n))?;
    let params = MontyParams::new_vartime(n_odd);
    let nm1 = n.wrapping_sub(&U256::ONE);
    let neg_one = MontyForm::new(&nm1, params);
    let exp = nm1.shr(2); // (n-1)/4

    // √(-1) mod n = b^((n-1)/4) for a quadratic non-residue b; scan small b.
    let mut root: Option<U256> = None;
    let mut b: u64 = 2;
    while b < 1000 {
        let t = MontyForm::new(&U256::from_u64(b), params).pow(&exp);
        if t * t == neg_one {
            root = Some(t.retrieve());
            break;
        }
        b += 1;
    }
    let mut x0 = root?;
    // Take the larger root (n/2 < x0 < n), per Cohen's formulation.
    let half = n.shr(1);
    if x0 <= half {
        x0 = n.wrapping_sub(&x0);
    }

    // Euclidean descent: a, b = n, x0; stop at the first b ≤ ⌊√n⌋.
    let limit = isqrt_u256(n);
    let mut a = *n;
    let mut bb = x0;
    while bb > limit {
        let (_, r) = a.div_rem(&NonZero::new(bb).expect("descent divisor is positive"));
        a = bb;
        bb = r;
    }
    // n = bb² + s²; recover s.
    let bb2 = bb.wrapping_mul(&bb);
    let c = n.wrapping_sub(&bb2);
    let s = isqrt_u256(&c);
    if s.wrapping_mul(&s) != c {
        return None;
    }
    // Canonicalise: a1 odd, a2 even.
    let bb_odd = bool::from(bb.is_odd());
    if bb_odd {
        Some((bb, s))
    } else {
        Some((s, bb))
    }
}

/// `N = 2^f - q`, then [`sum_of_two_squares`]. At Level 1 `f = 136`.
pub fn norm_equation_2f_minus_q(f: u32, q: &U256) -> Option<(U256, U256)> {
    let n = U256::ONE.shl(f).wrapping_sub(q);
    sum_of_two_squares(&n)
}

// Symplectic matrices over Z/2^f (f = 70 inside KaniEndoHalf).

/// Modulus exponent of the `KaniEndoHalf` symplectic matrices at Level 1
/// (`= r`). Every `matrix_f`/kernel entry lives in `[0, 2^70)`.
pub const F_MATRIX_L1: u32 = 70;
const MASK70: u128 = (1u128 << 70) - 1;

#[inline]
fn r70(x: u128) -> u128 {
    x & MASK70
}
#[inline]
fn neg70(x: u128) -> u128 {
    0u128.wrapping_sub(x & MASK70) & MASK70
}

/// `matrix_F` (sage `kani_base_change.matrix_F`): the matrix of `F(B1)` in
/// `B1`, an `8×8` matrix over `Z/2^f`. `a1, a2, q` must be reduced mod `2^f`.
pub fn matrix_f(a1: u128, a2: u128, q: u128) -> [[u128; 8]; 8] {
    let (a1, a2, q) = (r70(a1), r70(a2), r70(q));
    let n1 = neg70(1);
    [
        [a1, a2, q, 0, 0, 0, 0, 0],
        [neg70(a2), a1, 0, q, 0, 0, 0, 0],
        [n1, 0, a1, neg70(a2), 0, 0, 0, 0],
        [0, n1, a2, a1, 0, 0, 0, 0],
        [0, 0, 0, 0, a1, a2, 1, 0],
        [0, 0, 0, 0, neg70(a2), a1, 0, 1],
        [0, 0, 0, 0, neg70(q), 0, a1, neg70(a2)],
        [0, 0, 0, 0, 0, neg70(q), a2, a1],
    ]
}

/// `matrix_F_dual` (sage `kani_base_change.matrix_F_dual`).
pub fn matrix_f_dual(a1: u128, a2: u128, q: u128) -> [[u128; 8]; 8] {
    let (a1, a2, q) = (r70(a1), r70(a2), r70(q));
    let n1 = neg70(1);
    [
        [a1, neg70(a2), neg70(q), 0, 0, 0, 0, 0],
        [a2, a1, 0, neg70(q), 0, 0, 0, 0],
        [1, 0, a1, a2, 0, 0, 0, 0],
        [0, 1, neg70(a2), a1, 0, 0, 0, 0],
        [0, 0, 0, 0, a1, neg70(a2), n1, 0],
        [0, 0, 0, 0, a2, a1, 0, n1],
        [0, 0, 0, 0, q, 0, a1, a2],
        [0, 0, 0, 0, 0, q, neg70(a2), a1],
    ]
}

/// Kernel blocks `(C, D)` of `complete_kernel_matrix_F1` - the symplectic
/// generators of `B_Kp1` in `B1`, i.e. exactly the columns `kernel_basis`
/// applies to the torsion points. `4×4` each, over `Z/2^f`.
pub fn kernel_matrix_f1(a1: u128, a2: u128, q: u128) -> ([[u128; 4]; 4], [[u128; 4]; 4]) {
    let (a1, a2, q) = (r70(a1), r70(a2), r70(q));
    let c = [
        [a1, 0, neg70(a2), 0],
        [a2, 0, a1, 0],
        [1, 0, 0, 0],
        [0, 0, 1, 0],
    ];
    let d = [
        [0, a1, 0, neg70(a2)],
        [0, a2, 0, a1],
        [0, q, 0, 0],
        [0, 0, 0, q],
    ];
    (c, d)
}

/// Kernel blocks `(C, D)` of `complete_kernel_matrix_F2_dual`.
pub fn kernel_matrix_f2_dual(a1: u128, a2: u128, q: u128) -> ([[u128; 4]; 4], [[u128; 4]; 4]) {
    let (a1, a2, q) = (r70(a1), r70(a2), r70(q));
    let n1 = neg70(1);
    let c = [
        [a1, 0, a2, 0],
        [neg70(a2), 0, a1, 0],
        [n1, 0, 0, 0],
        [0, 0, n1, 0],
    ];
    let d = [
        [0, a1, 0, a2],
        [0, neg70(a2), 0, a1],
        [0, neg70(q), 0, 0],
        [0, 0, 0, neg70(q)],
    ];
    (c, d)
}

// Gluing matrices over Z/4.

#[inline]
fn m4(x: i64) -> u8 {
    x.rem_euclid(4) as u8
}

/// Inverse of an odd value mod 4.
#[inline]
fn inv4(a: u8) -> u8 {
    // a·a ≡ 1 (mod 4) for odd a (1·1=1, 3·3=9≡1).
    a & 3
}

/// `gluing_base_change_matrix_dim2_F1` - `4×4` symplectic matrix over `Z/4`.
/// `a1, a2, q` are reduced mod 4 internally (`a1` must be odd).
pub fn gluing_dim2_f1(a1: u128, a2: u128, q: u128) -> [[u8; 4]; 4] {
    let (a1, a2, q) = ((a1 % 4) as i64, (a2 % 4) as i64, (q % 4) as i64);
    let mu = inv4(a1 as u8) as i64;
    [
        [0, m4(mu), m4(a1), m4(a2)],
        [0, 0, 1, 0],
        [0, 0, m4(-a2), m4(a1)],
        [m4(-1), m4(-mu * a2), 0, m4(q)],
    ]
}

/// `gluing_base_change_matrix_dim2_F2` - `4×4` symplectic matrix over `Z/4`.
/// `a1, a2, q` are reduced mod 4 internally (`a1` must be odd).
pub fn gluing_dim2_f2(a1: u128, a2: u128, q: u128) -> [[u8; 4]; 4] {
    let (a1, a2, q) = ((a1 % 4) as i64, (a2 % 4) as i64, (q % 4) as i64);
    let mu = inv4(a1 as u8) as i64;
    [
        [0, m4(mu), m4(a1), m4(-a2)],
        [0, 0, m4(-1), 0],
        [0, 0, m4(a2), m4(a1)],
        [1, m4(-mu * a2), 0, m4(-q)],
    ]
}

// Symplectic completion over Z/2^k (general; validated by the symplectic
// property, not byte-matched to PARI's matsolvemod - see module docs).

/// Inverse of an odd `a` modulo `2^k` (`mask = 2^k - 1`). Public wrapper used by
/// the gluing chain (`inverse_mod(q, 2^{m+3})` etc.).
#[inline]
pub fn inverse_mod_pow2(a: u128, mask: u128) -> u128 {
    inv_pow2(a, mask)
}

/// Inverse of an odd `a` modulo `2^k` (`mask = 2^k - 1`), by Hensel doubling.
fn inv_pow2(a: u128, mask: u128) -> u128 {
    let a = a & mask;
    let mut inv = 1u128;
    for _ in 0..7 {
        inv = inv.wrapping_mul(2u128.wrapping_sub(a.wrapping_mul(inv))) & mask;
    }
    inv & mask
}

/// Solve `rows · x = rhs` over `Z/(mask+1)` (a power of two) with the
/// free-variable-zero particular solution, via unit-pivot Gaussian
/// elimination. `nr ≤ 7` rows, 8 columns. Returns `None` if inconsistent.
fn solve_pow2(rows: &[[u128; 8]], rhs: &[u128], nr: usize, mask: u128) -> Option<[u128; 8]> {
    let mut aug = [[0u128; 9]; 7];
    for i in 0..nr {
        for j in 0..8 {
            aug[i][j] = rows[i][j] & mask;
        }
        aug[i][8] = rhs[i] & mask;
    }
    let mut pivots = [(0usize, 0usize); 8];
    let mut np = 0;
    let mut pr = 0;
    for col in 0..8 {
        if pr >= nr {
            break;
        }
        let s = match (pr..nr).find(|&i| aug[i][col] & 1 == 1) {
            Some(s) => s,
            None => continue,
        };
        aug.swap(pr, s);
        let inv = inv_pow2(aug[pr][col], mask);
        for x in aug[pr].iter_mut() {
            *x = x.wrapping_mul(inv) & mask;
        }
        let piv = aug[pr]; // pivot row is fixed while eliminating this column
        for (i, row) in aug.iter_mut().enumerate().take(nr) {
            if i != pr && row[col] & mask != 0 {
                let f = row[col] & mask;
                for (x, &p) in row.iter_mut().zip(piv.iter()) {
                    *x = x.wrapping_sub(f.wrapping_mul(p)) & mask;
                }
            }
        }
        pivots[np] = (pr, col);
        np += 1;
        pr += 1;
    }
    let mut x = [0u128; 8];
    for &(p, c) in pivots.iter().take(np) {
        x[c] = aug[p][8] & mask;
    }
    // Verify (catches inconsistency / a missing unit pivot).
    for i in 0..nr {
        let mut acc = 0u128;
        for j in 0..8 {
            acc = acc.wrapping_add((rows[i][j] & mask).wrapping_mul(x[j])) & mask;
        }
        if acc != rhs[i] & mask {
            return None;
        }
    }
    Some(x)
}

/// Complete the symplectic basis: given the kernel blocks `(c, d)` (the right
/// half of an `8×8` symplectic matrix over `Z/(mask+1)`), return the full
/// matrix `[[A, C], [B, D]]` where `(A, B)` is a symplectic complement.
///
/// Mirrors `base_change_dim4.complete_symplectic_matrix_dim4`: solve
/// `[D^T | -C^T] x = e_i` for each completion column, accumulating the
/// orthogonality constraints `[B_t | -A_t]` against the previously chosen
/// columns. The particular solution is the free-variable-zero one (a *valid*
/// completion; not necessarily PARI's `matsolvemod` representative).
pub fn complete_symplectic_dim4(
    c: &[[u128; 4]; 4],
    d: &[[u128; 4]; 4],
    mask: u128,
) -> Option<[[u128; 8]; 8]> {
    // L_DC = [D^T | -C^T], 4×8.
    let mut l_dc = [[0u128; 8]; 7];
    for i in 0..4 {
        for k in 0..4 {
            l_dc[i][k] = d[k][i] & mask; // D^T
            l_dc[i][4 + k] = (0u128.wrapping_sub(c[k][i])) & mask; // -C^T
        }
    }

    let mut cols = [[0u128; 8]; 4]; // completion columns x_0..x_3
    for i in 0..4 {
        let mut rows = [[0u128; 8]; 7];
        rows[..4].copy_from_slice(&l_dc[..4]);
        // accumulated constraints [B_t_j | -A_t_j] for j < i
        for j in 0..i {
            for k in 0..4 {
                rows[4 + j][k] = cols[j][4 + k] & mask; // B_t_j
                rows[4 + j][4 + k] = (0u128.wrapping_sub(cols[j][k])) & mask; // -A_t_j
            }
        }
        let nr = 4 + i;
        let mut rhs = [0u128; 7];
        rhs[i] = 1;
        cols[i] = solve_pow2(&rows, &rhs, nr, mask)?;
    }

    // M = [[A, C], [B, D]]: columns 0..4 = completion x_i, columns 4..8 = (C; D).
    let mut m = [[0u128; 8]; 8];
    for i in 0..4 {
        for r in 0..8 {
            m[r][i] = cols[i][r] & mask;
        }
    }
    for j in 0..4 {
        for r in 0..4 {
            m[r][4 + j] = c[r][j] & mask;
            m[4 + r][4 + j] = d[r][j] & mask;
        }
    }
    Some(m)
}

/// `4×4` matrix product over `Z/(mask+1)`.
fn mat4_mul(x: &[[u128; 4]; 4], y: &[[u128; 4]; 4], mask: u128) -> [[u128; 4]; 4] {
    let mut z = [[0u128; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut acc = 0u128;
            for k in 0..4 {
                acc = acc.wrapping_add(x[i][k].wrapping_mul(y[k][j])) & mask;
            }
            z[i][j] = acc;
        }
    }
    z
}

#[inline]
fn mat4_t(x: &[[u128; 4]; 4]) -> [[u128; 4]; 4] {
    let mut t = [[0u128; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            t[i][j] = x[j][i];
        }
    }
    t
}

fn block(m: &[[u128; 8]; 8], r0: usize, c0: usize) -> [[u128; 4]; 4] {
    let mut b = [[0u128; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            b[i][j] = m[r0 + i][c0 + j];
        }
    }
    b
}

/// Test whether `m` is symplectic over `Z/(mask+1)`: with the block
/// decomposition `[[A,C],[B,D]]`, `BᵀA = AᵀB`, `CᵀD = DᵀC`, `AᵀD - BᵀC = I`.
pub fn is_symplectic_dim4(m: &[[u128; 8]; 8], mask: u128) -> bool {
    let a = block(m, 0, 0);
    let c = block(m, 0, 4);
    let b = block(m, 4, 0);
    let d = block(m, 4, 4);
    let eq = |x: &[[u128; 4]; 4], y: &[[u128; 4]; 4]| -> bool {
        (0..4).all(|i| (0..4).all(|j| x[i][j] & mask == y[i][j] & mask))
    };
    let bta = mat4_mul(&mat4_t(&b), &a, mask);
    let atb = mat4_mul(&mat4_t(&a), &b, mask);
    if !eq(&bta, &atb) {
        return false;
    }
    let ctd = mat4_mul(&mat4_t(&c), &d, mask);
    let dtc = mat4_mul(&mat4_t(&d), &c, mask);
    if !eq(&ctd, &dtc) {
        return false;
    }
    let atd = mat4_mul(&mat4_t(&a), &d, mask);
    let btc = mat4_mul(&mat4_t(&b), &c, mask);
    let mut id = [[0u128; 4]; 4];
    for (i, row) in id.iter_mut().enumerate() {
        row[i] = 1;
    }
    (0..4).all(|i| (0..4).all(|j| atd[i][j].wrapping_sub(btc[i][j]) & mask == id[i][j]))
}

// Self-derivation of the starting symplectic matrices M1, M2 and the dim-4
// gluing change of basis M_gluing_1/2 (`starting_two_symplectic_matrices` +
// `gluing_base_change_matrix_dim2_dim4_F1/F2`).
//
// M1 and M2 share a *single* symplectic completion (`complete_kernel_matrix_F1`'s)
// propagated through the closed-form `matrix_F`/`matrix_F_dual` - so any valid
// completion keeps F1 and F2 dual halves of the same F. The completion is
// self-consistent (5b.4); the gluing output is therefore completion-dependent
// but yields the invariant middle-codomain match.

/// `(M·N[:, 0..4])` for two `8×8` matrices over `Z/(mask+1)`: an `8×4` result.
fn mat8_mul_cols0to3(m: &[[u128; 8]; 8], n: &[[u128; 8]; 8], mask: u128) -> [[u128; 4]; 8] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut acc = 0u128;
            for k in 0..8 {
                acc = acc.wrapping_add(m[i][k].wrapping_mul(n[k][j])) & mask;
            }
            acc
        })
    })
}

/// `starting_two_symplectic_matrices(a1, a2, q, 2^f)`: the two starting
/// symplectic matrices `(M1, M2)` over `Z/2^f` (`mask = 2^f - 1`). Returns
/// `None` if a completion is not solvable.
#[allow(clippy::type_complexity)] // (M1, M2), two 8×8 matrices.
pub fn starting_two_symplectic_matrices(
    a1: u128,
    a2: u128,
    q: u128,
    mask: u128,
) -> Option<([[u128; 8]; 8], [[u128; 8]; 8])> {
    // M1_0 = complete_kernel_matrix_F1 (the single shared completion).
    let (c1, d1) = kernel_matrix_f1(a1, a2, q);
    let m1_0 = complete_symplectic_dim4(&c1, &d1, mask)?;

    // M2 kernel block = matrix_F · M1_0[:, 0..4]; complete it.
    let matf = matrix_f(a1, a2, q);
    let br2 = mat8_mul_cols0to3(&matf, &m1_0, mask);
    let c_br2: [[u128; 4]; 4] = core::array::from_fn(|i| br2[i]);
    let d_br2: [[u128; 4]; 4] = core::array::from_fn(|i| br2[4 + i]);
    let m2 = complete_symplectic_dim4(&c_br2, &d_br2, mask)?;

    // M1 = [ M1_0[:, 0..4] | -(matrix_F_dual · M2[:, 0..4]) ].
    let matf_dual = matrix_f_dual(a1, a2, q);
    let br1 = mat8_mul_cols0to3(&matf_dual, &m2, mask);
    let mut m1 = [[0u128; 8]; 8];
    for r in 0..8 {
        for j in 0..4 {
            m1[r][j] = m1_0[r][j] & mask;
            m1[r][4 + j] = (0u128.wrapping_sub(br1[r][j])) & mask;
        }
    }
    Some((m1, m2))
}

/// The four `4×4` blocks of an `8×8` matrix: `(A, B, C, D)` with
/// `A = M[0..4][0..4]`, `B = M[4..8][0..4]`, `C = M[0..4][4..8]`, `D = M[4..8][4..8]`.
#[allow(clippy::type_complexity)] // four 4×4 blocks (A, B, C, D).
fn blocks8(m: &[[u128; 8]; 8]) -> ([[u128; 4]; 4], [[u128; 4]; 4], [[u128; 4]; 4], [[u128; 4]; 4]) {
    (block(m, 0, 0), block(m, 4, 0), block(m, 0, 4), block(m, 4, 4))
}

/// `gluing_base_change_matrix_dim2_dim4_F1(a1, a2, q, m, M1)` - the dim-4 gluing
/// change of basis (`8×8` over `Z/4`, returned as `i64`). `f1 = true` selects the
/// F1 sign convention, `false` the F2_dual one.
fn gluing_bc_dim4(
    a1: u128,
    a2: u128,
    q: u128,
    m: usize,
    m1: &[[u128; 8]; 8],
    f1: bool,
) -> [[i64; 8]; 8] {
    let modp2 = 1i128 << (m + 2); // 2^(m+2)
    let maskp2 = (1u128 << (m + 2)) - 1;
    let two_mp1 = 1i128 << (m + 1); // 2^(m+1)
    let red = |x: i128| -> i128 { x.rem_euclid(modp2) };

    let a1i = (a1 & maskp2) as i128;
    let a2i = (a2 & maskp2) as i128;
    let qi = (q & maskp2) as i128;
    let inv_a1 = inv_pow2(a1, maskp2) as i128;
    let inv_q = inv_pow2(q, maskp2) as i128;

    let lamb = two_mp1;
    let cc = two_mp1;
    let (mu, bq) = if f1 {
        (red((1 - two_mp1 * qi) * inv_a1), red(-1 - two_mp1 * a1i))
    } else {
        (red((1 + two_mp1 * qi) * inv_a1), red(1 + two_mp1 * a1i))
    };
    let a = red(two_mp1 * a2i * inv_q);
    let dq = red(-mu * a2i);

    // Block entries of M1, reduced mod 2^(m+2) (signed-friendly).
    let (ba, bb, bc, bd) = blocks8(m1);
    let g = |blk: &[[u128; 4]; 4], i: usize, j: usize| -> i128 { (blk[i][j] & maskp2) as i128 };

    let mut out = [[0i64; 8]; 8];
    let m_shift = m as u32;
    for j in 0..4 {
        let (a0j, a1j, a2j, a3j) = (g(&ba, 0, j), g(&ba, 1, j), g(&ba, 2, j), g(&ba, 3, j));
        let (b0j, b1j, b2j, b3j) = (g(&bb, 0, j), g(&bb, 1, j), g(&bb, 2, j), g(&bb, 3, j));
        let (c0j, c1j, c2j, c3j) = (g(&bc, 0, j), g(&bc, 1, j), g(&bc, 2, j), g(&bc, 3, j));
        let (d0j, d1j, d2j, d3j) = (g(&bd, 0, j), g(&bd, 1, j), g(&bd, 2, j), g(&bd, 3, j));

        // Ap (mod 4).
        let ap = if f1 {
            [
                -b0j * a1i - a0j * a2i - b2j,
                a0j * a1i - b0j * a2i + a2j * qi,
                -b1j * a1i - a1j * a2i - b3j,
                a1j * a1i - b1j * a2i + a3j * qi,
            ]
        } else {
            [
                -b0j * a1i + a0j * a2i + b2j,
                a0j * a1i + b0j * a2i - a2j * qi,
                -b1j * a1i + a1j * a2i + b3j,
                a1j * a1i + b1j * a2i - a3j * qi,
            ]
        };
        // Bp (2^m · X mod 4) - identical formulas for F1 and F2.
        let bp_x = [
            b2j * a - a0j * lamb - a2j * bq,
            b2j * cc + b0j * mu - a2j * dq,
            b3j * a - a1j * lamb - a3j * bq,
            b3j * cc + b1j * mu - a3j * dq,
        ];
        // Cp ((num) // 2^m mod 4).
        let cp_num = if f1 {
            [
                -d0j * a1i - c0j * a2i - d2j,
                c0j * a1i - d0j * a2i + c2j * qi,
                -d1j * a1i - c1j * a2i - d3j,
                c1j * a1i - d1j * a2i + c3j * qi,
            ]
        } else {
            [
                -d0j * a1i + c0j * a2i + d2j,
                c0j * a1i + d0j * a2i - c2j * qi,
                -d1j * a1i + c1j * a2i + d3j,
                c1j * a1i + d1j * a2i - c3j * qi,
            ]
        };
        // Dp (mod 4) - identical formulas for F1 and F2.
        let dp = [
            d2j * a - c0j * lamb - c2j * bq,
            d0j * mu + d2j * cc - c2j * dq,
            d3j * a - c1j * lamb - c3j * bq,
            d1j * mu + d3j * cc - c3j * dq,
        ];

        for r in 0..4 {
            // Ap → rows 0..4, cols 0..4
            out[r][j] = ap[r].rem_euclid(4) as i64;
            // Cp → rows 0..4, cols 4..8
            out[r][4 + j] = ((cp_num[r].rem_euclid(modp2) >> m_shift) & 3) as i64;
            // Bp → rows 4..8, cols 0..4
            out[4 + r][j] = ((1i128 << m_shift) * bp_x[r]).rem_euclid(4) as i64;
            // Dp → rows 4..8, cols 4..8
            out[4 + r][4 + j] = dp[r].rem_euclid(4) as i64;
        }
    }
    out
}

/// `gluing_base_change_matrix_dim2_dim4_F1`.
pub fn gluing_bc_dim4_f1(a1: u128, a2: u128, q: u128, m: usize, m1: &[[u128; 8]; 8]) -> [[i64; 8]; 8] {
    gluing_bc_dim4(a1, a2, q, m, m1, true)
}

/// `gluing_base_change_matrix_dim2_dim4_F2`.
pub fn gluing_bc_dim4_f2(a1: u128, a2: u128, q: u128, m: usize, m2: &[[u128; 8]; 8]) -> [[i64; 8]; 8] {
    gluing_bc_dim4(a1, a2, q, m, m2, false)
}
