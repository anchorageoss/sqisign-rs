//! The quaternion layer of the SQIsign round-3 reference (`src/quaternion`,
//! spec Section 4.3): the algebra `B_{p,inf} = (-1, -p)`, its elements,
//! integer vectors and matrices, lattices, left `O0`-ideals in the inert
//! Hermite normal form, constant-time lattice reduction, and the norm
//! equations of Qlapoty.
//!
//! Everything is generic over the integer container size `N` and takes the
//! parameter set as a [`QuatAlg`] value, which carries `p`, `isqrt(p)`,
//! the maximal order `O0 = Z<1, i, (i+j)/2, (1+k)/2>` and the per-level
//! integer bounds the reference derives in `lll_config.h`.

// The code below is a limb-by-limb port of C; index loops, explicit range
// checks and long argument lists mirror the reference on purpose.
#![allow(
    clippy::needless_range_loop,
    clippy::manual_range_contains,
    clippy::manual_div_ceil,
    clippy::too_many_arguments,
    clippy::manual_clamp
)]

use crate::mp::Ibz;
use zeroize::Zeroize;

pub mod algebra;
pub mod dim2;
pub mod dim4;
pub mod ideal;
pub mod integers;
pub mod lattice;
pub mod lll;
pub mod params;
pub mod protocol;
pub mod qlapoty;
pub mod stdorder;

/// Vector of two integers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Vec2<const N: usize>(pub [Ibz<N>; 2]);

/// Vector of four integers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Vec4<const N: usize>(pub [Ibz<N>; 4]);

/// `2 x 2` integer matrix, row-major.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mat2x2<const N: usize>(pub [[Ibz<N>; 2]; 2]);

/// `4 x 4` integer matrix, row-major.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mat4x4<const N: usize>(pub [[Ibz<N>; 4]; 4]);

/// An element of the quaternion algebra in the basis `(1, i, j, ij)`:
/// `coord / denom`. Not necessarily normalised.
#[derive(Clone, Copy, Debug)]
pub struct QuatAlgElem<const N: usize> {
    /// Common denominator, non-zero.
    pub denom: Ibz<N>,
    /// Numerators of the four coordinates.
    pub coord: Vec4<N>,
}

/// A full-rank lattice: the columns of `basis` divided by `denom`.
#[derive(Clone, Copy, Debug)]
pub struct QuatLattice<const N: usize> {
    /// Denominator, non-zero.
    pub denom: Ibz<N>,
    /// Integer basis, columns are the generators.
    pub basis: Mat4x4<N>,
}

/// A left `O0`-ideal in inert Hermite normal form, on denominator 2:
///
/// ```text
/// 2n  0    2x        2(n-y)-1
/// 0   2n   2y+1      2x
/// 0   0    1         0
/// 0   0    0         1
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QuatIdeal<const N: usize> {
    /// `x` of the inert form.
    pub x: Ibz<N>,
    /// `y` of the inert form.
    pub y: Ibz<N>,
    /// The norm.
    pub norm: Ibz<N>,
}

impl<const N: usize> Zeroize for Vec2<N> {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<const N: usize> Zeroize for Vec4<N> {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<const N: usize> Zeroize for Mat2x2<N> {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<const N: usize> Zeroize for Mat4x4<N> {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<const N: usize> Zeroize for QuatAlgElem<N> {
    fn zeroize(&mut self) {
        self.denom.zeroize();
        self.coord.zeroize();
    }
}

impl<const N: usize> Zeroize for QuatLattice<N> {
    fn zeroize(&mut self) {
        self.denom.zeroize();
        self.basis.zeroize();
    }
}

impl<const N: usize> Zeroize for QuatIdeal<N> {
    fn zeroize(&mut self) {
        self.x.zeroize();
        self.y.zeroize();
        self.norm.zeroize();
    }
}

/// The per-level constants the lattice reduction is configured with
/// (the reference's `lll_config.h`), all public.
#[derive(Clone, Copy, Debug)]
pub struct LllConfig {
    /// Security parameter.
    pub lambda: i32,
    /// Bits of `p`.
    pub prime_bits: i32,
    /// Bits of `isqrt(p)`.
    pub sqrt_p_bits: i32,
    /// Bound on the double ideal norm `2N + 1`.
    pub norm: i32,
    /// Public bound on the D-profile gap: `2 norm - prime`.
    pub gap_pub: i32,
    /// Bound on the total unimodular transform.
    pub u_bits: i32,
    /// Fixed-point scale.
    pub p_bits: i32,
    /// Headroom for `L`.
    pub head: i32,
    /// Total bit size of `L`, `mu`, `eta`.
    pub l_bits: i32,
    /// Public down-scaling of the D entries.
    pub d_shift: i32,
    /// Bound on the D entries after the down-scaling.
    pub e_bits: i32,
    /// Per-block transform entries.
    pub block_u_bits: i32,
    /// Outer Lagrange-Gauss passes at round 0.
    pub lg_outer_its: i32,
    /// Truncated block Gram width.
    pub gram_tot_bits: i32,
    /// Working width of the Gram updates.
    pub gram_work_bits: i32,
    /// Newton-Raphson rounds of the reciprocal.
    pub nr_rounds: u32,
    /// LLL tours.
    pub tours: i32,
    /// Bound on the inert HNF entries.
    pub hnf_bits: i32,
    /// Bound on the inverse transform.
    pub ainv_bits: i32,
    /// Bound on the reduced basis.
    pub basis_bits: i32,
    /// Bound on the Gram diagonal.
    pub gram_bits: i32,
    /// Shift of the Gram diagonal.
    pub gram_shift: i32,
    /// Bound on the shifted Gram diagonal.
    pub gram_out_bits: i32,
}

const fn limb_align(bits: i32) -> i32 {
    ((bits + 63) / 64) * 64
}

const fn maxb(a: i32, b: i32) -> i32 {
    if a > b {
        a
    } else {
        b
    }
}

const fn log2(x: i32) -> i32 {
    31 - (x as u32).leading_zeros() as i32
}

impl LllConfig {
    /// Derive the configuration from `(lambda, log2 p, log2 isqrt(p))`, as
    /// `lll_config.h` does, including its static assertions.
    pub const fn derive(lambda: i32, prime_bits: i32, sqrt_p_bits: i32, ibz_max_bits: i32) -> Self {
        let norm = prime_bits + lambda;
        let gap_pub = 2 * norm - prime_bits;
        let u_bits = limb_align(gap_pub / 4 + lambda / 2);
        let p_bits = limb_align(2 * u_bits + lambda);
        let head = u_bits;
        let l_bits = limb_align(p_bits + head);
        let d_shift = prime_bits - lambda - lambda / 2;
        let e_bits = limb_align(2 * norm) - d_shift;
        let block_u_bits = limb_align(gap_pub / 4 + lambda / 2);
        let lg_outer_its = block_u_bits / 8;
        let gram_tot_bits = limb_align(maxb(gap_pub + lambda, 0));
        let gram_work_bits = gram_tot_bits + 32;
        let nr_rounds = 5;
        let tours = log2(gap_pub) + 4;
        let hnf_bits = norm + 2;
        let ainv_bits = u_bits + 32;
        let basis_bits = limb_align(hnf_bits + ainv_bits + 2);
        let gram_bits = e_bits + 40;
        let gram_shift = d_shift;
        let gram_out_bits = gram_bits + gram_shift;
        assert!(
            ((61i64 - 2) << nr_rounds) + 2 >= (p_bits as i64) + 4,
            "nr_rounds too small"
        );
        assert!(gram_work_bits <= ibz_max_bits, "bsum must fit a single Ibz");
        assert!(gap_pub > 0, "(2N)^2 must exceed p");
        assert!(l_bits >= p_bits + head, "l_bits under-sized");
        assert!(e_bits + d_shift >= 2 * norm, "D must hold (2N)^2");
        assert!(
            gram_tot_bits >= gap_pub + lambda,
            "gram_tot_bits: precision"
        );
        assert!(
            block_u_bits > gap_pub / 4 + 1,
            "block_u_bits below the worst case"
        );
        assert!(
            ibz_max_bits >= e_bits + block_u_bits,
            "container too small for the block update"
        );
        assert!(
            ibz_max_bits >= l_bits + block_u_bits,
            "container too small for the L-domain multiplies"
        );
        assert!(
            ibz_max_bits >= hnf_bits + ainv_bits,
            "container too small for hnf * Ainv"
        );
        assert!(
            ibz_max_bits >= basis_bits
                && ibz_max_bits >= gram_bits
                && ibz_max_bits >= gram_out_bits
        );
        Self {
            lambda,
            prime_bits,
            sqrt_p_bits,
            norm,
            gap_pub,
            u_bits,
            p_bits,
            head,
            l_bits,
            d_shift,
            e_bits,
            block_u_bits,
            lg_outer_its,
            gram_tot_bits,
            gram_work_bits,
            nr_rounds,
            tours,
            hnf_bits,
            ainv_bits,
            basis_bits,
            gram_bits,
            gram_shift,
            gram_out_bits,
        }
    }

    /// Outer Lagrange-Gauss passes for round `round`: halved every three
    /// rounds from round 1, floor 4. A function of the round index only.
    pub fn lg_outer_its_for_round(&self, round: i32) -> i32 {
        let shift = (if round > 1 { round - 1 } else { 0 }) / 3;
        let its = self.lg_outer_its >> (if shift > 24 { 24 } else { shift });
        if its > 4 {
            its
        } else {
            4
        }
    }
}

/// The quaternion algebra `(-1, -p)` with the per-level constants the
/// quaternion module needs.
#[derive(Clone, Copy, Debug)]
pub struct QuatAlg<const N: usize> {
    /// The prime, `p = 3 mod 4`.
    pub p: Ibz<N>,
    /// `isqrt(p)`.
    pub sqrt_p: Ibz<N>,
    /// Bits of `p`.
    pub p_bits: i32,
    /// Bits of `isqrt(p)`.
    pub sqrt_p_bits: i32,
    /// Security parameter.
    pub lambda: i32,
    /// Miller-Rabin rounds.
    pub primality_num_iter: u32,
    /// Sampling bound of EquivalentCoprimeIdeal.
    pub equiv_bound_coeff: i32,
    /// Power of two used by Qlapoty: `f - 2`.
    pub qlapoty_used_power_of_two: u32,
    /// The prime cofactor of RandomIdealGivenNorm.
    pub prime_cofactor: Ibz<N>,
    /// The maximal order `O0` as a lattice.
    pub o0: QuatLattice<N>,
    /// Lattice reduction configuration.
    pub lll: LllConfig,
}

impl<const N: usize> QuatAlg<N> {
    /// Build the parameter set from little-endian encodings of `p`,
    /// `isqrt(p)` and the prime cofactor, and the level constants.
    pub fn new(
        p_le: &[u8],
        sqrt_p_le: &[u8],
        prime_cofactor_le: &[u8],
        lambda: i32,
        primality_num_iter: u32,
        equiv_bound_coeff: i32,
        qlapoty_used_power_of_two: u32,
    ) -> Self {
        let p = Ibz::<N>::from_le_bytes(p_le);
        let p = p.with_bound(p.bitsize() + 1);
        let sqrt_p = Ibz::<N>::from_le_bytes(sqrt_p_le);
        let sqrt_p = sqrt_p.with_bound(sqrt_p.bitsize() + 1);
        let prime_cofactor = Ibz::<N>::from_le_bytes(prime_cofactor_le);
        let prime_cofactor = prime_cofactor.with_bound(prime_cofactor.bitsize() + 1);
        let p_bits = p.bitsize();
        let sqrt_p_bits = sqrt_p.bitsize();
        Self {
            p,
            sqrt_p,
            p_bits,
            sqrt_p_bits,
            lambda,
            primality_num_iter,
            equiv_bound_coeff,
            qlapoty_used_power_of_two,
            prime_cofactor,
            o0: QuatLattice::o0(),
            lll: LllConfig::derive(lambda, p_bits, sqrt_p_bits, Ibz::<N>::MAX_BITS),
        }
    }
}

impl<const N: usize> Vec2<N> {
    /// `(a0, a1)`.
    pub fn set(a0: i32, a1: i32) -> Self {
        Self([Ibz::set(a0 as i64, 32), Ibz::set(a1 as i64, 32)])
    }

    /// The zero vector.
    pub fn zero() -> Self {
        Self([Ibz::zero(); 2])
    }
}

impl<const N: usize> Vec4<N> {
    /// The zero vector.
    pub fn zero() -> Self {
        Self([Ibz::zero(); 4])
    }

    /// `(c0, c1, c2, c3)`.
    pub fn set(c0: i32, c1: i32, c2: i32, c3: i32) -> Self {
        Self([
            Ibz::set(c0 as i64, 32),
            Ibz::set(c1 as i64, 32),
            Ibz::set(c2 as i64, 32),
            Ibz::set(c3 as i64, 32),
        ])
    }
}

impl<const N: usize> Mat2x2<N> {
    /// The zero matrix.
    pub fn zero() -> Self {
        Self([[Ibz::zero(); 2]; 2])
    }

    /// `[[a00, a01], [a10, a11]]`.
    pub fn set(a00: i32, a01: i32, a10: i32, a11: i32) -> Self {
        Self([
            [Ibz::set(a00 as i64, 32), Ibz::set(a01 as i64, 32)],
            [Ibz::set(a10 as i64, 32), Ibz::set(a11 as i64, 32)],
        ])
    }
}

impl<const N: usize> Mat4x4<N> {
    /// The zero matrix.
    pub fn zero() -> Self {
        Self([[Ibz::zero(); 4]; 4])
    }

    /// The identity.
    pub fn identity() -> Self {
        let mut m = Self::zero();
        for i in 0..4 {
            m.0[i][i] = Ibz::set(1, 2);
        }
        m
    }
}

impl<const N: usize> QuatAlgElem<N> {
    /// Zero with denominator 1.
    pub fn zero() -> Self {
        Self {
            denom: Ibz::set(1, 2),
            coord: Vec4::zero(),
        }
    }

    /// `(c0 + c1 i + c2 j + c3 ij) / denom` from small integers.
    pub fn set(denom: i32, c0: i32, c1: i32, c2: i32, c3: i32) -> Self {
        debug_assert!(denom != 0);
        Self {
            denom: Ibz::set(denom as i64, 32),
            coord: Vec4::set(c0, c1, c2, c3),
        }
    }

    /// From integer coordinates and a denominator.
    pub fn from_ibz(denom: &Ibz<N>, c: [&Ibz<N>; 4]) -> Self {
        Self {
            denom: *denom,
            coord: Vec4([*c[0], *c[1], *c[2], *c[3]]),
        }
    }

    /// The scalar `numerator / denominator`.
    pub fn scalar(numerator: &Ibz<N>, denominator: &Ibz<N>) -> Self {
        let mut e = Self::zero();
        e.denom = *denominator;
        e.coord.0[0] = *numerator;
        e
    }
}

impl<const N: usize> QuatLattice<N> {
    /// The trivial lattice `Z^4` with denominator 1 (an initialised value).
    pub fn init() -> Self {
        Self {
            denom: Ibz::set(1, 2),
            basis: Mat4x4::zero(),
        }
    }

    /// The maximal order `O0 = Z<1, i, (i+j)/2, (1+ij)/2>`.
    pub fn o0() -> Self {
        let mut o = Self::init();
        o.denom = Ibz::set(2, 3);
        o.basis.0[0][0] = Ibz::set(2, 3);
        o.basis.0[1][1] = Ibz::set(2, 3);
        o.basis.0[2][2] = Ibz::set(1, 2);
        o.basis.0[1][2] = Ibz::set(1, 2);
        o.basis.0[3][3] = Ibz::set(1, 2);
        o.basis.0[0][3] = Ibz::set(1, 2);
        o
    }
}

impl<const N: usize> QuatIdeal<N> {
    /// The unit ideal `O0` (`x = 0`, `y = 0`, `norm = 1`).
    pub fn init() -> Self {
        Self {
            x: Ibz::zero(),
            y: Ibz::zero(),
            norm: Ibz::zero(),
        }
    }
}
