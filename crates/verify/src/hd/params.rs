//! The parameters of the compact format per prime (`docs/COMPACT_R3.md`):
//! the embedding exponent `e` (the response degree is `q < 2^e`), the
//! torsion `r = ⌈e/2⌉ + 2` the half-chains use and the response scalars
//! are transmitted modulo, the challenge length `λ`, and the library's
//! non-residue tables for the basis-from-hint convention. Level I only.

use crate::fp::FpBackend;
use crate::params::P324_3;
use crate::precomp::PrimePrecomp;

use super::nqr_tables_p324_3::{NQR_TABLE_P324_3, Z_NQR_TABLE_P324_3};

/// A prime with compact-format parameters.
pub trait HdLevel: FpBackend + PrimePrecomp {
    /// `λ`: the challenge isogeny has degree `2^LAMBDA`.
    const LAMBDA: u32;
    /// `e`: the response degree is below `2^E_EMBED`.
    const E_EMBED: u32;
    /// `r = ⌈e/2⌉ + 2`: the torsion of the half-chains, the scalars' modulus.
    const R: u32;
    /// `⌈e/2⌉`, the length of the first half-chain.
    const E1: u32;
    /// Bytes of one encoded `F_p^2` element (`2 FP`).
    const FP2_BYTES: usize;
    /// Bytes of the response degree on the wire, `⌈e/8⌉`.
    const Q_BYTES: usize;
    /// Bytes of one response scalar on the wire, `⌈r/8⌉`.
    const SCALAR_BYTES: usize;
    /// Bytes of the challenge, `λ/8`.
    const CHALLENGE_BYTES: usize;
    /// The twenty quadratic non-residues, canonical encodings.
    fn nqr_table() -> [&'static [u8]; 20];
    /// The twenty squares `x` with `x − 1` a non-residue, canonical encodings.
    fn z_nqr_table() -> [&'static [u8]; 20];
    /// The odd cofactor `(p + 1) / 2^f` as limbs, and its bit length.
    fn cofactor() -> ([u64; 1], usize) {
        ([Self::COFACTOR], Self::COFACTOR_BITLENGTH)
    }
}

/// The wire sizes of a level.
pub const fn sig_wire_bytes<L: HdLevel>() -> usize {
    L::FP2_BYTES + L::Q_BYTES + 3 * L::SCALAR_BYTES + 2
}

/// The public key: `A_pk` and two hint bytes.
pub const fn pk_wire_bytes<L: HdLevel>() -> usize {
    L::FP2_BYTES + 2
}

impl HdLevel for P324_3 {
    const LAMBDA: u32 = 128;
    const E_EMBED: u32 = 174;
    const R: u32 = 89;
    const E1: u32 = 87;
    const FP2_BYTES: usize = 82;
    const Q_BYTES: usize = 22;
    const SCALAR_BYTES: usize = 12;
    const CHALLENGE_BYTES: usize = 16;
    fn nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &NQR_TABLE_P324_3[i][..])
    }
    fn z_nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &Z_NQR_TABLE_P324_3[i][..])
    }
}

/// The largest signature of the implemented levels.
pub const MAX_SIG_WIRE_BYTES: usize = sig_wire_bytes::<P324_3>();
/// The largest public key of the implemented levels.
pub const MAX_PK_WIRE_BYTES: usize = pk_wire_bytes::<P324_3>();
