//! The parameters of the compact format per prime (`docs/COMPACT_R3.md`):
//! the embedding exponent `e` (the response degree is `q < 2^e`), the
//! torsion `r = ⌈e/2⌉ + 2` the half-chains use and the response scalars
//! are transmitted modulo, the challenge length `λ`, and the non-residue
//! tables of the basis-from-hint convention: the library's seeded tables
//! at level I, the rule of [`super::nqr_rule`] at levels III and V.

use crate::fp::{Fp2, FpBackend};
use crate::params::{P324_3, P500_27, P664_17};
use crate::precomp::PrimePrecomp;

use super::nqr_rule::{NQR_RULE_P500_27, NQR_RULE_P664_17, Z_NQR_RULE_P500_27, Z_NQR_RULE_P664_17};
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
    /// `⌈log2 p⌉`, the bits of one coordinate of `A_com` on the wire.
    const P_BITS: u32 = Self::TWO_ADIC_EXPONENT + Self::COFACTOR_BITLENGTH as u32;
    /// Bytes of the challenge, `λ/8`.
    const CHALLENGE_BYTES: usize;
    /// The twenty quadratic non-residues `x(P)` candidates are drawn from.
    fn nqr_table() -> [Fp2<Self>; 20];
    /// The twenty squares `x` with `x − 1` a non-residue, for `x(Q)`.
    fn z_nqr_table() -> [Fp2<Self>; 20];
    /// The odd cofactor `(p + 1) / 2^f` as limbs, and its bit length.
    fn cofactor() -> ([u64; 1], usize) {
        ([Self::COFACTOR], Self::COFACTOR_BITLENGTH)
    }
}

/// The bits of a signature: `2 ⌈log2 p⌉ + (e − 2) + 3 r + 16`.
pub const fn sig_wire_bits<L: HdLevel>() -> usize {
    (2 * L::P_BITS + (L::E_EMBED - 2) + 3 * L::R + 16) as usize
}

/// The bytes of a signature: its bits padded once to a whole byte.
pub const fn sig_wire_bytes<L: HdLevel>() -> usize {
    sig_wire_bits::<L>().div_ceil(8)
}

/// The public key: `A_pk` and two hint bytes.
pub const fn pk_wire_bytes<L: HdLevel>() -> usize {
    L::FP2_BYTES + 2
}

fn decode_table<L: FpBackend, const B: usize>(t: &'static [[u8; B]; 20]) -> [Fp2<L>; 20] {
    core::array::from_fn(|i| {
        Fp2::<L>::decode(&t[i][..]).expect("invariant: table entries are canonical field elements")
    })
}

/// `j + i` for the small integers `j` of a rule table.
fn small_table<L: FpBackend>(t: &'static [u16; 20]) -> [Fp2<L>; 20] {
    core::array::from_fn(|k| Fp2::<L>::from_small(t[k] as u64).add(&Fp2::<L>::i_element()))
}

impl HdLevel for P324_3 {
    const LAMBDA: u32 = 128;
    const E_EMBED: u32 = 174;
    const R: u32 = 89;
    const E1: u32 = 87;
    const FP2_BYTES: usize = 82;
    const CHALLENGE_BYTES: usize = 16;
    fn nqr_table() -> [Fp2<Self>; 20] {
        decode_table(&NQR_TABLE_P324_3)
    }
    fn z_nqr_table() -> [Fp2<Self>; 20] {
        decode_table(&Z_NQR_TABLE_P324_3)
    }
}

impl HdLevel for P500_27 {
    const LAMBDA: u32 = 192;
    const E_EMBED: u32 = 264;
    const R: u32 = 134;
    const E1: u32 = 132;
    const FP2_BYTES: usize = 128;
    const CHALLENGE_BYTES: usize = 24;
    fn nqr_table() -> [Fp2<Self>; 20] {
        small_table(&NQR_RULE_P500_27)
    }
    fn z_nqr_table() -> [Fp2<Self>; 20] {
        small_table(&Z_NQR_RULE_P500_27)
    }
}

impl HdLevel for P664_17 {
    const LAMBDA: u32 = 256;
    const E_EMBED: u32 = 346;
    const R: u32 = 175;
    const E1: u32 = 173;
    const FP2_BYTES: usize = 168;
    const CHALLENGE_BYTES: usize = 32;
    fn nqr_table() -> [Fp2<Self>; 20] {
        small_table(&NQR_RULE_P664_17)
    }
    fn z_nqr_table() -> [Fp2<Self>; 20] {
        small_table(&Z_NQR_RULE_P664_17)
    }
}

/// The largest encoded `F_p^2` element of the implemented levels.
pub const MAX_FP2_BYTES: usize = <P664_17 as HdLevel>::FP2_BYTES;
/// The largest signature of the implemented levels.
pub const MAX_SIG_WIRE_BYTES: usize = sig_wire_bytes::<P664_17>();
/// The largest public key of the implemented levels.
pub const MAX_PK_WIRE_BYTES: usize = pk_wire_bytes::<P664_17>();

const _: () = assert!(sig_wire_bits::<P324_3>() == 1107 && sig_wire_bytes::<P324_3>() == 139);
const _: () = assert!(sig_wire_bits::<P500_27>() == 1690 && sig_wire_bytes::<P500_27>() == 212);
const _: () = assert!(sig_wire_bits::<P664_17>() == 2223 && sig_wire_bytes::<P664_17>() == 278);
const _: () = assert!(pk_wire_bytes::<P324_3>() == 84);
const _: () = assert!(pk_wire_bytes::<P500_27>() == 130);
const _: () = assert!(pk_wire_bytes::<P664_17>() == 170);
