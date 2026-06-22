//! Phase 5b.1 - HD torsion-basis-from-hint recovery.
//!
//! Recovers the canonical `2^f`-torsion basis `(P, Q)` of a Montgomery curve
//! `E_A : y² = x³ + A x² + x` from the integer hints `(hP, hQ)`, in the
//! SQIsignHD convention (`utilities/basis_from_hints.py::
//! torsion_basis_2f_from_hint`).
//!
//! # How HD differs from the dim-2 basis-from-hint
//!
//! Only the **x-coordinate selection** differs. The dim-2
//! `ec_curve_to_basis_2f_from_hint` takes a bit-packed hint and sets
//! `xQ = -(A + xP)` (so `P` and `Q` sit above a shared 2-torsion point). HD
//! instead picks the two x-coordinates **independently** from precomputed
//! tables: `xP = NQR_TABLE[hP]` (a quadratic non-residue) and
//! `xQ = Z_NQR_TABLE[hQ] · α`, where `α = (-A + √(A²-4)) / 2` is a 2-torsion
//! x-coordinate. For a hint `≥ 20` the candidate is `h + i` rather than a
//! table entry (not exercised by the Level-1 oracle vectors, whose hints are
//! all `< 20`).
//!
//! # What is reused unchanged
//!
//! After picking `xP, xQ`, the recovery clears the odd cofactor by an x-only
//! scalar multiplication and then lifts to full points. The cofactor clear
//! ([`ec_mul`]), the difference point ([`difference_point`]) and the
//! Okeya-Sakurai y-recovery ([`lift_basis`]) are the **dim-2 primitives,
//! reused verbatim** - they are identical formulas (both SQIsign2D-West
//! derived) with the same canonical `√` sign convention as HD's
//! `sqrt_Fp2_det`. Normalising `KP, KQ` to affine (`z = 1`) before
//! `difference_point` makes its conjugate-normalisation trivial, so the result
//! matches the HD `difference_point` byte-for-byte.

use crate::ec::basis::{difference_point, is_on_curve, lift_basis};
use crate::ec::point::ec_mul;
use crate::ec::{EcBasis, EcCurve, EcPoint, JacPoint};
use crate::params::SecurityLevel;
use crate::precomp::LevelPrecomp;
use crate::{Fp2, FpBackend, Level1, Level3, Level5};

use crate::hd::nqr_tables_l1::{NQR_TABLE_L1, Z_NQR_TABLE_L1};
use crate::hd::nqr_tables_l3::{NQR_TABLE_L3, Z_NQR_TABLE_L3};
use crate::hd::nqr_tables_l5::{NQR_TABLE_L5, Z_NQR_TABLE_L5};

/// Per-level access to the HD non-residue tables - the only basis input that is
/// genuinely level-specific. The odd cofactor `c = (p+1)/2^f` and the `2^f`
/// power both come from [`LevelPrecomp`] (shared with the dim-2 basis recovery),
/// so they need no HD-specific duplication. Each table entry is a canonical
/// little-endian `Fp2` encoding (`re‖im`, `Fp2EncodedBytes` long).
pub trait HdNqr: SecurityLevel {
    /// 20 quadratic-non-residue x-candidates (`xP = NQR_TABLE[hP]`).
    fn nqr_table() -> [&'static [u8]; 20];
    /// 20 candidates for the `Q`-side (`xQ = Z_NQR_TABLE[hQ]·α`).
    fn z_nqr_table() -> [&'static [u8]; 20];
}
impl HdNqr for Level1 {
    fn nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &NQR_TABLE_L1[i][..])
    }
    fn z_nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &Z_NQR_TABLE_L1[i][..])
    }
}
impl HdNqr for Level3 {
    fn nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &NQR_TABLE_L3[i][..])
    }
    fn z_nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &Z_NQR_TABLE_L3[i][..])
    }
}
impl HdNqr for Level5 {
    fn nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &NQR_TABLE_L5[i][..])
    }
    fn z_nqr_table() -> [&'static [u8]; 20] {
        core::array::from_fn(|i| &Z_NQR_TABLE_L5[i][..])
    }
}

/// Decode a table entry (canonical little-endian `re‖im`) to `Fp2<L>`.
#[inline]
fn fp2_from_le<L: FpBackend>(bytes: &[u8]) -> Fp2<L> {
    Fp2::<L>::decode(bytes).expect("non-residue table entry must be a canonical field element")
}

/// Affine `(x, y)` of a Jacobian point `(X : Y : Z)`: `(X/Z², Y/Z³)`.
#[inline]
pub fn jac_to_affine<L: FpBackend>(j: &JacPoint<L>) -> (Fp2<L>, Fp2<L>) {
    let z2 = j.z.sqr();
    let z3 = z2.mul(&j.z);
    (j.x.mul(&z2.inv()), j.y.mul(&z3.inv()))
}

/// Recover the `2^f`-torsion basis `(P, Q)` of `E_A` from hints `(hp, hq)` in
/// the HD convention, using the supplied non-residue tables and odd cofactor.
///
/// Generic over the level; reuses the dim-2 lifting. Returns `None` if `A` is
/// not a valid Montgomery coefficient. Points are returned in Jacobian
/// coordinates (as from the dim-2 `lift_basis`); use [`jac_to_affine`] to
/// compare against affine references.
pub fn torsion_basis_2f_from_hint<L: FpBackend>(
    curve_a: &Fp2<L>,
    hp: u32,
    hq: u32,
    nqr_table: &[Fp2<L>],
    z_nqr_table: &[Fp2<L>],
    cofactor: &[u64],
    cofactor_bits: usize,
) -> Option<(JacPoint<L>, JacPoint<L>)> {
    let mut curve = EcCurve::from_a(curve_a)?;
    let i = Fp2::<L>::i_element();

    // xP candidate.
    let xp = if (hp as usize) < nqr_table.len() {
        nqr_table[hp as usize].clone()
    } else {
        Fp2::<L>::from_small(hp as u64).add(&i)
    };

    // α = (-A + √(A²-4)) / 2  (a 2-torsion x-coordinate).
    let disc = curve_a.sqr().sub(&Fp2::<L>::from_small(4));
    let alpha = disc.sqrt().sub(curve_a).half();

    // xQ candidate = (table entry or hq+i) · α.
    let xq_base = if (hq as usize) < z_nqr_table.len() {
        z_nqr_table[hq as usize].clone()
    } else {
        Fp2::<L>::from_small(hq as u64).add(&i)
    };
    let xq = xq_base.mul(&alpha);

    // Clear the odd cofactor (x-only), then normalise to affine x (z = 1).
    let mut kp = ec_mul(&EcPoint::new(xp, Fp2::one()), cofactor, cofactor_bits, &mut curve);
    let mut kq = ec_mul(&EcPoint::new(xq, Fp2::one()), cofactor, cofactor_bits, &mut curve);
    kp.normalize();
    kq.normalize();

    let kpq = difference_point(&kp, &kq, &curve);

    let mut basis = EcBasis::new(kp, kq, kpq);
    let (p, q, _on_curve) = lift_basis(&mut basis, &mut curve);
    Some((p, q))
}

/// Recompute the canonical `2^f`-torsion basis hints `(hP, hQ)` for a Level-1
/// Montgomery curve from its coefficient `A` alone.
///
/// The hints are a deterministic function of the curve: each is the smallest
/// table index whose x-candidate lies on the curve, exactly mirroring the C
/// reference's `ec_curve_to_point_2f_{not_,}above_montgomery`
/// (`is_on_curve(NQR_TABLE[hP])` for `P`; `is_on_curve(Z_NQR_TABLE[hQ]·α)` for
/// `Q`, with `α = (-A + √(A²-4))/2` a 2-torsion x-coordinate). Because the
/// signer also picks the smallest valid index, the transmitted hints are
/// redundant - the verifier reconstructs them here, so they need not be sent
/// on the wire.
///
/// Returns `None` if `A` is not a valid Montgomery coefficient, or (only with
/// probability `≈ 2^-20` per curve) if no table entry yields an on-curve
/// candidate; the 20-entry tables make the fallback path of the reference
/// unreachable for realistic inputs (all Level-1 oracle vectors use hints `< 4`).
pub fn canonical_hints<L: FpBackend + HdNqr>(curve_a: &Fp2<L>) -> Option<(u32, u32)> {
    let curve = EcCurve::from_a(curve_a)?;
    let nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::nqr_table()[k]));
    let z_nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::z_nqr_table()[k]));

    // hP: smallest table index whose NQR entry lies on the curve.
    let hp = (0..20).find(|&k| bool::from(is_on_curve(&nqr[k], &curve)))? as u32;

    // hQ: x = Z_NQR_TABLE[hQ] · α, where α is a 2-torsion x-coordinate. The √
    // sign matches the HD `torsion_basis_2f_from_hint` α above, so the
    // reconstructed basis is identical.
    let alpha = curve_a
        .sqr()
        .sub(&Fp2::<L>::from_small(4))
        .sqrt()
        .sub(curve_a)
        .half();
    let hq = (0..20).find(|&k| bool::from(is_on_curve(&z_nqr[k].mul(&alpha), &curve)))? as u32;

    Some((hp, hq))
}

/// Recover the canonical `2^f`-torsion basis from hints at any security level,
/// using the per-level NQR tables ([`HdNqr`]) and the odd cofactor / `2^f` power
/// from [`LevelPrecomp`] (shared with the dim-2 basis recovery).
pub fn hd_torsion_basis<L: FpBackend + LevelPrecomp + HdNqr>(
    curve_a: &Fp2<L>,
    hp: u32,
    hq: u32,
) -> Option<(JacPoint<L>, JacPoint<L>)> {
    let nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::nqr_table()[k]));
    let z_nqr: [Fp2<L>; 20] = core::array::from_fn(|k| fp2_from_le::<L>(L::z_nqr_table()[k]));
    torsion_basis_2f_from_hint(
        curve_a,
        hp,
        hq,
        &nqr,
        &z_nqr,
        L::p_cofactor_for_2f(),
        L::p_cofactor_for_2f_bitlength() as usize,
    )
}

/// Level-1 convenience wrappers (preserve the original API used by the Level-1
/// challenge/response/wire paths). They are exactly the generic functions at
/// `L = Level1`.
pub fn canonical_hints_l1(curve_a: &Fp2<Level1>) -> Option<(u32, u32)> {
    canonical_hints::<Level1>(curve_a)
}
pub fn hd_torsion_basis_l1(
    curve_a: &Fp2<Level1>,
    hp: u32,
    hq: u32,
) -> Option<(JacPoint<Level1>, JacPoint<Level1>)> {
    hd_torsion_basis::<Level1>(curve_a, hp, hq)
}
