//! Phase 9c - structural validation of the level-generic HD torsion-basis
//! recovery at Level 3 and Level 5 (no sage / oracle vectors required).
//!
//! For a real reference public-key curve at each level (from the SQIsignHD
//! reference `Verification/Data/Public_keys_lvl{3,5}.txt`), this checks:
//!
//! 1. the canonical hints recomputed from `A_pk` by [`canonical_hints`] match
//!    the reference public key's transmitted `(hint_pk_P, hint_pk_Q)` - which
//!    validates the transcribed L3/L5 NQR tables and the generic hint logic;
//! 2. [`hd_torsion_basis`] recovers a basis whose two generators lie on the
//!    curve and each have order exactly `2^f` (`f = torsion_even_power`).
//!
//! Together these exercise the generic basis path (NQR tables + the
//! `LevelPrecomp` odd cofactor) end-to-end at L3/L5 without needing the
//! per-step oracle vectors (which require sage to generate). Level 1 is already
//! covered by the oracle tests; the generic functions are the same code at
//! `L = Level1`.

use sqisign_verify::ec::basis::is_on_curve;
use sqisign_verify::ec::jacobian::jac_dbl;
use sqisign_verify::ec::{EcCurve, JacPoint};
use sqisign_verify::hd::{canonical_hints, hd_torsion_basis, jac_to_affine};
use sqisign_verify::precomp::LevelPrecomp;
use sqisign_verify::{Fp2, FpBackend, Level3, Level5};

// Reference public keys (curve coefficient A_pk, canonical little-endian re‖im)
// and their transmitted basis hints, from Verification/Data/Public_keys_lvl*.txt.
const A_PK_L3: [u8; 96] = [
    171, 112, 103, 7, 36, 218, 227, 138, 128, 0, 253, 125, 238, 198, 196, 249, 26, 83, 88, 40, 84,
    230, 23, 109, 50, 82, 73, 17, 137, 245, 123, 76, 210, 59, 119, 80, 155, 134, 82, 182, 200, 86,
    232, 16, 9, 168, 238, 20, 15, 158, 222, 43, 82, 245, 217, 69, 66, 39, 47, 225, 169, 0, 230,
    248, 37, 2, 105, 163, 103, 110, 73, 80, 76, 99, 167, 137, 17, 248, 36, 144, 225, 1, 202, 86,
    197, 161, 105, 126, 168, 50, 140, 102, 216, 80, 126, 4,
];
const HINTS_L3: (u32, u32) = (0, 4);

const A_PK_L5: [u8; 128] = [
    25, 7, 230, 107, 250, 146, 30, 17, 27, 172, 170, 211, 84, 200, 42, 35, 174, 245, 103, 137, 129,
    30, 172, 215, 253, 227, 96, 42, 80, 253, 107, 141, 86, 178, 121, 23, 200, 55, 98, 254, 184,
    142, 174, 149, 31, 68, 184, 83, 206, 21, 203, 211, 126, 209, 125, 2, 74, 50, 94, 78, 88, 27,
    46, 0, 234, 157, 105, 121, 105, 106, 72, 147, 114, 195, 31, 115, 201, 230, 189, 132, 125, 223,
    113, 61, 27, 151, 199, 153, 119, 195, 78, 234, 219, 108, 180, 87, 51, 147, 117, 107, 219, 182,
    240, 249, 208, 150, 113, 209, 208, 148, 139, 54, 98, 47, 28, 68, 184, 53, 6, 41, 41, 5, 204,
    134, 228, 28, 185, 0,
];
const HINTS_L5: (u32, u32) = (0, 0);

fn is_identity<L: FpBackend>(j: &JacPoint<L>) -> bool {
    bool::from(j.z.ct_is_zero())
}

fn check_level<L: FpBackend + LevelPrecomp + sqisign_verify::hd::HdNqr>(
    a_bytes: &[u8],
    ref_hints: (u32, u32),
) {
    let a = Fp2::<L>::decode(a_bytes).expect("A_pk must be a canonical Fp2 element");
    let curve = EcCurve::<L>::from_a(&a).expect("A_pk must be a valid Montgomery coefficient");

    // (1) Recomputed canonical hints match the reference public key's hints.
    let (hp, hq) =
        canonical_hints::<L>(&a).expect("canonical hints must be found in the NQR table");
    assert_eq!(
        (hp, hq),
        ref_hints,
        "recomputed hints must match the reference public key"
    );

    // (2) The basis recovers and both generators lie on the curve.
    let (p, q) = hd_torsion_basis::<L>(&a, hp, hq).expect("basis must recover");
    for pt in [&p, &q] {
        let (x, _y) = jac_to_affine(pt);
        assert!(
            bool::from(is_on_curve(&x, &curve)),
            "recovered basis generator must lie on the curve"
        );
    }

    // (3) Each generator has order exactly 2^f.
    let f = L::torsion_even_power() as usize;
    for pt in [&p, &q] {
        let mut acc = pt.clone();
        for _ in 0..f - 1 {
            acc = jac_dbl(&acc, &curve);
        }
        assert!(!is_identity(&acc), "order must exceed 2^(f-1)");
        acc = jac_dbl(&acc, &curve);
        assert!(
            is_identity(&acc),
            "order must divide 2^f (so it is exactly 2^f)"
        );
    }
}

#[test]
fn level3_basis_recovery() {
    check_level::<Level3>(&A_PK_L3, HINTS_L3);
}

#[test]
fn level5_basis_recovery() {
    check_level::<Level5>(&A_PK_L5, HINTS_L5);
}
