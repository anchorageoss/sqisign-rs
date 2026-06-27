//!
//! The [`SigningPrecomp`] struct aggregates all level-specific precomputed
//! data that the signing functions need. It is constructed once per level
//! and passed by reference to every signing function, keeping those
//! functions free of direct precomp imports.

extern crate alloc;

use crate::quaternion::intbig::Ibz;
use crate::quaternion::types::{
    IbzMat2x2, IbzMat4x4, IbzVec4, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal,
    QuatPExtremalMaximalOrder,
};
use alloc::vec::Vec;
use num_traits::{One, Zero};
use sqisign_verify::ec::{EcBasis, EcCurve, EcPoint};
use sqisign_verify::fp::{Fp2, FpBackend};

/// Precomputed 2×2 action matrices for one endomorphism curve.
pub struct ActionMatrices {
    pub i: IbzMat2x2,
    pub j: IbzMat2x2,
    pub gen2: IbzMat2x2,
    pub gen3: IbzMat2x2,
    pub gen4: IbzMat2x2,
}

impl ActionMatrices {
    /// Look up an action matrix by generator name.
    pub fn by_name(&self, gen: &str) -> &IbzMat2x2 {
        match gen {
            "I" => &self.i,
            "J" => &self.j,
            "GEN2" => &self.gen2,
            "GEN3" => &self.gen3,
            "GEN4" => &self.gen4,
            _ => unreachable!("unknown generator"),
        }
    }
}

/// Bundle of precomputed data needed by signing-side id2iso functions.
///
/// Constructed once per security level and passed by reference. This
/// keeps the signing functions free of direct `sqisign_verify::precomp::levelN`
/// imports.
pub struct SigningPrecomp<L: FpBackend> {
    pub extremal_orders: Vec<QuatPExtremalMaximalOrder>,
    pub connecting_ideals: Vec<QuatLeftIdeal>,
    pub endomorphism_curves: Vec<EcCurve<L>>,
    pub endomorphism_bases: Vec<EcBasis<L>>,
    pub action_matrices: Vec<ActionMatrices>,
    pub algebra: QuatAlg,
    /// `2^TORSION_EVEN_POWER` as a big integer.
    pub torsion_plus_2power: Ibz,
    /// Cofactor `(p+1)/2^TORSION_EVEN_POWER` as a digit array.
    pub p_cofactor_for_2f: &'static [u64],
    /// Bit length of the cofactor.
    pub p_cofactor_for_2f_bitlength: usize,
    pub finduv_box_size: u32,
    pub finduv_cube_size: u32,
    pub quat_repres_bound_input: u32,
    pub quat_primality_num_iter: u32,
    pub quat_equiv_bound_coeff: i32,
    pub num_alternate_extremal_orders: usize,
    /// Degree of the commitment isogeny.
    pub com_degree: Ibz,
    /// Prime cofactor for non-prime-norm ideal generation.
    pub quat_prime_cofactor: Ibz,
    /// Degree of the secret isogeny (norm of the keygen ideal).
    pub sec_degree: Ibz,
    /// Byte length of 2^TORSION_EVEN_POWER for encoding.
    pub torsion_2power_bytes: usize,
    /// X-coordinate bytes for canonical basis P on E0.
    pub basis_e0_px_bytes: &'static [u8],
    /// X-coordinate bytes for canonical basis Q on E0.
    pub basis_e0_qx_bytes: &'static [u8],
}

/// Generates a private module containing builder helpers that reference
/// the correct level's `endomorphism_action` and `quaternion_data`
/// modules. The inner macros (`get_mat!`, `build_order!`, etc.) use
/// `paste::paste!` to construct constant names from literal indices,
/// which requires compile-time module resolution, hence the per-level
/// module approach.
macro_rules! define_level_builders {
    ($mod_name:ident, $level:ident, [$($idx:literal),+]) => {
        mod $mod_name {
            use super::*;
            use crate::precomp_signing::$level::endomorphism_action;
            use crate::precomp_signing::$level::quaternion_data;

            pub fn get_action_matrix(curve_idx: usize, gen: &str) -> IbzMat2x2 {
                macro_rules! get_mat {
                    ($i:literal, $gen:ident) => {
                        paste::paste! {
                            IbzMat2x2([
                                [
                                    endomorphism_action::[<ENDOMORPHISM_ $i _ACTION_ $gen>]()[0].clone(),
                                    endomorphism_action::[<ENDOMORPHISM_ $i _ACTION_ $gen>]()[1].clone(),
                                ],
                                [
                                    endomorphism_action::[<ENDOMORPHISM_ $i _ACTION_ $gen>]()[2].clone(),
                                    endomorphism_action::[<ENDOMORPHISM_ $i _ACTION_ $gen>]()[3].clone(),
                                ],
                            ])
                        }
                    };
                }

                macro_rules! for_idx {
                    ($i:literal) => {
                        match gen {
                            "I" => get_mat!($i, I),
                            "J" => get_mat!($i, J),
                            "K" => get_mat!($i, K),
                            "GEN2" => get_mat!($i, GEN2),
                            "GEN3" => get_mat!($i, GEN3),
                            "GEN4" => get_mat!($i, GEN4),
                            _ => unreachable!("unknown generator"),
                        }
                    };
                }

                match curve_idx {
                    $( $idx => for_idx!($idx), )+
                    _ => unreachable!("curve index out of range"),
                }
            }

            pub fn build_extremal_order(idx: usize) -> QuatPExtremalMaximalOrder {
                macro_rules! build_order {
                    ($i:literal) => {
                        paste::paste! {
                            {
                                let mut basis = IbzMat4x4::default();
                                for r in 0..4 {
                                    for c in 0..4 {
                                        basis.0[r][c] = match (r, c) {
                                            (0, 0) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_00>]().clone(),
                                            (0, 1) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_01>]().clone(),
                                            (0, 2) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_02>]().clone(),
                                            (0, 3) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_03>]().clone(),
                                            (1, 0) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_10>]().clone(),
                                            (1, 1) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_11>]().clone(),
                                            (1, 2) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_12>]().clone(),
                                            (1, 3) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_13>]().clone(),
                                            (2, 0) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_20>]().clone(),
                                            (2, 1) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_21>]().clone(),
                                            (2, 2) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_22>]().clone(),
                                            (2, 3) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_23>]().clone(),
                                            (3, 0) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_30>]().clone(),
                                            (3, 1) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_31>]().clone(),
                                            (3, 2) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_32>]().clone(),
                                            (3, 3) => quaternion_data::[<EXTREMAL_ORDER_ $i _BASIS_33>]().clone(),
                                            _ => unreachable!(),
                                        };
                                    }
                                }
                                let order = QuatLattice {
                                    denom: quaternion_data::[<EXTREMAL_ORDER_ $i _DENOM>]().clone(),
                                    basis,
                                };
                                let z = QuatAlgElem {
                                    denom: quaternion_data::[<EXTREMAL_ORDER_ $i _I_DENOM>]().clone(),
                                    coord: IbzVec4([
                                        quaternion_data::[<EXTREMAL_ORDER_ $i _I_COORD_0>]().clone(),
                                        quaternion_data::[<EXTREMAL_ORDER_ $i _I_COORD_1>]().clone(),
                                        quaternion_data::[<EXTREMAL_ORDER_ $i _I_COORD_2>]().clone(),
                                        quaternion_data::[<EXTREMAL_ORDER_ $i _I_COORD_3>]().clone(),
                                    ]),
                                };
                                let t = QuatAlgElem {
                                    denom: Ibz::one(),
                                    coord: IbzVec4([Ibz::zero(), Ibz::zero(), Ibz::one(), Ibz::zero()]),
                                };
                                QuatPExtremalMaximalOrder {
                                    q: quaternion_data::[<EXTREMAL_ORDER_ $i _Q>],
                                    order,
                                    z,
                                    t,
                                }
                            }
                        }
                    };
                }

                match idx {
                    $( $idx => build_order!($idx), )+
                    _ => unreachable!("index out of range"),
                }
            }

            pub fn build_connecting_ideal(idx: usize, parent_order: &QuatLattice) -> QuatLeftIdeal {
                macro_rules! build_ideal {
                    ($i:literal) => {
                        paste::paste! {
                            {
                                let mut basis = IbzMat4x4::default();
                                for r in 0..4 {
                                    for c in 0..4 {
                                        basis.0[r][c] = match (r, c) {
                                            (0, 0) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_00>]().clone(),
                                            (0, 1) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_01>]().clone(),
                                            (0, 2) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_02>]().clone(),
                                            (0, 3) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_03>]().clone(),
                                            (1, 0) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_10>]().clone(),
                                            (1, 1) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_11>]().clone(),
                                            (1, 2) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_12>]().clone(),
                                            (1, 3) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_13>]().clone(),
                                            (2, 0) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_20>]().clone(),
                                            (2, 1) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_21>]().clone(),
                                            (2, 2) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_22>]().clone(),
                                            (2, 3) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_23>]().clone(),
                                            (3, 0) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_30>]().clone(),
                                            (3, 1) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_31>]().clone(),
                                            (3, 2) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_32>]().clone(),
                                            (3, 3) => quaternion_data::[<CONNECTING_IDEAL_ $i _BASIS_33>]().clone(),
                                            _ => unreachable!(),
                                        };
                                    }
                                }
                                let lattice = QuatLattice {
                                    denom: quaternion_data::[<CONNECTING_IDEAL_ $i _DENOM>]().clone(),
                                    basis,
                                };
                                QuatLeftIdeal {
                                    lattice,
                                    norm: quaternion_data::[<CONNECTING_IDEAL_ $i _NORM>]().clone(),
                                    parent_order: parent_order.clone(),
                                }
                            }
                        }
                    };
                }

                match idx {
                    $( $idx => build_ideal!($idx), )+
                    _ => unreachable!("index out of range"),
                }
            }

            pub fn build_endomorphism_basis<L: FpBackend>(idx: usize) -> EcBasis<L> {
                macro_rules! build_basis {
                    ($i:literal) => {
                        paste::paste! {
                            {
                                let p = EcPoint::new(
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_P_X_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_P_X_IM>],
                                    ),
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_P_Z_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_P_Z_IM>],
                                    ),
                                );
                                let q = EcPoint::new(
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_Q_X_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_Q_X_IM>],
                                    ),
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_Q_Z_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_Q_Z_IM>],
                                    ),
                                );
                                let pmq = EcPoint::new(
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_PMQ_X_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_PMQ_X_IM>],
                                    ),
                                    Fp2::from_limbs(
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_PMQ_Z_RE>],
                                        &endomorphism_action::[<ENDOMORPHISM_ $i _BASIS_PMQ_Z_IM>],
                                    ),
                                );
                                EcBasis { p, q, pmq }
                            }
                        }
                    };
                }

                match idx {
                    $( $idx => build_basis!($idx), )+
                    _ => unreachable!("curve index out of range"),
                }
            }

            pub fn build_endomorphism_curve<L: FpBackend>(idx: usize) -> EcCurve<L> {
                macro_rules! build_curve {
                    ($i:literal) => {
                        paste::paste! {
                            {
                                let a = Fp2::from_limbs(
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A_RE>],
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A_IM>],
                                );
                                let c = Fp2::from_limbs(
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_C_RE>],
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_C_IM>],
                                );
                                let a24_x = Fp2::from_limbs(
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A24_X_RE>],
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A24_X_IM>],
                                );
                                let a24_z = Fp2::from_limbs(
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A24_Z_RE>],
                                    &endomorphism_action::[<ENDOMORPHISM_ $i _CURVE_A24_Z_IM>],
                                );
                                EcCurve {
                                    a,
                                    c,
                                    a24: EcPoint::new(a24_x, a24_z),
                                    is_a24_computed_and_normalized: false,
                                }
                            }
                        }
                    };
                }

                match idx {
                    $( $idx => build_curve!($idx), )+
                    _ => unreachable!("curve index out of range"),
                }
            }
        }
    };
}

define_level_builders!(l1_builders, level1, [0, 1, 2, 3, 4, 5, 6]);
define_level_builders!(l3_builders, level3, [0, 1, 2, 3, 4, 5, 6, 7]);
define_level_builders!(l5_builders, level5, [0, 1, 2, 3, 4, 5, 6]);

/// Generates a constructor for `SigningPrecomp` using the specified level's
/// precomp modules and the generated builder module.
macro_rules! impl_signing_precomp_constructor {
    ($fn_name:ident, $level:ident, $builders:ident) => {
        impl<L: FpBackend> SigningPrecomp<L> {
            pub fn $fn_name() -> Self {
                use crate::precomp_signing::$level::quaternion_constants;
                use crate::precomp_signing::$level::quaternion_data;
                use crate::precomp_signing::$level::torsion_constants;

                let num_curves = quaternion_constants::NUM_ALTERNATE_EXTREMAL_ORDERS + 1;

                let mut extremal_orders = Vec::with_capacity(num_curves);
                for i in 0..num_curves {
                    extremal_orders.push($builders::build_extremal_order(i));
                }

                let mut connecting_ideals = Vec::with_capacity(num_curves);
                for i in 0..num_curves {
                    connecting_ideals.push($builders::build_connecting_ideal(
                        i,
                        &extremal_orders[0].order,
                    ));
                }

                let mut endomorphism_curves = Vec::with_capacity(num_curves);
                for i in 0..num_curves {
                    endomorphism_curves.push($builders::build_endomorphism_curve::<L>(i));
                }

                let mut endomorphism_bases = Vec::with_capacity(num_curves);
                for i in 0..num_curves {
                    endomorphism_bases.push($builders::build_endomorphism_basis::<L>(i));
                }

                let mut action_matrices = Vec::with_capacity(num_curves);
                for i in 0..num_curves {
                    action_matrices.push(ActionMatrices {
                        i: $builders::get_action_matrix(i, "I"),
                        j: $builders::get_action_matrix(i, "J"),
                        gen2: $builders::get_action_matrix(i, "GEN2"),
                        gen3: $builders::get_action_matrix(i, "GEN3"),
                        gen4: $builders::get_action_matrix(i, "GEN4"),
                    });
                }

                let algebra = QuatAlg::new(&quaternion_data::QUATALG_P());
                let torsion_plus_2power = torsion_constants::TORSION_PLUS_2POWER().clone();

                SigningPrecomp {
                    extremal_orders,
                    connecting_ideals,
                    endomorphism_curves,
                    endomorphism_bases,
                    action_matrices,
                    algebra,
                    torsion_plus_2power,
                    p_cofactor_for_2f: sqisign_verify::precomp::$level::P_COFACTOR_FOR_2F,
                    p_cofactor_for_2f_bitlength:
                        sqisign_verify::precomp::$level::P_COFACTOR_FOR_2F_BITLENGTH as usize,
                    finduv_box_size: quaternion_constants::FINDUV_BOX_SIZE,
                    finduv_cube_size: quaternion_constants::FINDUV_CUBE_SIZE,
                    quat_repres_bound_input: quaternion_constants::QUAT_REPRES_BOUND_INPUT,
                    quat_primality_num_iter: quaternion_constants::QUAT_PRIMALITY_NUM_ITER,
                    quat_equiv_bound_coeff: quaternion_constants::QUAT_EQUIV_BOUND_COEFF,
                    num_alternate_extremal_orders:
                        quaternion_constants::NUM_ALTERNATE_EXTREMAL_ORDERS,
                    com_degree: torsion_constants::COM_DEGREE().clone(),
                    quat_prime_cofactor: quaternion_data::QUAT_PRIME_COFACTOR().clone(),
                    sec_degree: torsion_constants::SEC_DEGREE().clone(),
                    torsion_2power_bytes: torsion_constants::TORSION_2POWER_BYTES,
                    basis_e0_px_bytes: &sqisign_verify::precomp::$level::BASIS_E0_PX_BYTES,
                    basis_e0_qx_bytes: &sqisign_verify::precomp::$level::BASIS_E0_QX_BYTES,
                }
            }
        }
    };
}

impl_signing_precomp_constructor!(level1, level1, l1_builders);
impl_signing_precomp_constructor!(level3, level3, l3_builders);
impl_signing_precomp_constructor!(level5, level5, l5_builders);

/// Trait for constructing [`SigningPrecomp`] generically over security levels.
///
/// Implemented by `Level1`, `Level3`, and `Level5`. Allows generic code
/// to obtain precomputed data without naming the level-specific constructor.
pub trait HasSigningPrecomp: FpBackend {
    fn signing_precomp() -> SigningPrecomp<Self>;
}

impl HasSigningPrecomp for sqisign_verify::params::Level1 {
    fn signing_precomp() -> SigningPrecomp<Self> {
        SigningPrecomp::level1()
    }
}

impl HasSigningPrecomp for sqisign_verify::params::Level3 {
    fn signing_precomp() -> SigningPrecomp<Self> {
        SigningPrecomp::level3()
    }
}

impl HasSigningPrecomp for sqisign_verify::params::Level5 {
    fn signing_precomp() -> SigningPrecomp<Self> {
        SigningPrecomp::level5()
    }
}
