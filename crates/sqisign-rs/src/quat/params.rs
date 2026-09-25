//! The quaternion parameter sets of the three SQIsign round-3 levels, built
//! from the generated constants.

use super::QuatAlg;
use crate::precomp;
use sqisign_verify::params::{sqisign_v3, Prime, P324_3, P500_27, P664_17};

/// Level I: `p324_3`, 30 limbs.
pub type QuatAlgLevel1 = QuatAlg<{ precomp::p324_3::IBZ_NLIMBS }>;
/// Level III: `p500_27`, 40 limbs.
pub type QuatAlgLevel3 = QuatAlg<{ precomp::p500_27::IBZ_NLIMBS }>;
/// Level V: `p664_17`, 55 limbs.
pub type QuatAlgLevel5 = QuatAlg<{ precomp::p664_17::IBZ_NLIMBS }>;

/// The level I parameter set.
pub fn level1() -> QuatAlgLevel1 {
    use sqisign_v3::p324_3 as l;
    QuatAlg::new(
        P324_3::prime_le_bytes(),
        &precomp::p324_3::SQRT_P_LE,
        &l::RI_COFACTOR_LE,
        l::LAMBDA as i32,
        l::PRIMALITY_NUM_ITER,
        l::EQUIV_BOUND_COEFF as i32,
        precomp::p324_3::QLAPOTY_USED_POWER_OF_TWO,
    )
}

/// The level III parameter set.
pub fn level3() -> QuatAlgLevel3 {
    use sqisign_v3::p500_27 as l;
    QuatAlg::new(
        P500_27::prime_le_bytes(),
        &precomp::p500_27::SQRT_P_LE,
        &l::RI_COFACTOR_LE,
        l::LAMBDA as i32,
        l::PRIMALITY_NUM_ITER,
        l::EQUIV_BOUND_COEFF as i32,
        precomp::p500_27::QLAPOTY_USED_POWER_OF_TWO,
    )
}

/// The level V parameter set.
pub fn level5() -> QuatAlgLevel5 {
    use sqisign_v3::p664_17 as l;
    QuatAlg::new(
        P664_17::prime_le_bytes(),
        &precomp::p664_17::SQRT_P_LE,
        &l::RI_COFACTOR_LE,
        l::LAMBDA as i32,
        l::PRIMALITY_NUM_ITER,
        l::EQUIV_BOUND_COEFF as i32,
        precomp::p664_17::QLAPOTY_USED_POWER_OF_TWO,
    )
}
