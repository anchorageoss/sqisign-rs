//! The `E0` endomorphism actions of the three SQIsign round-3 levels, from
//! the generated constants.

use super::E0Actions;
use crate::precomp;

/// Level I: `p324_3`.
pub fn level1() -> E0Actions<{ precomp::p324_3::IBZ_NLIMBS }> {
    use precomp::p324_3 as c;
    E0Actions::from_bytes(
        &c::ACTION_I,
        &c::ACTION_J,
        &c::ACTION_GEN2,
        &c::ACTION_GEN3,
        &c::ACTION_GEN4,
    )
}

/// Level III: `p500_27`.
pub fn level3() -> E0Actions<{ precomp::p500_27::IBZ_NLIMBS }> {
    use precomp::p500_27 as c;
    E0Actions::from_bytes(
        &c::ACTION_I,
        &c::ACTION_J,
        &c::ACTION_GEN2,
        &c::ACTION_GEN3,
        &c::ACTION_GEN4,
    )
}

/// Level V: `p664_17`.
pub fn level5() -> E0Actions<{ precomp::p664_17::IBZ_NLIMBS }> {
    use precomp::p664_17 as c;
    E0Actions::from_bytes(
        &c::ACTION_I,
        &c::ACTION_J,
        &c::ACTION_GEN2,
        &c::ACTION_GEN3,
        &c::ACTION_GEN4,
    )
}
