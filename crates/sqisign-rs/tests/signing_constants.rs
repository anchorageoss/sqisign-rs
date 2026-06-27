//! Validation tests for signing-layer precomputed constants.
//!
//! Verifies that all constants load successfully and have expected properties.
//! These tests require the `signing` feature.

use num_bigint::BigInt;
use num_traits::{One, Zero};

// ====================================================================
// Level 1
// ====================================================================
mod level1 {
    use super::*;
    use sqisign_rs::precomp_signing::level1::{
        endomorphism_action, quaternion_constants, quaternion_data, torsion_constants,
    };

    #[test]
    fn torsion_constants_load() {
        assert!(!torsion_constants::TWO_TO_SECURITY_BITS().is_zero());
        assert!(!torsion_constants::TORSION_PLUS_2POWER().is_zero());
        assert!(!torsion_constants::SEC_DEGREE().is_zero());
        assert!(!torsion_constants::COM_DEGREE().is_zero());

        assert_eq!(
            torsion_constants::TWO_TO_SECURITY_BITS(),
            BigInt::one() << 128
        );
        assert_ne!(torsion_constants::TORSION_2POWER_BYTES, 0);
    }

    #[test]
    fn quaternion_constants_values() {
        assert_eq!(quaternion_constants::NUM_EXTREMAL_ORDERS, 7);
        assert_eq!(quaternion_constants::NUM_ALTERNATE_EXTREMAL_ORDERS, 6);
        assert_eq!(quaternion_constants::QUAT_PRIMALITY_NUM_ITER, 32);
    }

    #[test]
    fn quaternion_data_loads() {
        assert!(!quaternion_data::QUATALG_P().is_zero());
        assert!(!quaternion_data::QUAT_PRIME_COFACTOR().is_zero());
        assert_eq!(quaternion_data::NUM_ORDERS, 7);

        assert_eq!(quaternion_data::EXTREMAL_ORDER_0_Q, 1);
        assert!(!quaternion_data::EXTREMAL_ORDER_0_DENOM().is_zero());

        assert!(
            !quaternion_data::EXTREMAL_ORDER_0_BASIS_00().is_zero(),
            "basis[0][0] should be nonzero"
        );

        assert!(!quaternion_data::CONNECTING_IDEAL_0_DENOM().is_zero());
        assert!(!quaternion_data::CONNECTING_IDEAL_0_NORM().is_zero());
        assert!(!quaternion_data::CONJUGATING_ELEM_0_DENOM().is_zero());
    }

    #[test]
    fn quaternion_data_p_matches_torsion() {
        let p = &quaternion_data::QUATALG_P();
        let tors2 = &torsion_constants::TORSION_PLUS_2POWER();
        assert!(tors2 > &BigInt::zero());
        assert!(
            (p + BigInt::one()) % tors2 == BigInt::zero(),
            "p+1 should be divisible by the 2-power torsion"
        );
    }

    #[test]
    fn endomorphism_action_curve0_is_e0() {
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_RE, [0u64; 5]);
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_IM, [0u64; 5]);
        assert_ne!(endomorphism_action::ENDOMORPHISM_0_CURVE_C_RE, [0u64; 5]);
        assert_eq!(endomorphism_action::NUM_CURVES_WITH_ENDOMORPHISMS, 7);
    }

    #[test]
    fn endomorphism_action_matrices_load() {
        for entry in endomorphism_action::ENDOMORPHISM_0_ACTION_I().iter() {
            let _val: &BigInt = entry;
        }
        for entry in endomorphism_action::ENDOMORPHISM_0_ACTION_J().iter() {
            let _val: &BigInt = entry;
        }
    }

    #[test]
    fn e0_basis_limbs_nonzero() {
        use sqisign_rs::precomp::level1::e0_basis;
        assert_ne!(e0_basis::BASIS_E0_P_X_RE, [0u64; 5]);
        assert_ne!(e0_basis::BASIS_E0_P_X_IM, [0u64; 5]);
        assert_ne!(e0_basis::BASIS_E0_Q_X_RE, [0u64; 5]);
        assert_ne!(e0_basis::BASIS_E0_Q_X_IM, [0u64; 5]);
    }

    #[test]
    fn ec_params_match_existing() {
        use sqisign_rs::precomp::level1::ec_params;
        assert_eq!(ec_params::TORSION_EVEN_POWER, 248);
        assert_eq!(ec_params::P_COFACTOR_FOR_2F_BITLENGTH, 3);
        assert_eq!(ec_params::P_COFACTOR_FOR_2F, &[0x5u64]);
    }
}

// ====================================================================
// Level 3
// ====================================================================
mod level3 {
    use super::*;
    use sqisign_rs::precomp_signing::level3::{
        endomorphism_action, quaternion_constants, quaternion_data, torsion_constants,
    };

    #[test]
    fn torsion_constants_load() {
        assert!(!torsion_constants::TWO_TO_SECURITY_BITS().is_zero());
        assert!(!torsion_constants::TORSION_PLUS_2POWER().is_zero());
        assert!(!torsion_constants::SEC_DEGREE().is_zero());
        assert!(!torsion_constants::COM_DEGREE().is_zero());

        assert_eq!(
            torsion_constants::TWO_TO_SECURITY_BITS(),
            BigInt::one() << 192
        );
    }

    #[test]
    fn quaternion_constants_values() {
        assert_eq!(quaternion_constants::NUM_EXTREMAL_ORDERS, 8);
        assert_eq!(quaternion_constants::NUM_ALTERNATE_EXTREMAL_ORDERS, 7);
    }

    #[test]
    fn quaternion_data_loads() {
        assert!(!quaternion_data::QUATALG_P().is_zero());
        assert_eq!(quaternion_data::NUM_ORDERS, 8);
        assert_eq!(quaternion_data::EXTREMAL_ORDER_0_Q, 1);
    }

    #[test]
    fn quaternion_data_p_matches_torsion() {
        let p = &quaternion_data::QUATALG_P();
        let tors2 = &torsion_constants::TORSION_PLUS_2POWER();
        assert!(
            (p + BigInt::one()) % tors2 == BigInt::zero(),
            "p+1 should be divisible by the 2-power torsion"
        );
    }

    #[test]
    fn endomorphism_action_curve0_is_e0() {
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_RE, [0u64; 7]);
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_IM, [0u64; 7]);
        assert_ne!(endomorphism_action::ENDOMORPHISM_0_CURVE_C_RE, [0u64; 7]);
        assert_eq!(endomorphism_action::NUM_CURVES_WITH_ENDOMORPHISMS, 8);
    }

    #[test]
    fn ec_params() {
        use sqisign_rs::precomp::level3::ec_params;
        assert_eq!(ec_params::TORSION_EVEN_POWER, 376);
        assert_eq!(ec_params::P_COFACTOR_FOR_2F, &[0x41u64]);
    }

    #[test]
    fn e0_basis_limbs_nonzero() {
        use sqisign_rs::precomp::level3::e0_basis;
        assert_ne!(e0_basis::BASIS_E0_P_X_RE, [0u64; 7]);
        assert_ne!(e0_basis::BASIS_E0_Q_X_RE, [0u64; 7]);
    }
}

// ====================================================================
// Level 5
// ====================================================================
mod level5 {
    use super::*;
    use sqisign_rs::precomp_signing::level5::{
        endomorphism_action, quaternion_constants, quaternion_data, torsion_constants,
    };

    #[test]
    fn torsion_constants_load() {
        assert!(!torsion_constants::TWO_TO_SECURITY_BITS().is_zero());
        assert!(!torsion_constants::TORSION_PLUS_2POWER().is_zero());
        assert!(!torsion_constants::SEC_DEGREE().is_zero());
        assert!(!torsion_constants::COM_DEGREE().is_zero());

        assert_eq!(
            torsion_constants::TWO_TO_SECURITY_BITS(),
            BigInt::one() << 256
        );
    }

    #[test]
    fn quaternion_constants_values() {
        assert_eq!(quaternion_constants::NUM_EXTREMAL_ORDERS, 7);
        assert_eq!(quaternion_constants::NUM_ALTERNATE_EXTREMAL_ORDERS, 6);
    }

    #[test]
    fn quaternion_data_loads() {
        assert!(!quaternion_data::QUATALG_P().is_zero());
        assert_eq!(quaternion_data::NUM_ORDERS, 7);
        assert_eq!(quaternion_data::EXTREMAL_ORDER_0_Q, 1);
    }

    #[test]
    fn quaternion_data_p_matches_torsion() {
        let p = &quaternion_data::QUATALG_P();
        let tors2 = &torsion_constants::TORSION_PLUS_2POWER();
        assert!(
            (p + BigInt::one()) % tors2 == BigInt::zero(),
            "p+1 should be divisible by the 2-power torsion"
        );
    }

    #[test]
    fn endomorphism_action_curve0_is_e0() {
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_RE, [0u64; 9]);
        assert_eq!(endomorphism_action::ENDOMORPHISM_0_CURVE_A_IM, [0u64; 9]);
        assert_ne!(endomorphism_action::ENDOMORPHISM_0_CURVE_C_RE, [0u64; 9]);
        assert_eq!(endomorphism_action::NUM_CURVES_WITH_ENDOMORPHISMS, 7);
    }

    #[test]
    fn ec_params() {
        use sqisign_rs::precomp::level5::ec_params;
        assert_eq!(ec_params::TORSION_EVEN_POWER, 500);
        assert_eq!(ec_params::P_COFACTOR_FOR_2F, &[0x1bu64]);
    }

    #[test]
    fn e0_basis_limbs_nonzero() {
        use sqisign_rs::precomp::level5::e0_basis;
        assert_ne!(e0_basis::BASIS_E0_P_X_RE, [0u64; 9]);
        assert_ne!(e0_basis::BASIS_E0_Q_X_RE, [0u64; 9]);
    }
}
