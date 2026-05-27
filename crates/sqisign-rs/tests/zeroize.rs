//! Compile-time and runtime tests for secret key zeroization.

use num_traits::Zero;
use sqisign_rs::keygen::keypair;
use sqisign_rs::keygen::SecretKey;
use sqisign_rs::params::{Level1, Level3, Level5};
use zeroize::{Zeroize, ZeroizeOnDrop};

#[test]
fn secret_key_implements_zeroize_on_drop() {
    fn assert_zod<T: ZeroizeOnDrop>() {}
    assert_zod::<SecretKey<Level1>>();
    assert_zod::<SecretKey<Level3>>();
    assert_zod::<SecretKey<Level5>>();
}

#[test]
fn secret_key_zeroize_clears_data() {
    let mut rng = rand::thread_rng();

    let (_pk, sk) = keypair::<Level1>(&mut rng);

    let sk_bytes = sk.to_bytes().expect("sk encoding must succeed");
    assert!(
        sk_bytes.iter().any(|&b| b != 0),
        "sk should contain non-zero data before zeroize"
    );

    let mut sk = sk;
    sk.zeroize();

    assert!(
        sk.secret_ideal.norm.is_zero(),
        "secret ideal norm should be zero after zeroize"
    );
    assert!(
        sk.secret_ideal.lattice.denom.is_zero(),
        "lattice denom should be zero after zeroize"
    );
    for row in sk.mat_ba_can_to_ba0_two.0.iter() {
        for v in row.iter() {
            assert!(v.is_zero(), "basis-change matrix entry should be zero");
        }
    }
}
