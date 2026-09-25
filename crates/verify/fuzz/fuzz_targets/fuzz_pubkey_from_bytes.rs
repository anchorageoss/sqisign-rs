#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::sqisign::{level1, public_key_from_bytes, PreparedPublicKey};
use sqisign_verify::P324_3;

// Arbitrary bytes into the round-3 public-key decoder and the key
// preparation (bounded basis search, fails closed on crafted curves).
fuzz_target!(|data: &[u8]| {
    let params = level1();
    if let Some(pk) = public_key_from_bytes::<P324_3>(&params, data) {
        let _ = PreparedPublicKey::new(&params, &pk);
    }
});
