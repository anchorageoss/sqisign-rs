#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::sqisign::{level1, public_key_from_bytes, signature_from_bytes, verify};
use sqisign_verify::P324_3;

const PK_BYTES: &[u8] = include_bytes!("l1_pk.bin");

// Arbitrary signature bytes against a fixed valid round-3 public key: the
// verifier must reject or accept without panicking.
fuzz_target!(|data: &[u8]| {
    let params = level1();
    let pk = public_key_from_bytes::<P324_3>(&params, PK_BYTES).expect("KAT key");
    if let Some(sig) = signature_from_bytes::<P324_3>(&params, data) {
        let _ = verify(&params, &pk, &sig, b"fuzz message");
    }
});
