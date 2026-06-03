#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{CompressedSignature, Level5, PublicKey, Verifier};

const PK_BYTES: &[u8] = include_bytes!("l5_pk.bin");

fuzz_target!(|data: &[u8]| {
    let pk = match PublicKey::<Level5>::from_bytes(PK_BYTES) {
        Ok(pk) => pk,
        Err(_) => return,
    };
    let sig = match CompressedSignature::<Level5>::from_bytes(data) {
        Ok(sig) => sig,
        Err(_) => return,
    };
    let _ = pk.verify(b"fuzz message", &sig);
});
