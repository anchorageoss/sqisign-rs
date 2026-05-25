#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{Level3, PublicKey, Signature};

const PK_BYTES: &[u8] = include_bytes!("l3_pk.bin");

fuzz_target!(|data: &[u8]| {
    let pk = match PublicKey::<Level3>::from_bytes(PK_BYTES) {
        Ok(pk) => pk,
        Err(_) => return,
    };
    let sig = match Signature::<Level3>::from_bytes(data) {
        Ok(sig) => sig,
        Err(_) => return,
    };
    let _ = sig.verify(&pk, b"fuzz message");
});
