#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{Level1, PublicKey};

fuzz_target!(|data: &[u8]| {
    let _ = PublicKey::<Level1>::from_bytes(data);
});
