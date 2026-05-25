#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{Level5, PublicKey};

fuzz_target!(|data: &[u8]| {
    let _ = PublicKey::<Level5>::from_bytes(data);
});
