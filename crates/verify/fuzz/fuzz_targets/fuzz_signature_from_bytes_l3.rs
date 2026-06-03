#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{Level3, Signature};

fuzz_target!(|data: &[u8]| {
    let _ = Signature::<Level3>::from_bytes(data);
});
