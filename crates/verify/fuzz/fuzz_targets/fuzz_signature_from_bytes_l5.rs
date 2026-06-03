#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::{Level5, Signature};

fuzz_target!(|data: &[u8]| {
    let _ = Signature::<Level5>::from_bytes(data);
});
