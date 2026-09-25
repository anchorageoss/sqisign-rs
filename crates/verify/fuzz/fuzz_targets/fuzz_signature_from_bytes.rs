#![no_main]
use libfuzzer_sys::fuzz_target;
use sqisign_verify::sqisign::{level1, signature_from_bytes, signature_to_bytes};
use sqisign_verify::P324_3;

// Arbitrary bytes into the round-3 signature decoder; what decodes must
// re-encode to the same bytes.
fuzz_target!(|data: &[u8]| {
    let params = level1();
    if let Some(sig) = signature_from_bytes::<P324_3>(&params, data) {
        let mut out = [0u8; 200];
        let n = signature_to_bytes(&params, &sig, &mut out).expect("fits");
        assert_eq!(&out[..n], data);
    }
});
