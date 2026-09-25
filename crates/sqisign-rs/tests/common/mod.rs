//! A deterministic `rand_core` generator over the crate's SHAKE stream,
//! shared by the integration tests.

#![allow(dead_code)]

use sqisign_rs::mp::ShakeRng;

pub struct Det(pub ShakeRng);

impl rand_core::RngCore for Det {
    fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.fill_bytes(&mut b);
        u32::from_le_bytes(b)
    }
    fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.fill_bytes(&mut b);
        u64::from_le_bytes(b)
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        use sqisign_verify::rng::Rng;
        assert!(self.0.fill(dest));
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}
impl rand_core::CryptoRng for Det {}

pub fn rng(label: &str, domain: &[u8]) -> Det {
    Det(ShakeRng::new(label.as_bytes(), domain))
}

/// Levels III and V sign slowly in a debug build; they run in release
/// builds and under `SQISIGN_ALL_LEVELS=1`.
pub fn all_levels() -> bool {
    !cfg!(debug_assertions) || std::env::var("SQISIGN_ALL_LEVELS").is_ok_and(|v| v == "1")
}
