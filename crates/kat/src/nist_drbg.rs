//! The NIST AES-256 CTR-DRBG of the KAT generator (`rng.c` of the NIST
//! submission package), so the `.rsp` files can be reproduced: seeded with
//! an entry's 48-byte `seed`, its first `randombytes(48)` is key
//! generation's root seed and its second is signing's. Not part of any
//! protocol.

use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;
use sqisign_verify::rng::Rng;

/// AES-256 CTR-DRBG state: `randombytes_init(seed, NULL, 256)`.
pub struct NistDrbg {
    key: [u8; 32],
    v: [u8; 16],
}

impl NistDrbg {
    /// `randombytes_init(seed, NULL, 256)`.
    pub fn new(seed: &[u8; 48]) -> Self {
        let mut d = Self {
            key: [0u8; 32],
            v: [0u8; 16],
        };
        d.update(Some(seed));
        d
    }

    fn update(&mut self, provided: Option<&[u8; 48]>) {
        let cipher = Aes256::new((&self.key).into());
        let mut temp = [0u8; 48];
        for chunk in temp.chunks_mut(16) {
            increment(&mut self.v);
            let mut block = aes::Block::clone_from_slice(&self.v);
            cipher.encrypt_block(&mut block);
            chunk.copy_from_slice(&block);
        }
        if let Some(data) = provided {
            for (t, d) in temp.iter_mut().zip(data.iter()) {
                *t ^= d;
            }
        }
        self.key.copy_from_slice(&temp[..32]);
        self.v.copy_from_slice(&temp[32..]);
    }

    /// `randombytes(out, len)`.
    pub fn generate(&mut self, out: &mut [u8]) {
        let cipher = Aes256::new((&self.key).into());
        for chunk in out.chunks_mut(16) {
            increment(&mut self.v);
            let mut block = aes::Block::clone_from_slice(&self.v);
            cipher.encrypt_block(&mut block);
            chunk.copy_from_slice(&block[..chunk.len()]);
        }
        self.update(None);
    }
}

fn increment(v: &mut [u8; 16]) {
    for b in v.iter_mut().rev() {
        if *b == 0xff {
            *b = 0;
        } else {
            *b += 1;
            break;
        }
    }
}

impl Rng for NistDrbg {
    fn fill(&mut self, out: &mut [u8]) -> bool {
        self.generate(out);
        true
    }
}
