//!
//! This is a deterministic PRNG used by the NIST KAT generation tool.
//! It is NOT part of the SQIsign protocol. It exists solely to reproduce
//! byte-identical output from the NIST KAT generator when seeded with the
//! same 48-byte entropy input from the `.rsp` files.

use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;

/// NIST AES-256-CTR-DRBG state.
pub struct NistDrbg {
    key: [u8; 32],
    v: [u8; 16],
    reseed_counter: u64,
}

impl NistDrbg {
    /// Initialize from a 48-byte seed (entropy_input).
    ///
    /// Equivalent to the NIST `randombytes_init(seed, NULL, 256)` API.
    pub fn new(seed: &[u8; 48]) -> Self {
        let mut drbg = NistDrbg {
            key: [0u8; 32],
            v: [0u8; 16],
            reseed_counter: 1,
        };
        Self::update(Some(seed), &mut drbg.key, &mut drbg.v);
        drbg
    }

    /// AES-256-CTR-DRBG Update function.
    fn update(provided_data: Option<&[u8; 48]>, key: &mut [u8; 32], v: &mut [u8; 16]) {
        let mut temp = [0u8; 48];

        for i in 0..3 {
            // Increment V (big-endian counter)
            increment_v(v);

            let cipher = Aes256::new((&*key).into());
            let mut block = aes::Block::clone_from_slice(v);
            cipher.encrypt_block(&mut block);
            temp[16 * i..16 * (i + 1)].copy_from_slice(&block);
        }

        if let Some(data) = provided_data {
            for i in 0..48 {
                temp[i] ^= data[i];
            }
        }

        key.copy_from_slice(&temp[..32]);
        v.copy_from_slice(&temp[32..48]);
    }

    /// Fill `buf` with deterministic random bytes.
    ///
    /// Equivalent to the NIST `randombytes(buf, len)` API.
    pub fn fill(&mut self, buf: &mut [u8]) {
        let mut offset = 0;
        let mut remaining = buf.len();

        while remaining > 0 {
            increment_v(&mut self.v);

            let cipher = Aes256::new((&self.key as &[u8; 32]).into());
            let mut block = aes::Block::clone_from_slice(&self.v);
            cipher.encrypt_block(&mut block);

            if remaining > 15 {
                buf[offset..offset + 16].copy_from_slice(&block);
                offset += 16;
                remaining -= 16;
            } else {
                buf[offset..offset + remaining].copy_from_slice(&block[..remaining]);
                remaining = 0;
            }
        }

        Self::update(None, &mut self.key, &mut self.v);
        self.reseed_counter += 1;
    }
}

/// Increment V as a big-endian 128-bit counter.
fn increment_v(v: &mut [u8; 16]) {
    for j in (0..16).rev() {
        if v[j] == 0xff {
            v[j] = 0x00;
        } else {
            v[j] += 1;
            break;
        }
    }
}

// Implement rand_core traits so NistDrbg can be passed to keygen.

impl rand_core::RngCore for NistDrbg {
    fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill(&mut buf);
        u32::from_le_bytes(buf)
    }

    fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill(&mut buf);
        u64::from_le_bytes(buf)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.fill(dest);
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill(dest);
        Ok(())
    }
}

impl rand_core::CryptoRng for NistDrbg {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_drbg_deterministic() {
        let seed = [0u8; 48];
        let mut drbg1 = NistDrbg::new(&seed);
        let mut drbg2 = NistDrbg::new(&seed);

        let mut buf1 = [0u8; 64];
        let mut buf2 = [0u8; 64];
        drbg1.fill(&mut buf1);
        drbg2.fill(&mut buf2);

        assert_eq!(buf1, buf2);
    }

    #[test]
    fn test_drbg_kat_seed_generation() {
        // The C KAT generator seeds with entropy_input = [0, 1, 2, ..., 47],
        // then calls randombytes(seed, 48) to get the first test case's seed.
        // The first seed in the Level 1 .rsp file is:
        // 061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1
        let mut entropy_input = [0u8; 48];
        for (i, byte) in entropy_input.iter_mut().enumerate() {
            *byte = i as u8;
        }
        let mut outer_drbg = NistDrbg::new(&entropy_input);

        let mut seed = [0u8; 48];
        outer_drbg.fill(&mut seed);

        let expected = hex::decode(
            "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1"
        ).unwrap();

        assert_eq!(&seed[..], &expected[..], "first KAT seed mismatch");
    }
}
