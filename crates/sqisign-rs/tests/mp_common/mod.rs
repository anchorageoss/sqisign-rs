//! Shared helpers for the integer-layer tests: the SHAKE256 input source
//! matching `tools/c-validate/gen_mp_vectors.c`, conversions to and from
//! `num-bigint`, and the vector text format.

#![allow(dead_code)]

use num_bigint::{BigInt, BigUint, Sign};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use sqisign_rs::mp::{Ibz, Rng, ShakeRng};

pub struct DetRng {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl DetRng {
    pub fn new(label: &str) -> Self {
        let mut h = Shake256::default();
        h.update(b"prism-mp-test-rng/");
        h.update(label.as_bytes());
        Self {
            reader: h.finalize_xof(),
        }
    }

    pub fn fill(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }

    pub fn u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.fill(&mut b);
        u32::from_le_bytes(b)
    }

    /// Same construction as the harness's `rng_ibz`.
    pub fn ibz<const N: usize>(&mut self, bits: i32, signed: bool) -> Ibz<N> {
        let nbytes = ((bits + 7) / 8) as usize;
        let mut buf = [0u8; 8 * 64];
        self.fill(&mut buf[..nbytes]);
        let mut d = [0u64; 64];
        for i in 0..nbytes {
            d[i / 8] |= (buf[i] as u64) << (8 * (i % 8));
        }
        if bits > 0 {
            let top = (bits - 1) as usize;
            d[top / 64] &= if top % 64 == 63 {
                u64::MAX
            } else {
                (1u64 << (top % 64 + 1)) - 1
            };
            d[top / 64] |= 1u64 << (top % 64);
        }
        let n = ((bits + 63) / 64) as usize;
        let mut x = Ibz::<N>::from_bits(&d[..n.max(1)], bits);
        if signed {
            let mut s = [0u8; 1];
            self.fill(&mut s);
            if s[0] & 1 == 1 {
                x = x.neg();
            }
        }
        x
    }
}

impl Rng for DetRng {
    fn fill(&mut self, out: &mut [u8]) -> bool {
        DetRng::fill(self, out);
        true
    }
}

/// The harness's `seed_default_domain`.
pub fn default_domain(label: &str) -> ShakeRng {
    let mut h = Shake256::default();
    h.update(b"prism-mp-prng-seed/");
    h.update(label.as_bytes());
    let mut seed = [0u8; 48];
    h.finalize_xof().read(&mut seed);
    ShakeRng::new(&seed, b"dom")
}

pub const SIZES: [i32; 17] = [
    1, 5, 31, 32, 63, 64, 65, 127, 128, 200, 326, 400, 511, 512, 513, 700, 900,
];

pub fn to_big<const N: usize>(x: &Ibz<N>) -> BigInt {
    let mut buf = [0u8; 2 * 64 * 16 + 2];
    let s = x.to_hex(&mut buf);
    let (neg, digits) = if let Some(rest) = s.strip_prefix('-') {
        (true, rest)
    } else {
        (false, s)
    };
    let mag = BigUint::parse_bytes(digits.as_bytes(), 16).expect("hex");
    BigInt::from_biguint(if neg { Sign::Minus } else { Sign::Plus }, mag)
}

pub fn from_big<const N: usize>(v: &BigInt) -> Ibz<N> {
    let s = v.to_str_radix(16);
    Ibz::<N>::from_str_radix(&s, 16).expect("fits")
}

/// Parse `[-]hex:bitlen`.
pub fn parse_ibz<const N: usize>(tok: &str) -> Ibz<N> {
    let (h, b) = tok.split_once(':').expect("hex:bitlen");
    let mut x =
        Ibz::<N>::from_str_radix(h, 16).unwrap_or_else(|| panic!("bad integer token `{tok}`"));
    let bound: i32 = b.parse().expect("bitlen");
    x.set_bound(bound);
    x
}

/// `[-]hex:bitlen` of an integer.
pub fn fmt_ibz<const N: usize>(x: &Ibz<N>) -> String {
    let mut buf = [0u8; 2 * 64 * 16 + 2];
    format!("{}:{}", x.to_hex(&mut buf), x.get_bound())
}

pub fn ibz_eq_exact<const N: usize>(a: &Ibz<N>, b: &Ibz<N>) -> bool {
    a == b && a.get_bound() == b.get_bound()
}
