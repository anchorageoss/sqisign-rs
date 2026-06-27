# sqisign-rs

[![crates.io](https://img.shields.io/crates/v/sqisign-rs.svg)](https://crates.io/crates/sqisign-rs)
[![docs.rs](https://docs.rs/sqisign-rs/badge.svg)](https://docs.rs/sqisign-rs)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

A pure Rust implementation of SQIsign v2.0, the post-quantum signature scheme with the smallest signatures of any NIST PQC candidate, down to 108 bytes. Fully `no_std` (uses `alloc`, requires no OS) and passes all 300 NIST KAT vectors across Levels 1, 3, and 5.

> **Not audited.** The verification path is designed to be constant-time; signing is inherently variable-time. See [SECURITY.md](SECURITY.md). Use at your own discretion.

## Quick start

```rust
use sqisign_rs::{generate, PublicKey, SigningKey, Verifier};

let mut rng = rand::rngs::OsRng;

// Level 1 (the default): generate, sign, and verify in three lines.
let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;

// Compress to a 129-byte wire format; it verifies with the same call.
pk.verify(b"hello world", &sig.compress())?;
```

The compact scheme produces the smallest signature (108 bytes), and higher security levels are a single type parameter:

```rust
use sqisign_rs::{generate, generate_compact, Level3, Verifier};

// Compact: 108-byte signatures, verified via a dimension-4 isogeny.
let (cpk, csk) = generate_compact(&mut rng);
let csig = csk.sign(b"hello world", &mut rng)?;
cpk.verify(b"hello world", &csig)?;

// Levels 3 and 5 (dimension-2 formats):
let (pk, sk) = generate::<Level3>(&mut rng); // or Level5
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;
```

Verification auto-detects the wire format from its byte length. Standard and compact public keys are separate, non-interchangeable schemes (chosen at keygen).

## Performance

A signature can be carried in four formats that trade wire size for verification speed. Level 1:

| Format             | Size  | Verify   | Levels |
|--------------------|-------|----------|--------|
| Compact (dim-4)    | 108 B | ~33 ms   | 1      |
| Compressed         | 129 B | ~1.4 ms  | 1/3/5  |
| Standard (default) | 148 B | ~1.4 ms  | 1/3/5  |
| Expanded           | 212 B | ~1.4 ms  | 1/3/5  |

Dimension-2 verification (the three larger formats) is at parity with the C reference on portable builds with fat LTO (~1.4 ms at L1, ~8 ms at L5 on an Apple M4 Pro; expanded is fastest, compressed slowest). The compact format trades a heavier dimension-4 verify for the smallest signature. Key generation and signing are ~3-4x the C reference (the pure-Rust `num-bigint` quaternion layer instead of GMP) and affect only the signer, never verification.

For a standalone constant-time, zero-allocation dim-2 verifier, depend on [`sqisign-verify`](https://crates.io/crates/sqisign-verify) directly. The compact format is documented in [COMPRESSION.md](COMPRESSION.md).

## About

SQIsign is advancing through NIST's post-quantum signature standardization and has **not** been standardized. It is the only isogeny-based candidate, and the only one with signatures this small. The 2022 attacks that broke the SIDH/SIKE key exchange do not apply to it; in fact the higher-dimensional isogeny techniques those attacks introduced are now used constructively to build and speed up SQIsign (the dimension-2 and dimension-4 constructions in this library), so the mathematics that ended SIKE strengthens SQIsign. The scheme and its implementations are also young and moving fast, with substantial engineering headroom still to capture (faster signing, optimized field and quaternion backends); this library tracks that progress.

## License

Apache-2.0 OR MIT
