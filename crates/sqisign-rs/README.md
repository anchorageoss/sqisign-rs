# sqisign-rs

[![crates.io](https://img.shields.io/crates/v/sqisign-rs.svg)](https://crates.io/crates/sqisign-rs)
[![docs.rs](https://docs.rs/sqisign-rs/badge.svg)](https://docs.rs/sqisign-rs)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

A pure Rust implementation of SQIsign v2.0.

SQIsign is a post-quantum digital signature scheme based on isogenies between supersingular elliptic curves. It was proposed for NIST PQC standardization and produces the smallest signatures of any post-quantum candidate at comparable security levels.

This implementation passes all 300 NIST KAT vectors (100 per level) across Level 1, Level 3, and Level 5. Verification is 29-41% faster than the C reference implementation.

> **This library has not been audited.** The verification path is designed to be constant-time, but no formal verification or third-party audit has been performed. The signing path is inherently variable-time due to SQIsign's algorithmic structure. Use at your own discretion.

For verify-only usage (`no_std`, no heap), depend on [`sqisign-verify`](https://crates.io/crates/sqisign-verify) directly.

## Quick start

Generate a keypair, sign, and verify:

```rust
use sqisign_rs::{generate, PublicKey, SigningKey, Verifier};

fn main() -> Result<(), sqisign_rs::Error> {
    let mut rng = rand::rngs::OsRng;
    let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
    let sig = sk.sign(b"hello world", &mut rng)?;
    pk.verify(b"hello world", &sig)?;
    Ok(())
}
```

Compress a signature for minimal wire size (129 bytes at Level 1):

```rust
use sqisign_rs::{Signature, PublicKey, Verifier};

fn example(sig: &Signature, pk: &PublicKey) -> Result<(), sqisign_rs::Error> {
    let compressed = sig.compress();
    pk.verify(b"hello world", &compressed)?;
    Ok(())
}
```

Use a higher security level:

```rust
use sqisign_rs::{generate, Level3, Verifier};

fn main() -> Result<(), sqisign_rs::Error> {
    let mut rng = rand::rngs::OsRng;
    let (pk, sk) = generate::<Level3>(&mut rng);
    let sig = sk.sign(b"hello world", &mut rng)?;
    pk.verify(b"hello world", &sig)?;
    Ok(())
}
```

## Signature formats

Three wire formats with different size/speed tradeoffs. The verifier determines the format from the byte length.

| Format | L1 | L3 | L5 |
|---|---|---|---|
| Compressed | 129 bytes | 196 bytes | 257 bytes |
| Standard | 148 bytes | 224 bytes | 292 bytes |
| Expanded | 212 bytes | 316 bytes | 420 bytes |

Public keys are 65 bytes (L1), 97 bytes (L3), 129 bytes (L5) across all formats.

## Performance

Intel Xeon @ 2.80 GHz, `--release`, `target-cpu=native`:

| Operation | L1 | vs C reference |
|---|---|---|
| Verify (expanded) | 3.82 ms | 41% faster |
| Verify (standard) | 4.65 ms | 29% faster |
| Verify (compressed) | 6.83 ms | comparable |
| Key generation | ~185 ms | ~2.5x slower |
| Signing | ~185 ms | ~2.5x slower |

Verification outperforms the C reference because LLVM aggressively inlines field arithmetic across crate boundaries. Signing and keygen are slower because the quaternion algebra layer uses `num-bigint` (heap-allocated arbitrary precision integers) instead of GMP. This is an explicit tradeoff to keep the library pure Rust with no FFI dependencies.

## Memory security

Secret key material is protected by a three-tier system:

1. `SecretKey` implements `ZeroizeOnDrop`. The key is zeroed when it goes out of scope.
2. Intermediate signing values are explicitly zeroed at the end of the signing protocol.
3. A `ZeroizingAllocator` clears all freed heap memory, catching residual copies left by `num-bigint`'s internal reallocations.

The allocator is enabled by default with zero measured overhead. Disable it if you use a custom allocator:

```toml
sqisign-rs = { version = "0.2", default-features = false }
```

## References

SQIsign was introduced by De Feo, Kohel, Leroux, Petit, and Wesolowski (ASIACRYPT 2020). The v2.0 construction uses (2,2)-isogenies in the theta model following the SQIsign2D-West approach (Basso, Dartois, De Feo, Leroux, Maino, Pope, Robert, Wesolowski, ASIACRYPT 2024).

## License

Apache-2.0 OR MIT
