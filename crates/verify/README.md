# sqisign-verify

[![crates.io](https://img.shields.io/crates/v/sqisign-verify.svg)](https://crates.io/crates/sqisign-verify)
[![docs.rs](https://docs.rs/sqisign-verify/badge.svg)](https://docs.rs/sqisign-verify)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

SQIsign signature verification in pure Rust. `no_std`-compatible, zero heap allocation, suitable for embedded targets.

SQIsign is a post-quantum digital signature scheme based on isogenies between supersingular elliptic curves. It was proposed for NIST PQC standardization and produces the smallest signatures of any post-quantum candidate at comparable security levels.

This crate provides only the verification path. It has no dependency on the quaternion algebra stack, `num-bigint`, or any heap allocator. Tested on `thumbv7em-none-eabihf`.

For keygen and signing, use [`sqisign-rs`](https://crates.io/crates/sqisign-rs).

## Usage

Pass raw signature bytes; the format (standard, expanded, or compressed) is auto-detected from the byte length:

```rust
use sqisign_verify::PublicKey;

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::from_bytes(pk_bytes)?;
    pk.verify_bytes(msg, sig_bytes)?;
    Ok(())
}
```

With typed signatures, use the `Verifier` trait:

```rust
use sqisign_verify::{PublicKey, Signature, Verifier};

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::from_bytes(pk_bytes)?;
    let sig = Signature::from_bytes(sig_bytes)?;
    pk.verify(msg, &sig).map_err(|_| sqisign_verify::Error::InvalidSignature)?;
    Ok(())
}
```

Level 1 is the default type parameter. For higher security levels, specify explicitly:

```rust
use sqisign_verify::{PublicKey, Level3};

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::<Level3>::from_bytes(pk_bytes)?;
    pk.verify_bytes(msg, sig_bytes)?;
    Ok(())
}
```

## Signature formats

Three wire formats with different size/speed tradeoffs. Format detection is purely length-based (each format has a unique byte count per level).

| Format | L1 | L3 | L5 | Verify (L1) |
|---|---|---|---|---|
| Compressed | 129 bytes | 196 bytes | 257 bytes | 6.83 ms |
| Standard | 148 bytes | 224 bytes | 292 bytes | 4.65 ms |
| Expanded | 212 bytes | 316 bytes | 420 bytes | 3.82 ms |

Public keys are 65 bytes (L1), 97 bytes (L3), 129 bytes (L5).

## Performance

Intel Xeon @ 2.80 GHz, `--release`, `target-cpu=native`:

| Operation | L1 | vs C reference |
|---|---|---|
| Verify (expanded) | 3.82 ms | 41% faster |
| Verify (standard) | 4.65 ms | 29% faster |
| Verify (compressed) | 6.83 ms | comparable |

Verification outperforms the C reference because LLVM aggressively inlines field arithmetic across crate boundaries.

## Types

- `PublicKey<L>`: Montgomery curve coefficient + torsion hint byte
- `Signature<L>`: standard NIST v2.0 format (2x2 matrix + hints)
- `ExpandedSignature<L>`: pre-evaluated kernel points (fastest verification)
- `CompressedSignature<L>`: 3-of-4 matrix entries, 4th recovered via Weil pairing
- `Scalar<L>`: fixed-width multi-precision integer for matrix entries and challenge

The default type parameter is `Level1`, so `PublicKey` and `PublicKey<Level1>` are equivalent. Use `Level3` or `Level5` for higher security levels.

## `no_std`

This crate is `no_std` by default with zero heap allocation. Add the `std` feature if you need `std::error::Error` impls:

```toml
sqisign-verify = { version = "0.2", features = ["std"] }
```

## References

SQIsign was introduced by De Feo, Kohel, Leroux, Petit, and Wesolowski (ASIACRYPT 2020). The v2.0 construction uses (2,2)-isogenies in the theta model following the SQIsign2D-West approach (Basso, Dartois, De Feo, Leroux, Maino, Pope, Robert, Wesolowski, ASIACRYPT 2024).

## License

Apache-2.0 OR MIT
