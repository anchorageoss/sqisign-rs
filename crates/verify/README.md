# sqisign-verify

[![crates.io](https://img.shields.io/crates/v/sqisign-verify.svg)](https://crates.io/crates/sqisign-verify)
[![docs.rs](https://docs.rs/sqisign-verify/badge.svg)](https://docs.rs/sqisign-verify)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

SQIsign signature verification in pure Rust. `no_std`-compatible and `#![forbid(unsafe_code)]`, suitable for embedded targets.

SQIsign is a post-quantum digital signature scheme based on isogenies between supersingular elliptic curves. It was proposed for NIST PQC standardization and produces the smallest signatures of any post-quantum candidate at comparable security levels.

This crate provides only the verification path, with no dependency on the quaternion algebra stack or `num-bigint`. It verifies both the **dimension-2** formats (standard / expanded / compressed, at Levels 1/3/5) and the **compact** 108-byte format (dimension-4, Level 1). The dim-2 verification path is heap-free and tested on `thumbv7em-none-eabihf`; the compact verifier keeps a small, bounded heap allocation (the optimal-strategy chain loop) off the constant-time path, so it requires `alloc`.

All formats - including the compact 108-byte one - are auto-detected from byte length by `AnySignature::from_bytes` and verified through the same RustCrypto `Verifier` trait: **108 = compact, 129 = compressed, 148 = standard, 212 = expanded** (Level 1), a collision-free mapping.

For keygen and signing, use [`sqisign-rs`](https://crates.io/crates/sqisign-rs).

## Usage

All verification goes through `pk.verify(msg, &sig)` via the RustCrypto `Verifier` trait:

```rust
use sqisign_verify::{PublicKey, Signature, Verifier};

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::from_bytes(pk_bytes)?;
    let sig = Signature::from_bytes(sig_bytes)?;
    pk.verify(msg, &sig)?;
    Ok(())
}
```

For raw bytes where the format is unknown, parse into `AnySignature` first (auto-detects from byte length):

```rust
use sqisign_verify::{formats::AnySignature, PublicKey, Verifier};

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::from_bytes(pk_bytes)?;
    let sig = AnySignature::from_bytes(sig_bytes)?;
    pk.verify(msg, &sig)?;
    Ok(())
}
```

Compact (108-byte) signatures verify with a `CompactPublicKey` (a distinct
scheme from the dim-2 keys, Level 1 only):

```rust
use sqisign_verify::{CompactPublicKey, CompactSignature, Verifier};

fn verify_compact(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = CompactPublicKey::from_bytes(pk_bytes)?;   // 64 bytes
    let sig = CompactSignature::from_bytes(sig_bytes)?; // 108 bytes
    pk.verify(msg, &sig)?;
    Ok(())
}
```

Level 1 is the default type parameter. For higher security levels, specify explicitly:

```rust
use sqisign_verify::{PublicKey, Signature, Level3, Verifier};

fn verify(pk_bytes: &[u8], sig_bytes: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::<Level3>::from_bytes(pk_bytes)?;
    let sig = Signature::<Level3>::from_bytes(sig_bytes)?;
    pk.verify(msg, &sig)?;
    Ok(())
}
```

## Signature formats

Four wire formats, with different size/speed tradeoffs. Format **and dimension**
are detected purely by length - each combination has a unique byte count, so
`AnySignature::from_bytes` selects the right one with no tag byte. (Compact is
implemented at Level 1; its L3/L5 sizes are listed for the collision-free table.)

| Format | L1 | L3 | L5 | Relative verify cost |
|---|---|---|---|---|
| Compact (dim-4) | 108 bytes | 161 bytes | 213 bytes | dim-4 (see Performance) |
| Compressed | 129 bytes | 196 bytes | 257 bytes | slowest dim-2 |
| Standard | 148 bytes | 224 bytes | 292 bytes | baseline dim-2 |
| Expanded | 212 bytes | 316 bytes | 420 bytes | fastest dim-2 |

Compact signatures are the smallest (108 bytes at L1). They verify through the
same `Verifier` trait, but with a **compact public key** (`CompactPublicKey`):
the compact and dim-2 schemes use different torsion-basis conventions and their
keys are not interchangeable. `AnySignature` autodetects the format, but the
public key must match the scheme - a dim-2 `PublicKey` rejects a compact
signature and vice versa. The dim-4 verifier itself lives in the `hd` module.

Public keys are 65 bytes (L1 dim-2), 97 bytes (L3), 129 bytes (L5); a compact
public key is 64 bytes (L1).

## Performance

Dim-2 (standard) verification, Apple M4 Pro, single thread, portable build,
maximal (fat) LTO:

| Level | Rust | C reference |
|---|---:|---:|
| L1 | 1.37 ms | 1.53 ms |
| L3 | 4.13 ms | 4.13 ms |
| L5 | 8.06 ms | 8.66 ms |

Verification is at parity with the C reference on portable ARM builds with
maximal LTO; differences under ~4% are within thermal noise. (An earlier
faster-than-C result came from comparing against a thin-LTO C build.)

Compact (dim-4, Level 1) verification is **~33 ms serial**, or **~20.5 ms** with
the `parallel` feature (the two independent half-chains run on separate threads;
the result is bit-identical to the serial path). These compact numbers are from
the development environment, not the M4 Pro above, and are unoptimized portable
Rust.

## Features

- `std` - adds `std::error::Error` impls (off by default).
- `parallel` - runs the two dim-4 half-chains on separate threads (`std`-only,
  off by default). Affects only compact-verification latency, not the total work
  or any result; the serial `no_std` path is the canonical implementation.
- `vartime` - routes the dominant dim-4 field inversion through a variable-time
  binary-GCD path. Verification is on public data, so this is sound. Off by
  default; currently no measurable speedup (the profile is multiply-bound, not
  inversion-bound), kept for future backends.

## Types

- `PublicKey<L>`: Montgomery curve coefficient + torsion hint byte (dim-2)
- `Signature<L>`: standard NIST v2.0 format (2×2 matrix + hints)
- `ExpandedSignature<L>`: pre-evaluated kernel points (fastest verification)
- `CompressedSignature<L>`: 3-of-4 matrix entries, 4th recovered via Weil pairing
- `CompactPublicKey<L>` / `CompactSignature<L>`: the compact 108-byte scheme (Level 1)
- `AnySignature<L>`: any of the above, auto-detected from byte length
- `Scalar<L>`: fixed-width multi-precision integer for matrix entries and challenge

The default type parameter is `Level1`, so `PublicKey` and `PublicKey<Level1>` are equivalent. Use `Level3` or `Level5` for higher security levels (dim-2 only; the compact scheme is Level 1).

## `no_std`

This crate is `no_std` by default. The dim-2 verification path is heap-free
(tested on `thumbv7em-none-eabihf`); the compact (dim-4) verifier keeps a small,
bounded heap allocation off the constant-time path, so it requires `alloc`. Add
the `std` feature if you need `std::error::Error` impls:

```toml
sqisign-verify = { version = "0.4", features = ["std"] }
```

## References

SQIsign was introduced by De Feo, Kohel, Leroux, Petit, and Wesolowski (ASIACRYPT 2020). The v2.0 construction uses (2,2)-isogenies in the theta model following the SQIsign2D-West approach (Basso, Dartois, De Feo, Leroux, Maino, Pope, Robert, Wesolowski, ASIACRYPT 2024).

## License

Apache-2.0 OR MIT
