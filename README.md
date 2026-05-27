# sqisign-rs

[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

A pure Rust implementation of SQIsign v2.0.

SQIsign is a post-quantum digital signature scheme based on isogenies between supersingular elliptic curves. It was proposed for NIST PQC standardization and produces the smallest signatures of any post-quantum candidate at comparable security levels.

This implementation passes all 300 NIST KAT vectors (100 per level) across Level 1, Level 3, and Level 5. Verification is observed to be 29% faster than the C reference implementation. See performance numbers below for more details.

> **This library has not been audited.** The verification path is designed to be constant-time, but no formal verification or third-party audit has been performed. The signing path is inherently variable-time due to SQIsign's algorithmic structure. See [SECURITY.md](SECURITY.md) for details. Use at your own discretion.

## Quick start

Generate a keypair, sign, and verify:

```rust
use sqisign_rs::{generate, PublicKey, SigningKey, Verifier};

let mut rng = rand::rngs::OsRng;
let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;
```

Compress a signature for minimal wire size (129 bytes at Level 1):

```rust
let compressed = sig.compress();
pk.verify(b"hello world", &compressed)?;
```

Use a higher security level:

```rust
use sqisign_rs::{generate, Level3, Verifier};

// Level5 is available as well:
let (pk, sk) = generate::<Level3>(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;
```

## Signature formats

There are three wire formats available with different size and speed tradeoffs. The verifier determines the format from the byte length.

### Sizes

| Format | L1 | L3 | L5 |
|---|---|---|---|
| Compressed | 129 bytes | 196 bytes | 257 bytes |
| Standard | 148 bytes | 224 bytes | 292 bytes |
| Expanded | 212 bytes | 316 bytes | 420 bytes |

Public keys are 65 bytes (L1), 97 bytes (L3), 129 bytes (L5) across all formats.

### Standard (148 bytes at L1)

Contains the auxiliary curve, a challenge scalar, the full 2x2 basis-change matrix, and hint bytes for canonical basis reconstruction.

| Field | Size (L1) | Description |
|---|---|---|
| e_aux_a | 64 bytes | Auxiliary curve A-coefficient (Fp2) |
| backtracking | 1 byte | Backtracking amount n_bt (0-3) |
| two_resp_length | 1 byte | Dim-1 response length r' (0-8) |
| M[0][0], M[0][1], M[1][0], M[1][1] | 4 x 16 bytes | Basis-change matrix coefficients |
| challenge | 16 bytes | Challenge scalar |
| hint_aux, hint_chall | 2 bytes | Canonical basis selection hints |

### Expanded (212 bytes at L1)

Stores pre-evaluated kernel point x-coordinates instead of the matrix. The verifier skips the biscalar multiplication during verification, trading 64 extra bytes for ~17% faster verification.

| Field | Size (L1) | Description |
|---|---|---|
| e_aux_a | 64 bytes | Auxiliary curve A-coefficient (Fp2) |
| backtracking + flags | 1 byte | Backtracking, kernel_is_q flag, pmq sign hint |
| two_resp_length | 1 byte | Dim-1 response length r' |
| challenge | 16 bytes | Challenge scalar |
| P_chl_x | 64 bytes | Challenge kernel point P x-coordinate (Fp2) |
| Q_chl_x | 64 bytes | Challenge kernel point Q x-coordinate (Fp2) |
| hint_aux, hint_chall | 2 bytes | Canonical basis selection hints |

### Compressed (129 bytes at L1)

Drops one matrix coefficient and recovers it during verification using a Weil pairing determinant formula. Canonical basis hints are recomputed from the curves. The three remaining metadata values (backtracking, response length, and a 2-bit determinant hint) are packed into a single byte.

| Field | Size (L1) | Description |
|---|---|---|
| e_aux_a | 64 bytes | Auxiliary curve A-coefficient (Fp2) |
| packed_meta | 1 byte | [trl:4 \| det_hint:2 \| bt:2] packed LSB-first |
| M[0][0] | 16 bytes | Matrix coefficient (always stored) |
| M[0][1] | 16 bytes | Matrix coefficient (always stored) |
| M[var] | 16 bytes | M[1][0] if M[0][0] is odd, M[1][1] if even |
| challenge | 16 bytes | Challenge scalar |

The dropped coefficient is recovered via `det(M) = dlog(omega_aux^{-1}, omega_f^{-1})` where omega_f and omega_aux are Weil pairings of the canonical bases of E_chall and E_aux. See [COMPRESSION.md](COMPRESSION.md) for the full algorithm.

## Performance

Intel Xeon @ 2.80 GHz, `--release`, `target-cpu=native`:

| Operation | L1 | vs C reference |
|---|---|---|
| Verify (expanded) | 3.82 ms | 41% faster |
| Verify (standard) | 4.65 ms | 29% faster |
| Verify (compressed) | 6.83 ms | comparable |
| Key generation | ~185 ms | ~2.5x slower |
| Signing | ~185 ms | ~2.5x slower |

We observe that verification in this library outperforms the C reference, likely because LLVM aggressively inlines field arithmetic across crate boundaries. Signing and keygen are slower than the C reference because the quaternion algebra layer uses `num-bigint` (heap-allocated arbitrary precision integers) instead of GMP. We explicitly pay this performance tax in order to keep the library pure-rust and avoid an FFI dependency through `rug`. `crypto-bigint` was evaluated as a replacement but would have degraded performance further and requires signed operations which are currently unsupported. Future work targets improvements to the big-integer arithmetic performance in the quaternion crate and elsewhere.

## Memory security

Secret key material is protected by a three-tier system:

1. `SecretKey` implements `ZeroizeOnDrop`. The key is zeroed when it goes out of scope.
2. Intermediate signing values are explicitly zeroed at the end of `protocols_sign`.
3. A `ZeroizingAllocator` clears all freed heap memory. This catches residual copies left by `num-bigint`'s internal reallocations.

The allocator is enabled by default through the `sqisign-core` facade crate with zero measured overhead. If you use a custom allocator (jemalloc, mimalloc, etc.), disable it:

```toml
sqisign-core = { version = "0.1", default-features = false }
```

Disabling the allocator means heap memory freed by `num-bigint` during signing will NOT be zeroed. Ghost copies of secret intermediate values may persist in freed heap pages until the memory is reused. Tier 1 (SecretKey zeroization) and Tier 2 (explicit intermediate zeroization) remain active regardless of allocator choice.

## Building

```
cargo build --release
cargo test --workspace --release
cargo bench -p sqisign-verify
cargo bench -p sqisign-sign
```

For best performance:

```toml
# .cargo/config.toml
[target.x86_64-unknown-linux-gnu]
rustflags = ["-C", "target-cpu=native"]

[profile.release]
lto = "thin"
codegen-units = 1
```

## Workspace

```
params        Security level trait and per-level constants
fp            Fp and Fp2 field arithmetic (Montgomery form)
ec            Elliptic curves, isogenies, pairings
precomp       Precomputed constants for all 3 security levels
theta         (2,2)-isogenies in the theta model
quaternion    Quaternion orders, ideals, lattices
id2iso        Ideal-to-isogeny translation (Deuring correspondence)
keygen        Key generation and SecretKey type
sign          Signing protocol
verify        Verification protocol and signature formats (no_std)
core          Facade crate with re-exports and ZeroizingAllocator
alloc         ZeroizingAllocator implementation
kat           KAT test infrastructure (dev-only)
```

The verification path (verify and its dependencies: ec, fp, theta, precomp, params) compiles for bare-metal `no_std` targets without an allocator. It has zero dependency on the quaternion algebra layer.

## References

SQIsign was introduced by De Feo, Kohel, Leroux, Petit, and Wesolowski (ASIACRYPT 2020). The v2.0 construction uses (2,2)-isogenies in the theta model following the SQIsign2D-West approach (Basso, Dartois, De Feo, Leroux, Maino, Pope, Robert, Wesolowski, ASIACRYPT 2024).

## Self-Signature

This library signs its own source code using SQIsign compressed signatures.
The signature below covers the SHA-256 hash of all Rust source files in the
`crates/` directory, sorted lexicographically.

**Public key** (65 bytes):
```
1a5ec6377182ec36e4ab8d29ca525bbfb8dac4d547e985
7774df20afd726690006ddca551230c4b1b416c53016ec
4008becfda012d4afcfd198dae9f55e30b020c
```

**Source hash** (SHA-256):
```
cf122e19f05cfa13668cbe6808050717567ef0220d1c0fdf83fc98cacd88ca21
```

**Signature** (129 bytes, SQIsign Level 1 compressed):
```
47b243c0c560f3b88c4a0bf2c65afe11019e660deeacbc
26c4d6239067b5af03abb51bc0fe5c48523181b1084490
57cd9cc0492ce2355257de9a228d858bde0108a78eed72
39da7df4b2d9ec60139e9461ece6e6c569978eb57c322d
232d73e89226d01042bcfb14dc55e8d0ac30fe3a9c5545
53c2ec6a63e15079ee7aa4e46b00
```

To verify:
```bash
cargo run --release -p project-sign -- hash
# Compare with the source hash above, then:
cargo run --release -p project-sign -- verify \
  --public-key $(cat PROJECT_KEY) \
  --message <hash> \
  --signature <sig>
```

For comparison, an ML-DSA (Dilithium) signature for the same message would be
2,420 bytes. An SLH-DSA (SPHINCS+) signature would be 7,856 bytes.
This SQIsign signature is **129 bytes**.

## License

Apache-2.0 OR MIT
