# sqisign-rs

[![crates.io](https://img.shields.io/crates/v/sqisign-rs.svg)](https://crates.io/crates/sqisign-rs)
[![docs.rs](https://docs.rs/sqisign-rs/badge.svg)](https://docs.rs/sqisign-rs)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

A pure Rust implementation of SQIsign v2.0.

SQIsign is a post-quantum digital signature scheme based on isogenies between supersingular elliptic curves. It was proposed for NIST PQC standardization and produces the smallest signatures of any post-quantum candidate at comparable security levels.

This implementation passes all 300 NIST KAT vectors (100 per level) across Level 1, Level 3, and Level 5. Verification is observed to be 29% faster than the C reference implementation. See performance numbers below for more details.

> **This library has not been audited.** The verification path is designed to be constant-time, but no formal verification or third-party audit has been performed. The signing path is inherently variable-time due to SQIsign's algorithmic structure. See [SECURITY.md](SECURITY.md) for details. Use at your own discretion.

For verify-only usage (`no_std`), depend on [`sqisign-verify`](https://crates.io/crates/sqisign-verify) directly. Its dim-2 verification path is heap-free; the compact (dim-4) verifier needs only `alloc`.

## Quick start

Generate a keypair, sign, and verify - choose the scheme at keygen time:

```rust
use sqisign_rs::{generate, generate_compact, AnySignature, PublicKey, SigningKey, Verifier};

let mut rng = rand::rngs::OsRng;

// Standard (dimension-2, 148 bytes at Level 1):
let (pk, sk): (PublicKey, SigningKey) = generate(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;

// Compact - the smallest post-quantum signature (108 bytes at Level 1):
let (cpk, csk) = generate_compact(&mut rng);
let csig = csk.sign(b"hello world", &mut rng)?;
cpk.verify(b"hello world", &csig)?;
```

Verification autodetects the format from its byte length. The public key must
match the scheme: standard keys verify standard/compressed/expanded signatures,
compact keys verify compact signatures (they are not interchangeable).

```rust
let any = AnySignature::from_bytes(&sig_bytes)?;
pk.verify(b"hello world", &any)?;       // 129/148/212-byte dim-2 formats
// cpk.verify(b"hello world", &any)?;   // 108-byte compact format
```

Compress a standard signature for a smaller wire size (129 bytes at Level 1):

```rust
let compressed = sig.compress();
pk.verify(b"hello world", &compressed)?;
```

Use a higher security level (dim-2; the compact scheme is Level 1):

```rust
use sqisign_rs::{generate, Level3, Verifier};

// Level5 is available as well:
let (pk, sk) = generate::<Level3>(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;
```

## Signature formats

Four wire formats are available across two dimensions, with different size and
speed tradeoffs. The verifier determines both the format and the dimension from
the byte length alone (no tag byte) - `AnySignature::from_bytes` autodetects,
and the same `pk.verify(msg, &any)` call handles all of them.

### Sizes

| Format | L1 | L3 | L5 |
|---|---|---|---|
| Compact (dim-4) | 108 bytes | 161 bytes | 213 bytes |
| Compressed | 129 bytes | 196 bytes | 257 bytes |
| Standard | 148 bytes | 224 bytes | 292 bytes |
| Expanded | 212 bytes | 316 bytes | 420 bytes |

The compact (dimension-4) format produces the smallest signatures (108 bytes at
Level 1, autodetected and verified through the same API) and is implemented at
Level 1 today; the dim-2 formats are available at Levels 1/3/5. Public keys are
65 bytes (L1 dim-2), 97 bytes (L3), 129 bytes (L5); a compact public key is 64
bytes (L1).

### Two distinct key schemes

Standard and compact keys are **separate, non-interchangeable** cryptographic
objects (different torsion-basis conventions and precomputed data). You choose
the scheme at keygen time - `generate` for standard keys, `generate_compact` for
compact keys - and verify with the matching public-key type:

- a `PublicKey` verifies `Signature` / `CompressedSignature` / `ExpandedSignature`;
- a `CompactPublicKey` verifies `CompactSignature`.

`AnySignature::from_bytes` autodetects the wire format by length, but the public
key must belong to the right scheme; a mismatched key is rejected.

### Compact (108 bytes at L1)

The smallest format. The signer transmits the commitment curve and the response
isogeny's action as three rescaled scalars plus the response degree; the verifier
reconstructs and checks the response via a dimension-4 isogeny (the Kani
embedding). Implemented at Level 1.

| Field | Size (L1) | Description |
|---|---|---|
| A_com | 64 bytes | Commitment curve A-coefficient (𝔽p²); basis hints packed into spare bits |
| q | 17 bytes | Response degree (q < 2¹³⁶) |
| a, b, c_or_d | 3 × 9 bytes | Response scalars (the 4th is determinant-recovered) |

The challenge is not transmitted - the verifier recomputes it as a Fiat-Shamir
hash of the curves and message.

### Standard (148 bytes at L1)

Contains the auxiliary curve, a challenge scalar, the full 2×2 basis-change matrix, and hint bytes for canonical basis reconstruction.

| Field | Size (L1) | Description |
|---|---|---|
| E_aux(A) | 64 bytes | Auxiliary curve A-coefficient (𝔽p²) |
| backtracking | 1 byte | Backtracking amount n_bt (0-3) |
| two_resp_length | 1 byte | Dim-1 response length r′ (0-8) |
| M₀₀, M₀₁, M₁₀, M₁₁ | 4 × 16 bytes | Basis-change matrix coefficients |
| challenge | 16 bytes | Challenge scalar |
| hint_aux, hint_chall | 2 bytes | Canonical basis selection hints |

### Expanded (212 bytes at L1)

Stores pre-evaluated kernel point x-coordinates instead of the matrix. The verifier skips the bi-scalar multiplication during verification, trading 64 extra bytes for ~17% faster verification.

| Field | Size (L1) | Description |
|---|---|---|
| E_aux(A) | 64 bytes | Auxiliary curve A-coefficient (𝔽p²) |
| backtracking + flags | 1 byte | Backtracking, kernel_is_Q flag, P-Q sign hint |
| two_resp_length | 1 byte | Dim-1 response length r′ |
| challenge | 16 bytes | Challenge scalar |
| P_chl(x) | 64 bytes | Challenge kernel point P x-coordinate (𝔽p²) |
| Q_chl(x) | 64 bytes | Challenge kernel point Q x-coordinate (𝔽p²) |
| hint_aux, hint_chall | 2 bytes | Canonical basis selection hints |

### Compressed (129 bytes at L1)

Drops one matrix coefficient and recovers it during verification using a Weil pairing determinant formula. Canonical basis hints are recomputed from the curves. The three remaining metadata values (backtracking, response length, and a 2-bit determinant hint) are packed into a single byte.

| Field | Size (L1) | Description |
|---|---|---|
| E_aux(A) | 64 bytes | Auxiliary curve A-coefficient (𝔽p²) |
| packed_meta | 1 byte | [trl:4 \| det_hint:2 \| bt:2] packed LSB-first |
| M₀₀ | 16 bytes | Matrix coefficient (always stored) |
| M₀₁ | 16 bytes | Matrix coefficient (always stored) |
| M_var | 16 bytes | M₁₀ if M₀₀ is odd, M₁₁ if even |
| challenge | 16 bytes | Challenge scalar |

The dropped coefficient is recovered via det(M) = dlog(ω_aux⁻¹, ω_f⁻¹) where ω_f and ω_aux are Weil pairings of the canonical bases of E_chall and E_aux. See [COMPRESSION.md](COMPRESSION.md) for the full algorithm.

## Performance

Dimension-2 operations, Intel Xeon @ 2.80 GHz, `--release`, `target-cpu=native`:

| Operation | L1 | vs C reference |
|---|---|---|
| Verify (expanded) | 3.82 ms | 41% faster |
| Verify (standard) | 4.65 ms | 29% faster |
| Verify (compressed) | 6.83 ms | comparable |
| Key generation | ~185 ms | ~2.5x slower |
| Signing | ~185 ms | ~2.5x slower |

Compact (dimension-4, Level 1), measured on the development machine:

| Operation | Compact |
|---|---|
| Key generation | ~47 ms |
| Signing | ~51 ms |
| Verify (serial) | ~33 ms |
| Verify (`parallel` feature) | ~20.5 ms |

The two tables use different reference machines, so compare within a table, not
across. The compact `parallel` feature runs the two independent dim-4
half-chains on separate threads; the result is bit-identical to the serial path.

We observe that dim-2 verification in this library outperforms the C reference, likely because LLVM aggressively inlines field arithmetic across crate boundaries. Signing and keygen are slower than the C reference because the quaternion algebra layer uses `num-bigint` (heap-allocated arbitrary precision integers) instead of GMP. We explicitly pay this performance tax in order to keep the library pure-rust and avoid an FFI dependency through `rug`. `crypto-bigint` was evaluated as a replacement but would have degraded performance further and requires signed operations which are currently unsupported. Future work targets improvements to the big-integer arithmetic performance in the quaternion crate and elsewhere.

## Memory security

Secret key material is protected by a three-tier system:

1. `SecretKey` implements `ZeroizeOnDrop`. The key is zeroed when it goes out of scope.
2. Intermediate signing values are explicitly zeroed at the end of `protocols_sign`.
3. A `ZeroizingAllocator` clears all freed heap memory. This catches residual copies left by `num-bigint`'s internal reallocations.

The allocator is enabled by default via the `zeroize-alloc` feature of `sqisign-rs` with zero measured overhead. If you use a custom allocator (jemalloc, mimalloc, etc.), disable it:

```toml
sqisign-rs = { version = "0.4", default-features = false }
```

Disabling the allocator means heap memory freed by `num-bigint` during signing will NOT be zeroed. Ghost copies of secret intermediate values may persist in freed heap pages until the memory is reused. Tier 1 (SecretKey zeroization) and Tier 2 (explicit intermediate zeroization) remain active regardless of allocator choice.

## Building

```
cargo build --release
cargo test --workspace --release
cargo bench -p sqisign-verify
cargo bench -p sqisign-rs
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

Two published crates:

- **`sqisign-verify`** - verification only (`no_std`, `#![forbid(unsafe_code)]`):
  field arithmetic (`fp`, `params`), elliptic curves and isogenies (`ec`),
  precomputed constants (`precomp`), the dim-2 (2,2)-theta model (`theta`), the
  dimension-4 verifier (`hd`), and the verification protocol plus signature
  formats (`verify`, `formats`, `compact`).
- **`sqisign-rs`** - key generation and signing, re-exporting everything from
  `sqisign-verify`: the quaternion layer (`quaternion`), ideal-to-isogeny
  translation (`id2iso`, Deuring correspondence), key generation (`keygen`), the
  signing protocols (`sign`, including the compact signer), and the optional
  `ZeroizingAllocator` (`alloc`).

The verification path (`sqisign-verify` and its dependencies) has zero dependency
on the quaternion algebra layer. The dim-2 path compiles for bare-metal `no_std`
without an allocator; the dimension-4 (compact) verifier additionally needs
`alloc`.

## References

SQIsign was introduced by De Feo, Kohel, Leroux, Petit, and Wesolowski (ASIACRYPT 2020). The v2.0 construction uses (2,2)-isogenies in the theta model following the SQIsign2D-West approach (Basso, Dartois, De Feo, Leroux, Maino, Pope, Robert, Wesolowski, ASIACRYPT 2024).

## License

Apache-2.0 OR MIT
