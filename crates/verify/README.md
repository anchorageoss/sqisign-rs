# sqisign-verify

SQIsign signature verification. `no_std`-compatible, heap-free, and independent of the quaternion algebra stack.

## Usage

```rust
use sqisign_verify::{PublicKey, Signature};

let pk = PublicKey::from_bytes(&pk_bytes).expect("valid public key");
let sig = Signature::from_bytes(&sig_bytes).expect("valid signature");
sig.verify(&pk, msg)?;
```

Auto-detect format from wire length:

```rust
use sqisign_verify::AnySignature;

let sig = AnySignature::from_bytes(&wire_bytes).expect("valid signature");
sig.verify(&pk, msg)?;
```

## Supported formats

| Format | L1 | L3 | L5 | Verify (L1) |
|---|---|---|---|---|
| Compressed | 129 B | 196 B | 257 B | 6.83 ms |
| Standard | 148 B | 224 B | 292 B | 4.65 ms |
| Expanded | 212 B | 316 B | 420 B | 3.82 ms |

Public keys are 65 B (L1), 97 B (L3), 129 B (L5).

Format detection is purely length-based (each format has a unique byte count per level).

## Key types

- `PublicKey<L>`: Montgomery curve coefficient + torsion hint byte
- `Signature<L>`: standard NIST v2.0 format (2x2 matrix + hints)
- `ExpandedSignature<L>`: pre-evaluated kernel points (fastest verify)
- `CompressedSignature<L>`: 3-of-4 matrix entries, 4th recovered via Weil pairing
- `AnySignature<L>`: auto-detecting wrapper over all three formats
- `Scalar<L>`: fixed-width multi-precision integer for matrix entries and challenge

## Design constraints

This crate deliberately excludes `sqisign-quaternion` and `num-bigint` from its dependency tree. Verification uses only elliptic curve and field arithmetic, suitable for embedded targets (`thumbv7em-none-eabihf` tested).

For the full API with keygen and signing, use `sqisign-core`.
