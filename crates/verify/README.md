# sqisign-verify

The verification half of [sqisign-rs](https://crates.io/crates/sqisign-rs):
SQIsign as submitted to the third round of NIST's post-quantum signature
process (the-sqisign at tag `nist-v3`), levels I, III and V, in pure Rust.
`no_std`, heap-free, no quaternion algebra and no big integers: field and
curve arithmetic, the theta-model isogeny chain, the protocol's
verification, the byte codecs and a compressed signature format of this
crate's own. The 300 known-answer vectors of the reference are reproduced
byte for byte by the workspace's KAT crate.

Round 2 is not supported from 0.6 on; round-2 keys and signatures do not
decode and there is no conversion.

## Usage

```rust
use sqisign_verify::{Level1, PublicKey, Signature, Verifier};

fn verify(pk: &[u8], sig: &[u8], msg: &[u8]) -> Result<(), sqisign_verify::Error> {
    let pk = PublicKey::<Level1>::from_bytes(pk)?;        // 83 bytes
    let sig = Signature::<Level1>::from_bytes(sig)?;      // 200 bytes
    pk.verify(msg, &sig)                                  // or Verifier::verify(&pk, msg, &sig)
}
```

`Level3` and `Level5` are the other two type parameters (129 / 306 and
169 / 406 bytes). A `PreparedPublicKey` (`pk.prepare()`) does the per-key
work once for a verifier that sees many signatures under one key, and
`verify_batch` takes a slice of `(key, message, signature)` items. The
compressed format (`CompressedSignature`, 176 / 269 / 353 bytes) verifies
with `pk.verify_compressed` or the same `Verifier` trait; it is not in the
specification and needs this crate to open.

The slice API of `sqisign_verify::sqisign` (`verify`, `verify_prepared`,
`public_key_from_bytes`, `signature_from_bytes`, ...) and of
`sqisign_verify::compressed` takes `&VerifyParams` (`level1()`,
`level3()`, `level5()`) and works without `alloc`.

## Formats

| format | I | III | V | in the specification |
|---|---|---|---|---|
| public key | 83 | 129 | 169 | yes |
| signature | 200 | 306 | 406 | yes |
| compressed signature | 176 | 269 | 353 | no (this crate) |

Decoders are strict about lengths and field elements; the verifier applies
the specification's range and sign check on the basis-change matrix.
SQIsign is not strongly unforgeable, and this crate does not try to make
its formats canonical beyond that check (see the workspace's
`SECURITY.md`).

## Arithmetic

On x86-64 with BMI2 and ADX the field arithmetic is inline assembly of the
reference's `broadwell` schedule, chosen at run time from CPUID; without
them, portable kernels of the same layout; on other architectures, the
generated radix backends. No compiler flags are needed;
`SQISIGN_FORCE_PORTABLE=1` selects the portable kernels for testing
(`fp::force_portable_arithmetic`).

## Performance

Session of 2026-09-25 (Intel Xeon 2.8 GHz, `rdtsc` cycles, criterion
means; the C reference is the-sqisign `nist-v3` built `broadwell`, its
medians; the workspace's `BENCH.md` has the whole table):

| level | verify from bytes | prepared key | compressed, prepared | C `broadwell` | C `ref` |
|---|---|---|---|---|---|
| I | 12.3 M | 11.9 M | 14.1 M | 11.3 M | 18.9 M |
| III | 30.9 M | 29.9 M | 35.1 M | 28.7 M | 53.3 M |
| V | 88.5 M | 86.2 M | 97.2 M | 79.3 M | 101.8 M |

## Features

- `std` (default): `std::error::Error` for `Error`. Off, the crate is
  `no_std` and heap-free; it is checked on `thumbv7em-none-eabihf`.

```toml
sqisign-verify = { version = "0.6", default-features = false }
```

## Types

- `PublicKey<L>`, `Signature<L>`, `CompressedSignature<L>` with fixed-length
  byte arrays (`hybrid_array`), `TryFrom<&[u8]>`, `SignatureEncoding`.
- `PreparedPublicKey<L>`, `BatchItem`, `verify_batch`.
- `Level`: `Level1` (`P324_3`), `Level3` (`P500_27`), `Level5` (`P664_17`).
- `Error`: `InvalidSignature`, `MalformedInput`, `InvalidLength`, `InternalError`.

## References

- SQIsign, round-3 submission (2026-09-01), the-sqisign at `nist-v3`.
- Wesolowski, ePrint 2026/1486, the security rationale of the round-3 parameters.
- "On the Auxiliary Isogeny Freedom in SQIsign's Two-Dimensional
  Representation", ePrint 2026/1305, and the specification's Section 8.4.2,
  on unforgeability.

## License

Apache-2.0 OR MIT
