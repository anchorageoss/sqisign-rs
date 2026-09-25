# sqisign-rs

[![crates.io](https://img.shields.io/crates/v/sqisign-rs.svg)](https://crates.io/crates/sqisign-rs)
[![docs.rs](https://docs.rs/sqisign-rs/badge.svg)](https://docs.rs/sqisign-rs)
[![KAT](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/kat.yml)
[![Tests](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/anchorageoss/sqisign-rs/actions/workflows/tests.yml)

A pure Rust implementation of SQIsign as submitted to the third round of
NIST's post-quantum signature process (the submission of 2026-09-01,
the-sqisign at tag `nist-v3`): key generation, signing and verification
at levels I, III and V, with 83 / 129 / 169-byte public keys and
200 / 306 / 406-byte signatures, and a compressed format of
176 / 269 / 353 bytes of this crate's own. Round 3 changes the primes
(`3·2^324 − 1`, `27·2^500 − 1`, `17·2^664 − 1`) and the parameters of
the two-dimensional verification; for the security rationale of the new
parameters see Wesolowski, [ePrint 2026/1486](https://eprint.iacr.org/2026/1486).
The verifier is `no_std` and heap-free; the 300 known-answer vectors of
the reference are reproduced byte for byte. Round-2 keys and signatures
are not supported from 0.6 on, and there is no conversion
([CHANGELOG.md](CHANGELOG.md)). MSRV 1.80.

> **Not audited.** Verification is designed to be constant-time; signing
> is variable-time by design. See [SECURITY.md](SECURITY.md).

## Quick start

```rust
use sqisign_rs::{generate, Level1, Level3, PublicKey, Signature, SigningKey, Verifier};

let mut rng = rand::rngs::OsRng;

// generate, sign, verify
let (pk, sk): (PublicKey<Level1>, SigningKey<Level1>) = generate(&mut rng);
let sig = sk.sign(b"hello world", &mut rng)?;
pk.verify(b"hello world", &sig)?;

// the wire: fixed-length arrays, strict decoders
let pk_bytes = pk.to_bytes();                       // 83 bytes
let sig_bytes = sig.to_bytes();                     // 200 bytes
let pk = PublicKey::<Level1>::from_bytes(&pk_bytes)?;
let sig = Signature::<Level1>::from_bytes(&sig_bytes)?;
pk.verify(b"hello world", &sig)?;

// the compressed format, 176 bytes; needs the key to open
pk.verify_compressed(b"hello world", &sig.compress())?;

// levels III and V are a type parameter
let (pk3, sk3) = generate::<Level3>(&mut rng);
pk3.verify(b"hello world", &sk3.sign(b"hello world", &mut rng)?)?;
```

`SigningKey::to_bytes` / `from_bytes` carry the secret key (270 / 417 /
549 bytes). The RustCrypto `signature` traits are implemented
(`Verifier`, `SignatureEncoding`, `RandomizedSigner`), and a
`PreparedPublicKey` (`pk.prepare()`) does the per-key work once for a
verifier that sees many signatures under one key. For a verifier only,
depend on [`sqisign-verify`](https://crates.io/crates/sqisign-verify):
`no_std`, no heap, the same types and a slice API
(`sqisign_verify::sqisign`).

## Which arithmetic a build gets

| target | field arithmetic |
|---|---|
| x86-64 with BMI2 and ADX (Intel Broadwell and later, AMD Zen and later) | inline assembly (`mulx`/`adcx`/`adox`, the reference's `broadwell` schedule), chosen at run time from CPUID |
| x86-64 without them | the same saturated layout's portable kernels, about 2.8x slower at verification |
| other architectures | the generated radix backends (unsaturated limbs), the portable code path |

No compiler flags are needed or wanted: `-C target-cpu=native` made this
build slower. `SQISIGN_FORCE_PORTABLE=1` runs the portable kernels for
testing ([PORT.md](PORT.md)).

## Performance

Benchmarks are of the default build (`cargo bench -p sqisign-rs --bench
bench`, no features, no `RUSTFLAGS`), in `rdtsc` cycles, one session with
the C reference's `broadwell` and `ref` builds on the same machine; the
full table with the C numbers is in [BENCH.md](BENCH.md).

Session of 2026-09-25 (Intel(R) Xeon(R) CPU @ 2.80GHz), millions of cycles:

| level | verify (Rust, cold) | verify (C `broadwell`) | verify (C `ref`) | verify compressed (Rust, prepared key) | sign (Rust) | sign (C `broadwell`) |
|---|---|---|---|---|---|---|
| I | 12.3 | 11.3 | 18.9 | 14.1 | 119.3 | 103.8 |
| III | 30.9 | 28.7 | 53.3 | 35.1 | 498.2 | 341.9 |
| V | 88.5 | 79.3 | 101.8 | 97.2 | 937.3 | 754.3 |

Verification is within 7 to 12 % of the reference's assembly-backed
build and 1.5 to 1.7x faster than its portable `ref` build; a prepared
key saves another 3 %. The compressed format's recovery adds 13 to 18 %
to a verification. Signing is 1.1 to 1.5x the reference's.

The numbers published for 0.4 were taken with `-C target-cpu=native` in
`.cargo/config.toml` and against the C reference's `ref` build; neither
holds here.

## Formats

| format | I | III | V | verifies with |
|---|---|---|---|---|
| specification | 200 B | 306 B | 406 B | any round-3 verifier |
| compressed ([COMPRESSION.md](COMPRESSION.md)) | 176 B | 269 B | 353 B | this crate |

The compressed format drops one entry of the basis-change matrix and
recovers it from two Weil pairings the verifier computes on data it
already has. It is not in the specification.

### Compact (dimension 4), experimental

Behind the `compact` feature, level I only: our instantiation of SQIsignHD
(Dartois, Leroux, Robert, Wesolowski, [ePrint 2023/436](https://eprint.iacr.org/2023/436))
on the round-3 prime `3·2^324 − 1`, with the round-3 key generation,
commitment and challenge and SQIsignHD's response, a degree-`q` isogeny
embedded in a `(2^174)^4`-isogeny the verifier computes in the theta
model. The SQIsignHD authors have published parameters only for their own
primes; ours are derived in [docs/COMPACT_R3.md](docs/COMPACT_R3.md).
Signatures are **142 bytes** and public keys 84 bytes; a compact key is
not a round-3 key (its basis convention is the SQIsignHD library's), and
its signatures verify only with `CompactPublicKey`.

| level I, own run 2026-09-25 | compact | round-3 standard |
|---|---|---|
| signature | 142 B | 200 B |
| verify (Mcycles) | 122 | 12.3 |
| sign (Mcycles, mean) | 538 | 119 |
| keygen (Mcycles) | 43 | 39 |

Verification costs about ten times the two-dimensional verifier; signing
is a rejection loop over response samples until `2^174 − q` is a prime
`≡ 1 mod 4` (about one in 240). The format is **experimental and not for
production**: it has not been verified against an independent
implementation beyond the check that the SQIsignHD library's Sage
verifier, re-parameterised to this prime, accepts our vectors and rejects
tampered ones (`sqisignhd-harness/round3`, workflow `compact-oracle`).

```rust
use sqisign_rs::{generate_compact, CompactPublicKey, CompactSignature, Verifier};
const N: usize = sqisign_rs::precomp::p324_3::IBZ_NLIMBS;
let (pk, sk) = generate_compact::<sqisign_rs::P324_3, N>(sqisign_rs::sqisign::level1(), &mut rng)?;
let sig = sk.sign(b"hello world", &mut rng)?;      // 142 bytes
pk.verify(b"hello world", &sig)?;
```

SQIsign is not, and cannot be, SUF-CMA: a valid signature has other
byte strings the verifier accepts, in the specification's format and in
this one, and this crate does not try to change that. See "On the
Auxiliary Isogeny Freedom in SQIsign's Two-Dimensional Representation",
[ePrint 2026/1305](https://eprint.iacr.org/2026/1305), and the round-3
specification's own statement in its Section 8.4.2.

## Not in this release

SQIsign-RK (rerandomizable keys) is not planned for this crate; it lives
in prism-rs on ECLIPSE keys. The compact format at levels III and V is
derived in [docs/COMPACT_R3.md](docs/COMPACT_R3.md) but not built.

## About

SQIsign is advancing through NIST's post-quantum signature
standardization and has **not** been standardized. It is the only
isogeny-based candidate and the one with the smallest signatures. The
2022 attacks that broke the SIDH/SIKE key exchange do not apply to it;
the higher-dimensional isogeny techniques those attacks introduced are
now used constructively to build and speed up SQIsign, including the
dimension-2 verification implemented here. The scheme and its
implementations are young and moving fast; this library tracks that
progress.

## License

Apache-2.0 OR MIT
