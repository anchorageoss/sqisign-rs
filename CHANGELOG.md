# Changelog

## 0.6 (unreleased)

**The version number.** crates.io has `sqisign-rs` and `sqisign-verify` at
0.5.0 (a round-2 release with the `kat-compat` and `sqisign-rk` features),
while the manifest on `main` says 0.4.0 and the release tags stop at
v0.4.31; the local branches assumed 0.5.0 for this work. This release is
the next breaking bump from what is published: the 0.6 line, both crates
together. The manifest carries 0.6.0 and the release workflow mints the
patch number from the commit height on `main` (`scripts/auto-version.sh`),
so the published version is 0.6.x. The workflow tags and makes the
GitHub release but does not publish to crates.io; publishing is
`cargo publish -p sqisign-verify` then `cargo publish -p sqisign-rs` from
the tagged release commit, by a maintainer with the credentials.

SQIsign round 3 is the crate. The round-2 implementation, its formats and
its extensions are removed; there is no compatibility layer and no
conversion of round-2 keys or signatures.

### Added

- SQIsign as submitted to NIST's third round (the-sqisign at tag
  `nist-v3`, commit `6d017708`): primes `3·2^324 − 1`, `27·2^500 − 1`,
  `17·2^664 − 1`; public keys of 83 / 129 / 169 bytes, signatures of
  200 / 306 / 406 bytes, secret keys of 270 / 417 / 549 bytes. The 300
  known-answer vectors of the reference are reproduced byte for byte
  (keys, signed messages), on the assembly kernels and on the portable
  fallback (`crates/kat`).
- Typed API: `Level1` / `Level3` / `Level5` (`sqisign_verify::Level`),
  `PublicKey<L>`, `Signature<L>`, `SigningKey<L>`, `generate`,
  `SigningKey::sign`, `PublicKey::{verify, prepare}`, byte codecs with
  `hybrid_array` lengths, the RustCrypto `signature` traits
  (`Verifier`, `SignatureEncoding`, `RandomizedSigner`), `zeroize` on
  the signing key. The verify crate is `no_std` and heap-free; signing
  uses fixed-precision integers (no `num-bigint`).
- The compressed format for round 3 ([COMPRESSION.md](COMPRESSION.md)):
  176 / 269 / 353 bytes, three matrix entries and four bits, the fourth
  entry recovered from two Weil pairings. `Signature::compress`,
  `CompressedSignature::{from_bytes, to_bytes, decompress}`,
  `PublicKey::verify_compressed`, and the slice API in
  `sqisign_verify::compressed`. Not in the specification. Recovery costs
  13 to 18 % of a verification (2.2 / 5.2 / 11.0 Mcycles).
- Field arithmetic chosen at run time: on x86-64 with BMI2 and ADX,
  inline-assembly Montgomery kernels of the reference's `broadwell`
  schedule (`mulx`/`adcx`/`adox`) with results in registers and an
  in-place API through the curve and theta layers; portable kernels of
  the same limb layout otherwise; the generated radix backends on other
  architectures. `SQISIGN_FORCE_PORTABLE=1` selects the portable kernels
  for testing (`sqisign_verify::fp::{force_portable_arithmetic,
  arithmetic_backend}`).
- `PreparedPublicKey` (the per-key work once), `verify_batch`,
  `verify_bytes_prepared`.
- Benchmarks in `rdtsc` cycles (`cargo bench -p sqisign-rs --bench
  bench`), a one-session table in [BENCH.md](BENCH.md) with the C
  reference's `broadwell` and `ref` builds named.
- The compact (dimension-4) format at level I behind the `compact`
  feature, experimental ([docs/COMPACT_R3.md](docs/COMPACT_R3.md)): our
  instantiation of SQIsignHD (ePrint 2023/436) on `3·2^324 − 1` with
  `e = 174`, `r = 89`; the round-2 dimension-4 layer re-parameterised over
  the prime, without `crypto-bigint`; 142-byte signatures, 84-byte keys;
  verification 122 Mcycles (ten times the two-dimensional verifier),
  signing 538 Mcycles on average (a rejection loop on the response
  degree). `CompactPublicKey`, `CompactSignature`, `CompactSigningKey`,
  `generate_compact`. Self-generated vectors are accepted by the SQIsignHD
  library's unmodified Sage verifier at these parameters and the tampered
  ones rejected (`sqisignhd-harness/round3`). The default build carries
  no dimension-4 code. The 0.4.x compact verifier's final image check
  (the paper's) is present here as it was there.

### Removed

- Everything round 2: keys, signatures, the standard / expanded /
  compressed / compact (dimension-4) formats, the `sqisign-rk`
  rerandomizable keys, the `kat-compat` feature and the signer-side
  canonical matrix, the `num-bigint` quaternion layer, the round-2 KATs.
- `-C target-cpu=native` from `.cargo/config.toml`: it made the
  saturated-limb build slower (level V verification 97.5 against 89.5
  Mcycles without it). Builds are of the default target; the kernels
  dispatch on CPUID.

### Not in this release

- SQIsign-RK on round 3: not planned; it lives in prism-rs on ECLIPSE keys.
- The compact format at levels III and V (derived, not built).

### Migration

0.5.0 (and 0.4) round-2 keys and signatures do not decode: `PublicKey::from_bytes`
and `Signature::from_bytes` take the round-3 lengths (83 / 129 / 169 and
200 / 306 / 406 bytes). Generate round-3 keys with `generate::<Level1>`.
The default level is no longer implied: `PublicKey<Level1>` where 0.5 had
`PublicKey`. `SecretKey` is `SigningKey<L>`; `keypair` is `generate`;
`sqisign_rs::sign::sign` is `SigningKey::sign`. The compressed format is
`Signature::compress()` and `PublicKey::verify_compressed` (or the
`Verifier` impl); verification no longer detects the format from the
length.

## 0.5.0 and 0.4.x

SQIsign round 2 (v2.0): the standard, expanded, compressed and compact
formats, SQIsign-RK. See the git history before `75671c1`.
