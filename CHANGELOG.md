# Changelog

## 0.7.0 (unreleased)

The compressed formats changed; nothing else did. A signature compressed
by any earlier version must be recompressed from the standard format. No
verification result of a signature that did decompress was affected.

### Changed

- Round 3: the compressed format is one bit string, each field at its
  own width: the two coordinates of `A_aux` at `ceil(log2 p)` bits (326 /
  505 / 669), the low `RESPONSE_BITS + 1` bits of three matrix entries
  (197 / 303 / 401 each), three top bits (bit `RESPONSE_BITS + 1` of the
  dropped entry and of the two transmitted entries other than the first
  entry in wire order with a non-zero low part, whose top bit the
  signer's normalisation fixes at zero), the challenge at `lambda` bits,
  the two hints at 8 bits each, zero-padded to a byte: 174 / 267 / 352
  bytes (`compressed_bits`, `compressed_bytes`; `entry_bits` replaces
  `entry_bytes`, and `VerifyParams` gains `fp_bits`). The fourth entry
  is recovered modulo `2^(RESPONSE_BITS + 1)` from the determinant, which
  the two Weil pairings give to that modulus, and its top bit from the
  string, so decompression returns the signer's signature exactly;
  `Signature::compress` returns `Err` for a matrix the verifier's checks
  reject. The decoder rejects a non-zero padding bit. Signatures
  compressed by `sqisign-verify` and `sqisign-rs` 0.6.0 and 0.6.1 on
  crates.io, and by the GitHub releases v0.6.34 to v0.6.38 (176 / 269 /
  353 bytes, four bits in a byte), do not decode.
- Round 2, on the `fix/compressed-nbt-round2` branch from tag v0.4.31
  (the round-2 code is not in this tree): the metadata byte is
  `[r:4 | h:1 | n_bt:3]`, the hint is one bit (the top bit of the dropped
  entry; the pairings give the determinant modulo `2^(N - 1)`), and
  `compress()` rejects `n_bt >= 8` and `r >= 16` instead of truncating.
  Every published round-2 version, `sqisign-verify` and `sqisign-rs`
  0.1.0 to 0.5.0, silently truncated `n_bt` to two bits and could not
  decompress signatures with `n_bt >= 4`: one of the 300 round-2 KAT
  vectors, level III entry 97. Those versions are deprecated.

### Migration

Recompress from the standard format; the standard and expanded formats
and every key encoding are unchanged.

## 0.6.0 and 0.6.1

**The version number.** 0.6.0 and 0.6.1 are the round-3 release on
crates.io, both crates together. The release workflow mints its tags from
the commit height on `main` (`scripts/auto-version.sh`), so the GitHub
releases of the same line are v0.6.34 to v0.6.38. The workflow tags and
makes the GitHub release but does not publish to crates.io; publishing is
`cargo publish -p sqisign-verify` then `cargo publish -p sqisign-rs`, by
a maintainer with the credentials.

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
  (`Verifier`, `SignatureEncoding`, `RandomizedSigner`). The signing
  key, the protocol-level secret keys and the secret-key encoding
  (`Zeroizing<Vec<u8>>`) zeroize on drop, and the signing path's secret
  intermediates are held in `Zeroizing` wrappers (SECURITY.md lists what
  is and is not covered). The verify crate is `no_std` and heap-free;
  signing uses fixed-precision integers (no `num-bigint`).
- `paper-dim4/`: a short paper on the dimension-4 format at the round-3
  parameters (one `main.tex`, `refs.bib`, figures and tables generated
  from `paper-dim4/BENCH.md` and `docs/COMPACT_R3.md`), with its
  one-session benchmark harness `paper-dim4/bench` and the session it
  reads in `paper-dim4/BENCH.md`.
  `compact_sign_stats` returns the signer's rejection-loop counts next to
  the signature length.
- The compressed format for round 3: 176 / 269 / 353 bytes, three matrix
  entries and four bits, the fourth entry recovered from two Weil
  pairings (changed in 0.7.0 to a bit-packed 174 / 267 / 352 bytes). `Signature::compress`,
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
- `ZeroizingAllocator`, the `enable_secure_allocator!` macro and the
  `crates/alloc` audit crate: they scrubbed the round-2 signer's heap
  big-integer temporaries, which the round-3 signer does not have (its
  integers are fixed-precision and live on the stack or in the key).
- `-C target-cpu=native` from `.cargo/config.toml`: it made the
  saturated-limb build slower (level V verification 97.5 against 89.5
  Mcycles without it). Builds are of the default target; the kernels
  dispatch on CPUID.

### Not in this release

- SQIsign-RK on round 3: not planned; it lives in eclipse-rs on ECLIPSE keys.
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
