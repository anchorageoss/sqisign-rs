# The round-3 port

sqisign-rs 0.6 implements SQIsign as submitted to NIST's third round
(the-sqisign at tag `nist-v3`, commit `6d017708`), ported from the
prism-rs substrate that the ECLIPSE work built and measured against the
same reference. This file records what the modules are, what gates them,
and what is not here.

## Layout

| module | what | origin in prism-rs |
|---|---|---|
| `sqisign_verify::params` | the three primes as types (`P324_3`, `P500_27`, `P664_17`), the level constants | `prism-verify/src/params` |
| `sqisign_verify::fp` | `Fp`, `Fp2`; per-prime backends: on x86-64 saturated 64-bit limbs with inline-assembly Montgomery multiplication, fused `Fp2` product, addition and subtraction (`mulx`/`adcx`/`adox`, the reference's `broadwell` schedule), results in registers and an in-place API, when CPUID reports BMI2 and ADX; portable kernels of the same layout otherwise (`fp::dispatch`, one cached atomic); elsewhere the generated unsaturated radix backends | `prism-verify/src/fp`, `tools/gen-fp-asm`, `tools/gen-fp` |
| `sqisign_verify::ec` | Kummer-line arithmetic, the deterministic torsion basis with a bounded entangled-basis search that fails closed on crafted curves, Weil and Tate pairings, the two-adic discrete logarithms | `prism-verify/src/ec` |
| `sqisign_verify::theta` | `(2^n, 2^n)`-isogeny chains in the theta model, round-3 formulas, gluing and splitting, in-place evaluation | `prism-verify/src/theta` |
| `sqisign_verify::precomp` | the deterministic basis of `E0[2^f]` and the small sizes per prime | `prism-verify/src/precomp`, `tools/gen-precomp` |
| `sqisign_verify::sqisign` | `verify`, `PreparedPublicKey` / `verify_prepared`, `verify_batch`, `hash_to_challenge`, the public-key and signature encodings on slices, `VerifyParams` (`level1/3/5`) | `prism-substrate-check/src/{protocol,encode,params}.rs`, made `no_std` |
| `sqisign_verify::compressed` | the compressed format ([COMPRESSION.md](COMPRESSION.md)) | this repository (the method of prism-rs's ECLIPSE encoding) |
| `sqisign_verify::types` | `Level`, `PublicKey<L>`, `Signature<L>`, `CompressedSignature<L>`, the `signature` traits | this repository |
| `sqisign_rs::mp` | fixed-precision integers (`Ibz<N>`), the constant-time layer, the SHAKE-based `Rng` | `prism/src/mp` |
| `sqisign_rs::quat` | the quaternion algebra in inert representation, ideals, constant-time lattice reduction, Qlapoty | `prism/src/quat` |
| `sqisign_rs::id2iso` | ideal-to-isogeny translation with basis images, endomorphism application | `prism/src/id2iso` |
| `sqisign_rs::precomp` | the action of `O0` on the basis of `E0[2^f]`, quaternion-layer sizes | `prism/src/precomp` |
| `sqisign_rs::sqisign` | `keygen`, `sign`, `Params` (`level1/3/5`), the secret-key encoding | `prism-substrate-check/src/protocol.rs` |
| `sqisign_rs` (root) | `SigningLevel`, `SigningKey<L>`, `generate`, the `RandomizedSigner` impl | this repository |

## Gates

- `crates/kat/tests/kat.rs`: the 300 known-answer vectors of the
  reference (`crates/kat/kat/*.rsp`), byte for byte: the DRBG seeded with
  the entry's seed, key generation's keys, signing's `sm = sig || msg`
  from the re-decoded secret key, cold and prepared verification, a
  modified message rejected. CI runs it on the detected kernels and under
  `SQISIGN_FORCE_PORTABLE=1` (the portable kernels). Debug builds run two
  entries per level unless `SQISIGN_KAT_ENTRIES` says otherwise.
- `crates/verify/tests/fp_backends_differential.rs`: every field
  operation against big integers on random inputs, both kernel paths, and
  the two paths against each other limb for limb.
- `crates/verify/tests/{fp,ec,theta}_props.rs`,
  `crates/sqisign-rs/tests/{mp,quat,id2iso}_props.rs`, `theta_chain.rs`:
  the property tests of each layer with the reference's vectors.
- `crates/sqisign-rs/tests/api_roundtrip.rs`: the typed API end to end,
  the codecs, the traits, the rejections; `compressed.rs`: the compressed
  format's round trip at each level, the decoder, tampering, and the
  acceptance survey on generated signatures.
- Generators: `python3 tools/gen-fp-asm/gen_fp_asm.py --check
  --without-mont crates/verify/src/fp`, `cargo run -p gen-fp -- --check
  crates/verify/src/fp`, `cargo run --release -p gen-precomp -- --check .`
  must find every generated file current (CI job `generated`).
- Fuzz (`cargo +nightly fuzz run <target>` in `crates/verify/fuzz`):
  `fuzz_pubkey_from_bytes` (decode and prepare arbitrary key bytes),
  `fuzz_signature_from_bytes` (decode, re-encode, compare), `fuzz_verify`
  (arbitrary signature bytes against a KAT key).
- Whole-workspace gates: `fmt`, `clippy -D warnings` on all targets, the
  `thumbv7em-none-eabihf` check of `sqisign-verify`, `cargo doc` with
  warnings denied, `cargo audit`. `#![forbid(unsafe_code)]` is `deny` on
  x86-64 with the `asm!` blocks and the `cpuid` read as the only opt-outs,
  and stays `forbid` elsewhere.

## Bounds and timing

Every entry point checks lengths before arithmetic
(`public_key_from_bytes`, `signature_from_bytes`, `compressed_from_bytes`,
`secret_key_from_bytes` return `None` on a wrong length or a
non-canonical field element) and the verifier rejects a non-canonical
basis-change matrix before any curve operation. The entangled-basis
search is bounded and fails closed. The field kernels are straight-line;
the two kernel paths compute identical limbs. [SECURITY.md](SECURITY.md)
states the timing status.

## Benchmarks

`cargo bench -p sqisign-rs --bench bench` measures key generation, signing
and verification (from bytes, cold and with a prepared key, standard and
compressed) in `rdtsc` cycles, criterion means. Benchmarks are of the
default build (no features, no `RUSTFLAGS`); every C reference number
names its build. [BENCH.md](BENCH.md) holds the session table.

`-C target-cpu=native` is not set: with it the saturated-limb build was
slower (level V verification 97.5 against 89.5 Mcycles).

## Not here

- The compact (dimension-4) format at levels III and V: derived in
  `docs/COMPACT_R3.md`, not built (level I is, behind `compact`).
- SQIsign-RK on round 3: a point release.
- Everything PRISM and ECLIPSE: the hash to prime degrees and its
  Miller-Rabin kernels, the ECLIPSE verifier, the canonical encoding.
- Round 2: removed entirely in 0.6.
