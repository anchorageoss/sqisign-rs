# Benchmarks

One session, all three levels, the default build (`cargo bench -p
sqisign-rs --bench bench`, no features, no `RUSTFLAGS`), the same binary
under `SQISIGN_FORCE_PORTABLE=1` for the portable column, and the C
reference (the-sqisign at `nist-v3`, commit `6d017708`) built twice on the
same machine: `broadwell` (`-DSQISIGN_BUILD_TYPE=broadwell`, the
assembly-backed field arithmetic) and `ref`. Rust numbers are criterion
means in `rdtsc` cycles; C numbers are the reference benchmark's medians
(`apps/benchmark_<prime> --iterations=300/100/100`). Every ratio names
the C build it is against.

## Session 2026-09-25

Machine: Intel(R) Xeon(R) CPU @ 2.80GHz; `rustc 1.92.0 (ded5c06cf 2025-12-08)`; sqisign-rs at `aaadbae` (the 0.6 line, this branch);
the-sqisign at `nist-v3` (`6d01770`) built with CMake `Release` twice,
`-DSQISIGN_BUILD_TYPE=broadwell` (`libsqisign_gf_<prime>` with the
`mulx`/`adcx`/`adox` field arithmetic; `-O2 -flto=auto -funroll-loops`,
the project's defaults) and `ref`. Rust: `cargo bench -p sqisign-rs
--bench bench`, criterion means of `rdtsc` cycles (20 samples of 20 s for
key generation and signing, 100 of 10 s for verification), then the same
binary under `SQISIGN_FORCE_PORTABLE=1`. C: `apps/benchmark_<prime>
--iterations=300` (level I), `100` (III, V), medians. The compressed row
has no C counterpart (the format is this crate's); its ratio is against
the C `broadwell` cold verification like the prepared row's.

| level | operation | default build (Mcycles) | portable, forced | C `broadwell` (median) | C `ref` (median) | default / `broadwell` |
|---|---|---|---|---|---|---|
| I | keygen | 39.2 | 86.4 | 34.6 | 49.8 | 1.13 |
| I | sign | 119.3 | 247.7 | 103.8 | 142.5 | 1.15 |
| I | verify from bytes | 12.3 | 35.2 | 11.3 | 18.9 | 1.09 |
| I | verify, prepared key | 11.9 | 33.9 |  |  | 1.05 (prepared / `broadwell` cold) |
| I | verify compressed, prepared key | 14.1 | 40.3 |  |  | 1.24 (prepared / `broadwell` cold) |
| III | keygen | 185.2 | 295.4 | 113.6 | 170.8 | 1.63 |
| III | sign | 498.2 | 799.6 | 341.9 | 506.5 | 1.46 |
| III | verify from bytes | 30.9 | 85.2 | 28.7 | 53.3 | 1.07 |
| III | verify, prepared key | 29.9 | 82.3 |  |  | 1.04 (prepared / `broadwell` cold) |
| III | verify compressed, prepared key | 35.1 | 97.3 |  |  | 1.22 (prepared / `broadwell` cold) |
| V | keygen | 306.7 | 510.3 | 239.6 | 285.6 | 1.28 |
| V | sign | 937.3 | 1451.9 | 754.3 | 868.1 | 1.24 |
| V | verify from bytes | 88.5 | 182.6 | 79.3 | 101.8 | 1.12 |
| V | verify, prepared key | 86.2 | 176.3 |  |  | 1.09 (prepared / `broadwell` cold) |
| V | verify compressed, prepared key | 97.2 | 208.1 |  |  | 1.22 (prepared / `broadwell` cold) |

Reading the table:

- Verification from bytes is within 7 to 12 % of the reference's
  `broadwell` build at the three levels and 1.15 to 1.7x faster than
  its `ref` build (1.5 / 1.7 / 1.15: level V's `ref` build is close to
  its assembly build). A prepared key saves 3 % (the per-key basis).
- The compressed format's recovery (two Weil pairings and one discrete
  logarithm at `2^(RESPONSE_BITS + 2)`) costs 2.2 / 5.2 / 11.0 Mcycles
  over the prepared verification: 18 / 17 / 13 % of it, for 12 % fewer
  bytes.
- Key generation and signing are 1.1 to 1.6x the reference's `broadwell`
  build (the quaternion layer on fixed-precision integers, variable-time
  by design); level III is the slowest ratio.
- The portable kernels (x86-64 without BMI2 and ADX) verify 2.1 to 2.9x
  slower than the assembly kernels; they are the correctness fallback.

### Field arithmetic, same session

| level | operation (dependent chain of 1000, cycles per operation) | default build | portable, forced |
|---|---|---|---|
| I | Fp mul | 47.2 | 207.8 |
| I | Fp2 mul | 178.1 | 557.9 |
| I | Fp2 sqr | 114.9 | 418.9 |
| I | Fp2 add | 24.7 | 24.7 |
| III | Fp mul | 78.7 | 310.4 |
| III | Fp2 mul | 308.1 | 845.0 |
| III | Fp2 sqr | 184.1 | 631.8 |
| III | Fp2 add | 32.9 | 32.7 |
| V | Fp mul | 165.8 | 477.5 |
| V | Fp2 mul | 671.6 | 1512.4 |
| V | Fp2 sqr | 391.7 | 993.6 |
| V | Fp2 add | 48.1 | 48.0 |

The `Fp` multiplication chain matches the reference's `broadwell` field
arithmetic (47 cycles at level I in prism-rs's C harness, the same
kernels); see prism-rs `docs/BENCH.md` Sections 11 and 12 for the
kernel-level comparison and the in-place API's level V number.

The compressed row is the compressed format's verification with a
prepared key ([COMPRESSION.md](COMPRESSION.md)): the standard
verification plus two Weil pairings and one discrete logarithm.

### Compact format, level I (separate run, same machine, feature `compact`, 2026-09-25)

`cargo bench -p sqisign-rs --features compact --bench bench -- compact`;
no C reference exists at these parameters (the SQIsignHD authors published
none for round-3 primes), so the comparison is with the round-3
two-dimensional rows of the session above.

| operation | Mcycles | against the two-dimensional format |
|---|---|---|
| keygen | 42.8 | 1.1x (39.2) |
| sign (mean; 461 to 619 over the samples, the rejection loop) | 538 | 4.5x (119.3) |
| verify from bytes, 142 B | 122.2 | 9.9x (12.3) |

The round-2 compact verifier cost 7.7x its two-dimensional verifier at
`5·2^248 − 1` (P18); at the round-3 prime the half-chains have 87 steps
instead of 68 on a 6-limb field instead of 4.

### Throughput (separate run, same machine and build, 2026-09-25)

`verify_batch` of 16 signatures under one prepared key, `cargo bench -p
sqisign-rs --bench bench -- verify_batch`, run after the session above
(the bench binary gained this row afterwards, so it is not in the
session's table):

| level | batch of 16 (Mcycles) | per signature (Mcycles) | signatures per second at 2.8 GHz, one thread |
|---|---|---|---|
| I | 192.6 | 12.0 | 233 |
| III | 482.6 | 30.2 | 93 |
| V | 1308.3 | 81.8 | 34 |

The batch is the prepared verification in a loop with the decoding
inside; it saves nothing beyond the prepared key, and the per-signature
numbers agree with the session's prepared rows within their spread.

The rows are not comparable across sessions; a new table replaces this one
whole.
