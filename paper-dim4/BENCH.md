# The dimension-4 paper's benchmark session

This file is the data source of `paper-dim4/main.tex`: its figure scripts
parse the section below. It was recorded on 2026-09-26 by
`paper-dim4/bench/session.sh` and lived in the repository root's BENCH.md,
which left the repository with the 0.6 working notes; it was restored here
on 2026-10-03 from the session's output, numbers unchanged. The crate's
own current numbers are in the README.

## Session 2026-09-26 (the dimension-4 paper)

One binary, `paper-dim4/bench`, measures the Rust code of both rounds in
one process: this crate at the 0.6 merge (`8718c97` on `main`; round 3,
dimension 2 and the compact dimension-4 format at level I) and the last
round-2 release, tag `v0.4.31` (`a03b834`), pulled from GitHub as a
dependency (round 2, dimension 2 at three levels and the 108-byte compact
format at level I). Both are compiled by the same `rustc` in release mode
with `lto = "fat"`, `codegen-units = 1` and no target-CPU flag; the
round-3 code picks its assembly field kernels at run time from CPUID, the
round-2 crate has none and signs on `num-bigint`. Each call is timed with
`rdtsc` and the wall clock; medians over the stated runs, seeded
randomness. The C reference is the-sqisign at `nist-v3` (`6d01770`) and
`nist-v2` (`91e9e46`), each built with CMake `Release` and
`-DSQISIGN_BUILD_TYPE=broadwell` (the `mulx`/`adcx`/`adox` field
arithmetic, the project's default flags), measured by its own
`apps/benchmark_<variant> --iterations=300` (level I) or `100` (III, V);
the C `verify` row is that benchmark's verification of the encoded
signature. The compact signer's rejection loop is counted by
`compact_sign_stats`; the memory row is a counting allocator around one
verification. Machine: Intel(R) Xeon(R) CPU @ 2.80GHz (8 vCPUs, 31 GiB);
`rustc 1.92.0 (ded5c06cf 2025-12-08)`; Linux 6.12.105+; gcc as installed;
`paper-dim4/bench/session.sh`.

#### Sizes observed (bytes)

| implementation | round | dimension | level | public key | signature |
|---|---|---|---|---|---|
| SQIsign (Rust, this crate) | 3 | 2 | I | 83 | 200 |
| SQIsign compressed (Rust, this crate) | 3 | 2 | I | 83 | 176 |
| SQIsign (Rust, this crate) | 3 | 2 | III | 129 | 306 |
| SQIsign compressed (Rust, this crate) | 3 | 2 | III | 129 | 269 |
| SQIsign (Rust, this crate) | 3 | 2 | V | 169 | 406 |
| SQIsign compressed (Rust, this crate) | 3 | 2 | V | 169 | 353 |
| SQIsign (Rust, this crate) | 3 | 4 | I | 84 | 142 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | 65 | 148 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | III | 97 | 224 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | V | 129 | 292 |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | 64 | 108 |

#### Key generation, signing, verification (Rust; megacycles over `runs`, and wall-clock milliseconds)

| implementation | round | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|---|
| SQIsign (Rust, this crate) | 3 | 2 | I | keygen | 50 | 39.0 | 39.5 | 33.4 | 50.8 | 13.93 |
| SQIsign (Rust, this crate) | 3 | 2 | I | sign | 50 | 113.5 | 116.8 | 97.5 | 173.2 | 40.55 |
| SQIsign (Rust, this crate) | 3 | 2 | I | verify from bytes | 200 | 12.2 | 12.3 | 12.0 | 14.6 | 4.37 |
| SQIsign (Rust, this crate) | 3 | 2 | I | verify compressed from bytes | 200 | 14.4 | 14.4 | 14.2 | 16.3 | 5.14 |
| SQIsign (Rust, this crate) | 3 | 2 | III | keygen | 50 | 137.8 | 204.6 | 86.9 | 882.8 | 49.22 |
| SQIsign (Rust, this crate) | 3 | 2 | III | sign | 50 | 448.1 | 490.7 | 287.7 | 844.9 | 160.03 |
| SQIsign (Rust, this crate) | 3 | 2 | III | verify from bytes | 200 | 32.7 | 32.8 | 32.3 | 35.8 | 11.69 |
| SQIsign (Rust, this crate) | 3 | 2 | III | verify compressed from bytes | 200 | 35.9 | 36.1 | 35.5 | 42.7 | 12.83 |
| SQIsign (Rust, this crate) | 3 | 2 | V | keygen | 50 | 280.1 | 298.2 | 227.9 | 461.6 | 100.04 |
| SQIsign (Rust, this crate) | 3 | 2 | V | sign | 50 | 922.6 | 973.4 | 683.7 | 1639.2 | 329.52 |
| SQIsign (Rust, this crate) | 3 | 2 | V | verify from bytes | 200 | 83.0 | 83.1 | 82.3 | 87.4 | 29.63 |
| SQIsign (Rust, this crate) | 3 | 2 | V | verify compressed from bytes | 200 | 96.7 | 97.3 | 95.9 | 130.4 | 34.54 |
| SQIsign (Rust, this crate) | 3 | 4 | I | keygen | 50 | 41.7 | 41.6 | 35.4 | 52.4 | 14.88 |
| SQIsign (Rust, this crate) | 3 | 4 | I | sign | 100 | 336.6 | 546.2 | 49.8 | 3712.3 | 120.21 |
| SQIsign (Rust, this crate) | 3 | 4 | I | verify from bytes | 100 | 120.3 | 120.6 | 119.1 | 124.2 | 42.98 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | keygen | 50 | 175.8 | 179.0 | 127.9 | 253.1 | 62.77 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | sign | 50 | 429.3 | 422.7 | 312.5 | 628.2 | 153.33 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | verify from bytes | 200 | 13.6 | 13.6 | 13.5 | 15.8 | 4.85 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | III | keygen | 50 | 412.2 | 490.0 | 298.0 | 2536.0 | 147.22 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | III | sign | 50 | 1653.7 | 1998.8 | 933.1 | 6812.1 | 590.61 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | III | verify from bytes | 200 | 38.2 | 38.8 | 37.9 | 73.2 | 13.65 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | V | keygen | 50 | 1169.8 | 1421.1 | 564.7 | 4611.3 | 417.80 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | V | sign | 50 | 3036.1 | 3153.6 | 1902.6 | 6639.8 | 1084.34 |
| SQIsign (Rust, v0.4.31) | 2 | 2 | V | verify from bytes | 200 | 77.1 | 77.2 | 76.4 | 80.0 | 27.54 |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | keygen | 50 | 173.8 | 185.2 | 132.6 | 361.9 | 62.07 |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | sign | 100 | 204.6 | 221.5 | 156.4 | 363.2 | 73.09 |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | verify from bytes | 100 | 108.6 | 108.7 | 107.6 | 111.6 | 38.80 |

#### The compact signer's rejection loop (round 3, level I, the `sign` rows' signatures)

| signatures | response samples per signature (mean) | median | min | max | primality tests per signature (mean) |
|---|---|---|---|---|---|
| 100 | 113.9 | 68.0 | 1 | 835 | 57.0 |

#### Memory of one compact verification (round 3, level I)

| peak heap during the call (bytes) | heap live before the call (bytes) | process VmHWM after (KiB) |
|---|---|---|
| 1714576 | 3596 | 9184 |

#### The C reference (`broadwell` builds; its own benchmark, medians over iterations)

| implementation | round | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|---|
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | I | keygen | 300 | 35.472 | 36.371 | 28.816 | 61.927 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | I | sign | 300 | 101.428 | 104.667 | 86.431 | 169.930 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | I | verify from bytes | 300 | 11.328 | 11.567 | 11.136 | 19.937 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | III | keygen | 100 | 111.415 | 131.260 | 76.273 | 323.587 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | III | sign | 100 | 370.970 | 395.354 | 232.108 | 725.127 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | III | verify from bytes | 100 | 28.701 | 29.396 | 28.242 | 47.440 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | V | keygen | 100 | 241.455 | 260.506 | 205.267 | 533.992 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | V | sign | 100 | 780.634 | 823.224 | 610.295 | 1450.498 | |
| SQIsign (C reference, nist-v3, broadwell) | 3 | 2 | V | verify from bytes | 100 | 79.778 | 81.213 | 78.898 | 108.082 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | I | keygen | 300 | 48.405 | 49.259 | 41.868 | 68.156 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | I | sign | 300 | 112.398 | 113.681 | 97.246 | 192.396 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | I | verify from bytes | 300 | 5.310 | 5.583 | 5.082 | 10.404 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | III | keygen | 100 | 136.905 | 141.313 | 118.232 | 295.231 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | III | sign | 100 | 314.365 | 326.838 | 269.434 | 687.346 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | III | verify from bytes | 100 | 17.022 | 17.070 | 16.773 | 19.392 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | V | keygen | 100 | 210.627 | 219.721 | 180.199 | 325.786 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | V | sign | 100 | 527.441 | 540.963 | 419.725 | 804.966 | |
| SQIsign (C reference, nist-v2, broadwell) | 2 | 2 | V | verify from bytes | 100 | 32.580 | 32.669 | 32.158 | 34.084 | |

#### M4 Pro (wall-clock milliseconds; the author's run of `paper-dim4/bench` on that machine, pending)

| implementation | round | dimension | level | operation | runs | median (ms) |
|---|---|---|---|---|---|---|
| SQIsign (Rust, this crate) | 3 | 2 | I | keygen | 50 | pending |
| SQIsign (Rust, this crate) | 3 | 2 | I | sign | 50 | pending |
| SQIsign (Rust, this crate) | 3 | 2 | I | verify from bytes | 200 | pending |
| SQIsign (Rust, this crate) | 3 | 2 | I | verify compressed from bytes | 200 | pending |
| SQIsign (Rust, this crate) | 3 | 4 | I | keygen | 50 | pending |
| SQIsign (Rust, this crate) | 3 | 4 | I | sign | 100 | pending |
| SQIsign (Rust, this crate) | 3 | 4 | I | verify from bytes | 100 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | keygen | 50 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | sign | 50 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 2 | I | verify from bytes | 200 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | keygen | 50 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | sign | 100 | pending |
| SQIsign (Rust, v0.4.31) | 2 | 4 | I | verify from bytes | 100 | pending |
