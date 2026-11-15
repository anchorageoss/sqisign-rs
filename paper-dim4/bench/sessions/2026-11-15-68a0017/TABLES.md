# TABLES: paper (Nov 4) / last session (2026-11-14, noisy) / this session (2026-11-15, quiet)

This session: 3 passes, value of a row: median of the three passes' values. Rows that moved by more than 1.5 percent between passes 1 and 2: c-III-keygen, c-III-sign, c-V-keygen, c-V-sign. All four are key generation and signing of the C reference, whose benchmark draws a fresh random seed on every run. Every Rust row moved by at most 1.1 percent across the three passes.
Megacycles differ between the paper's machine and this one; the ratios are comparable.

## Table 2: public key + signature, bytes

Equal in all three columns on every cell: standard 83+200, 129+306, 169+406; compressed 83+176, 129+269, 169+353; compact 84+139, 130+212, 170+278.

## Table 3: megacycles, medians (paper / last / this)

| implementation | level | keygen | sign | verify | x dim. 2 | x dim. 2 without the pre-chain checks (this session) |
|---|---|---|---|---|---|---|
| C reference, dimension 2 | I | 19.2 / 31.4 / 31.2 | 56.2 / 92.7 / 92.0 | 6.9 / 11.0 / 11.0 | | |
| C reference, dimension 2 | III | 62.1 / 101.9 / 106.2 | 181.6 / 312.4 / 323.1 | 18.0 / 29.1 / 29.0 | | |
| C reference, dimension 2 | V | 133.0 / 218.4 / 216.9 | 384.2 / 690.5 / 679.6 | 44.4 / 75.0 / 75.3 | | |
| Rust, dimension 2 | I | 23.7 / 32.9 / 32.9 | 67.2 / 99.9 / 100.0 | 7.1 / 11.2 / 11.2 | | |
| Rust, dimension 2 | III | 72.6 / 108.2 / 109.0 | 219.2 / 377.9 / 375.4 | 18.6 / 29.9 / 29.8 | | |
| Rust, dimension 2 | V | 143.6 / 249.1 / 248.6 | 453.6 / 772.4 / 769.7 | 50.0 / 80.3 / 80.0 | | |
| Rust, dimension 2, compressed verify | I | | | 8.4 / 13.3 / 13.3 | | |
| Rust, dimension 4 (compact) | I | 22.6 / 35.3 / 35.3 | 219.5 / 355.4 / 360.5 | 45.6 / 76.8 / 76.5 | 6.4 / 6.9 / 6.8 | 6.7 |
| two points | I | | | 28.5 / 48.4 / 48.0 | 4.0 / 4.3 / 4.3 | 4.2 |
| two points, no check at the gluing | I | | | not in the paper / 47.5 / 47.3 | | |
| Rust, dimension 4 (compact) | III | 75.1 / 126.9 / 126.3 | 610.7 / 916.8 / 913.9 | 118.9 / 195.3 / 194.6 | 6.4 / 6.5 / 6.5 | 6.4 |
| two points | III | | | 72.8 / 120.6 / 120.1 | 3.9 / 4.0 / 4.0 | 3.9 |
| two points, no check at the gluing | III | | | not in the paper / 119.2 / 118.6 | | |
| Rust, dimension 4 (compact) | V | 163.3 / 250.2 / 249.8 | 921.7 / 1711.4 / 1714.5 | 288.1 / 484.6 / 482.2 | 5.8 / 6.0 / 6.0 | 6.0 |
| two points | V | | | 174.1 / 292.5 / 291.6 | 3.5 / 3.6 / 3.6 | 3.6 |
| two points, no check at the gluing | V | | | not in the paper / 290.2 / 289.2 | | |

## Two-point over four-point verification, from the one-decimal cells

| level | paper | last session | this session | this session without the pre-chain checks |
|---|---|---|---|---|
| I | 0.63 | 0.63 | 48.0 / 76.5 = 0.63 | 0.62 |
| III | 0.61 | 0.62 | 120.1 / 194.6 = 0.62 | 0.61 |
| V | 0.60 | 0.60 | 291.6 / 482.2 = 0.60 | 0.60 |

## The pre-chain checks on their own, megacycles per verification (this session, medians of the profile runs of pass 1)

| level | decode (key, signature, challenge hash) | hints of both curves and the order of e(P, Q) | goodness test on q | checks together | share of a four-point verification | share of a two-point verification |
|---|---|---|---|---|---|---|
| I | 0.020 | 0.28 | 0.94 | 1.22 | 1.6 % | 2.6 % |
| III | 0.023 | 0.78 | 2.42 | 3.20 | 1.6 % | 2.7 % |
| V | 0.031 | 1.74 | 3.76 | 5.50 | 1.1 % | 1.9 % |

The decode was part of the paper's timings too (verification from bytes); the two checks were not. The check that k is odd costs one test and is not timed.

## Table 4: wall-clock milliseconds at level I (paper / last / this)

| format | operation | ms |
|---|---|---|
| dimension 2 (standard) | key generation | 10.3 / 11.8 / 11.8 |
|  | signing | 29.2 / 35.7 / 35.7 |
|  | verification | 3.1 / 4.0 / 4.0 |
|  | verification, compressed | 3.6 / 4.7 / 4.7 |
| dimension 4 (compact) | key generation | 9.8 / 12.6 / 12.6 |
|  | signing | 95.5 / 126.9 / 128.8 |
|  | verification | 19.8 / 27.4 / 27.3 |
|  | verification, two points | 12.4 / 17.3 / 17.2 |

## Table 5: the good-degree loop at level I (paper / last / this); the 100 signatures are the same in both of our sessions

| | value | expected (paper) |
|---|---|---|
| signatures measured | 100 / 100 / 100 | |
| response samples per signature, mean | 107 / 109 / 109 | 121 |
| response samples per signature, median (min, max) | 80 (1, 566) / 95 (1, 446) / 95 (1, 446) | |
| primality tests per signature, mean | 53 / 55 / 55 | 60 |
| signing, Mcycles, median (min, max) | 220 (27, 1400) / 359 (46, 1530) / 361 (46, 1546) | |
| signing, dimension 2, Mcycles, median | 67 / 100 / 100 | |

## Heap peak of one verification, KiB (paper / last / this; mean over five signatures)

| row | KiB |
|---|---|
| dim4-I-four | 1673 / 1655 / 1655 |
| dim4-I-two | 1637 / 1631 / 1631 |
| dim4-III-four | 3474 / 3443 / 3443 |
| dim4-III-two | not given / 3379 / 3379 |
| dim4-V-four | 6271 / 6160 / 6160 |
| dim4-V-two | not given / 6072 / 6072 |
| dim2-I-verify | 0 / 0 / 0 |

## Other ratios the prose quotes (paper / last / this)

| quantity | I | III | V |
|---|---|---|---|
| compact signing over dimension-2 signing | 3.3 / 3.6 / 3.6 | 2.8 / 2.4 / 2.4 | 2.0 / 2.2 / 2.2 |
| Rust dimension 2 over C: verify | 1.03 / 1.02 / 1.02 | | |
| Rust dimension 2 over C: keygen | 1.23 / 1.05 / 1.05 | | |
| Rust dimension 2 over C: sign | 1.20 / 1.08 / 1.09 | | |
| four-point verifier over the C reference's verify | 6.6 / 7.0 / 7.0 | | |

## The same harness on an Apple M4 Pro (wall-clock milliseconds, portable field backend, 2026-11-14, head 68a0017)

Unpacked in /home/user/work/results/20261114T225808Z-SECENG-2636-mac/m4/. Two passes within 3 percent on every row. The C reference was not run there.

| implementation | level | keygen | sign | verify | x dim. 2 |
|---|---|---:|---:|---:|---:|
| Rust, dimension 2 | I | 6.68 | 19.23 | 2.19 | |
| Rust, dimension 2 | III | 23.12 | 80.64 | 6.41 | |
| Rust, dimension 2 | V | 42.18 | 128.14 | 12.36 | |
| Rust, dimension 2, compressed verify | I | | | 2.60 | |
| Rust, dimension 2, compressed verify | III | | | 7.54 | |
| Rust, dimension 2, compressed verify | V | | | 14.45 | |
| Rust, dimension 4 (compact) | I | 6.79 | 65.88 | 14.71 | 6.7 |
| two points | I | | | 9.07 | 4.1 |
| two points, no check at the gluing | I | | | 8.89 | 4.1 |
| Rust, dimension 4 (compact) | III | 28.02 | 166.94 | 42.17 | 6.6 |
| two points | III | | | 25.43 | 4.0 |
| two points, no check at the gluing | III | | | 25.10 | 3.9 |
| Rust, dimension 4 (compact) | V | 41.47 | 307.46 | 77.84 | 6.3 |
| two points | V | | | 46.39 | 3.8 |
| two points, no check at the gluing | V | | | 45.97 | 3.7 |

Two-point over four-point on the M4 Pro: 0.62 / 0.60 / 0.60 at levels I / III / V.
