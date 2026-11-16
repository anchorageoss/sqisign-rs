# TABLES: date: 2026-11-14T23:12:21Z

Field kernels: portable (aarch64). Timer: the monotonic wall clock (std::time::Instant), nanoseconds.
Rows that moved by more than 3 percent between the two passes are listed in numbers.json.

# Wall-clock milliseconds, medians (the paper's Table 3 layout)

Passes: 2; session value of a row: pass 1 (two passes within 3 percent on every row).
Compact verification at levels III and V: the mean over six signatures of the median of 15 runs.

| implementation | level | keygen | sign | verify | x dim. 2 |
|---|---|---:|---:|---:|---:|
| C reference, dimension 2 (portable build) | I | not run | not run | not run | |
| C reference, dimension 2 (portable build) | III | not run | not run | not run | |
| C reference, dimension 2 (portable build) | V | not run | not run | not run | |
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

## Two-point over four-point verification

| level | two / four | ratio |
|---|---|---:|
| I | 9.07 / 14.71 | 0.62 |
| III | 25.43 / 42.17 | 0.60 |
| V | 46.39 / 77.84 | 0.60 |

## The good-degree loop at level I

| | value |
|---|---|
| signatures measured | 100 |
| response samples per signature, mean | 109 |
| response samples per signature, median (min, max) | 95 (1, 446) |
| primality tests per signature, mean | 55 |
| signing, ms, median (min, max) | 65.9 (8.9, 284.3) |
| commitments abandoned for a fresh one | 0 |

## Heap peak of one verification, KiB (counting allocator around the call)

| row | mean over 5 signatures (max) |
|---|---|
| dim2-I-verify | 0 (0) |
| dim2-III-verify | 0 (0) |
| dim2-V-verify | 0 (0) |
| dim4-I-four | 1655 (1680) |
| dim4-I-two | 1631 (1656) |
| dim4-III-four | 3873 (3882) |
| dim4-III-two | 3801 (3810) |
| dim4-V-four | 6160 (6183) |
| dim4-V-two | 6072 (6095) |
