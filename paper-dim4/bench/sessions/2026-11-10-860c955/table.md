#### Rust rows (this crate; megacycles over `runs`, medians; the signing rows with mean and maximum; one pinned core, round-robin over the session's rounds)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|
| Rust, dimension 2 | 2 | I | keygen | 52 | 38.5 | 39.4 | 33.3 | 58.4 | 12.63 |
| Rust, dimension 2 | 2 | I | sign | 52 | 112.3 | 114.5 | 97.4 | 165.1 | 36.82 |
| Rust, dimension 2 | 2 | I | verify from bytes | 200 | 13.3 | 13.2 | 12.8 | 13.8 | 4.36 |
| Rust, dimension 2 | 2 | I | verify compressed from bytes | 200 | 15.3 | 15.4 | 15.2 | 16.1 | 5.01 |
| Rust, dimension 2 | 2 | III | keygen | 52 | 141.3 | 167.7 | 89.1 | 794.2 | 46.32 |
| Rust, dimension 2 | 2 | III | sign | 52 | 392.1 | 437.6 | 260.7 | 1437.6 | 128.55 |
| Rust, dimension 2 | 2 | III | verify from bytes | 200 | 33.3 | 33.0 | 32.3 | 34.6 | 10.91 |
| Rust, dimension 2 | 2 | III | verify compressed from bytes | 200 | 39.0 | 38.8 | 38.1 | 39.6 | 12.78 |
| Rust, dimension 2 | 2 | V | keygen | 52 | 255.3 | 268.8 | 214.0 | 421.8 | 83.71 |
| Rust, dimension 2 | 2 | V | sign | 52 | 782.2 | 832.8 | 624.1 | 1524.8 | 256.46 |
| Rust, dimension 2 | 2 | V | verify from bytes | 200 | 80.2 | 80.2 | 79.7 | 80.7 | 26.30 |
| Rust, dimension 2 | 2 | V | verify compressed from bytes | 200 | 95.8 | 95.5 | 93.5 | 97.5 | 31.40 |
| Rust, dimension 4 (compact) | 4 | I | keygen | 52 | 40.5 | 41.9 | 36.7 | 56.2 | 13.28 |
| Rust, dimension 4 (compact) | 4 | I | sign | 100 | 326.8 | 436.4 | 48.4 | 1880.9 | 107.14 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, four points | 200 | 84.7 | 84.5 | 82.3 | 86.6 | 27.76 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points with the image check | 200 | 59.9 | 59.9 | 58.9 | 61.0 | 19.65 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points without the image check | 200 | 59.0 | 58.9 | 57.8 | 60.6 | 19.34 |

#### The C reference (the-sqisign nist-v3, Release, broadwell; its own benchmark and cycle counter, medians over iterations)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max |
|---|---|---|---|---|---|---|---|---|
| C reference, dimension 2 | 2 | I | keygen | 300 | 34.825 | 35.605 | 30.262 | 57.754 |
| C reference, dimension 2 | 2 | I | sign | 300 | 101.572 | 104.959 | 89.496 | 188.076 |
| C reference, dimension 2 | 2 | I | verify from bytes | 300 | 12.540 | 12.674 | 12.431 | 13.441 |
| C reference, dimension 2 | 2 | III | keygen | 100 | 107.255 | 132.613 | 80.526 | 741.099 |
| C reference, dimension 2 | 2 | III | sign | 100 | 373.414 | 398.064 | 240.575 | 1108.083 |
| C reference, dimension 2 | 2 | III | verify from bytes | 100 | 31.028 | 31.356 | 30.884 | 33.249 |
| C reference, dimension 2 | 2 | V | keygen | 100 | 229.433 | 250.747 | 194.210 | 491.839 |
| C reference, dimension 2 | 2 | V | sign | 100 | 729.884 | 778.408 | 568.435 | 1389.495 |
| C reference, dimension 2 | 2 | V | verify from bytes | 100 | 77.161 | 77.981 | 76.719 | 82.304 |

#### Ratios from the printed cells

- two_point_with_check_over_four_point: 0.71
- two_point_without_check_over_four_point: 0.7
- four_point_over_c_reference_verify_I: 6.8
- two_point_with_check_over_c_reference_verify_I: 4.8
- four_point_over_dim2_standard_verify_I: 6.4
- two_point_with_check_over_dim2_standard_verify_I: 4.5
