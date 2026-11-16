#### Rust rows (this crate; megacycles over `runs`, medians; the signing rows with mean and maximum; one pinned core, round-robin over the session's rounds)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|
| Rust, dimension 2 | 2 | I | keygen | 52 | 39.0 | 39.9 | 33.9 | 58.9 | 12.80 |
| Rust, dimension 2 | 2 | I | sign | 52 | 113.5 | 115.7 | 99.0 | 163.5 | 37.23 |
| Rust, dimension 2 | 2 | I | verify from bytes | 200 | 14.2 | 14.1 | 13.8 | 14.4 | 4.65 |
| Rust, dimension 2 | 2 | I | verify compressed from bytes | 200 | 16.2 | 16.4 | 16.1 | 21.1 | 5.32 |
| Rust, dimension 2 | 2 | III | keygen | 52 | 143.7 | 170.7 | 90.0 | 810.1 | 47.12 |
| Rust, dimension 2 | 2 | III | sign | 52 | 394.1 | 441.4 | 264.7 | 1450.6 | 129.20 |
| Rust, dimension 2 | 2 | III | verify from bytes | 200 | 35.1 | 34.6 | 33.9 | 35.9 | 11.50 |
| Rust, dimension 2 | 2 | III | verify compressed from bytes | 200 | 40.5 | 40.4 | 39.8 | 41.5 | 13.27 |
| Rust, dimension 2 | 2 | V | keygen | 52 | 257.0 | 270.7 | 215.9 | 425.8 | 84.25 |
| Rust, dimension 2 | 2 | V | sign | 52 | 787.0 | 835.8 | 628.7 | 1535.0 | 258.03 |
| Rust, dimension 2 | 2 | V | verify from bytes | 200 | 83.8 | 83.8 | 83.4 | 85.6 | 27.47 |
| Rust, dimension 2 | 2 | V | verify compressed from bytes | 200 | 99.9 | 99.2 | 97.3 | 102.5 | 32.75 |
| Rust, dimension 4 (compact) | 4 | I | keygen | 52 | 41.2 | 42.5 | 37.3 | 56.4 | 13.51 |
| Rust, dimension 4 (compact) | 4 | I | sign | 100 | 331.9 | 442.4 | 48.9 | 1905.4 | 108.82 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, four points | 200 | 132.7 | 133.7 | 130.8 | 142.6 | 43.52 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points with the image check | 200 | 103.8 | 104.9 | 101.6 | 114.1 | 34.04 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points without the image check | 200 | 100.7 | 102.1 | 98.6 | 111.1 | 33.03 |

#### The C reference (the-sqisign nist-v3, Release, broadwell; its own benchmark and cycle counter, medians over iterations)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max |
|---|---|---|---|---|---|---|---|---|
| C reference, dimension 2 | 2 | I | keygen | 300 | 34.931 | 36.003 | 30.797 | 62.749 |
| C reference, dimension 2 | 2 | I | sign | 300 | 102.574 | 106.089 | 90.196 | 190.744 |
| C reference, dimension 2 | 2 | I | verify from bytes | 300 | 12.646 | 12.783 | 12.540 | 13.523 |
| C reference, dimension 2 | 2 | III | keygen | 100 | 110.113 | 124.346 | 80.424 | 359.301 |
| C reference, dimension 2 | 2 | III | sign | 100 | 327.448 | 347.971 | 240.079 | 791.134 |
| C reference, dimension 2 | 2 | III | verify from bytes | 100 | 31.248 | 31.603 | 31.025 | 33.349 |
| C reference, dimension 2 | 2 | V | keygen | 100 | 230.537 | 250.734 | 196.414 | 503.859 |
| C reference, dimension 2 | 2 | V | sign | 100 | 705.940 | 758.025 | 578.618 | 1286.058 |
| C reference, dimension 2 | 2 | V | verify from bytes | 100 | 77.825 | 79.032 | 77.397 | 83.363 |

#### Ratios from the printed cells

- two_point_with_check_over_four_point: 0.78
- two_point_without_check_over_four_point: 0.76
- four_point_over_c_reference_verify_I: 10.5
- two_point_with_check_over_c_reference_verify_I: 8.2
- four_point_over_dim2_standard_verify_I: 9.3
- two_point_with_check_over_dim2_standard_verify_I: 7.3
