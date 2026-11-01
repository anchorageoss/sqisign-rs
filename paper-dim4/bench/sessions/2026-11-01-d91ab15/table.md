#### Rust rows (this crate; megacycles over `runs`, medians; the signing rows with mean and maximum; one pinned core, round-robin over the session's rounds)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|
| Rust, dimension 2 | 2 | I | keygen | 52 | 39.0 | 39.9 | 33.9 | 58.7 | 12.77 |
| Rust, dimension 2 | 2 | I | sign | 52 | 113.8 | 116.0 | 99.0 | 166.2 | 37.30 |
| Rust, dimension 2 | 2 | I | verify from bytes | 200 | 14.3 | 14.1 | 13.8 | 15.5 | 4.68 |
| Rust, dimension 2 | 2 | I | verify compressed from bytes | 200 | 16.3 | 16.4 | 16.2 | 17.0 | 5.33 |
| Rust, dimension 2 | 2 | III | keygen | 52 | 142.5 | 169.8 | 90.2 | 815.8 | 46.72 |
| Rust, dimension 2 | 2 | III | sign | 52 | 395.4 | 441.2 | 263.3 | 1457.5 | 129.64 |
| Rust, dimension 2 | 2 | III | verify from bytes | 200 | 34.6 | 34.5 | 33.8 | 35.3 | 11.35 |
| Rust, dimension 2 | 2 | III | verify compressed from bytes | 200 | 40.4 | 40.3 | 39.5 | 41.1 | 13.23 |
| Rust, dimension 2 | 2 | V | keygen | 52 | 256.7 | 270.7 | 215.2 | 427.6 | 84.17 |
| Rust, dimension 2 | 2 | V | sign | 52 | 787.3 | 836.0 | 627.9 | 1536.7 | 258.12 |
| Rust, dimension 2 | 2 | V | verify from bytes | 200 | 84.2 | 84.0 | 83.3 | 85.5 | 27.59 |
| Rust, dimension 2 | 2 | V | verify compressed from bytes | 200 | 100.8 | 99.4 | 97.6 | 103.9 | 33.05 |
| Rust, dimension 4 (compact) | 4 | I | keygen | 52 | 41.2 | 42.6 | 37.6 | 56.5 | 13.50 |
| Rust, dimension 4 (compact) | 4 | I | sign | 100 | 330.6 | 443.7 | 49.5 | 1902.2 | 108.38 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, four points | 200 | 102.2 | 102.1 | 100.8 | 103.8 | 33.50 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points with the image check | 200 | 75.2 | 75.2 | 74.4 | 76.2 | 24.67 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points without the image check | 200 | 72.5 | 72.5 | 71.5 | 73.4 | 23.76 |

#### The C reference (the-sqisign nist-v3, Release, broadwell; its own benchmark and cycle counter, medians over iterations)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max |
|---|---|---|---|---|---|---|---|---|
| C reference, dimension 2 | 2 | I | keygen | 300 | 35.156 | 35.836 | 30.609 | 49.961 |
| C reference, dimension 2 | 2 | I | sign | 300 | 102.579 | 106.899 | 89.543 | 223.018 |
| C reference, dimension 2 | 2 | I | verify from bytes | 300 | 12.620 | 12.771 | 12.527 | 13.621 |
| C reference, dimension 2 | 2 | III | keygen | 100 | 110.323 | 132.335 | 80.994 | 439.546 |
| C reference, dimension 2 | 2 | III | sign | 100 | 340.671 | 366.524 | 238.236 | 778.206 |
| C reference, dimension 2 | 2 | III | verify from bytes | 100 | 31.243 | 31.571 | 31.060 | 33.371 |
| C reference, dimension 2 | 2 | V | keygen | 100 | 228.152 | 248.095 | 197.046 | 512.307 |
| C reference, dimension 2 | 2 | V | sign | 100 | 716.683 | 766.894 | 571.732 | 1302.315 |
| C reference, dimension 2 | 2 | V | verify from bytes | 100 | 77.663 | 78.565 | 77.382 | 82.827 |

#### Ratios from the printed cells

- two_point_with_check_over_four_point: 0.74
- two_point_without_check_over_four_point: 0.71
- four_point_over_c_reference_verify_I: 8.1
- two_point_with_check_over_c_reference_verify_I: 6.0
- four_point_over_dim2_standard_verify_I: 7.1
- two_point_with_check_over_dim2_standard_verify_I: 5.3
