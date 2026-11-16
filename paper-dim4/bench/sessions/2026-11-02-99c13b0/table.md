#### Rust rows (this crate; megacycles over `runs`, medians; the signing rows with mean and maximum; one pinned core, round-robin over the session's rounds)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |
|---|---|---|---|---|---|---|---|---|---|
| Rust, dimension 2 | 2 | I | keygen | 52 | 38.5 | 39.6 | 33.8 | 58.9 | 12.62 |
| Rust, dimension 2 | 2 | I | sign | 52 | 112.5 | 115.2 | 98.4 | 163.9 | 36.89 |
| Rust, dimension 2 | 2 | I | verify from bytes | 200 | 13.3 | 13.2 | 12.9 | 13.6 | 4.35 |
| Rust, dimension 2 | 2 | I | verify compressed from bytes | 200 | 15.4 | 15.5 | 15.3 | 16.1 | 5.04 |
| Rust, dimension 2 | 2 | III | keygen | 52 | 142.2 | 169.7 | 89.6 | 806.1 | 46.62 |
| Rust, dimension 2 | 2 | III | sign | 52 | 392.9 | 439.4 | 262.9 | 1441.5 | 128.83 |
| Rust, dimension 2 | 2 | III | verify from bytes | 200 | 33.7 | 33.3 | 32.5 | 34.4 | 11.04 |
| Rust, dimension 2 | 2 | III | verify compressed from bytes | 200 | 39.2 | 39.1 | 38.4 | 39.9 | 12.86 |
| Rust, dimension 2 | 2 | V | keygen | 52 | 257.0 | 270.4 | 215.5 | 426.0 | 84.25 |
| Rust, dimension 2 | 2 | V | sign | 52 | 787.5 | 837.5 | 628.2 | 1539.4 | 258.20 |
| Rust, dimension 2 | 2 | V | verify from bytes | 200 | 80.7 | 80.7 | 80.3 | 81.5 | 26.44 |
| Rust, dimension 2 | 2 | V | verify compressed from bytes | 200 | 96.3 | 96.2 | 94.2 | 99.0 | 31.56 |
| Rust, dimension 4 (compact) | 4 | I | keygen | 52 | 40.8 | 42.1 | 36.9 | 56.1 | 13.37 |
| Rust, dimension 4 (compact) | 4 | I | sign | 100 | 329.0 | 438.2 | 48.8 | 1900.9 | 107.88 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, four points | 200 | 85.4 | 85.2 | 82.6 | 86.8 | 28.01 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points with the image check | 200 | 60.6 | 60.6 | 59.6 | 62.7 | 19.86 |
| Rust, dimension 4 (compact) | 4 | I | verify from bytes, two points without the image check | 200 | 59.6 | 59.6 | 58.6 | 60.6 | 19.54 |

#### The C reference (the-sqisign nist-v3, Release, broadwell; its own benchmark and cycle counter, medians over iterations)

| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max |
|---|---|---|---|---|---|---|---|---|
| C reference, dimension 2 | 2 | I | keygen | 300 | 35.116 | 36.264 | 30.603 | 65.473 |
| C reference, dimension 2 | 2 | I | sign | 300 | 103.753 | 106.190 | 90.243 | 172.688 |
| C reference, dimension 2 | 2 | I | verify from bytes | 300 | 12.642 | 12.782 | 12.528 | 13.604 |
| C reference, dimension 2 | 2 | III | keygen | 100 | 104.951 | 139.942 | 82.512 | 799.648 |
| C reference, dimension 2 | 2 | III | sign | 100 | 359.671 | 374.923 | 243.619 | 623.526 |
| C reference, dimension 2 | 2 | III | verify from bytes | 100 | 31.254 | 31.600 | 31.020 | 33.422 |
| C reference, dimension 2 | 2 | V | keygen | 100 | 218.470 | 232.698 | 196.393 | 483.795 |
| C reference, dimension 2 | 2 | V | sign | 100 | 730.877 | 766.652 | 567.541 | 1409.017 |
| C reference, dimension 2 | 2 | V | verify from bytes | 100 | 77.577 | 78.381 | 77.247 | 82.672 |

#### Ratios from the printed cells

- two_point_with_check_over_four_point: 0.71
- two_point_without_check_over_four_point: 0.7
- four_point_over_c_reference_verify_I: 6.8
- two_point_with_check_over_c_reference_verify_I: 4.8
- four_point_over_dim2_standard_verify_I: 6.4
- two_point_with_check_over_dim2_standard_verify_I: 4.6
