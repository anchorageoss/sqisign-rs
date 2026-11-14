# TABLES: the session's numbers beside the paper's (SECENG-2636)

Session: /home/user/work/results/20261114T213508Z-SECENG-2636/session. Passes: 3; session value of a row: median of the three passes' values.
Rows whose median moved by more than 3 percent between passes 1 and 2: c-III-keygen (95.61, 107.59), c-III-sign (297.58, 312.39), c-V-keygen (218.42, 227.38).
Every column "paper" is the number in main_standalone.tex (sha256 ce4b629c...) with its line; "diff" is ours minus the paper's.
The paper's numbers are from another machine (one vCPU of a Xeon at 2.30 GHz, 2026-11-04); ours are from this workstation (machine.txt). Megacycle counts are not comparable across machines; the ratios are.

## Table 2 (tex 1579 to 1581): public key + signature = total, bytes

| format | level | ours | paper | diff |
|---|---|---|---|---|
| standard (specification) | I | 83 + 200 = 283 | 83 + 200 = 283 | +0 |
| standard (specification) | III | 129 + 306 = 435 | 129 + 306 = 435 | +0 |
| standard (specification) | V | 169 + 406 = 575 | 169 + 406 = 575 | +0 |
| compressed (this crate) | I | 83 + 176 = 259 | 83 + 176 = 259 | +0 |
| compressed (this crate) | III | 129 + 269 = 398 | 129 + 269 = 398 | +0 |
| compressed (this crate) | V | 169 + 353 = 522 | 169 + 353 = 522 | +0 |
| compact, dim. 4 (this crate) | I | 84 + 139 = 223 | 84 + 139 = 223 | +0 |
| compact, dim. 4 (this crate) | III | 130 + 212 = 342 | 130 + 212 = 342 | +0 |
| compact, dim. 4 (this crate) | V | 170 + 278 = 448 | 170 + 278 = 448 | +0 |

## Table 3 (tex 1618 to 1630): megacycles, medians; last column the dimension-4 verifier over the dimension-2 verifier of the same code

| implementation | level | keygen: ours / paper / diff | sign: ours / paper / diff | verify: ours / paper / diff | x dim. 2: ours / paper |
|---|---|---|---|---|---|
| C reference, dimension 2 | I | 31.4 / 19.2 / +12.2 (+64 %) | 92.7 / 56.2 / +36.5 (+65 %) | 11.0 / 6.9 / +4.1 (+59 %) | |
| C reference, dimension 2 | III | 101.9 / 62.1 / +39.8 (+64 %) | 312.4 / 181.6 / +130.8 (+72 %) | 29.1 / 18.0 / +11.1 (+62 %) | |
| C reference, dimension 2 | V | 218.4 / 133.0 / +85.4 (+64 %) | 690.5 / 384.2 / +306.3 (+80 %) | 75.0 / 44.4 / +30.6 (+69 %) | |
| Rust, dimension 2 | I | 32.9 / 23.7 / +9.2 (+39 %) | 99.9 / 67.2 / +32.7 (+49 %) | 11.2 / 7.1 / +4.1 (+58 %) | |
| Rust, dimension 2 | III | 108.2 / 72.6 / +35.6 (+49 %) | 377.9 / 219.2 / +158.7 (+72 %) | 29.9 / 18.6 / +11.3 (+61 %) | |
| Rust, dimension 2 | V | 249.1 / 143.6 / +105.5 (+73 %) | 772.4 / 453.6 / +318.8 (+70 %) | 80.3 / 50.0 / +30.3 (+61 %) | |
| Rust, dimension 2, compressed verify | I | | | 13.3 / 8.4 / +4.9 (+58 %) | |
| Rust, dimension 4 (compact) | I | 35.3 / 22.6 / +12.7 (+56 %) | 355.4 / 219.5 / +135.9 (+62 %) | 76.8 / 45.6 / +31.2 (+68 %) | 6.9 / 6.4 |
| two points | I | | | 48.4 / 28.5 / +19.9 (+70 %) | 4.3 / 4.0 |
| two points, no check at the gluing (not in the paper) | I | | | 47.5 | |
| Rust, dimension 4 (compact) | III | 126.9 / 75.1 / +51.8 (+69 %) | 916.8 / 610.7 / +306.1 (+50 %) | 195.3 / 118.9 / +76.4 (+64 %) | 6.5 / 6.4 |
| two points | III | | | 120.6 / 72.8 / +47.8 (+66 %) | 4.0 / 3.9 |
| two points, no check at the gluing (not in the paper) | III | | | 119.2 | |
| Rust, dimension 4 (compact) | V | 250.2 / 163.3 / +86.9 (+53 %) | 1711.4 / 921.7 / +789.7 (+86 %) | 484.6 / 288.1 / +196.5 (+68 %) | 6.0 / 5.8 |
| two points | V | | | 292.5 / 174.1 / +118.4 (+68 %) | 3.6 / 3.5 |
| two points, no check at the gluing (not in the paper) | V | | | 290.2 | |

## The two-point ratios (tex 109, 1180, 1779): verify two-point over four-point

| level | ours, from the cells | paper | diff |
|---|---|---|---|
| I | 48.4 / 76.8 = 0.63 | 0.63 | +0.00 |
| III | 120.6 / 195.3 = 0.62 | 0.61 | +0.01 |
| V | 292.5 / 484.6 = 0.60 | 0.60 | +0.00 |

## Table 4 (tex 1679 to 1687): wall-clock milliseconds at level I, medians

| format | operation | ours / paper / diff |
|---|---|---|
| dimension 2 (standard) | key generation | 11.8 / 10.3 / +1.5 (+15 %) |
|  | signing | 35.7 / 29.2 / +6.5 (+22 %) |
|  | verification | 4.0 / 3.1 / +0.9 (+29 %) |
|  | verification, compressed | 4.7 / 3.6 / +1.1 (+31 %) |
| dimension 4 (compact) | key generation | 12.6 / 9.8 / +2.8 (+29 %) |
|  | signing | 126.9 / 95.5 / +31.4 (+33 %) |
|  | verification | 27.4 / 19.8 / +7.6 (+38 %) |
|  | verification, two points | 17.3 / 12.4 / +4.9 (+40 %) |

## Table 5 (tex 1711 to 1717): the compact signer's good-degree loop at level I

The loop's counts and its signing times are those of pass 1 (the same 100 signatures in every pass); the Table 3 cell is the median of the passes' medians.

| | ours | paper | expected (paper) |
|---|---|---|---|
| signatures measured | 100 | 100 | |
| response samples per signature, mean | 109 | 107 | 121 |
| response samples per signature, median (min, max) | 95 (1, 446) | 80 (1, 566) | |
| primality tests per signature, mean | 55 | 53 | 60 |
| signing, Mcycles, median (min, max) | 359 (46, 1530) | 220 (27, 1400) | |
| signing, dimension 2, Mcycles, median | 100 | 67 | |
| commitments abandoned for a fresh one | 0 | not reported | |

## The heap paragraph (tex 1722 to 1726): peak of one verification, KiB, counting allocator around the call

| row | ours, mean over 5 signatures (max) | paper |
|---|---|---|
| dim4-I-four | 1655 (1680) | 1673 |
| dim4-I-two | 1631 (1656) | 1637 |
| dim4-III-four | 3443 (3451) | 3474 |
| dim4-III-two | 3379 (3387) | not given |
| dim4-V-four | 6160 (6183) | 6271 |
| dim4-V-two | 6072 (6095) | not given |
| dim2-I-verify | 0 (0) | 0 ("allocates nothing") |
| dim2-III-verify | 0 (0) | 0 |
| dim2-V-verify | 0 (0) | 0 |

## Other numbers the prose quotes, recomputed from the cells

| quantity (tex line) | ours | paper |
|---|---|---|
| dimension-4 signing over dimension-2 signing, I / III / V (1590 to 1593, 1777) | 3.6 / 2.4 / 2.2 | 3.3 / 2.8 / 2.0 |
| dimension-4 verify over dimension-2, four points, I / III / V (1589 to 1592, 1775) | 6.9 / 6.5 / 6.0 | 6.4 / 6.4 / 5.8 |
| the same with two points (1592, 1774) | 4.3 / 4.0 / 3.6 | 4.0 / 3.9 / 3.5 |
| four-point verifier over the C reference's verify, level I (529) | 7.0 | 6.6 |
| Rust dimension 2 over the C reference at level I: verify / keygen / sign (1599 to 1601) | 1.02 / 1.05 / 1.08 | 1.03 / 1.23 / 1.20 |
| compressed recovery adds to a verification (1601) | 19 % | 18 % |
| key generation, compact against standard, level I (1598) | 35.3 against 32.9 | 22.6 against 23.7 |
| middle half of the compact signing times, III (1595) | 397 to 1456 | 300 to 1289 |
| middle half of the compact signing times, V (1595) | 867 to 3383 | 567 to 1596 |
| megacycles per step, two points, level I (1798) | 0.28 | 0.16 |
| cut of the verifier's cost by two points, I / III / V (abstract, 53) | 37 / 38 / 40 % | 37 to 40 % |
| profile of the four-point verifier at level I: kernel points doubled and pushed / codomains (531 to 534) | 59 % / 15 % | 59 % / 16 % |

### The phase profile at level I (share of one verification from bytes)

| phase | four points | two points |
|---|---|---|
| decode_share | 0.0 % | 0.0 % |
| basis generation: recover_challenge, recover_response, hd_torsion_basis | 8.9 % | 14.2 % |
| Cornacchia and kernel setup: norm_equation_2f_minus_q, make_canonical, product null, starting_two_symplectic_matrices, gluing base changes, b_kpp | 6.2 % | 9.9 % |
| dimension-2 steps: IsogenyChainDim2::new, both halves | 0.5 % | 0.8 % |
| gluing: KaniGluingChainHalf::new after the dimension-2 chain, and the gluing images of the generators, both halves | 4.2 % | 6.7 % |
| plain dimension-4 steps, point pushes: double_iter and image of the saved bases, both halves | 59.5 % | 51.4 % |
| plain dimension-4 steps, codomain: IsogenyDim4::from_kernel and the next structure's precompute, both halves | 15.0 % | 8.2 % |
| final comparisons: middle_codomain_matches | 0.0 % | 0.0 % |
| image check: hd_image_check | 4.5 % | 7.1 % |
| total, megacycles | 76.9 | 48.3 |
