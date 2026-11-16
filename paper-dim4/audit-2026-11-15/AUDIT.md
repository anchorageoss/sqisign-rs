# AUDIT: every number of paper-dim4/main_standalone.tex against the session and the code

Second audit, 2026-11-15, run from scratch on the edited paper after the open items of the first were closed. File: paper-dim4/main_standalone.tex (2539 lines), the one current paper file. Its parent, the text with the session of 2026-11-04, is kept as paper-dim4/main_standalone.2026-11-04.orig.tex (sha256 ce4b629c8449064bc4ac35921952077651b4c5d26852b255607cf05356b87661). Line numbers are those of main_standalone.tex.

## Sources

- Session: the quiet three-pass session of 2026-11-15 at 68a0017 (paper-dim4/bench/sessions/2026-11-15-68a0017/; ~/work/results/20261114T224715Z-SECENG-2636-clean/session/). recompute.py reads only its raw files and prints every median, ratio, share, quartile, loop statistic and heap mean (recompute.out). A session value is the median of the three passes' values; ratios come from the unrounded medians.
- The pass after the divisor tests: one quiet pinned pass at bf2153b (pass-compare.txt). No verification row moved by more than 0.54 percent, so the session's tables stand.
- Instrumented measurements at bf2153b (feature instrument, paper-dim4/bench/src/paper_instrument.rs): instrument/gram.csv and instrument/hints.csv, with gram_minima.py's output in instrument/minima-all.txt (4 180 commitments) and minima-first100.txt (the loop table's 100).
- MIKE: session/mike-step.csv, paper-dim4/bench/mike-step against the copy of mike_rs (upstream 0d79121) on this machine, CPU 13.
- The 0.07 count: paper-dim4/bench/short_vector_count.py.
- Tests at the head of agent/dim4-two-point with the divisor tests: tests-head.log (142 passed, 0 failed, 3 ignored), the known-answer tests 3 of 3, the SQIsignHD-lib verifier on the five level-I vectors and the five tampered ones (sage-crosscheck/).
- Code: file and line under crates/ at the branch head.

The files named without a path are in this folder or in ~/work/results/20261115T081239Z-SECENG-2636-close/.

Status is one of: matches / derived (recomputed) / cited (not ours, left alone) / needs an author. No row is stale and none is without a source.


## Abstract (plain-text copy and typeset copy)

| line | text | number | source | status |
|---|---|---|---|---|
| 48 | plain-text: automorphism of order 4 | 4 | mathematics of the paper (Lemma J(i)); not a measurement | matches |
| 48 | plain-text: ten coordinates instead of sixteen | 10, 16 | code: ORBIT_REPS has 10 entries, THETA_DIM4_N = 16 (crates/verify/src/hd/two_point.rs:62) | matches |
| 48 | plain-text: cuts the verifier's cost by 37 to 40 % | 37, 40 | 1 - two/four from the unrounded session medians: 37.2, 38.3, 39.5 % at I, III, V (recompute.py) | derived (recomputed) |
| 48 | plain-text: 139 bytes against 200; key and signature 223 | 139, 200, 223 | raw headers: "84 + 139 bytes, packed layout" (pass1/raw/dim4-I-verify-four.csv), "83 + 200 bytes" (dim2-I-verify.csv); 84 + 139 = 223 | matches |
| 48 | plain-text: verification costs 4.3 times that of the dimension-2 scheme | 4.3 | two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (recompute.py). The paper of Nov 4 had 4.0, which "four" described; 4.3 is now written out | derived (recomputed) |
| 53 | typeset: automorphism of order 4 | 4 | mathematics of the paper (Lemma J(i)); not a measurement | matches |
| 53 | typeset: ten coordinates instead of sixteen | 10, 16 | code: ORBIT_REPS has 10 entries, THETA_DIM4_N = 16 (crates/verify/src/hd/two_point.rs:62) | matches |
| 53 | typeset: cuts the verifier's cost by 37 to 40 % | 37, 40 | 1 - two/four from the unrounded session medians: 37.2, 38.3, 39.5 % at I, III, V (recompute.py) | derived (recomputed) |
| 55 | typeset: 139 bytes against 200; key and signature 223 | 139, 200, 223 | raw headers: "84 + 139 bytes, packed layout" (pass1/raw/dim4-I-verify-four.csv), "83 + 200 bytes" (dim2-I-verify.csv); 84 + 139 = 223 | matches |
| 55 | typeset: verification costs 4.3 times that of the dimension-2 scheme | 4.3 | two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (recompute.py). The paper of Nov 4 had 4.0, which "four" described; 4.3 is now written out | derived (recomputed) |

## Introduction

| line | text | number | source | status |
|---|---|---|---|---|
| 68 | 2022 (the SIDH attacks) | 2022 | cited (CD23, MMP+23, Rob23) | cited (not ours, left alone) |
| 86 | 2026; p^(1/3+o(1)); September 2026 | 2026, 1/3 | cited (Wes26, AAA+26) | cited (not ours, left alone) |
| 89 | 3*2^324 - 1, 27*2^500 - 1, 17*2^664 - 1 | 3 primes | cited (AAA+26); the same primes in the code (crates/verify/src/params) | cited (not ours, left alone) |
| 109 | contribution list: 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | two/four from the unrounded medians: 0.6278, 0.6169, 0.6046; from the one-decimal cells: 48.0/76.5, 120.1/194.6, 291.6/482.2 give the same | derived (recomputed) |
| 124 | about once in 2^44 at level I | 2^44 | instrument/minima-all.txt: tail constant 1.53 over 4 180 commitments; 1.53 * (2^11.21)^-4 = 2^-44.2 (gram_minima.py) | derived (recomputed) |
| 134 | 139 bytes against 200 | 139, 200 | as in the abstract | matches |
| 136 | both verifiers at all three levels | 3 | tests at head: conformance_level1, _level3, _level5 pass (tests-head-8747fa9.log) | matches |

## Section 2, Preliminaries

| line | text | number | source | status |
|---|---|---|---|---|
| 216 | every S and star costs 16 squarings or multiplications; every H costs 64 additions | 16, 64 | 16 coordinates in dimension 4; code: hadamard_assign, 64 additions (doc at two_point.rs:80-81) | matches |
| 242 | glossary: dimensions 2 and 4; 2^g coordinates; four and sixteen theta coordinates | 2, 4, 16 | definitions | matches |
| 263 | megacycle: one million cycles read with rdtsc | 10^6 | harness: cycles / 1e6 (paper-dim4/bench/paper_report.py); timer.txt of each pass: "rdtsc and the monotonic wall clock" | matches |

## Section 3, Related work

| line | text | number | source | status |
|---|---|---|---|---|
| 297 | library rule e = lambda + ceil(log2 2 lambda), r = ceil(e/2) + 2, primes of about 2 lambda bits | formulas | cited (DLRW24, DLRW25) | cited (not ours, left alone) |
| 319 | MIKE: 228 steps in 16.0 megacycles | 228, 16.0 | cited (BDD+26, Table 1) | cited (not ours, left alone) |
| 321 | MIKE proposes log2 p = 317, 484, 651 instead of 325, 504, 668 | 317, 484, 651, 325, 504, 668 | cited (BDD+26, Section 4.5); the second triple is floor(log2 p) of the round-3 primes: 325.58, 504.75, 668.09 | cited (not ours, left alone) |

## Section 4.1, What is sent

| line | text | number | source | status |
|---|---|---|---|---|
| 359 | 2 * 41 = 82 bytes of curve, 22 bytes of q, three 12-byte scalars, two one-byte hints: 142 bytes | 41, 82, 22, 12, 2, 142 | ceil(326/8) = 41; ceil(174/8) = 22; ceil(89/8) = 12; 82 + 22 + 36 + 2 = 142 | derived (recomputed) |
| 360 | against 200 for the specification; matrix modulo 2^198 | 200, 198 | 200: raw header of dim2-I-verify.csv; 2^198: cited (AAA+26, Chapter 6) | matches |
| 365 | 82 bytes of the auxiliary curve replaced by 82 bytes of the commitment curve | 82, 82 | 2 * 41 | derived (recomputed) |
| 367 | four matrix entries take 100 bytes and the challenge 16, against 36 for three entries and 22 for q, for 58 bytes in all | 100, 16, 36, 22, 58 | 4 * ceil(198/8) = 100; 128/8 = 16; 3 * 12 = 36; (100 + 16) - (36 + 22) = 58 = 200 - 142 | derived (recomputed) |
| 372 | bit packing saves only 3 more, for 139 bytes | 3, 139 | 142 - 139 = 3; 1107 bits = 139 bytes | derived (recomputed) |

## Section 4.2, The verifier

| line | text | number | source | status |
|---|---|---|---|---|
| 539 | 87 steps from each end at level I; sixteen coordinates instead of four | 87, 16, 4 | ceil(174/2) = 87 | derived (recomputed) |
| 540 | costs 6.8 times the dimension-2 verifier (76.5 against 11.2 megacycles), and 7.0 times the reference | 6.8, 76.5, 11.2, 7.0 | medians of three passes: dim4-I-verify-four 76.540, dim2-I-verify 11.185, c-I-verify 10.982 (pass*/raw/*.csv, pass*/c-ref-I.txt); 76.540/11.185 = 6.84; 76.540/10.982 = 6.97 | derived (recomputed) |
| 544 | profile: 60 % to doubling and pushing the kernel points | 60 % | pass*/raw/dim4-I-profile-four.csv, slot4 / total: 59.7, 59.6, 59.8 % in the three passes | derived (recomputed) |
| 546 | 15 % to the codomains and their constants | 15 % | same files, slot5 / total: 15.0, 15.1, 15.1 % | derived (recomputed) |
| 438 | multiples by 2^(f - lambda - r); degree 2^lambda | formula | code: crates/verify/src/hd/challenge.rs:163-222 | matches |
| 468 | m > 40 in the list of rejects | 40 | MAX_GLUING_M = 40 | matches |

## Section 5, Two-point verification

| line | text | number | source | status |
|---|---|---|---|---|
| 592 | pi exchanges index bits 0 and 2 and bits 1 and 3 of the sixteen coordinates | 0, 2, 1, 3, 16 | code: PI = [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15] (two_point.rs:44) | matches |
| 1027 | (e, r) = (174, 89), (264, 134), (346, 175) | 174, 89, 264, 134, 346, 175 | e = ceil(log2 p / 2) + 11 with log2 p = 325.58, 504.75, 668.09; r = ceil(e/2) + 2; the same in crates/verify/src/hd/params.rs | derived (recomputed) |
| 1137 | 10-dimensional space of vectors fixed by pi; ten orbits | 10 | 4 fixed indices + 6 two-cycles (PI_FIXED, PI_PAIRS, two_point.rs:47-58) | matches |
| 1193 | 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | as in the introduction | derived (recomputed) |
| 1207 | tree on the ten orbits, 25 multiplications | 10, 9, 1, 25 | code: along_tree, TREE_ORBITS (two_point.rs:140-171): 8 + 8 + 7 + 2 = 25 multiplications per vector, called once for the dual null point and once for its inverse (two_point.rs:199-200) | derived (recomputed) |
| 1208 | Hadamard transform of a pi-fixed vector: 32 additions and 4 doublings instead of 64 additions | 32, 4, 64 | code: hadamard_pi_assign, doc and body (two_point.rs:74-82) | derived (recomputed) |
| 1216 | doubling constants needed after about 43 % of the steps at level I | 43 % | the strategy of the code (optimised_strategy, mul_c = 1.4) replayed in strat.py: 36 of 85 codomains (42.4 %) when m = 1, 36 of 84 (42.9 %) when m = 2, 36 of 83 (43.4 %) when m = 3; 360 of the 844 plain steps of the five level-I test vectors = 42.7 % | derived (recomputed) |

## Section 6, The change of basis at the gluing

| line | text | number | source | status |
|---|---|---|---|---|
| 1293 | of the six pairings, four are units or zero, the other two reduce to the norm equation | 6, 4, 2 | mathematics of the proof; not a measurement | matches |

## Section 7.1, Parameters

| line | text | number | source | status |
|---|---|---|---|---|
| 1326 | every class contains an ideal of norm at most about 0.90 sqrt(p) | 0.90 | cited (DLRW24, Lemma 12) | cited (not ours, left alone) |
| 1331 | the library rule amounts to log2 sqrt(p) plus 10 to 13 bits on its primes | 10, 13 | on the round-2 primes 5*2^248-1, 65*2^376-1, 27*2^500-1: e - log2(p)/2 = 10.8, 10.0, 12.6 (computed here; the tex does not name those primes) | derived (recomputed) |
| 1332 | 70 bits more than 2^(2 lambda); the formula gives 2^136, below 2^162.6 | 70, 136, 162.6 | 325.58 - 256 = 69.6; 128 + ceil(log2 256) = 136; 325.58/2 + log2 0.90 = 162.64 | derived (recomputed) |
| 1335 | sqrt(p) plus about ten bits | 10 | as above | derived (recomputed) |
| 1337 | e = ceil(log2 p / 2) + 11 = 174/264/346 | 11, 174, 264, 346 | 163 + 11, 253 + 11, 335 + 11; hd/params.rs | derived (recomputed) |
| 1339 | r = 89/134/175; margin of 11 bits | 89, 134, 175, 11 | 87 + 2, 132 + 2, 173 + 2 | derived (recomputed) |
| 1341 | SQIsignHD's own prime uses 15.3 | 15.3 | cited parameters: 142 - log2(13 * 2^126 * 3^78 - 1)/2 = 15.34 if the cited prime and e = 142 are as recalled; not checked against DLRW24 here | cited (not ours, left alone) |
| 1344 | about 2 pi^2 2^(2e)/p = 2^26.7 vectors at level I; one good norm in a few hundred | 26.7 | log2(2 pi^2) + 348 - 325.585 = 26.72; one good norm in 121 odd ones (241 of random parity) | derived (recomputed) |
| 1349 | 4 lambda/3 = 170.7/256.0/341.3 against 174/264/346; a margin of 7 would fail at level I | 170.7, 256.0, 341.3, 7 | 4 * 128/3, 4 * 192/3, 4 * 256/3; 163 + 7 = 170 < 170.7. The 4 lambda/3 rule itself: cited (NO25) | derived (recomputed) |
| 1356 | 40 Miller-Rabin rounds | 40 | code: crates/verify/src/hd/goodness.rs (40 rounds on the first 40 primes); the signer calls the same test (crates/sqisign-rs/src/compact.rs) | matches |
| 1360 | density 1/(e ln 2): one candidate in 121 at level I, of which 60 reach a primality test | 121, 60 | 174 ln 2 = 120.6; half of it = 60.3 | derived (recomputed) |
| 1370 | a nonzero c or d gives a norm of at least p/4, which exceeds 2^e at every level | p/4 | log2(p/4) = 323.6, 502.8, 666.1 against e = 174, 264, 346 | derived (recomputed) |
| 1377 | norm form of determinant p^2/16; Minkowski: N1 N2 N3 N4 <= 64 p^2 / pi^4 | p^2/16, 64/pi^4 | (2^4 * (p/4) / (pi^2/2))^2 = 64 p^2 / pi^4 | derived (recomputed) |
| 1382 | at least 0.13 X / sqrt(N1 N2) vectors per translate; covering radius at most sqrt(X/4) | 0.13 | (pi/2) (1 - 1/sqrt 2)^2 = 0.135, a disc of radius sqrt(X/2) shrunk by the covering radius | derived (recomputed) |
| 1386 | at least about 0.07 X^2 / p vectors outside L2 | 0.07 | paper-dim4/bench/short_vector_count.py: 4 c1^2 = 0.0726 when N4 <= X/2 and c1 pi^2/16 = 0.0831 when N4 > X/2, with c1 = (pi/2)(1 - 1/sqrt 2)^2 = 0.1348 | derived (recomputed) |
| 1387 | 2^18.6, 2^19.4, 2^20.1 | 18.6, 19.4, 20.1 | log2 0.07 + 2e - log2 p = 18.58, 19.41, 20.08 (given the 0.07) | derived (recomputed) |
| 1390 | about 2^9.7 good at level I, v and -v counted once | 9.7 | short_vector_count.py: 18.58 - 1 (odd) - 1 (v and -v once) - log2 120.6 = 9.66 | derived (recomputed) |
| 1393 | for O_0: N1 = N2 = 1 and N3 >= p/4 | 1, p/4 | mathematics of the example | matches |
| 1397 | 4 180 commitments: the 100 signatures of Table 5 and 4 080 more | 4180, 100, 4080 | instrument/gram.csv: 4 180 rows, signatures 0 to 4179, no restart; the first 100 have the samples and tests of the loop table (129;66, 140;65, 5;3, ...) | matches |
| 1401 | third minimum: median 2^0.39, maximum 2^3.65, that is 2^-8.6 X | 0.39, 3.65, -8.6 | instrument/minima-all.txt (gram_minima.py): nu3 median 2^0.39, max 2^3.65, N3 = 2^-8.6 X | matches |
| 1402 | tail 1.5 y^-4 over 1.4 <= y <= 4.8; between 1 588 and 6 points; median 1.53; between 1.3 and 1.6 with at least 50 points; 0.8 at the last 6 | 1.5, 1.4, 4.8, 1588, 6, 1.53, 1.3, 1.6, 50, 0.8 | minima-all.txt: range 1.41 to 4.84; y^4 Pr min 0.79, median 1.53, max 1.60; with at least 50 points 1.34 to 1.60; fitted exponent -4.15 | matches |
| 1410 | Siegel: Pr[nu1 <= x] <= pi^2 x^2 / 4 = 2.47 x^2 | 2.47 | pi^2/4 = 2.467; the theorem: cited (Sie45) | derived (recomputed) |
| 1411 | our sample gives between 1.8 x^2 and 3.0 x^2 over 2^-5 <= x <= 2^-1 | 1.8, 3.0 | minima-all.txt: Pr[nu1 <= x]/x^2 from 1.82 to 3.00 at the nine half-powers of 2 in the range | matches |
| 1413 | nu1 nu2 < 1024/(pi^4 y^2); determinant less than 32/(pi^2 y) | 1024, 32 | 64 p^2/pi^4 divided by p^2/16 = 1024/pi^4; its square root over y | derived (recomputed) |
| 1421 | a(4, 2) = 3/2 * V(4)/(V(1) zeta(4)) * V(3) zeta(2)/(V(2) zeta(3)) = 6.24 | 6.24 | the formula: cited (Kim22, Theorem 3); its value recomputed: 6.239 | derived (recomputed) |
| 1423 | Markov: at most 690 y^-4 | 690 | 6.24 * (32/pi^2)^4 = 689.6 | derived (recomputed) |
| 1429 | y = X/sqrt(p) = 2^11.2, 2^11.6, 2^12.0 | 11.2, 11.6, 12.0 | e - log2(p)/2 = 11.21, 11.62, 11.96 | derived (recomputed) |
| 1430 | probability about 2^-44, 2^-46, 2^-47 | -44, -46, -47 | minima-all.txt: 2^-44.2, 2^-45.9, 2^-47.2 from the measured 1.53 | derived (recomputed) |
| 1432 | random-lattice bound: at most 2^-35, 2^-37, 2^-38 | -35, -37, -38 | log2 690 - 4 * (11.21, 11.62, 11.96) = -35.4, -37.1, -38.4 | derived (recomputed) |
| 1437 | gives up after 4096 samples and restarts on a fresh commitment | 4096 | code: MAX_RESPONSE_TRIES = 4096 (crates/sqisign-rs/src/compact.rs:124); no restart in the 204 signatures of each pass (extra column of dim4-*-sign.csv) | matches |
| 1439 | stops with an error after 64 commitments | 64 | code: MAX_COMMITMENTS = 64 (crates/sqisign-rs/src/compact.rs) | matches |
| 1442 | (1 - 1/109)^4096 = 2^-54 | 109, 4096, -54 | 108.94 samples per good degree (dim4-I-sign.csv); 4096 log2(1 - 1/108.94) = -54.5 | derived (recomputed) |
| 1444 | a restart in fewer than one signature in 2^40 | 2^40 | 2^-44.2 + 2^-54.5 < 2^-40 | derived (recomputed) |
| 1446 | 64 in a row: about 2^-2800 | -2800 | 64 * -44.2 = -2830 (minima-all.txt) | derived (recomputed) |

## Table 1, parameters

| line | text | number | source | status |
|---|---|---|---|---|
| 1460 | p | 3 primes | cited (AAA+26); code | cited (not ours, left alone) |
| 1461 | log2 p = 325.58, 504.75, 668.09 | 325.58, 504.75, 668.09 | log2 3 + 324 = 325.585; log2 27 + 500 = 504.755; log2 17 + 664 = 668.087 | derived (recomputed) |
| 1462 | lambda = 128, 192, 256 | 128, 192, 256 | cited (AAA+26) | cited (not ours, left alone) |
| 1463 | e = 174, 264, 346 | 174, 264, 346 | as above; hd/params.rs | derived (recomputed) |
| 1464 | r = 89, 134, 175 | 89, 134, 175 | as above | derived (recomputed) |
| 1465 | half-chain steps 87, 132, 173 | 87, 132, 173 | ceil(e/2) | derived (recomputed) |
| 1467 | public key predicted 84, 130, 170 | 84, 130, 170 | 2 * ceil(ceil(log2 p)/8) + 2 = 2*41+2, 2*64+2, 2*84+2 | derived (recomputed) |
| 1468 | signature predicted 139, 212, 278 | 139, 212, 278 | ceil(1107/8), ceil(1690/8), ceil(2223/8) | derived (recomputed) |
| 1470 | public key measured 84, 130, 170 | 84, 130, 170 | raw headers of pass1/raw/dim4-{I,III,V}-verify-four.csv: "84 + 139", "130 + 212", "170 + 278 bytes, packed layout" | matches |
| 1471 | signature measured 139, 212, 278 | 139, 212, 278 | same headers | matches |
| 1472 | round-3 standard signature 200, 306, 406 | 200, 306, 406 | raw headers of dim2-{I,III,V}-verify.csv: "83 + 200", "129 + 306", "169 + 406 bytes" | matches |
| 1479 | note d: 2 ceil(log2 p) + (e - 2) + 3r + 16 bits; 137 / 210 / 276 without the hints | 16, 137, 210, 276 | (1107, 1690, 2223) - 16 = 1091, 1674, 2207 bits = 137, 210, 276 bytes | derived (recomputed) |

## Section 7.2, The wire format

| line | text | number | source | status |
|---|---|---|---|---|
| 1509 | each scalar 89 bits, not 12 bytes; q needs 172 bits | 89, 12, 172 | r = 89; e - 2 = 172 | derived (recomputed) |
| 1513 | hints in 8 bits each | 8 | code: crates/verify/src/hd/wire.rs | matches |
| 1515 | 1107, 1690, 2223 bits; 139, 212, 278 bytes; 142, 214, 280 with whole bytes | 1107, 1690, 2223, 139, 212, 278, 142, 214, 280 | 2*326+172+267+16, 2*505+262+402+16, 2*669+344+525+16; whole bytes: 82+22+36+2, 128+33+51+2, 168+44+66+2 | derived (recomputed) |
| 1524 | 97 of the 399 other pairs per signature on average over 100 signatures; median 95, between 39 and 224; all rejected with the check on | 97, 399, 100, 95, 39, 224 | instrument/hints.csv (paper_instrument, feature instrument): 20 x 20 - 1 = 399; valid_other mean 96.7, median 95, min 39, max 224; rejected_by_verifier equals valid_other on all 100 rows | matches |
| 1529 | unpacking costs a few hundred cycles | a few hundred | measured at 8747fa9 for the first audit: parse_signature at level I, median 390 cycles over 20 001 calls on CPU 13; the decoder is unchanged since | matches |
| 1532 | 137, 210 and 276 bytes without the hints | 137, 210, 276 | as in note d | derived (recomputed) |
| 1548 | flipping the top bit of a or b gave a valid signature in 2 of 12 trials; no other bit did | 2, 12 | test at head: "top bit of a or b: 2 of 12 flipped signatures verify; no other bit of a scalar does (264 flips on one signature)" (tests-head-8747fa9.log, compact_rejections). The 24-signature test of two_point.rs sees 5 valid of 24 such flips | matches |

## Section 8.1, Setup

| line | text | number | source | status |
|---|---|---|---|---|
| 1566 | 32 vCPUs, Intel Xeon at 2.80 GHz, 126 GiB, n2-standard-32 | 32, 2.80, 126 | session/machine.txt: lscpu "CPU(s): 32", "Intel(R) Xeon(R) CPU @ 2.80GHz", MemTotal 132002464 kB = 125.9 GiB, "machine type: n2-standard-32" | matches |
| 1568 | one core whose sibling was kept idle; 2026-11-15 | 2026-11-15 | machine.txt: "pinned core: 13 (siblings: 13,29)"; commands.txt: every process moved off CPUs 13 and 29; passes at 2026-11-15T00:00Z, 00:36Z, 00:51Z (the build started 2026-11-14T23:47Z) | matches |
| 1569 | rustc 1.100 | 1.100 | machine.txt: "rustc 1.100.0 (0ef1bc386 2026-11-09)" | matches |
| 1572 | dimension 2: 52 key generations, 52 signatures, 200 verifications at each level | 52, 52, 200 | raw file headers and line counts (dim2-*-keygen.csv 52 rows, -sign 52, -verify 200); src/paper.rs | matches |
| 1574 | compact, level I: 52 key generations, 100 signatures, 200 verifications of five signatures in each mode | 52, 100, 200, 5 | dim4-I-keygen.csv 52 rows, dim4-I-sign.csv 100, dim4-I-verify-*.csv "5 signatures x 40 runs" | matches |
| 1575 | levels III and V: 20 key generations, 52 signatures, mean over six signatures of the median of 15 runs | 20, 52, 6, 15 | dim4-{III,V}-keygen.csv 20 rows, -sign 52, -verify-* "6 signatures x 15 runs"; paper_report.py takes the mean of the six medians (reproduced by recompute.py) | matches |
| 1578 | the-sqisign at tag nist-v3, Release, broadwell | nist-v3 | machine.txt lines 62-66: "nist-v3", CMAKE_BUILD_TYPE=Release, SQISIGN_BUILD_TYPE=broadwell | matches |
| 1581 | iterations 300 at level I, 100 at III and V | 300, 100 | session.log: benchmark_p324_3 --iterations=300, _p500_27 and _p664_17 --iterations=100, in each pass | matches |
| 1581 | three passes; C reference keygen and sign at III and V differed by more than 1.5 % | 3, 1.5 % | numbers.json: passes 3; pass 1 against pass 2: c-III-keygen 4.7 %, c-III-sign 7.7 %, c-V-keygen 2.7 %, c-V-sign 2.9 %; every other row under 1.5 % | matches |
| 1581 | no Rust row moved by more than 1.2 % | 1.2 % | largest spread of a Rust row over the three passes: dim4-III-sign, 1.13 % | derived (recomputed) |
| 1590 | one further pass with the divisor tests: no verification row moved by more than 0.6 % | 0.6 % | pass of 2026-11-15T08:51Z at bf2153b against the session medians: largest move of a Rust verification row 0.54 % (dim4-III-verify-two-nocheck); C reference verify rows within 0.75 % (audit/pass-compare.txt) | derived (recomputed) |
| 1583 | timings taken with the bit-packed layout | packed | raw headers: "(84 + 139 bytes, packed layout)"; the harness calls hd_verify_bytes_with on the packed bytes | matches |
| 1586 | pre-chain checks cost 1.2, 3.2 and 5.5 megacycles | 1.2, 3.2, 5.5 | profile columns checks_hints_pairing + checks_goodness, medians, median of three passes: 0.286 + 0.952 = 1.24; 0.778 + 2.399 = 3.18; 1.745 + 3.767 = 5.51 (pass*/raw/dim4-*-profile-four.csv) | derived (recomputed) |
| 1587 | under 2 % of a four-point verification and under 3 % of a two-point one | 2 %, 3 % | 1.24/76.54 = 1.6 %, 3.18/194.63 = 1.6 %, 5.51/482.24 = 1.1 %; 1.24/48.05 = 2.6 %, 3.18/120.06 = 2.6 %, 5.51/291.58 = 1.9 % | derived (recomputed) |

## Section 8.2, Results: sizes

| line | text | number | source | status |
|---|---|---|---|---|
| 1596 | 223 bytes against 283; a saving of 30.5 % on the signature, of which bit packing accounts for 1.5 points | 223, 283, 30.5, 1.5 | 84 + 139; 83 + 200; (200 - 139)/200 = 30.5 %; 3/200 = 1.5 | derived (recomputed) |
| 1599 | 176-byte compressed signature | 176 | raw header of dim2-I-verify-compressed.csv: "(176 bytes)" | matches |
| 1613 | Table 2, standard: 83 + 200 = 283, 129 + 306 = 435, 169 + 406 = 575 | 9 numbers | raw headers of dim2-*-verify.csv; sums rechecked | derived (recomputed) |
| 1614 | Table 2, compressed: 83 + 176 = 259, 129 + 269 = 398, 169 + 353 = 522 | 9 numbers | raw headers of dim2-*-verify-compressed.csv (176, 269, 353); sums rechecked | derived (recomputed) |
| 1615 | Table 2, compact: 84 + 139 = 223, 130 + 212 = 342, 170 + 278 = 448 | 9 numbers | raw headers of dim4-*-verify-four.csv; sums rechecked | derived (recomputed) |

## Section 8.2, Results: time

| line | text | number | source | status |
|---|---|---|---|---|
| 1623 | level I: verifier 6.8 times dimension 2; signing 3.6 times | 6.8, 3.6 | 76.540/11.185 = 6.84; 360.544/100.023 = 3.60 | derived (recomputed) |
| 1625 | levels III and V: verifier 6.5 and 6.0 times | 6.5, 6.0 | 194.632/29.798 = 6.53; 482.237/80.017 = 6.03 | derived (recomputed) |
| 1626 | two-point verifier 4.0 and 3.6 times | 4.0, 3.6 | 120.059/29.798 = 4.03; 291.583/80.017 = 3.64 | derived (recomputed) |
| 1627 | signing 2.4 and 2.2 times | 2.4, 2.2 | 913.887/375.368 = 2.43; 1714.477/769.668 = 2.23 | derived (recomputed) |
| 1629 | middle half of the signatures: 398 to 1471 at level III, 869 to 3392 at level V | 398, 1471, 869, 3392 | quartiles of the 52 signatures of each pass (exclusive method), median of the three passes: III 398.3 and 1470.9; V 868.7 and 3391.95. The 3391 was pass 1 alone | derived (recomputed) |
| 1632 | key generation 35.3 against 32.9 megacycles | 35.3, 32.9 | dim4-I-keygen 35.269; dim2-I-keygen 32.938 | matches |
| 1633 | Rust dimension 2 over C at level I: 1.02 for verification | 1.02 | 11.185/10.982 = 1.019 | derived (recomputed) |
| 1634 | 1.06 for key generation, 1.09 for signing | 1.06, 1.09 | 32.938/31.185 = 1.056; 100.023/91.960 = 1.088. The one-decimal cells 32.9/31.2 give 1.054, which is where 1.05 came from; the unrounded medians are used for every ratio | derived (recomputed) |
| 1635 | compressed recovery adds 19 % to a verification | 19 % | 13.266/11.185 = 1.186 | derived (recomputed) |
| 1652 | Table 3, C reference I: 31.2, 92.0, 11.0 | 31.2, 92.0, 11.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 31.185, 91.960, 10.982 | matches |
| 1653 | Table 3, C reference III: 106.2, 323.1, 29.0 | 106.2, 323.1, 29.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 106.248, 323.050, 29.048 | matches |
| 1654 | Table 3, C reference V: 216.9, 679.6, 75.3 | 216.9, 679.6, 75.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 216.893, 679.626, 75.253 | matches |
| 1655 | Table 3, Rust dimension 2, I: 32.9, 100.0, 11.2 | 32.9, 100.0, 11.2 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 32.938, 100.023, 11.185 | matches |
| 1656 | Table 3, Rust dimension 2, III: 109.0, 375.4, 29.8 | 109.0, 375.4, 29.8 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 108.955, 375.368, 29.798 | matches |
| 1657 | Table 3, Rust dimension 2, V: 248.6, 769.7, 80.0 | 248.6, 769.7, 80.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 248.618, 769.668, 80.017 | matches |
| 1658 | Table 3, compressed verify, I: 13.3 | 13.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 13.266 | matches |
| 1659 | Table 3, compact I: 35.3, 360.5, 76.5; x dim. 2: 6.8 | 35.3, 360.5, 76.5; x dim. 2: 6.8 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 35.269, 360.544, 76.540; 6.84 | derived (recomputed) |
| 1660 | Table 3, two points I: 48.0; 4.3 | 48.0; 4.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 48.048; 4.30 | derived (recomputed) |
| 1661 | Table 3, compact III: 126.3, 913.9, 194.6; 6.5 | 126.3, 913.9, 194.6; 6.5 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 126.298, 913.887, 194.632; 6.53 | derived (recomputed) |
| 1662 | Table 3, two points III: 120.1; 4.0 | 120.1; 4.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 120.059; 4.03 | derived (recomputed) |
| 1663 | Table 3, compact V: 249.8, 1714.5, 482.2; 6.0 | 249.8, 1714.5, 482.2; 6.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 249.781, 1714.477, 482.237; 6.03 | derived (recomputed) |
| 1664 | Table 3, two points V: 291.6; 3.6 | 291.6; 3.6 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 291.583; 3.64 | derived (recomputed) |

## Figure 2, the trade-off

| line | text | number | source | status |
|---|---|---|---|---|
| 1674 | axes: x 120 to 240, y 0 to 90; ticks every 20 and every 30 | 120, 240, 0, 90 | layout; every point lies inside (largest y is 76.5) | matches |
| 1682 | compact: (139, 76.5) | 139, 76.5 | Table 2 and Table 3: dim4-I-verify-four 76.540 | matches |
| 1683 | compressed (176, 13.3); standard (200, 11.2) | 176, 13.3, 200, 11.2 | dim2-I-verify-compressed 13.266; dim2-I-verify 11.185 | matches |
| 1684 | two points: (139, 48.0) | 139, 48.0 | dim4-I-verify-two 48.048 | matches |
| 1685 | label positions (145,78), (145,49), (168,17), (203,15) | 8 numbers | layout only; each label sits next to its point | matches |

## Table 4, wall-clock milliseconds at level I

| line | text | number | source | status |
|---|---|---|---|---|
| 1713 | dimension 2 key generation 11.8 | 11.8 | ns column of the same raw files, median of three passes: 11.764 | matches |
| 1714 | dimension 2 signing 35.7 | 35.7 | ns column of the same raw files, median of three passes: 35.724 | matches |
| 1715 | dimension 2 verification 4.0 | 4.0 | ns column of the same raw files, median of three passes: 3.995 | matches |
| 1716 | compressed verification 4.7 | 4.7 | ns column of the same raw files, median of three passes: 4.738 | matches |
| 1718 | compact key generation 12.6 | 12.6 | ns column of the same raw files, median of three passes: 12.597 | matches |
| 1719 | compact signing 128.8 | 128.8 | ns column of the same raw files, median of three passes: 128.768 | matches |
| 1720 | compact verification 27.3 | 27.3 | ns column of the same raw files, median of three passes: 27.337 | matches |
| 1721 | two points 17.2 | 17.2 | ns column of the same raw files, median of three passes: 17.161 | matches |

## Section 8.2, the good-degree loop and Table 5

| line | text | number | source | status |
|---|---|---|---|---|
| 1729 | one candidate in two reaches a Miller-Rabin test | 1/2 | extra column of dim4-I-sign.csv: 54.87 tests for 108.94 samples per signature = 0.504 | derived (recomputed) |
| 1731 | median 95 samples, mean 109, range 1 to 446 | 95, 109, 1, 446 | extra column of dim4-I-sign.csv (the same 100 signatures in every pass): mean 108.94, median 95, min 1, max 446 | matches |
| 1745 | signatures measured: 100 | 100 | dim4-I-sign.csv, 100 rows | matches |
| 1746 | samples mean 109; expected 121 | 109, 121 | 108.94; 174 ln 2 = 120.6 | matches |
| 1747 | random parity: 241 | 241 | 2 * 174 ln 2 = 241.2 | derived (recomputed) |
| 1748 | median (min, max): 95 (1, 446) | 95, 1, 446 | as above | matches |
| 1749 | primality tests mean 55; expected 60 | 55, 60 | 54.87; 60.3 | matches |
| 1750 | signing 361 (46, 1529) megacycles | 361, 46, 1529 | median of the three passes of each statistic: median 360.54, min 45.51, max 1528.9. The 1546 was the maximum of pass 1 alone (passes: 1546.3, 1528.9, 1524.1) | derived (recomputed) |
| 1751 | dimension-2 signing 100 | 100 | dim2-I-sign 100.023 | matches |

## Section 8.2, heap

| line | text | number | source | status |
|---|---|---|---|---|
| 1757 | 1655 KiB (1.62 MiB) with four points, 1631 KiB with two | 1655, 1.62, 1631 | session/heap.csv, mean of five signatures: 1 694 442 bytes = 1654.7 KiB = 1.616 MiB; 1 669 866 bytes = 1630.7 KiB | derived (recomputed) |
| 1758 | 3443 and 6160 KiB at levels III and V; dimension 2 allocates nothing | 3443, 6160, 0 | heap.csv: 3 525 398 bytes = 3442.8 KiB; 6 307 749 bytes = 6159.9 KiB; dim2-*-verify 0 bytes, 0 allocations | derived (recomputed) |

## Section 8.3, Scope

| line | text | number | source | status |
|---|---|---|---|---|
| 1773 | 24 signatures at level I and 5 each at levels III and V | 24, 5, 5 | tests at head: "level I, 24 signatures: 24 accepted, 336 tampered variants rejected, ... 4022 plain steps"; "level III: 5 vectors"; "level V: 5 vectors" (tests-head-8747fa9.log) | matches |

## Section 9, Conclusion

| line | text | number | source | status |
|---|---|---|---|---|
| 1803 | 139 bytes against 200; 223 bytes | 139, 200, 223 | as in the abstract | matches |
| 1808 | two points: 4.3, 4.0 and 3.6 times dimension 2 | 4.3, 4.0, 3.6 | 4.30, 4.03, 3.64 (as in Table 3) | derived (recomputed) |
| 1809 | four points: 6.8, 6.5 and 6.0 times | 6.8, 6.5, 6.0 | 6.84, 6.53, 6.03 | derived (recomputed) |
| 1811 | signing 3.6, 2.4 and 2.2 times | 3.6, 2.4, 2.2 | 3.60, 2.43, 2.23 | derived (recomputed) |
| 1813 | 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | as in the introduction | derived (recomputed) |
| 1815 | a few tens of milliseconds at level I; 43 to 172 milliseconds at levels III and V | tens, 43, 172 | ns columns: level I 17.2 and 27.3 ms; level III 42.9 and 69.5; level V 104.1 and 172.2 (recompute.out) | matches |
| 1824 | probability about 2^-44 at level I | 2^-44 | as in Section 7.1: 2^-44.2 from the measured tail | derived (recomputed) |
| 1831 | 61 bytes saved at level I | 61 | 200 - 139 | derived (recomputed) |
| 1834 | MIKE on the session core: 228 steps in 20.9 megacycles, median of 200 runs, about 0.09 per step | 228, 20.9, 200, 0.09 | session/mike-step.csv in the results folder of 2026-11-15T08:12Z (paper-dim4/bench/mike-step, MIKE-I of mike_rs 0d79121, CPU 13): shared_secret median 20.932 megacycles; 20.932/228 = 0.0918. E = 228 in mike_rs public_parameters | matches |
| 1836 | about 0.28 for ours (48.0 megacycles over 2 * 87 steps) | 0.28, 48.0, 87 | 48.048/174 = 0.276, same core, same build settings | derived (recomputed) |
| 1838 | MIKE-I prime of 318 bits | 318 | 633 * 2^308 - 1 has 318 bits (mike_rs public_parameters) | derived (recomputed) |

## AI disclosure

| line | text | number | source | status |
|---|---|---|---|---|
| 1853 | the disclosure names Claude Opus 5.5 | 5.5 | an authors' statement; marked with the one TODO(author) that remains | needs an author (item 11) |

## Appendix A, proofs

| line | text | number | source | status |
|---|---|---|---|---|
| 1916 | f = 324, 500, 664 against e + 2 = 176, 266, 348 | 324, 500, 664, 176, 266, 348 | the exponents of 2 in p + 1; 174 + 2, 264 + 2, 346 + 2 | derived (recomputed) |
| 2330 | pi has trace 4; fixes 0, 5, 10, 15 | 4, 0, 5, 10, 15 | code: PI_FIXED = [0, 5, 10, 15] (two_point.rs:58) | matches |
| 2332 | the remaining numbers of Appendix A (orders 2^(4e), 16 fixed points, eigenvalues, denominators 4, sixteen terms) | several | mathematics of the proofs; not measurements and not checked by the audit beyond the two rows above | matches |

## Appendix B, Conformance of the implementation

| line | text | number | source | status |
|---|---|---|---|---|
| 2397 | m exceeds 40 when a2 = 0; every m > 40 rejected; 1 <= m <= 40 | 40 | code: MAX_GLUING_M = 40 (crates/verify/src/hd/self_contained.rs:67, test at :414) | matches |
| 2402 | 40 Miller-Rabin rounds on the first 40 primes | 40, 40 | code: crates/verify/src/hd/goodness.rs | matches |
| 2488 | four evaluations through the gluing, two projective comparisons of points, six comparisons of coordinates | 4, 2, 6 | code: gluing_images_paired (two projective_eq), null_fixed_by_pi (six), two_point.rs; four evaluations in self_contained.rs (evaluate4) | matches |
| 2500 | the same codomain at all 844 plain steps of five level-I vectors; 70 tampered variants | 844, 5, 70 | test at head: "level I: 5 vectors accepted in the three modes with no fall-back; 70 tampered variants rejected with the same reason in the three modes; ... 844 plain steps with the same codomain in both two-point modes" | matches |
| 2504 | levels III and V: 1 274 and 1 706 plain steps of five vectors each; 75 and 70 tampered variants | 1274, 1706, 5, 75, 70 | tests at head: "level III: 5 vectors ... 75 tampered variants ... 1274 plain steps"; "level V: 5 vectors ... 70 tampered variants ... 1706 plain steps" | matches |

## Appendix C, Beyond two-point evaluation

| line | text | number | source | status |
|---|---|---|---|---|
| 2525 | MIKE: the null point vanishes at the six coordinates that [i] negates; ten coordinates kept | 6, 10 | cited (BDD+26, Definition G.1, Proposition G.3, Remark E.2) | cited (not ours, left alone) |
| 2529 | doubling about 1.6 times faster on the quotient | 1.6 | cited (BDD+26, Remark E.2) | cited (not ours, left alone) |
| 2535 | up to about 1 - 1/1.6 = 38 % | 38 % | 1 - 1/1.6 = 37.5 % | derived (recomputed) |

184 rows: 90 derived, 81 matches, 12 cited, 1 needs an author.

## Statements checked against the code

Every sentence the first audit listed as not exact was either made true in the code or rewritten to say what the code does.

- The four divisions with no zero test now test their divisors and reject (commit 6302b7a): the dimension-2 chain (doubling constants, gluing codomain, gluing image, plain step), the special image of the dimension-4 gluing, the weights the image check applies, the inversion of e_4(T1, T2). One test per site (zero_divisor_tests). Item (j) of the appendix is true as written and now says what it covers.
- Section 4.2 and Algorithm 1: no basis of E_chl[2^f]; the order of the checks is the code's (e(P, Q), then a and b both even, then e(R, S)); the middle comparison is described; the list of rejects is complete.
- Section 5: a vanishing coordinate falls back for the codomain and rejects for the doubling constants.
- Section 7.1: the signer's full test; the stop after 64 commitments with its probability.
- Section 7.2: the hints leave a single encoding, the top bit aside.
- Appendix B: the order in (e); the order of i_0 in (h); the middle comparison as a new item; the counts of the check at the gluing.

## Cross-machine comparisons

None is left. The conclusion sets MIKE's step cost against ours with both measured on the same core: 20.9 megacycles for 228 steps, 0.09 per step, against 0.28. Line 319 quotes MIKE's published 16.0 megacycles on its own and compares it with nothing.

## Sanity

- 99 old values and phrases (those of the first audit, and 3.66, the 1.6 tail, "over two hundred", "always ends", 17 790, "0.07 megacycles", "different processor", 1.9 x^2, 2.8 x^2, "not re-measured", the E_chl[2^f] basis, 1 600) were searched outside comment lines (oldvalues.out). Eight raw hits, all inside other numbers or unrelated ("3.3" in 13.3, "6.9" in 216.9, "3.1" in 323.1, "9.8" in 29.8 and 249.8, "107" in 1107, "8.4" in "Section 8.4.2", "four times a kernel lift").
- latexmk -pdf: one error, the missing db.png; no undefined reference or citation; 34 pages.
- One TODO(author) in the file, at the AI disclosure.

## Measured values that differ from what the text said

- Hint pairs: the text said "over two hundred hint pairs per signature". Measured: 97 on average of the 399 other pairs, median 95, between 39 and 224. The conclusion drawn from it (the verifier must check the hints) does not change.
- The tail constant: 1.5 in place of 1.6 (median 1.53 over 4 180 commitments). The probabilities 2^-44, 2^-46, 2^-47 and the 2^40 do not change.
- MIKE per step on our machine: 0.09, where its published figure on its own machine gave 0.07. The sentence's conclusion does not change.

## Stale

None.

## Needs an author

1. The model name in the AI disclosure (line 1852 marks it). Nothing else.

## Cannot source

None.
