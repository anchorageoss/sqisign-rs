# AUDIT: every number of paper-2026-11-14-standalone.tex against the quiet session and the code

Audited 2026-11-15. File: paper-dim4/paper-2026-11-14-standalone.tex (2496 lines after this audit). The original main_standalone.tex is untouched (sha256 ce4b629c8449064bc4ac35921952077651b4c5d26852b255607cf05356b87661). Line numbers are those of the dated tex after the audit.

## Sources

- Session: ~/work/results/20261114T224715Z-SECENG-2636-clean/session/ (three passes, 2026-11-15, CPU 13 of an n2-standard-32, code at 68a0017). Raw files: pass1, pass2, pass3 /raw/*.csv (one line per call: cycles, nanoseconds, extra), pass*/c-ref-{I,III,V}.txt (the C reference's own benchmark output), heap.csv, machine.txt, session.log. The same records are in the branch under paper-dim4/bench/sessions/2026-11-15-68a0017/.
- Recomputation: audit/recompute.py in that results folder reads only the raw files. It reproduces all 36 rows of session/numbers.json and prints every ratio, share, quartile, loop statistic and heap mean used below (audit/recompute.out). A session value is the median of the three passes' values, as the Setup paragraph says. Every ratio is taken from the unrounded medians; where the one-decimal cells of Table 3 give a different last digit, the row says so (one case, key generation over the C reference).
- Tests: cargo test --release -p sqisign-verify -p sqisign-rs --features sqisign-rs/compact at the head of agent/dim4-two-point (8747fa9), run for this audit: every test passes; output in audit/tests-head-8747fa9.log.
- Code: file and line at 8747fa9 under crates/. The crate code at 8747fa9 is the code the session measured at 68a0017 (the commit between them adds session records and a harness script only).
- audit/strat.py replays the code's strategy to count the steps after which a doubling follows.

Status is one of: matches / stale (fixed) / derived (recomputed) / cited (not ours, left alone) / cannot source. A row that has to go to an author says "needs an author" and why; those rows are repeated in the list at the end.


## Abstract (plain-text copy and typeset copy)

| line | text | number | source | status |
|---|---|---|---|---|
| 48 | plain-text: automorphism of order 4 | 4 | mathematics of the paper (Lemma J(i)); not a measurement | matches |
| 48 | plain-text: ten coordinates instead of sixteen | 10, 16 | code: ORBIT_REPS has 10 entries, THETA_DIM4_N = 16 (crates/verify/src/hd/two_point.rs:62) | matches |
| 48 | plain-text: cuts the verifier's cost by 37 to 40 % | 37, 40 | 1 - two/four from the unrounded session medians: 37.2, 38.3, 39.5 % at I, III, V (recompute.py) | derived (recomputed) |
| 48 | plain-text: 139 bytes against 200; key and signature 223 | 139, 200, 223 | raw headers: "84 + 139 bytes, packed layout" (pass1/raw/dim4-I-verify-four.csv), "83 + 200 bytes" (dim2-I-verify.csv); 84 + 139 = 223 | matches |
| 48 | plain-text: verification costs 4.3 times that of the dimension-2 scheme (was "four times") | 4.3 | two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (recompute.py). The paper of Nov 4 had 4.0, which "four" described; 4.3 is now written out | stale (fixed) |
| 53 | typeset: automorphism of order 4 | 4 | mathematics of the paper (Lemma J(i)); not a measurement | matches |
| 53 | typeset: ten coordinates instead of sixteen | 10, 16 | code: ORBIT_REPS has 10 entries, THETA_DIM4_N = 16 (crates/verify/src/hd/two_point.rs:62) | matches |
| 53 | typeset: cuts the verifier's cost by 37 to 40 % | 37, 40 | 1 - two/four from the unrounded session medians: 37.2, 38.3, 39.5 % at I, III, V (recompute.py) | derived (recomputed) |
| 55 | typeset: 139 bytes against 200; key and signature 223 | 139, 200, 223 | raw headers: "84 + 139 bytes, packed layout" (pass1/raw/dim4-I-verify-four.csv), "83 + 200 bytes" (dim2-I-verify.csv); 84 + 139 = 223 | matches |
| 55 | typeset: verification costs 4.3 times that of the dimension-2 scheme (was "four times") | 4.3 | two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (recompute.py). The paper of Nov 4 had 4.0, which "four" described; 4.3 is now written out | stale (fixed) |

## Introduction

| line | text | number | source | status |
|---|---|---|---|---|
| 68 | 2022 (the SIDH attacks) | 2022 | cited (CD23, MMP+23, Rob23) | cited (not ours, left alone) |
| 86 | 2026; p^(1/3+o(1)); September 2026 | 2026, 1/3 | cited (Wes26, AAA+26) | cited (not ours, left alone) |
| 89 | 3*2^324 - 1, 27*2^500 - 1, 17*2^664 - 1 | 3 primes | cited (AAA+26); the same primes in the code (crates/verify/src/params) | cited (not ours, left alone) |
| 109 | contribution list: 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | two/four from the unrounded medians: 0.6278, 0.6169, 0.6046; from the one-decimal cells: 48.0/76.5, 120.1/194.6, 291.6/482.2 give the same | derived (recomputed) |
| 124 | about once in 2^44 at level I | 2^44 | computed from the tail constant 1.6 of the Gram-matrix measurement, which was not re-measured; arithmetic rechecked: 1.6 * (2^11.21)^-4 = 2^-44.2 | needs an author (rests on a figure not re-measured) |
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
| 526 | 87 steps from each end at level I; sixteen coordinates instead of four | 87, 16, 4 | ceil(174/2) = 87 | derived (recomputed) |
| 527 | costs 6.8 times the dimension-2 verifier (76.5 against 11.2 megacycles), and 7.0 times the reference | 6.8, 76.5, 11.2, 7.0 | medians of three passes: dim4-I-verify-four 76.540, dim2-I-verify 11.185, c-I-verify 10.982 (pass*/raw/*.csv, pass*/c-ref-I.txt); 76.540/11.185 = 6.84; 76.540/10.982 = 6.97 | derived (recomputed) |
| 531 | profile: 60 % to doubling and pushing the kernel points | 60 % | pass*/raw/dim4-I-profile-four.csv, slot4 / total: 59.7, 59.6, 59.8 % in the three passes | derived (recomputed) |
| 533 | 15 % to the codomains and their constants | 15 % | same files, slot5 / total: 15.0, 15.1, 15.1 % | derived (recomputed) |

## Section 5, Two-point verification

| line | text | number | source | status |
|---|---|---|---|---|
| 579 | pi exchanges index bits 0 and 2 and bits 1 and 3 of the sixteen coordinates | 0, 2, 1, 3, 16 | code: PI = [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15] (two_point.rs:44) | matches |
| 1014 | (e, r) = (174, 89), (264, 134), (346, 175) | 174, 89, 264, 134, 346, 175 | e = ceil(log2 p / 2) + 11 with log2 p = 325.58, 504.75, 668.09; r = ceil(e/2) + 2; the same in crates/verify/src/hd/params.rs | derived (recomputed) |
| 1124 | 10-dimensional space of vectors fixed by pi; ten orbits | 10 | 4 fixed indices + 6 two-cycles (PI_FIXED, PI_PAIRS, two_point.rs:47-58) | matches |
| 1180 | 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | as in the introduction | derived (recomputed) |
| 1194 | tree on the ten orbits, 25 multiplications (was: a fixed path, about 26) | 10, 9, 1, 25 | code: along_tree, TREE_ORBITS (two_point.rs:140-171): 8 + 8 + 7 + 2 = 25 multiplications per vector, called once for the dual null point and once for its inverse (two_point.rs:199-200) | stale (fixed) |
| 1195 | Hadamard transform of a pi-fixed vector: 32 additions and 4 doublings instead of 64 additions (was: 34 additions instead of 64) | 32, 4, 64 | code: hadamard_pi_assign, doc and body (two_point.rs:74-82) | stale (fixed) |
| 1201 | doubling constants needed after about 43 % of the steps at level I (was 42) | 43 % | the strategy of the code (optimised_strategy, mul_c = 1.4) replayed in strat.py: 36 of 85 codomains (42.4 %) when m = 1, 36 of 84 (42.9 %) when m = 2, 36 of 83 (43.4 %) when m = 3; 360 of the 844 plain steps of the five level-I test vectors = 42.7 % | stale (fixed) |

## Section 6, The change of basis at the gluing

| line | text | number | source | status |
|---|---|---|---|---|
| 1278 | of the six pairings, four are units or zero, the other two reduce to the norm equation | 6, 4, 2 | mathematics of the proof; not a measurement | matches |

## Section 7.1, Parameters

| line | text | number | source | status |
|---|---|---|---|---|
| 1311 | every class contains an ideal of norm at most about 0.90 sqrt(p) | 0.90 | cited (DLRW24, Lemma 12) | cited (not ours, left alone) |
| 1316 | the library rule amounts to log2 sqrt(p) plus 10 to 13 bits on its primes | 10, 13 | on the round-2 primes 5*2^248-1, 65*2^376-1, 27*2^500-1: e - log2(p)/2 = 10.8, 10.0, 12.6 (computed here; the tex does not name those primes) | derived (recomputed) |
| 1317 | 70 bits more than 2^(2 lambda); the formula gives 2^136, below 2^162.6 | 70, 136, 162.6 | 325.58 - 256 = 69.6; 128 + ceil(log2 256) = 136; 325.58/2 + log2 0.90 = 162.64 | derived (recomputed) |
| 1320 | sqrt(p) plus about ten bits | 10 | as above | derived (recomputed) |
| 1322 | e = ceil(log2 p / 2) + 11 = 174/264/346 | 11, 174, 264, 346 | 163 + 11, 253 + 11, 335 + 11; hd/params.rs | derived (recomputed) |
| 1324 | r = 89/134/175; margin of 11 bits | 89, 134, 175, 11 | 87 + 2, 132 + 2, 173 + 2 | derived (recomputed) |
| 1326 | SQIsignHD's own prime uses 15.3 | 15.3 | cited parameters: 142 - log2(13 * 2^126 * 3^78 - 1)/2 = 15.34 if the cited prime and e = 142 are as recalled; not checked against DLRW24 here | cited (not ours, left alone) |
| 1329 | about 2 pi^2 2^(2e)/p = 2^26.7 vectors at level I; one good norm in a few hundred | 26.7 | log2(2 pi^2) + 348 - 325.585 = 26.72; one good norm in 121 odd ones (241 of random parity) | derived (recomputed) |
| 1334 | 4 lambda/3 = 170.7/256.0/341.3 against 174/264/346; a margin of 7 would fail at level I | 170.7, 256.0, 341.3, 7 | 4 * 128/3, 4 * 192/3, 4 * 256/3; 163 + 7 = 170 < 170.7. The 4 lambda/3 rule itself: cited (NO25) | derived (recomputed) |
| 1341 | 40 Miller-Rabin rounds | 40 | code: crates/verify/src/hd/goodness.rs (40 rounds on the first 40 primes); the signer calls the same test (crates/sqisign-rs/src/compact.rs) | matches |
| 1344 | density 1/(e ln 2): one candidate in 121 at level I, of which 60 reach a primality test | 121, 60 | 174 ln 2 = 120.6; half of it = 60.3 | derived (recomputed) |
| 1354 | a nonzero c or d gives a norm of at least p/4, which exceeds 2^e at every level | p/4 | log2(p/4) = 323.6, 502.8, 666.1 against e = 174, 264, 346 | derived (recomputed) |
| 1361 | norm form of determinant p^2/16; Minkowski: N1 N2 N3 N4 <= 64 p^2 / pi^4 | p^2/16, 64/pi^4 | (2^4 * (p/4) / (pi^2/2))^2 = 64 p^2 / pi^4 | derived (recomputed) |
| 1366 | at least 0.13 X / sqrt(N1 N2) vectors per translate; covering radius at most sqrt(X/4) | 0.13 | (pi/2) (1 - 1/sqrt 2)^2 = 0.135, a disc of radius sqrt(X/2) shrunk by the covering radius | derived (recomputed) |
| 1370 | at least about 0.07 X^2 / p vectors outside L2 | 0.07 | the count of translates is not written out in the tex and was not re-derived in the audit | cannot source |
| 1371 | 2^18.6, 2^19.4, 2^20.1 | 18.6, 19.4, 20.1 | log2 0.07 + 2e - log2 p = 18.58, 19.41, 20.08 (given the 0.07) | derived (recomputed) |
| 1374 | about 2^9.7 good at level I; chance that none is: exp(-2^9.7) | 9.7 | 18.58 - 1 (odd) - 1 (v and -v counted once) - log2 120.6 = 9.67. Without the identification of v and -v the count is 2^10.7; the tex does not say which it means | derived (recomputed) |
| 1377 | for O_0: N1 = N2 = 1 and N3 >= p/4 | 1, p/4 | mathematics of the example | matches |
| 1382 | 4 180 commitments at level I | 4180 | earlier version of the crate; no instrumentation and no log at head | needs an author (not re-measured) |
| 1385 | third minimum: median 2^0.39, maximum 2^3.66, that is 2^-8.6 X | 0.39, 3.66, -8.6 | not re-measured. The conversion is consistent: 3.66 - 1 - 11.21 = -8.55 | needs an author (not re-measured) |
| 1386 | tail 1.6 y^-4 over 1.4 <= y <= 5.7; between 1 600 and 6 points; y^4 Pr between 1.5 and 1.8 | 1.6, 1.4, 5.7, 1600, 6, 1.5, 1.8 | not re-measured | needs an author (not re-measured) |
| 1393 | Siegel: Pr[nu1 <= x] <= pi^2 x^2 / 4 = 2.47 x^2 | 2.47 | pi^2/4 = 2.467; the theorem: cited (Sie45) | derived (recomputed) |
| 1394 | our sample gives between 1.9 x^2 and 2.8 x^2 over 2^-5 <= x <= 2^-1 | 1.9, 2.8, 2^-5, 2^-1 | not re-measured | needs an author (not re-measured) |
| 1396 | nu1 nu2 < 1024/(pi^4 y^2); determinant less than 32/(pi^2 y) | 1024, 32 | 64 p^2/pi^4 divided by p^2/16 = 1024/pi^4; its square root over y | derived (recomputed) |
| 1404 | a(4, 2) = 3/2 * V(4)/(V(1) zeta(4)) * V(3) zeta(2)/(V(2) zeta(3)) = 6.24 | 6.24 | the formula: cited (Kim22, Theorem 3); its value recomputed: 6.239 | derived (recomputed) |
| 1406 | Markov: at most 690 y^-4 | 690 | 6.24 * (32/pi^2)^4 = 689.6 | derived (recomputed) |
| 1412 | y = X/sqrt(p) = 2^11.2, 2^11.6, 2^12.0 | 11.2, 11.6, 12.0 | e - log2(p)/2 = 11.21, 11.62, 11.96 | derived (recomputed) |
| 1413 | probability about 2^-44, 2^-46, 2^-47 | -44, -46, -47 | arithmetic rechecked: log2 1.6 - 4 * (11.21, 11.62, 11.96) = -44.2, -45.8, -47.2; the 1.6 is the figure not re-measured | needs an author (rests on a figure not re-measured) |
| 1415 | random-lattice bound: at most 2^-35, 2^-37, 2^-38 | -35, -37, -38 | log2 690 - 4 * (11.21, 11.62, 11.96) = -35.4, -37.1, -38.4 | derived (recomputed) |
| 1420 | gives up after 4096 samples and restarts on a fresh commitment | 4096 | code: MAX_RESPONSE_TRIES = 4096 (crates/sqisign-rs/src/compact.rs:124); no restart in the 204 signatures of each pass (extra column of dim4-*-sign.csv) | matches |
| 1423 | fewer than one signature in 2^40 | 2^40 | follows from the 2^-44 above, with margin | needs an author (rests on a figure not re-measured) |

## Table 1, parameters

| line | text | number | source | status |
|---|---|---|---|---|
| 1437 | p | 3 primes | cited (AAA+26); code | cited (not ours, left alone) |
| 1438 | log2 p = 325.58, 504.75, 668.09 | 325.58, 504.75, 668.09 | log2 3 + 324 = 325.585; log2 27 + 500 = 504.755; log2 17 + 664 = 668.087 | derived (recomputed) |
| 1439 | lambda = 128, 192, 256 | 128, 192, 256 | cited (AAA+26) | cited (not ours, left alone) |
| 1440 | e = 174, 264, 346 | 174, 264, 346 | as above; hd/params.rs | derived (recomputed) |
| 1441 | r = 89, 134, 175 | 89, 134, 175 | as above | derived (recomputed) |
| 1442 | half-chain steps 87, 132, 173 | 87, 132, 173 | ceil(e/2) | derived (recomputed) |
| 1444 | public key predicted 84, 130, 170 | 84, 130, 170 | 2 * ceil(ceil(log2 p)/8) + 2 = 2*41+2, 2*64+2, 2*84+2 | derived (recomputed) |
| 1445 | signature predicted 139, 212, 278 | 139, 212, 278 | ceil(1107/8), ceil(1690/8), ceil(2223/8) | derived (recomputed) |
| 1447 | public key measured 84, 130, 170 | 84, 130, 170 | raw headers of pass1/raw/dim4-{I,III,V}-verify-four.csv: "84 + 139", "130 + 212", "170 + 278 bytes, packed layout" | matches |
| 1448 | signature measured 139, 212, 278 | 139, 212, 278 | same headers | matches |
| 1449 | round-3 standard signature 200, 306, 406 | 200, 306, 406 | raw headers of dim2-{I,III,V}-verify.csv: "83 + 200", "129 + 306", "169 + 406 bytes" | matches |
| 1456 | note d: 2 ceil(log2 p) + (e - 2) + 3r + 16 bits; 137 / 210 / 276 without the hints | 16, 137, 210, 276 | (1107, 1690, 2223) - 16 = 1091, 1674, 2207 bits = 137, 210, 276 bytes | derived (recomputed) |

## Section 7.2, The wire format

| line | text | number | source | status |
|---|---|---|---|---|
| 1486 | each scalar 89 bits, not 12 bytes; q needs 172 bits | 89, 12, 172 | r = 89; e - 2 = 172 | derived (recomputed) |
| 1490 | hints in 8 bits each | 8 | code: crates/verify/src/hd/wire.rs; see the code check below | matches |
| 1492 | 1107, 1690, 2223 bits; 139, 212, 278 bytes; 142, 214, 280 with whole bytes | 1107, 1690, 2223, 139, 212, 278, 142, 214, 280 | 2*326+172+267+16, 2*505+262+402+16, 2*669+344+525+16; whole bytes: 82+22+36+2, 128+33+51+2, 168+44+66+2 | derived (recomputed) |
| 1502 | over two hundred hint pairs per signature | 200+ | earlier version of the crate; no test at head counts them | needs an author (not re-measured) |
| 1504 | unpacking costs a few hundred cycles | a few hundred | measured for the audit at head: parse_signature at level I, median 390 cycles over 20 001 calls on CPU 13 (scratch example, not committed; the box was not quieted). The profile's "decode" of 0.020 megacycles also decodes the key and hashes the challenge | matches |
| 1507 | 137, 210 and 276 bytes without the hints | 137, 210, 276 | as in note d | derived (recomputed) |
| 1523 | flipping the top bit of a or b gave a valid signature in 2 of 12 trials; no other bit did | 2, 12 | test at head: "top bit of a or b: 2 of 12 flipped signatures verify; no other bit of a scalar does (264 flips on one signature)" (tests-head-8747fa9.log, compact_rejections). The 24-signature test of two_point.rs sees 5 valid of 24 such flips | matches |

## Section 8.1, Setup

| line | text | number | source | status |
|---|---|---|---|---|
| 1541 | 32 vCPUs, Intel Xeon at 2.80 GHz, 126 GiB, n2-standard-32 | 32, 2.80, 126 | session/machine.txt: lscpu "CPU(s): 32", "Intel(R) Xeon(R) CPU @ 2.80GHz", MemTotal 132002464 kB = 125.9 GiB, "machine type: n2-standard-32" | matches |
| 1543 | one core whose sibling was kept idle; 2026-11-15 | 2026-11-15 | machine.txt: "pinned core: 13 (siblings: 13,29)"; commands.txt: every process moved off CPUs 13 and 29; passes at 2026-11-15T00:00Z, 00:36Z, 00:51Z (the build started 2026-11-14T23:47Z) | matches |
| 1544 | rustc 1.100 | 1.100 | machine.txt: "rustc 1.100.0 (0ef1bc386 2026-11-09)" | matches |
| 1547 | dimension 2: 52 key generations, 52 signatures, 200 verifications at each level | 52, 52, 200 | raw file headers and line counts (dim2-*-keygen.csv 52 rows, -sign 52, -verify 200); src/paper.rs | matches |
| 1549 | compact, level I: 52 key generations, 100 signatures, 200 verifications of five signatures in each mode | 52, 100, 200, 5 | dim4-I-keygen.csv 52 rows, dim4-I-sign.csv 100, dim4-I-verify-*.csv "5 signatures x 40 runs" | matches |
| 1550 | levels III and V: 20 key generations, 52 signatures, mean over six signatures of the median of 15 runs | 20, 52, 6, 15 | dim4-{III,V}-keygen.csv 20 rows, -sign 52, -verify-* "6 signatures x 15 runs"; paper_report.py takes the mean of the six medians (reproduced by recompute.py) | matches |
| 1553 | the-sqisign at tag nist-v3, Release, broadwell | nist-v3 | machine.txt lines 62-66: "nist-v3", CMAKE_BUILD_TYPE=Release, SQISIGN_BUILD_TYPE=broadwell | matches |
| 1556 | iterations 300 at level I, 100 at III and V | 300, 100 | session.log: benchmark_p324_3 --iterations=300, _p500_27 and _p664_17 --iterations=100, in each pass | matches |
| 1556 | three passes; C reference keygen and sign at III and V differed by more than 1.5 % | 3, 1.5 % | numbers.json: passes 3; pass 1 against pass 2: c-III-keygen 4.7 %, c-III-sign 7.7 %, c-V-keygen 2.7 %, c-V-sign 2.9 %; every other row under 1.5 % | matches |
| 1556 | no Rust row moved by more than 1.2 % (was 1.1) | 1.2 % | largest spread of a Rust row over the three passes: dim4-III-sign 911.94 to 922.25, 1.13 %; next dim2-III-keygen 1.10 % | stale (fixed) |
| 1558 | timings taken with the bit-packed layout | packed | raw headers: "(84 + 139 bytes, packed layout)"; the harness calls hd_verify_bytes_with on the packed bytes | matches |
| 1561 | pre-chain checks cost 1.2, 3.2 and 5.5 megacycles | 1.2, 3.2, 5.5 | profile columns checks_hints_pairing + checks_goodness, medians, median of three passes: 0.286 + 0.952 = 1.24; 0.778 + 2.399 = 3.18; 1.745 + 3.767 = 5.51 (pass*/raw/dim4-*-profile-four.csv) | derived (recomputed) |
| 1562 | under 2 % of a four-point verification and under 3 % of a two-point one (added) | 2 %, 3 % | 1.24/76.54 = 1.6 %, 3.18/194.63 = 1.6 %, 5.51/482.24 = 1.1 %; 1.24/48.05 = 2.6 %, 3.18/120.06 = 2.6 %, 5.51/291.58 = 1.9 % | derived (recomputed) |

## Section 8.2, Results: sizes

| line | text | number | source | status |
|---|---|---|---|---|
| 1567 | 223 bytes against 283; a saving of 30.5 % on the signature, of which bit packing accounts for 1.5 points | 223, 283, 30.5, 1.5 | 84 + 139; 83 + 200; (200 - 139)/200 = 30.5 %; 3/200 = 1.5 | derived (recomputed) |
| 1570 | 176-byte compressed signature | 176 | raw header of dim2-I-verify-compressed.csv: "(176 bytes)" | matches |
| 1584 | Table 2, standard: 83 + 200 = 283, 129 + 306 = 435, 169 + 406 = 575 | 9 numbers | raw headers of dim2-*-verify.csv; sums rechecked | derived (recomputed) |
| 1585 | Table 2, compressed: 83 + 176 = 259, 129 + 269 = 398, 169 + 353 = 522 | 9 numbers | raw headers of dim2-*-verify-compressed.csv (176, 269, 353); sums rechecked | derived (recomputed) |
| 1586 | Table 2, compact: 84 + 139 = 223, 130 + 212 = 342, 170 + 278 = 448 | 9 numbers | raw headers of dim4-*-verify-four.csv; sums rechecked | derived (recomputed) |

## Section 8.2, Results: time

| line | text | number | source | status |
|---|---|---|---|---|
| 1594 | level I: verifier 6.8 times dimension 2; signing 3.6 times | 6.8, 3.6 | 76.540/11.185 = 6.84; 360.544/100.023 = 3.60 | derived (recomputed) |
| 1596 | levels III and V: verifier 6.5 and 6.0 times | 6.5, 6.0 | 194.632/29.798 = 6.53; 482.237/80.017 = 6.03 | derived (recomputed) |
| 1597 | two-point verifier 4.0 and 3.6 times | 4.0, 3.6 | 120.059/29.798 = 4.03; 291.583/80.017 = 3.64 | derived (recomputed) |
| 1598 | signing 2.4 and 2.2 times | 2.4, 2.2 | 913.887/375.368 = 2.43; 1714.477/769.668 = 2.23 | derived (recomputed) |
| 1600 | middle half of the signatures: 398 to 1471 at level III, 869 to 3392 at level V (was 3391) | 398, 1471, 869, 3392 | quartiles of the 52 signatures of each pass (exclusive method), median of the three passes: III 398.3 and 1470.9; V 868.7 and 3391.95. The 3391 was pass 1 alone | stale (fixed) |
| 1603 | key generation 35.3 against 32.9 megacycles | 35.3, 32.9 | dim4-I-keygen 35.269; dim2-I-keygen 32.938 | matches |
| 1604 | Rust dimension 2 over C at level I: 1.02 for verification | 1.02 | 11.185/10.982 = 1.019 | derived (recomputed) |
| 1605 | 1.06 for key generation (was 1.05), 1.09 for signing | 1.06, 1.09 | 32.938/31.185 = 1.056; 100.023/91.960 = 1.088. The one-decimal cells 32.9/31.2 give 1.054, which is where 1.05 came from; the unrounded medians are used for every ratio | stale (fixed) |
| 1606 | compressed recovery adds 19 % to a verification | 19 % | 13.266/11.185 = 1.186 | derived (recomputed) |
| 1623 | Table 3, C reference I: 31.2, 92.0, 11.0 | 31.2, 92.0, 11.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 31.185, 91.960, 10.982 | matches |
| 1624 | Table 3, C reference III: 106.2, 323.1, 29.0 | 106.2, 323.1, 29.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 106.248, 323.050, 29.048 | matches |
| 1625 | Table 3, C reference V: 216.9, 679.6, 75.3 | 216.9, 679.6, 75.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 216.893, 679.626, 75.253 | matches |
| 1626 | Table 3, Rust dimension 2, I: 32.9, 100.0, 11.2 | 32.9, 100.0, 11.2 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 32.938, 100.023, 11.185 | matches |
| 1627 | Table 3, Rust dimension 2, III: 109.0, 375.4, 29.8 | 109.0, 375.4, 29.8 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 108.955, 375.368, 29.798 | matches |
| 1628 | Table 3, Rust dimension 2, V: 248.6, 769.7, 80.0 | 248.6, 769.7, 80.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 248.618, 769.668, 80.017 | matches |
| 1629 | Table 3, compressed verify, I: 13.3 | 13.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 13.266 | matches |
| 1630 | Table 3, compact I: 35.3, 360.5, 76.5; x dim. 2: 6.8 | 35.3, 360.5, 76.5; x dim. 2: 6.8 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 35.269, 360.544, 76.540; 6.84 | derived (recomputed) |
| 1631 | Table 3, two points I: 48.0; 4.3 | 48.0; 4.3 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 48.048; 4.30 | derived (recomputed) |
| 1632 | Table 3, compact III: 126.3, 913.9, 194.6; 6.5 | 126.3, 913.9, 194.6; 6.5 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 126.298, 913.887, 194.632; 6.53 | derived (recomputed) |
| 1633 | Table 3, two points III: 120.1; 4.0 | 120.1; 4.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 120.059; 4.03 | derived (recomputed) |
| 1634 | Table 3, compact V: 249.8, 1714.5, 482.2; 6.0 | 249.8, 1714.5, 482.2; 6.0 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 249.781, 1714.477, 482.237; 6.03 | derived (recomputed) |
| 1635 | Table 3, two points V: 291.6; 3.6 | 291.6; 3.6 | medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes: 291.583; 3.64 | derived (recomputed) |

## Figure 2, the trade-off

| line | text | number | source | status |
|---|---|---|---|---|
| 1645 | axes: x 120 to 240, y 0 to 90; ticks every 20 and every 30 | 120, 240, 0, 90 | layout; every point lies inside (largest y is 76.5) | matches |
| 1653 | compact: (139, 76.5) | 139, 76.5 | Table 2 and Table 3: dim4-I-verify-four 76.540 | matches |
| 1654 | compressed (176, 13.3); standard (200, 11.2) | 176, 13.3, 200, 11.2 | dim2-I-verify-compressed 13.266; dim2-I-verify 11.185 | matches |
| 1655 | two points: (139, 48.0) | 139, 48.0 | dim4-I-verify-two 48.048 | matches |
| 1656 | label positions (145,78), (145,49), (168,17), (203,15) | 8 numbers | layout only; each label sits next to its point | matches |

## Table 4, wall-clock milliseconds at level I

| line | text | number | source | status |
|---|---|---|---|---|
| 1684 | dimension 2 key generation 11.8 | 11.8 | ns column of the same raw files, median of three passes: 11.764 | matches |
| 1685 | dimension 2 signing 35.7 | 35.7 | ns column of the same raw files, median of three passes: 35.724 | matches |
| 1686 | dimension 2 verification 4.0 | 4.0 | ns column of the same raw files, median of three passes: 3.995 | matches |
| 1687 | compressed verification 4.7 | 4.7 | ns column of the same raw files, median of three passes: 4.738 | matches |
| 1689 | compact key generation 12.6 | 12.6 | ns column of the same raw files, median of three passes: 12.597 | matches |
| 1690 | compact signing 128.8 | 128.8 | ns column of the same raw files, median of three passes: 128.768 | matches |
| 1691 | compact verification 27.3 | 27.3 | ns column of the same raw files, median of three passes: 27.337 | matches |
| 1692 | two points 17.2 | 17.2 | ns column of the same raw files, median of three passes: 17.161 | matches |

## Section 8.2, the good-degree loop and Table 5

| line | text | number | source | status |
|---|---|---|---|---|
| 1700 | one candidate in two reaches a Miller-Rabin test | 1/2 | extra column of dim4-I-sign.csv: 54.87 tests for 108.94 samples per signature = 0.504 | derived (recomputed) |
| 1702 | median 95 samples, mean 109, range 1 to 446 | 95, 109, 1, 446 | extra column of dim4-I-sign.csv (the same 100 signatures in every pass): mean 108.94, median 95, min 1, max 446 | matches |
| 1716 | signatures measured: 100 | 100 | dim4-I-sign.csv, 100 rows | matches |
| 1717 | samples mean 109; expected 121 | 109, 121 | 108.94; 174 ln 2 = 120.6 | matches |
| 1718 | random parity: 241 | 241 | 2 * 174 ln 2 = 241.2 | derived (recomputed) |
| 1719 | median (min, max): 95 (1, 446) | 95, 1, 446 | as above | matches |
| 1720 | primality tests mean 55; expected 60 | 55, 60 | 54.87; 60.3 | matches |
| 1721 | signing 361 (46, 1529) megacycles (max was 1546) | 361, 46, 1529 | median of the three passes of each statistic: median 360.54, min 45.51, max 1528.9. The 1546 was the maximum of pass 1 alone (passes: 1546.3, 1528.9, 1524.1) | stale (fixed) |
| 1722 | dimension-2 signing 100 | 100 | dim2-I-sign 100.023 | matches |

## Section 8.2, heap

| line | text | number | source | status |
|---|---|---|---|---|
| 1728 | 1655 KiB (1.62 MiB) with four points, 1631 KiB with two | 1655, 1.62, 1631 | session/heap.csv, mean of five signatures: 1 694 442 bytes = 1654.7 KiB = 1.616 MiB; 1 669 866 bytes = 1630.7 KiB | derived (recomputed) |
| 1729 | 3443 and 6160 KiB at levels III and V; dimension 2 allocates nothing | 3443, 6160, 0 | heap.csv: 3 525 398 bytes = 3442.8 KiB; 6 307 749 bytes = 6159.9 KiB; dim2-*-verify 0 bytes, 0 allocations | derived (recomputed) |

## Section 8.3, Scope

| line | text | number | source | status |
|---|---|---|---|---|
| 1744 | 24 signatures at level I and 5 each at levels III and V | 24, 5, 5 | tests at head: "level I, 24 signatures: 24 accepted, 336 tampered variants rejected, ... 4022 plain steps"; "level III: 5 vectors"; "level V: 5 vectors" (tests-head-8747fa9.log) | matches |

## Section 9, Conclusion

| line | text | number | source | status |
|---|---|---|---|---|
| 1774 | 139 bytes against 200; 223 bytes | 139, 200, 223 | as in the abstract | matches |
| 1779 | two points: 4.3, 4.0 and 3.6 times dimension 2 | 4.3, 4.0, 3.6 | 4.30, 4.03, 3.64 (as in Table 3) | derived (recomputed) |
| 1780 | four points: 6.8, 6.5 and 6.0 times | 6.8, 6.5, 6.0 | 6.84, 6.53, 6.03 | derived (recomputed) |
| 1782 | signing 3.6, 2.4 and 2.2 times | 3.6, 2.4, 2.2 | 3.60, 2.43, 2.23 | derived (recomputed) |
| 1784 | 0.63, 0.62 and 0.60 | 0.63, 0.62, 0.60 | as in the introduction | derived (recomputed) |
| 1786 | a verification of a few tens of milliseconds | tens of ms | true at level I (17.2 and 27.3 ms, Table 4). This session at level III: 42.9 and 69.5 ms; at level V: 104.1 and 172.2 ms (two and four points) | needs an author (true at level I only) |
| 1794 | probability about 2^-44 at level I | 2^-44 | as in Section 7.1 | needs an author (rests on a figure not re-measured) |
| 1801 | 61 bytes saved at level I | 61 | 200 - 139 | derived (recomputed) |
| 1803 | MIKE's chain costs about 0.07 megacycles per step at level I | 0.07 | 16.0/228 = 0.070, both cited (BDD+26, Table 1) | cited (not ours, left alone) |
| 1804 | against about 0.28 for ours (48.0 megacycles over 2 * 87 steps, with two points) | 0.28, 48.0, 87 | 48.048/174 = 0.276 on this session's machine. Set against a figure from another machine. The original read 0.16 and 28.5; an earlier pass had already put this session's figures here; this audit left the sentence as it found it and marked it | needs an author (cross-machine comparison) |

## AI disclosure

| line | text | number | source | status |
|---|---|---|---|---|
| 1820 | the disclosure names Claude Opus 5.5 | 5.5 | an authors' statement about how the work was done; nothing in the session or the code sources it | cannot source |

## Appendix A, proofs

| line | text | number | source | status |
|---|---|---|---|---|
| 1883 | f = 324, 500, 664 against e + 2 = 176, 266, 348 | 324, 500, 664, 176, 266, 348 | the exponents of 2 in p + 1; 174 + 2, 264 + 2, 346 + 2 | derived (recomputed) |
| 2297 | pi has trace 4; fixes 0, 5, 10, 15 | 4, 0, 5, 10, 15 | code: PI_FIXED = [0, 5, 10, 15] (two_point.rs:58) | matches |
| 2299 | the remaining numbers of Appendix A (orders 2^(4e), 16 fixed points, eigenvalues, denominators 4, sixteen terms) | several | mathematics of the proofs; not measurements and not checked by the audit beyond the two rows above | matches |

## Appendix B, Conformance of the implementation

| line | text | number | source | status |
|---|---|---|---|---|
| 2364 | m exceeds 40 when a2 = 0; every m > 40 rejected; 1 <= m <= 40 | 40 | code: MAX_GLUING_M = 40 (crates/verify/src/hd/self_contained.rs:67, test at :414); see the code check below | matches |
| 2369 | 40 Miller-Rabin rounds on the first 40 primes | 40, 40 | code: crates/verify/src/hd/goodness.rs; see the code check below | matches |
| 2441 | two-point check at the gluing: four evaluations through the gluing and six comparisons | 4, 6 | code: four evaluations as stated (gluing_chain.rs:496-504); six comparisons is the null-point test alone (null_fixed_by_pi, two_point.rs:214-233, a seventh when the four fixed coordinates all vanish); the pairing of the images adds two projective equalities of 16 cross-products each (gluing_images_paired, two_point.rs:238-240) | needs an author (count not exact) |
| 2453 | the same codomain at all 844 plain steps of five level-I vectors; 70 tampered variants | 844, 5, 70 | test at head: "level I: 5 vectors accepted in the three modes with no fall-back; 70 tampered variants rejected with the same reason in the three modes; ... 844 plain steps with the same codomain in both two-point modes" | matches |
| 2458 | earlier version: 32 honest and 72 forged signatures, 17 790 steps | 32, 72, 17790 | no log, test or record of this run in the branch or in ~/work/results; the sentence says it was an earlier version of the crate | cannot source |
| 2460 | levels III and V: 1 274 and 1 706 plain steps of five vectors each; 75 and 70 tampered variants | 1274, 1706, 5, 75, 70 | tests at head: "level III: 5 vectors ... 75 tampered variants ... 1274 plain steps"; "level V: 5 vectors ... 70 tampered variants ... 1706 plain steps" | matches |

## Appendix C, Beyond two-point evaluation

| line | text | number | source | status |
|---|---|---|---|---|
| 2481 | MIKE: the null point vanishes at the six coordinates that [i] negates; ten coordinates kept | 6, 10 | cited (BDD+26, Definition G.1, Proposition G.3, Remark E.2) | cited (not ours, left alone) |
| 2485 | doubling about 1.6 times faster on the quotient | 1.6 | cited (BDD+26, Remark E.2) | cited (not ours, left alone) |
| 2491 | up to about 1 - 1/1.6 = 38 % | 38 % | 1 - 1/1.6 = 37.5 % | derived (recomputed) |

178 rows: 71 derived, 70 matches, 13 cited, 12 needs an author, 9 stale, 3 cannot source.

## Statements checked against the code at 8747fa9

Fixed in the tex, because the code differs:

- Line 1192 to 1196. The orbits are reached along a tree (a path through nine orbits and one leaf), at 25 multiplications per vector, and both the dual null point and its inverse are computed that way. The text said a fixed path and about 26 multiplications.
- Line 1195 to 1196. The Hadamard transform of a vector fixed by pi is 32 additions and 4 doublings. The text said 34 additions.
- Line 2448. The equations are solved along a tree, not a single path.
- Line 1558 to 1562. The timings include the pre-chain checks; their cost and share are given.

Sentences of the form "our crate does X", read against the code. Each item of Remark (a) to (k), the verifier of Section 4.2, the wire format, the signer's loop and the tests paragraph was read; file and line for each are in audit/code-check.md. Exact as written: (a), (b), (c), (d), (e) except the order noted below, (f), (i), (k), the wire format (field order, widths, least significant bit first, padding and range rejects, canonical hints), the 4096-sample restart, the m bound, the 40 rounds on the first 40 primes. (g) was checked for the labelling of index bits only; its match to the numbered algorithms of the cited paper was not checked.

No longer exact, or never exact. None of these was edited; they are for the authors.

1. Lines 434 to 438 and Algorithm 1 step 4: "a basis of E_chl[2^f]". The code builds no such basis. It takes the hinted basis of the public-key curve, pushes multiples of it through the challenge isogeny to get P and Q of order 2^r, and computes e(P, Q) from a pairing on the public-key curve (crates/verify/src/hd/challenge.rs:163-222, self_contained.rs:365-373).
2. Lines 438 to 440, Algorithm 1 steps 4 and 5, and item (e) "Before that": the reject of a and b both even comes before the commitment basis and the check on e(R, S), not after both pairings (response.rs:158-160, 182-189).
3. Line 445: "the two middle codomains have the same theta null point". The code compares the first half-chain's last null point with the Hadamard transform of the second's (chain.rs:58-67).
4. Lines 455 to 459: the list of what a signature fails on leaves out a and b both even, a failed Cornacchia step, m above 40, and the failures of the canonical-basis and curve constructions (hd_verify.rs:40-77).
5. Item (h): "computed as in (a) and (b)". The dimension-2 change of coordinates scans i_0 in the order 0, 2, 1, 3, not 0, 1, 2, 3 (dim2.rs:196-215).
6. Item (j): "wherever these algorithms would divide by a vanishing value, our crate rejects". The code check found inversions with no zero test in the dimension-2 chain (dim2.rs:137, 544-547, 695), in the special image of the gluing (gluing.rs:128-139, 160-168), in the image check's use of the last weights (self_contained.rs:316-318) and at weil4(..).inv() (self_contained.rs:453). I confirmed the lines exist as described; I did not establish whether an earlier test makes a zero impossible there.
7. Lines 1198 to 1199: "a vanishing coordinate sends the verifier back to the full computation". True for the codomain (two_point.rs:193-199, 378-381). A zero met by the doubling constants on the orbits makes the chain reject instead (structure.rs:123-129, two_point.rs:347-355). Nothing divides by zero in either case.
8. Line 2441 to 2442: "six comparisons". Six is the null-point test alone; the pairing of the images adds two projective equalities (two_point.rs:214-240).
9. Lines 1339 to 1341: the signer keeps a candidate only if, besides the stated test, Cornacchia succeeds and m is at most 40 (crates/sqisign-rs/src/compact.rs:418-424). Item (c) of the appendix says so; Section 7.1 does not.
10. Lines 1420 to 1422: "so that signing always ends with a valid signature". The signer gives up with an error after 64 commitments (MAX_COMMITMENTS, compact.rs:128, 290-303).
11. Lines 1503 to 1504: "each signature has a single encoding". The same section says, and the tests at head show, that flipping the top bit of a or b gave a second valid signature in 2 of 12 trials.

## Cross-machine comparisons

One sentence group sets a megacycle figure of ours against a figure measured on another machine: lines 1803 to 1808, MIKE's 0.07 megacycles per step against "about 0.28 for ours". It is marked in the tex with a TODO(author) comment at line 1802 and its wording is left as the audit found it. An earlier pass had already replaced the original's 0.16 and 28.5 by this session's 0.28 and 48.0; the audit did not undo that, because 28.5 appears nowhere else in the dated copy any more.

Every other ratio in the paper is between rows of this session on this machine, the C reference included (line 529, lines 1604 to 1606). Line 319 (228 steps in 16.0 megacycles) and Appendix C quote MIKE's figures on their own and compare them with nothing of ours.

The two figures marked "not re-measured" have a TODO(author) comment each: the Gram-matrix paragraph (comment at line 1379) and the hint pairs (comment at line 1500).

## The M4 column of TABLES.md

Kept. The figures were produced by this branch's harness on Apple hardware:

- Machine: Apple M4 Pro, 12 cores, 24 GiB, macOS 27.2.1, on mains power (MACHINE.md of the archive).
- Revision: 68a0017a5133e011bb1aeb25c01a1f72ea4c013a, "git status: 0 changed files", rustc 1.99.0, release build with the bench crate's profile.
- Harness: scripts/bench-mac.sh of this branch; field kernels "portable (aarch64)"; timer "the monotonic wall clock"; no megacycles; the C reference not run.
- Who ran it: the author, on their own Mac, 2026-11-14 23:12Z. The archive 20261114T231221Z-m4.tar.gz (sha256 add2441c09fbc2282095a46712b592caf5ed9dcfad4dcbf72ea1c18138838849) was copied to the workstation by the author and unpacked under ~/work/results/20261114T225808Z-SECENG-2636-mac/m4/.
- Rechecked against the archive's numbers.json: every cell of the column is the first pass's value, as that file's rule says. One sentence of TABLES.md was wrong and is corrected: the two passes are not within 3 percent on every row; dimension-2 key generation at level I differs by 4.2 percent, every other row by at most 2.5.

None of the M4 figures is in the tex, so nothing was removed from it.

TABLES.md was also brought in line with the tex on two audit corrections (1.06 and the maximum 1529).

## Sanity

- Every old value replaced, by the earlier passes and by this audit (86 patterns, listed in audit/old-values.txt), was searched in the dated tex outside comment lines. Eight raw hits, all of them inside other numbers or unrelated: "3.3" in 13.3, "6.9" in 216.9, "3.1" in 323.1, "9.8" in 29.8 and 249.8, "107" in 1107, "8.4" in the citation "Section 8.4.2", "four times" in "four times a kernel lift" (line 1284). No old value remains.
- latexmk -pdf: one error, the missing db.png; no undefined reference or citation; 33 pages (audit/latexmk-audit.log, audit/build/).

## Fixed

1. Abstract, both copies: "four times" the dimension-2 scheme is now 4.3 times.
2. Section 5: tree and 25 multiplications in place of a path and about 26; 32 additions and 4 doublings in place of 34 additions.
3. Section 5: doubling constants needed after about 43 % of the steps, not 42.
4. Setup: no Rust row moved by more than 1.2 %, not 1.1.
5. Setup: the pre-chain checks are under 2 % of a four-point verification and under 3 % of a two-point one.
6. Results: 1.06 times the reference for key generation, not 1.05.
7. Results: the upper quartile of signing at level V is 3392, not 3391.
8. Table 5: the maximum signing time is 1529 megacycles, not 1546.
9. Appendix B: "along a tree", not "along a single path".
10. TABLES.md: the M4 passes' agreement, the 1.06 and the 1529.

Everything else the earlier passes changed was recomputed from the raw files and stands: Table 3 (every cell and ratio), Table 4, Table 5 otherwise, the heap figures, the figure's coordinates, the Setup paragraph, the ratios in the introduction, Sections 4, 5, 8 and the conclusion, and the appendix's test counts.

## Needs an author

1. The MIKE comparison in the conclusion (lines 1803 to 1808): megacycles from two machines.
2. The Gram-matrix figures (lines 1382 to 1397): 4 180 commitments, 2^0.39, 2^3.66, 2^-8.6 X, the tail 1.6 y^-4 and its range, 1.9 to 2.8. Not re-measured; no instrumentation at head.
3. What rests on that 1.6: about 2^-44, 2^-46, 2^-47 (line 1414 to 1415), "once in 2^44" (line 124), 2^-44 (line 1795), "one signature in 2^40" (line 1422). The arithmetic from 1.6 is right.
4. "Over two hundred hint pairs per signature" (line 1502). Not re-measured; no test at head.
5. "A verification of a few tens of milliseconds" (line 1786): true at level I (17.2 and 27.3 ms); at level V this session measured 104 and 172 ms.
6. "Four evaluations through the gluing and six comparisons" (line 2441): the count of comparisons.
7. The eleven statements listed above under "No longer exact".
8. "About 2^9.7 of them are good" (line 1374): right if v and -v are counted once, 2^10.7 otherwise; the text does not say.

## Cannot source

1. "At least about 0.07 X^2/p vectors" (line 1370): the count of translates is not written out and was not re-derived. The 2^18.6, 2^19.4, 2^20.1 that follow are right given the 0.07.
2. "32 honest and 72 forged signatures (17 790 steps)" (lines 2458 to 2459): an earlier version of the crate; no log, test or record of it in the branch or in the results folders.
3. The AI disclosure's model name (line 1820): an authors' statement.
