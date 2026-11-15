# Code check for the audit of 2026-11-15: "our crate does X" against the code at 8747fa9

Read-only reading of crates/verify/src/hd/ (H/) and crates/sqisign-rs/ (S/). The items marked "not exact" were looked at a second time at the lines given; AUDIT.md lists them.

## Order of the verifier's checks

parse signature and padding (H/wire.rs:108-147), parse key (:160-169), challenge (:235); then build_setup_prof in H/self_contained.rs: hints of both curves (:348-355), challenge isogeny (:357), order of e(P, Q) (:365-373), recover_response (:376: both-even reject at H/response.rs:158-160, commitment basis, e(R, S), k exists and is odd :182-189, completion :193), goodness (:400), Cornacchia (:405), m bound (:414); chains; middle match (:538); image check (:546).

## Remark (a) to (k)

| item | claim | verdict | evidence |
|---|---|---|---|
| (a) | index i0 + 2 i1 + 4 i2 + 8 i3 | exact | H/dim4.rs:18-20, H/point.rs:10 |
| (a) | component l times column l; blocks (U V; W U') | exact | H/dim4.rs:24-26, 47-53, 69-103 |
| (a) | generators from columns 4 to 7; copies 0, 1 on E_com, 2, 3 on E_chl | exact | H/gluing_chain.rs:169-188 |
| (a) | theta_pi(a) in position a | exact | H/two_point.rs:44, 52-55 |
| (b) | first i_0 in 0..15 with nonzero T, tested on sections; i_0 = 0 when V = 0 | exact | H/dim4.rs:69-85 |
| (c) | odd value as a1; rejects even or not 1 mod 4 | exact | H/kani.rs:55, 86-90 |
| (c) | m = trailing zeros of a2, above 40 at a2 = 0; m > 40 rejected | exact | H/self_contained.rs:67, 409, 414-416; H/wide.rs:108-115 (192 for zero) |
| (c) | goodness before Cornacchia; 40 rounds; first 40 primes | exact | H/goodness.rs:11-14, 38-53, 59-68 |
| (c) | signer: same test, same bound | exact | S/src/compact.rs:418-424 |
| (d) | kernel point only through ratios of its own coordinates | exact | H/isogeny.rs:238-251, 325-329; H/gluing.rs:64-71; H/two_point.rs:185-201 |
| (e) | tau, base points, lambda = 1/q; dropped entry by the rule; hints first | exact | H/response.rs:196-197; H/hd_verify.rs:84-98; H/self_contained.rs:348-355, 433-441 |
| (e) | "Before that": both-even reject after both pairings | not exact | H/response.rs:158-160 comes before the commitment basis and e(R, S) |
| (f) | u1..u4 and the assembly | exact | H/kani.rs:492-586 |
| (g) | kernel point l labels index bit l | exact | H/isogeny.rs:168-176, 200-229 |
| (g) | match to the numbered algorithms of Dar24 | not checked | the cited paper was not consulted |
| (h) | one chain for both pairs; kernel points; matrices modulo 4; one translate | exact | H/gluing_chain.rs:369-396, 521-526; H/kani.rs:199-229; H/dim2.rs:593-598 |
| (h) | "computed as in (a) and (b)" | not exact | H/dim2.rs:196-215 scans i_0 as 0, 2, 1, 3 |
| (i) | four kernel points and the sum of the first two; two translates | exact | H/gluing_chain.rs:412-424, 496-504; H/gluing.rs:38 |
| (j) | a vanishing divisor rejects | not exact | no zero test at H/dim2.rs:137, 544-547, 695; H/gluing.rs:128-139, 160-168; H/self_contained.rs:316-318, 453 |
| (k) | (R0, 0, 0, 0); x-coordinates of the first two; fourth at infinity | exact | H/self_contained.rs:264-267, 309-327; H/response.rs:164, 169 |

## Elsewhere

| tex | claim | verdict | evidence |
|---|---|---|---|
| 434-438 | a basis of E_chl[2^f] from the hints | not exact | H/challenge.rs:163-222: basis of E_pk[2^f], pushed through the challenge isogeny |
| 445 | same theta null point | not exact (literal) | H/chain.rs:58-67: against the Hadamard transform of the other |
| 455-459 | what a signature fails on | incomplete | H/hd_verify.rs:40-77 |
| 1198-1199 | vanishing coordinate, back to the full computation | codomain exact; constants reject | H/two_point.rs:193-199, 347-355, 378-381; H/structure.rs:123-129 |
| 1200-1201 | doubling constants only when a doubling follows | exact (two-point loop) | H/two_point.rs:347-355 |
| 2440-2442 | images paired, null fixed, four evaluations, fall-back | exact | H/two_point.rs:238-271 |
| 2441 | six comparisons | not exact | six for the null point (H/two_point.rs:214-233); two projective equalities more for the images |
| 2450-2460 | codomain at every plain step in both two-point modes; verdicts in three modes; five vectors | exact | S/tests/two_point.rs:35, 172-197, 229-290 |
| 1482-1506 | field order, widths, (q - 3)/4 in e - 2 bits, hints 8 bits, LSB first, padding and range rejects | exact | H/wire.rs:56-81, 108-147; H/params.rs:42-44, 119-124 |
| 1493-1497 | canonical hints: first table entries on the curve | exact | H/basis.rs:165-178 |
| 1339-1341 | signer's test | incomplete | S/src/compact.rs:394-427 also needs Cornacchia and m <= 40 |
| 1420-1422 | 4096 samples, fresh commitment | exact | S/src/compact.rs:124, 292-302 |
| 1421 | signing always ends with a valid signature | not exact | MAX_COMMITMENTS = 64 (S/src/compact.rs:128, 290-303) |

Tampered variants (S/tests/two_point.rs:119-166): per vector, a wrong message, bit 0 of A_com, the low bit of q, the low bit of a, a signature one byte short, a wrong key, and the low bit of nine bytes among the last 60: 15 per vector. A flip that lands on a scalar's top bit is set aside and compared across modes only: 5 at levels I and V, none at III, which gives 70, 75, 70.
