# CHANGES: every difference between main_standalone.2026-11-04.orig.tex and main_standalone.tex

One entry per changed block of lines, in file order: the original's lines, the old text, the new text, and the source of the new text. Session files: paper-dim4/bench/sessions/2026-11-15-68a0017/. Instrumented measurements, the MIKE run and the pass after the divisor tests: ~/work/results/20261115T081239Z-SECENG-2636-close/.


## 1. original line 48 (dated copy 48)

old:

    SQIsignHD verifies a dimension-4 isogeny where SQIsign verifies a dimension-2 one, which removes the auxiliary curve from the signature but makes verification much costlier. We show that the verifier's chain of dimension-4 isogenies respects an automorphism of order 4 that acts on its theta coordinates as one fixed permutation, for every input, honest or forged. After the gluing, two of the four kernel points can therefore be recovered from the other two by permuting coordinates, and the constants of each step are determined by ten coordinates instead of sixteen, which cuts the verifier's cost by 37 to 40%. Along the way we derive the change of basis that the verifier applies at the gluing, which Dartois states with only a brief justification, prove that the verifier's image check remains sound with the signs that its implementation leaves free, and show that the verifier must validate the basis hints, which a forger controls. We implement the scheme at the round-3 SQIsign primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs four times that of the dimension-2 scheme in the same code.

new:

    SQIsignHD verifies a dimension-4 isogeny where SQIsign verifies a dimension-2 one, which removes the auxiliary curve from the signature but makes verification much costlier. We show that the verifier's chain of dimension-4 isogenies respects an automorphism of order 4 that acts on its theta coordinates as one fixed permutation, for every input, honest or forged. After the gluing, two of the four kernel points can therefore be recovered from the other two by permuting coordinates, and the constants of each step are determined by ten coordinates instead of sixteen, which cuts the verifier's cost by 37 to 40%. Along the way we derive the change of basis that the verifier applies at the gluing, which Dartois states with only a brief justification, prove that the verifier's image check remains sound with the signs that its implementation leaves free, and show that the verifier must validate the basis hints, which a forger controls. We implement the scheme at the round-3 SQIsign primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs 4.3 times that of the dimension-2 scheme in the same code.

source: Two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (session of 2026-11-15, medians of three passes).

## 2. original line 55 (dated copy 55)

old:

    We implement the scheme at the round-3 {\SQIsign} primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs four times that of the dimension-2 scheme in the same code.

new:

    We implement the scheme at the round-3 {\SQIsign} primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs 4.3 times that of the dimension-2 scheme in the same code.

source: The same sentence in the typeset abstract.

## 3. original line 109 (dated copy 109)

old:

    0.63, 0.61 and 0.60 times the cycles of the four-point one at levels I,

new:

    0.63, 0.62 and 0.60 times the cycles of the four-point one at levels I,

source: two/four at level III: 120.059 / 194.632 = 0.617.

## 4. original line 434 (dated copy 434)

old:

    and checks that the basis hints are the canonical ones of their curves. From the hints it builds a basis $(R_0, S_0)$ of $\Ecom[2^f]$,

new:

    and checks that the basis hints are the canonical ones of their curves. From the hints it builds a basis of the $2^f$-torsion of the public-key curve,

source: Code: the hinted basis is that of the public-key curve (crates/verify/src/hd/challenge.rs:163-222).

## 5. original lines 436 to 438 (dated copy 436 to 444)

old:

    specification, and a basis of $\Echl[2^f]$; their multiples by $2^{f-r}$
    are the bases of the $2^r$-torsion in which the action is written, and the
    verifier checks that both pair to roots of unity of exact order $2^r$. It

new:

    specification, and a basis $(R_0, S_0)$ of $\Ecom[2^f]$. No basis of
    $\Echl[2^f]$ is built: the challenge isogeny, of degree $2^\lambda$, maps
    the multiples by $2^{f - \lambda - r}$ of the public-key basis to a basis
    $(P, Q)$ of $\Echl[2^r]$, and $(R, S) = ([2^{f-r}]R_0, [2^{f-r}]S_0)$ is
    the basis of $\Ecom[2^r]$; the action is written in these two bases. The
    verifier checks, in this order, that $e(P, Q)$ has exact order $2^r$, that
    the transmitted $a$ and $b$ are not both even, and that $e(R, S)$ has
    exact order $2^r$, which it does by computing the integer $k$ with $e(R,
    S) = e(P, Q)^k$ and requiring that it exist and be odd. It

source: Code: no basis of E_chl[2^f]; (P, Q) is the image of the public-key basis under the challenge isogeny; the order of the checks is e(P, Q), a and b both even, e(R, S) (self_contained.rs build_setup_prof, response.rs recover_response).

## 6. original line 446 (dated copy 452 to 454)

old:

    point, projectively, and the chain maps $(R_0, 0, 0, 0)$ to $(\pm[a_1]R_0,

new:

    point, projectively (our code compares the null point of the first with
    the Hadamard transform of that of the second, since it evaluates the
    second half-chain through its dual), and the chain maps $(R_0, 0, 0, 0)$ to $(\pm[a_1]R_0,

source: Code: middle_codomain_matches compares with the Hadamard transform of the second null point (chain.rs:58-67).

## 7. original lines 456 to 459 (dated copy 464 to 470)

old:

    fails on any decoding error, on a hint that is not canonical or a basis
    that does not pair to a root of unity of exact order $2^r$, on a degree
    that is not good, on a chain that lands on a product of curves where it
    should not, or on either condition.

new:

    fails on any decoding error, on a curve coefficient that is not that of a
    curve, on a hint that is not canonical, on a pairing that does not have
    exact order $2^r$, on transmitted $a$ and $b$ that are both even, on a
    degree that is not good, on a failure of Cornacchia's algorithm or a
    value $m > 40$ (Appendix~\ref{app:remains}), on a step that would divide
    by a vanishing value, as when the chain lands on a product of curves
    where it should not, or on either condition.

source: Code: the reject reasons of HdReject (hd_verify.rs:40-77) and the divisor tests added at 6302b7a.

## 8. original lines 483 to 485 (dated copy 494)

old:

    $\Echl \gets$ the challenge curve, recomputed from the public key,
    $\Ecom$ and $\mu$
    \State reject unless the hints of the public-key curve and of $\Ecom$

new:

    reject unless the hints of the public-key curve and of $\Ecom$

source: Code: the order of the verifier's steps (hints, then the challenge isogeny and e(P, Q)).

## 9. original lines 487 to 489 (dated copy 496 to 503)

old:

    \State bases of $\Ecom[2^f]$ and $\Echl[2^f]$ from the hints, and of
    the $2^r$-torsion from them; reject unless both pair to roots of unity
    of exact order $2^r$

new:

    \State $\Echl \gets$ the challenge curve, recomputed from the public
    key, $\Ecom$ and $\mu$; $(P, Q) \gets$ the basis of $\Echl[2^r]$ that
    is the image, under the challenge isogeny, of the hinted basis of
    the public-key curve; reject unless $e(P, Q)$ has exact order $2^r$
    \State reject if $a$ and $b$ are both even
    \State $(R, S) \gets$ the basis of $\Ecom[2^r]$ from the hinted
    basis $(R_0, S_0)$ of $\Ecom[2^f]$; reject unless $e(R, S)$ has
    exact order $2^r$

source: Code: as entry 5, in the algorithm.

## 10. original lines 491 to 492 (dated copy 505)

old:

    \end{smallmatrix}\right)$ by the rule of Section~\ref{sec:wire};
    reject if $a$ and $b$ are both even

new:

    \end{smallmatrix}\right)$ by the rule of Section~\ref{sec:wire}

source: Code: the both-even reject is its own step, before e(R, S) (response.rs:158-160).

## 11. original lines 527 to 529 (dated copy 540 to 542)

old:

    coordinates instead of the four of a dimension-2 step, and costs 6.4 times
    the dimension-2 verifier of the same crate at the same prime (45.6 against
    7.1 megacycles), and 6.6 times the reference's assembly build

new:

    coordinates instead of the four of a dimension-2 step, and costs 6.8 times
    the dimension-2 verifier of the same crate at the same prime (76.5 against
    11.2 megacycles), and 7.0 times the reference's assembly build

source: dim4-I-verify-four 76.540, dim2-I-verify 11.185, c-I-verify 10.982 megacycles; 6.84 and 6.97.

## 12. original line 531 (dated copy 544)

old:

    after the gluing. A profile of the verifier at level I attributes 59\,\%

new:

    after the gluing. A profile of the verifier at level I attributes 60\,\%

source: dim4-I-profile-four.csv, slot4 over total: 59.7, 59.6, 59.8 %.

## 13. original line 533 (dated copy 546)

old:

    them, with the other saved points, through the chain, and 16\,\% to

new:

    them, with the other saved points, through the chain, and 15\,\% to

source: dim4-I-profile-four.csv, slot5 over total: 15.0, 15.1, 15.1 %.

## 14. original line 1180 (dated copy 1193)

old:

    the two-point verifier took 0.63, 0.61 and 0.60 times the cycles of the four-point one at levels I, III and V, on

new:

    the two-point verifier took 0.63, 0.62 and 0.60 times the cycles of the four-point one at levels I, III and V, on

source: as entry 3.

## 15. original lines 1192 to 1195 (dated copy 1205 to 1209)

old:

    of its techniques: the inverse of the dual null point is read off the two
    pushed kernel points as products along a fixed path through the ten orbits,
    in about 26 multiplications, and the Hadamard transform of a vector fixed by
    $\pi$ takes 34 additions instead of 64. Like MIKE's, this computation uses

new:

    of its techniques: the dual null point and its inverse are each read off the
    two pushed kernel points as products along a fixed tree on the ten orbits,
    a path through nine of them and one leaf, in 25 multiplications, and the
    Hadamard transform of a vector fixed by $\pi$ takes 32 additions and 4
    doublings instead of 64 additions. Like MIKE's, this computation uses

source: Code: along_tree, 25 multiplications per vector; hadamard_pi_assign, 32 additions and 4 doublings (two_point.rs).

## 16. original lines 1197 to 1198 (dated copy 1211 to 1214)

old:

    applies to our verifier by Appendix~\ref{app:remains}. It never divides by zero, since a vanishing
    coordinate sends the verifier back to the full computation. Finally, as in

new:

    applies to our verifier by Appendix~\ref{app:remains}. It never divides by zero: a vanishing
    coordinate among those the codomain formula reads sends the step back to
    the full computation, and a vanishing coordinate among the doubling
    constants makes the verifier reject, as it does with four points. Finally, as in

source: Code: a zero among the tree's coordinates falls back to from_kernel; a zero in precompute_pi rejects (two_point.rs, structure.rs).

## 17. original line 1200 (dated copy 1216)

old:

    when a doubling follows, which happens after about 42\,\% of the steps of an

new:

    when a doubling follows, which happens after about 43\,\% of the steps of an

source: The code's strategy replayed (audit strat.py): 36 of 83 to 85 codomains; 360 of the 844 plain steps of the five level-I vectors = 42.7 %.

## 18. original line 1340 (dated copy 1356 to 1357)

old:

    probable prime (40 Miller--Rabin rounds). The sampled norms are odd, so one

new:

    probable prime (40 Miller--Rabin rounds), and for which the verifier's
    Cornacchia step succeeds with $m \le 40$ (Appendix~\ref{app:remains}). The sampled norms are odd, so one

source: Code: the signer also requires Cornacchia to succeed with m <= 40 (crates/sqisign-rs/src/compact.rs).

## 19. original line 1373 (dated copy 1390)

old:

    about $2^{9.7}$ of them are good at level I, and the chance that none is is about $\exp(-2^{9.7})$. A degenerate class

new:

    about $2^{9.7}$ of them are good at level I, counting $v$ and $-v$, which have the same norm, once, and the chance that none is is about $\exp(-2^{9.7})$. A degenerate class

source: paper-dim4/bench/short_vector_count.py: 2^9.66 with v and -v counted once.

## 20. original line 1380 (dated copy 1397 to 1399)

old:

    and computed its successive minima for 4\,180 commitments at level I.

new:

    and computed its successive minima for 4\,180 commitments at level I,
    those of the 100 signatures of Table~\ref{tab:loop} and of 4\,080 further
    signatures of the same key.

source: instrument/gram.csv: 4 180 commitments, the 100 signatures of the loop table and 4 080 more (paper_instrument at bf2153b).

## 21. original lines 1382 to 1385 (dated copy 1401 to 1405)

old:

    minimum has median $2^{0.39}$ and maximum $2^{3.66}$, that is, $2^{-8.6}X$,
    and its tail follows $\Pr[N_3 > y\sqrt{p}/2] \approx 1.6\,y^{-4}$ closely
    over $1.4 \le y \le 5.7$, the range in which the sample has between 1\,600
    and 6 points: the product $y^4 \Pr$ stays between 1.5 and 1.8. This is

new:

    minimum has median $2^{0.39}$ and maximum $2^{3.65}$, that is, $2^{-8.6}X$,
    and its tail follows $\Pr[N_3 > y\sqrt{p}/2] \approx 1.5\,y^{-4}$
    over $1.4 \le y \le 4.8$, the range in which the sample has between 1\,588
    and 6 points: the product $y^4 \Pr$ has median 1.53 there and stays between
    1.3 and 1.6 while at least 50 points remain, falling to 0.8 at the last 6. This is

source: instrument/minima-all.txt (gram_minima.py): nu3 median 2^0.39, maximum 2^3.65; tail constant median 1.53 over 1.41 <= y <= 4.84 (1 588 to 6 points), 1.34 to 1.60 with at least 50 points, 0.79 at the last 6.

## 22. original lines 1391 to 1392 (dated copy 1411 to 1412)

old:

    equality for small $x$; our sample gives between $1.9\,x^2$ and
    $2.8\,x^2$ over $2^{-5} \le x \le 2^{-1}$, within sampling error. For the third

new:

    equality for small $x$; our sample gives between $1.8\,x^2$ and
    $3.0\,x^2$ over $2^{-5} \le x \le 2^{-1}$, within sampling error. For the third

source: minima-all.txt: Pr[nu1 <= x]/x^2 between 1.82 and 3.00 over the range.

## 23. original lines 1418 to 1420 (dated copy 1438 to 1446)

old:

    degree and starts over with a fresh one, so that signing always ends with
    a valid signature; by the estimate above, this happens in fewer than one
    signature in $2^{40}$. The restart changes the distribution of the

new:

    degree and starts over with a fresh one, and it stops with an error after
    64 commitments. A commitment is given up when its class is degenerate,
    with probability about $2^{-44}$ at level I by the estimate above, or when
    4096 samples of a class that is not degenerate contain no good degree,
    with probability $(1 - 1/109)^{4096} \approx 2^{-54}$ at the measured
    density of one good degree in 109 samples (Table~\ref{tab:loop}). A
    restart thus happens in fewer than one signature in $2^{40}$, and the
    error, which takes 64 such commitments in a row, with probability about
    $2^{-2800}$. The restart changes the distribution of the

source: Code: MAX_COMMITMENTS = 64; (1 - 1/108.94)^4096 = 2^-54.5 from the loop table's density; 64 * -44.2 = -2830.

## 24. original lines 1497 to 1499 (dated copy 1523 to 1529)

old:

    valid encoding of the same signature, which in our tests happened for
    over two hundred hint pairs per signature at level I. With it, each
    signature has a single encoding, and unpacking costs a few hundred

new:

    valid encoding of the same signature. With the check switched off, 97 of
    the 399 other pairs of table entries gave such an encoding per signature
    on average at level I, over the 100 signatures of
    Table~\ref{tab:loop} (median 95, between 39 and 224), and the verifier
    rejected every one of them with the check on. With it, the hints
    leave a signature a single encoding (the top bit of a scalar, discussed
    below, is a separate matter), and unpacking costs a few hundred

source: instrument/hints.csv: mean 96.7, median 95, min 39, max 224 other valid encodings among 399 pairs, over the 100 signatures; every one rejected by the verifier with the check on. The sentence on a single encoding now sets the top bit aside.

## 25. original lines 1536 to 1538 (dated copy 1566 to 1569)

old:

    with one vCPU of an Intel Xeon at 2.30~GHz and 3~GiB of memory, on
    2026-11-04. One binary measures the Rust code, compiled by \texttt{rustc}
    1.89 in release mode with fat link-time optimization and no target-CPU

new:

    with 32 vCPUs of an Intel Xeon at 2.80~GHz and 126~GiB of memory (Google
    Compute Engine, \texttt{n2-standard-32}), with the measuring process
    pinned to one core whose sibling was kept idle, on 2026-11-15. One binary measures the Rust code, compiled by \texttt{rustc}
    1.100 in release mode with fat link-time optimization and no target-CPU

source: session/machine.txt: 32 CPUs, Intel Xeon 2.80 GHz, 125.9 GiB, n2-standard-32, core 13 with sibling 29, rustc 1.100.0, passes on 2026-11-15.

## 26. original line 1550 (dated copy 1581)

old:

    (\texttt{-{}-iterations=300} at level I, 100 at III and V), medians.

new:

    (\texttt{-{}-iterations=300} at level I, 100 at III and V), medians. Every row was measured in three passes of the session: key generation and signing of the C reference at levels III and V, whose benchmark draws a fresh seed on every run, differed by more than 1.5\,\% between the first two passes, so a third was run, and the tables report the median of the three passes' values; no row of the Rust code moved by more than 1.2\,\%.

source: numbers.json: three passes; the four C reference rows above 1.5 %; largest Rust spread 1.13 %.

## 27. original lines 1552 to 1557 (dated copy 1583 to 1591)

old:

    verify. The timings were taken with the byte-aligned layout that preceded
    the bit-packed one; decoding the packed layout costs a few hundred cycles
    more, far below the precision of the tables. They also precede the verifier's checks of the hints and of the
    orders of the pairings (Section~\ref{sec:wire}), which in a separate
    measurement at level I cost about 0.4\,\% of a four-point verification and
    0.6\,\% of a two-point one.

new:

    verify. The timings were taken with the bit-packed layout, and include the
    verifier's checks of the hints and of the orders of the pairings
    (Section~\ref{sec:wire}) and its goodness test, which together cost
    1.2, 3.2 and 5.5 megacycles per verification at levels I, III and V, under
    2\,\% of a four-point verification and under 3\,\% of a two-point one.
    The verifier's tests of its divisors (Remark~\ref{rem:conformance}(j))
    were added after this session; one further pass on the same core with
    them moved no verification row by more than 0.6\,\%, and the tables are
    those of the session.

source: raw headers "packed layout"; pre-chain checks 1.24, 3.18, 5.51 megacycles and their shares; the pass of 2026-11-15T08:51Z with the divisor tests: largest move of a verification row 0.54 %.

## 28. original lines 1589 to 1593 (dated copy 1623 to 1627)

old:

    megacycles. At level I the dimension-4 verifier costs 6.4 times the
    dimension-2 verifier of the same crate, and dimension-4 signing costs 3.3
    times dimension-2 signing, because the good degree is found by rejection. At levels III and V the dimension-4 verifier costs 6.4 and 5.8
    times the dimension-2 one, the two-point verifier 3.9 and 3.5 times, and
    signing 2.8 and 2.0 times as much in median. Signing times vary widely at

new:

    megacycles. At level I the dimension-4 verifier costs 6.8 times the
    dimension-2 verifier of the same crate, and dimension-4 signing costs 3.6
    times dimension-2 signing, because the good degree is found by rejection. At levels III and V the dimension-4 verifier costs 6.5 and 6.0
    times the dimension-2 one, the two-point verifier 4.0 and 3.6 times, and
    signing 2.4 and 2.2 times as much in median. Signing times vary widely at

source: Table 3's ratios from the unrounded medians: 6.84, 3.60; 6.53, 6.03; 4.03, 3.64; 2.43, 2.23.

## 29. original line 1595 (dated copy 1629)

old:

    the signatures took 300 to 1289 megacycles at level III and 567 to 1596 at

new:

    the signatures took 398 to 1471 megacycles at level III and 869 to 3392 at

source: quartiles of the 52 signatures per pass, median of three passes: 398.3 and 1470.9; 868.7 and 3391.95.

## 30. original lines 1598 to 1601 (dated copy 1632 to 1635)

old:

    (22.6 against 23.7 megacycles). As a check on the implementation, the
    crate's dimension-2 code runs at 1.03 times the reference's assembly build
    for verification, 1.23 times for key generation and 1.20 times for
    signing, and the compressed format's recovery adds 18\,\% to a

new:

    (35.3 against 32.9 megacycles). As a check on the implementation, the
    crate's dimension-2 code runs at 1.02 times the reference's assembly build
    for verification, 1.06 times for key generation and 1.09 times for
    signing, and the compressed format's recovery adds 19\,\% to a

source: 35.269 against 32.938; Rust over C at level I 1.019, 1.056, 1.088; compressed over standard 1.186.

## 31. original lines 1618 to 1630 (dated copy 1652 to 1664)

old:

    C reference, dimension 2 & I   & 19.2 & 56.2 & 6.9 & \\
    C reference, dimension 2 & III & 62.1 & 181.6 & 18.0 & \\
    C reference, dimension 2 & V   & 133.0 & 384.2 & 44.4 & \\
    Rust, dimension 2 & I   & 23.7 & 67.2 & 7.1 & \\
    Rust, dimension 2 & III & 72.6 & 219.2 & 18.6 & \\
    Rust, dimension 2 & V   & 143.6 & 453.6 & 50.0 & \\
    Rust, dimension 2, compressed verify & I & & & 8.4 & \\
    Rust, dimension 4 (compact) & I & 22.6 & 219.5 & 45.6 & 6.4 \\
    \quad two points & I & & & 28.5 & 4.0 \\
    Rust, dimension 4 (compact) & III & 75.1 & 610.7 & 118.9 & 6.4 \\
    \quad two points & III & & & 72.8 & 3.9 \\
    Rust, dimension 4 (compact) & V & 163.3 & 921.7 & 288.1 & 5.8 \\
    \quad two points & V & & & 174.1 & 3.5 \\

new:

    C reference, dimension 2 & I   & 31.2 & 92.0 & 11.0 & \\
    C reference, dimension 2 & III & 106.2 & 323.1 & 29.0 & \\
    C reference, dimension 2 & V   & 216.9 & 679.6 & 75.3 & \\
    Rust, dimension 2 & I   & 32.9 & 100.0 & 11.2 & \\
    Rust, dimension 2 & III & 109.0 & 375.4 & 29.8 & \\
    Rust, dimension 2 & V   & 248.6 & 769.7 & 80.0 & \\
    Rust, dimension 2, compressed verify & I & & & 13.3 & \\
    Rust, dimension 4 (compact) & I & 35.3 & 360.5 & 76.5 & 6.8 \\
    \quad two points & I & & & 48.0 & 4.3 \\
    Rust, dimension 4 (compact) & III & 126.3 & 913.9 & 194.6 & 6.5 \\
    \quad two points & III & & & 120.1 & 4.0 \\
    Rust, dimension 4 (compact) & V & 249.8 & 1714.5 & 482.2 & 6.0 \\
    \quad two points & V & & & 291.6 & 3.6 \\

source: every cell: medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt (audit recompute.out).

## 32. original line 1640 (dated copy 1674)

old:

    xmin=120, xmax=240, ymin=0, ymax=60,

new:

    xmin=120, xmax=240, ymin=0, ymax=90,

source: figure layout: the y axis has to reach 76.5.

## 33. original line 1644 (dated copy 1678)

old:

    xtick={120,140,...,240}, ytick={0,20,40,60},

new:

    xtick={120,140,...,240}, ytick={0,30,60,90},

source: figure layout: ticks for the new axis.

## 34. original lines 1648 to 1654 (dated copy 1682 to 1688)

old:

    \addplot[only marks, mark=*, mark size=2.5pt, color=teal] coordinates {(139,45.6)};
    \addplot[only marks, mark=*, mark size=2.5pt, color=blue!45] coordinates {(176,8.4) (200,7.1)};
    \addplot[only marks, mark=*, mark size=2.5pt, color=red] coordinates {(139,28.5)};
    \node[teal, anchor=west, font=\small] at (axis cs:145,47) {compact (dimension 4)};
    \node[red, anchor=west, font=\small] at (axis cs:145,29) {two points};
    \node[blue!55, anchor=south, font=\small] at (axis cs:168,11) {compressed};
    \node[blue!55, anchor=south west, font=\small] at (axis cs:203,10) {standard};

new:

    \addplot[only marks, mark=*, mark size=2.5pt, color=teal] coordinates {(139,76.5)};
    \addplot[only marks, mark=*, mark size=2.5pt, color=blue!45] coordinates {(176,13.3) (200,11.2)};
    \addplot[only marks, mark=*, mark size=2.5pt, color=red] coordinates {(139,48.0)};
    \node[teal, anchor=west, font=\small] at (axis cs:145,78) {compact (dimension 4)};
    \node[red, anchor=west, font=\small] at (axis cs:145,49) {two points};
    \node[blue!55, anchor=south, font=\small] at (axis cs:168,17) {compressed};
    \node[blue!55, anchor=south west, font=\small] at (axis cs:203,15) {standard};

source: the four points are Table 2's bytes against Table 3's level-I cells; labels moved next to their points.

## 35. original lines 1679 to 1682 (dated copy 1713 to 1716)

old:

    dimension 2 (standard) & key generation & 10.3 \\
    & signing & 29.2 \\
    & verification & 3.1 \\
    & verification, compressed & 3.6 \\

new:

    dimension 2 (standard) & key generation & 11.8 \\
    & signing & 35.7 \\
    & verification & 4.0 \\
    & verification, compressed & 4.7 \\

source: ns column of dim2-I-*.csv: 11.764, 35.724, 3.995, 4.738 ms.

## 36. original lines 1684 to 1687 (dated copy 1718 to 1721)

old:

    dimension 4 (compact) & key generation & 9.8 \\
    & signing & 95.5 \\
    & verification & 19.8 \\
    & verification, two points & 12.4 \\

new:

    dimension 4 (compact) & key generation & 12.6 \\
    & signing & 128.8 \\
    & verification & 27.3 \\
    & verification, two points & 17.2 \\

source: ns column of dim4-I-*.csv: 12.597, 128.768, 27.337, 17.161 ms.

## 37. original lines 1697 to 1698 (dated copy 1731 to 1732)

old:

    the spread of signing times: the median signature needs 80 samples, the
    mean 107, and the range is 1 to 566.

new:

    the spread of signing times: the median signature needs 95 samples, the
    mean 109, and the range is 1 to 446.

source: extra column of dim4-I-sign.csv: median 95, mean 108.94, min 1, max 446 samples.

## 38. original line 1712 (dated copy 1746)

old:

    response samples per signature, mean & 107 & 121 \\

new:

    response samples per signature, mean & 109 & 121 \\

source: as entry 37.

## 39. original lines 1714 to 1717 (dated copy 1748 to 1751)

old:

    response samples per signature, median (min, max) & 80 (1, 566) & \\
    primality tests per signature, mean & 53 & 60 \\
    signing, Mcycles, median (min, max) & 220 (27, 1400) & \\
    signing, dimension 2, Mcycles, median & 67 & \\

new:

    response samples per signature, median (min, max) & 95 (1, 446) & \\
    primality tests per signature, mean & 55 & 60 \\
    signing, Mcycles, median (min, max) & 361 (46, 1529) & \\
    signing, dimension 2, Mcycles, median & 100 & \\

source: dim4-I-sign.csv: 95 (1, 446); 54.87 tests; signing 360.54 (45.51, 1528.9), medians of the three passes; dim2-I-sign 100.02.

## 40. original lines 1723 to 1724 (dated copy 1757 to 1758)

old:

    1673~KiB (1.63~MiB) of heap with four points and 1637~KiB with two,
    measured with a counting allocator around the call, and 3474 and 6271~KiB

new:

    1655~KiB (1.62~MiB) of heap with four points and 1631~KiB with two,
    measured with a counting allocator around the call, and 3443 and 6160~KiB

source: session/heap.csv, means of five signatures: 1654.7, 1630.7, 3442.8, 6159.9 KiB.

## 41. original lines 1774 to 1777 (dated copy 1808 to 1811)

old:

    verification costs 4.0, 3.9 and 3.5 times a dimension-2 one in
    the same code at levels I, III and V with two points, and 6.4, 6.4 and
    5.8 times with four, and
    signing costs 3.3, 2.8 and 2.0 times as much in median, because of the search

new:

    verification costs 4.3, 4.0 and 3.6 times a dimension-2 one in
    the same code at levels I, III and V with two points, and 6.8, 6.5 and
    6.0 times with four, and
    signing costs 3.6, 2.4 and 2.2 times as much in median, because of the search

source: Table 3's ratios, with 4.30 and 6.84 at level I.

## 42. original line 1779 (dated copy 1813)

old:

    the two-point verifier takes 0.63, 0.61 and 0.60 times the cycles of

new:

    the two-point verifier takes 0.63, 0.62 and 0.60 times the cycles of

source: as entry 3.

## 43. original line 1782 (dated copy 1816 to 1817)

old:

    milliseconds is acceptable, while the specification's dimension-2 format

new:

    milliseconds at level I (Table~\ref{tab:ms}), and of 43 to 172
    milliseconds at levels III and V on the same machine, is acceptable, while the specification's dimension-2 format

source: ns columns: 17.2 and 27.3 ms at level I; 42.9 to 172.2 ms at levels III and V.

## 44. original lines 1797 to 1798 (dated copy 1832 to 1836)

old:

    having; there is reason to expect room. MIKE's chain costs about 0.07 megacycles per step
    at level I, against about 0.16 for ours (28.5 megacycles over $2 \cdot

new:

    having; there is reason to expect room. On the core of
    Section~\ref{sec:setup}, MIKE's Rust implementation~\cite{BDD+26} computes a
    level-I shared secret, a chain of 228 steps, in 20.9 megacycles (median of
    200 runs), about 0.09 megacycles per step,
    against about 0.28 for ours (48.0 megacycles over $2 \cdot

source: session/mike-step.csv (paper-dim4/bench/mike-step, MIKE-I of mike_rs 0d79121, CPU 13, 200 runs): shared secret 20.932 megacycles, 0.0918 per step of 228; ours 48.048/174 = 0.276. Both on one core of one machine.

## 45. original line 1800 (dated copy 1838)

old:

    on a different processor and our verifier also evaluates a point for its

new:

    at a prime of 318 bits, in constant time, and our verifier also evaluates a point for its

source: MIKE-I's prime 633 * 2^308 - 1 has 318 bits; mike_rs is constant-time by its README. "On a different processor" no longer applies.

## 46. inserted before original line 1814 (dated copy 1852)

old:

    (nothing)

new:

    % TODO(author): confirm the model name and version in the next sentence.

source: The one marker left for the authors: the model name in the disclosure.

## 47. original lines 2380 to 2384 (dated copy 2419 to 2427)

old:

    where the rule is not defined. Before that, it rejects hints that
    are not the canonical ones of their curves, and it rejects unless
    $e(P, Q)$ has exact order $2^r$ and the discrete logarithm $k$ of
    $e(R, S)$ to that base exists and is odd, that is, unless both
    pairings have exact order $2^r$. Every input that reaches our

new:

    where the rule is not defined. In order, it rejects hints that
    are not the canonical ones of their curves, then rejects unless
    $e(P, Q)$ has exact order $2^r$, then rejects $a$ and $b$ both
    even, then rejects unless the discrete logarithm $k$ of
    $e(R, S)$ to the base $e(P, Q)$ exists and is odd, that is, unless both
    pairings have exact order $2^r$; the entry is recovered last.
    The basis $(P, Q)$ is the image under the challenge isogeny of
    the hinted basis of the public-key curve, and $e(P, Q)$ is
    computed from the pairing of that basis. Every input that reaches our

source: Code: the order of the rejects, as entry 5.

## 48. original line 2405 (dated copy 2448 to 2450)

old:

    of~\cite[Theorem 12]{Dar24}, computed as in (a) and (b), for the

new:

    of~\cite[Theorem 12]{Dar24}, computed as in (a) and (b) except that
    its four candidates for $\mathbf i_0$ are tried in the order $0, 2, 1,
    3$, for the

source: Code: the dimension-2 change of coordinates tries i_0 in the order 0, 2, 1, 3 (dim2.rs).

## 49. original line 2416 (dated copy 2461 to 2464)

old:

    a vanishing value, our crate rejects the signature.

new:

    a vanishing value, our crate rejects the signature: it tests every
    divisor of the dimension-2 chain, of the two gluings, of the plain
    steps and of the image check before using it, and a test per site
    exercises the reject.

source: Code at 6302b7a: every divisor tested, one test per site (zero_divisor_tests in dim2.rs, gluing.rs, isogeny.rs, self_contained.rs).

## 50. inserted before original line 2423 (dated copy 2471 to 2474)

old:

    (nothing)

new:

    \item \emph{The middle comparison.} Our crate evaluates the second
    half-chain through its dual, and compares the null point of the
    first half-chain's last codomain with the Hadamard transform of
    that of the second's, projectively.

source: Code: chain.rs:58-67.

## 51. original lines 2435 to 2436 (dated copy 2487 to 2489)

old:

    evaluations through the gluing and six comparisons, and falls back to four
    points otherwise. These checks would make the two-point chain correct even

new:

    evaluations through the gluing, two projective comparisons of points for
    the images and six comparisons of coordinates for the null point, and
    falls back to four points otherwise. These checks would make the two-point chain correct even

source: Code: gluing_images_paired (two projective comparisons), null_fixed_by_pi (six comparisons), two_point.rs.

## 52. original line 2442 (dated copy 2495)

old:

    equations~\cite[Eq.~(21)]{Dar24} it solves along a single path, which the

new:

    equations~\cite[Eq.~(21)]{Dar24} it solves along a tree, which the

source: Code: the equations are solved along a tree (TREE_EDGES, along_tree).

## 53. original lines 2450 to 2454 (dated copy 2503 to 2504)

old:

    verdict on those vectors and on 70 tampered variants. An earlier version of
    the crate was checked more directly: $\theta(Dx) = \pi(\theta(x))$ held at
    every step of both half-chains of 32 honest and 72 forged signatures
    (17\,790 steps). At levels III and V, the same tests give the same codomain at all
    1\,282 and 1\,708 plain steps of five vectors each, and the same verdict on 75 and 70 tampered variants.

new:

    verdict on those vectors and on 70 tampered variants. At levels III and V, the same tests give the same codomain at all
    1\,274 and 1\,706 plain steps of five vectors each, and the same verdict on 75 and 70 tampered variants.

source: The sentence on 32 honest and 72 forged signatures (17 790 steps) had no record and is cut; tests at head: 1274 and 1706 plain steps at levels III and V.
