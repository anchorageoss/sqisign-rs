# CHANGES: every difference between main_standalone.tex and paper-2026-11-14-standalone.tex

Original: paper-dim4/main_standalone.tex (untouched, sha256 ce4b629c8449064bc4ac35921952077651b4c5d26852b255607cf05356b87661). Dated copy: paper-dim4/paper-2026-11-14-standalone.tex. One entry per changed block of lines, in file order: the original's lines, the old text, the new text, and the source of the new text. Entries marked "audit" were made or changed by the audit of 2026-11-15; the others were made on 2026-11-14 and 2026-11-15 from the session and were rechecked by the audit (AUDIT.md). Session files are under ~/work/results/20261114T224715Z-SECENG-2636-clean/session/ and, in the branch, paper-dim4/bench/sessions/2026-11-15-68a0017/.

## 1. original line 48 (dated copy 48)

old:

    SQIsignHD verifies a dimension-4 isogeny where SQIsign verifies a dimension-2 one, which removes the auxiliary curve from the signature but makes verification much costlier. We show that the verifier's chain of dimension-4 isogenies respects an automorphism of order 4 that acts on its theta coordinates as one fixed permutation, for every input, honest or forged. After the gluing, two of the four kernel points can therefore be recovered from the other two by permuting coordinates, and the constants of each step are determined by ten coordinates instead of sixteen, which cuts the verifier's cost by 37 to 40%. Along the way we derive the change of basis that the verifier applies at the gluing, which Dartois states with only a brief justification, prove that the verifier's image check remains sound with the signs that its implementation leaves free, and show that the verifier must validate the basis hints, which a forger controls. We implement the scheme at the round-3 SQIsign primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs four times that of the dimension-2 scheme in the same code.

new:

    SQIsignHD verifies a dimension-4 isogeny where SQIsign verifies a dimension-2 one, which removes the auxiliary curve from the signature but makes verification much costlier. We show that the verifier's chain of dimension-4 isogenies respects an automorphism of order 4 that acts on its theta coordinates as one fixed permutation, for every input, honest or forged. After the gluing, two of the four kernel points can therefore be recovered from the other two by permuting coordinates, and the constants of each step are determined by ten coordinates instead of sixteen, which cuts the verifier's cost by 37 to 40%. Along the way we derive the change of basis that the verifier applies at the gluing, which Dartois states with only a brief justification, prove that the verifier's image check remains sound with the signs that its implementation leaves free, and show that the verifier must validate the basis hints, which a forger controls. We implement the scheme at the round-3 SQIsign primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs 4.3 times that of the dimension-2 scheme in the same code.

source: audit. Two points over dimension 2 at level I: 48.048 / 11.185 = 4.30 (medians of three passes, dim4-I-verify-two.csv and dim2-I-verify.csv). "four times" described the 4.0 of the paper of Nov 4.

## 2. original line 55 (dated copy 55)

old:

    We implement the scheme at the round-3 {\SQIsign} primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs four times that of the dimension-2 scheme in the same code.

new:

    We implement the scheme at the round-3 {\SQIsign} primes, adopted after Wesolowski's attack. At level I a signature takes 139 bytes, against 200 for the round-3 specification, and key and signature together take 223 bytes, the smallest level-I combination of any post-quantum signature we are aware of at parameters that stand after that attack; verification costs 4.3 times that of the dimension-2 scheme in the same code.

source: audit. The same sentence in the typeset abstract.

## 3. original line 109 (dated copy 109)

old:

    0.63, 0.61 and 0.60 times the cycles of the four-point one at levels I,

new:

    0.63, 0.62 and 0.60 times the cycles of the four-point one at levels I,

source: two/four at level III: 120.059 / 194.632 = 0.617 (dim4-III-verify-two.csv, dim4-III-verify-four.csv).

## 4. original lines 527 to 529 (dated copy 527 to 529)

old:

    coordinates instead of the four of a dimension-2 step, and costs 6.4 times
    the dimension-2 verifier of the same crate at the same prime (45.6 against
    7.1 megacycles), and 6.6 times the reference's assembly build

new:

    coordinates instead of the four of a dimension-2 step, and costs 6.8 times
    the dimension-2 verifier of the same crate at the same prime (76.5 against
    11.2 megacycles), and 7.0 times the reference's assembly build

source: dim4-I-verify-four 76.540, dim2-I-verify 11.185, c-I-verify 10.982 megacycles; 76.540/11.185 = 6.84, 76.540/10.982 = 6.97.

## 5. original line 531 (dated copy 531)

old:

    after the gluing. A profile of the verifier at level I attributes 59\,\%

new:

    after the gluing. A profile of the verifier at level I attributes 60\,\%

source: dim4-I-profile-four.csv, slot4 over total: 59.7, 59.6, 59.8 % in the three passes.

## 6. original line 533 (dated copy 533)

old:

    them, with the other saved points, through the chain, and 16\,\% to

new:

    them, with the other saved points, through the chain, and 15\,\% to

source: dim4-I-profile-four.csv, slot5 over total: 15.0, 15.1, 15.1 %.

## 7. original line 1180 (dated copy 1180)

old:

    the two-point verifier took 0.63, 0.61 and 0.60 times the cycles of the four-point one at levels I, III and V, on

new:

    the two-point verifier took 0.63, 0.62 and 0.60 times the cycles of the four-point one at levels I, III and V, on

source: as entry 3.

## 8. original lines 1192 to 1195 (dated copy 1192 to 1196)

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

source: audit. Code: along_tree (crates/verify/src/hd/two_point.rs:140-171), 25 multiplications, called for the dual null point and for its inverse; hadamard_pi_assign (two_point.rs:74-82), 32 additions and 4 doublings.

## 9. original line 1200 (dated copy 1201)

old:

    when a doubling follows, which happens after about 42\,\% of the steps of an

new:

    when a doubling follows, which happens after about 43\,\% of the steps of an

source: audit. The code's strategy replayed (audit/strat.py): 36 of 85, 84 or 83 codomains for m = 1, 2, 3; 360 of the 844 plain steps of the five level-I test vectors = 42.7 %.

## 10. inserted before original line 1378 (dated copy 1379)

old:

    (nothing)

new:

    % TODO(author): every measured figure of this paragraph (4 180 commitments, the medians and maxima of the minima, the tail constant 1.6 and its range, the 1.9 to 2.8 of the first minimum) comes from an earlier version of the crate and was not re-measured; the instrumentation is not in the crate at the head of agent/dim4-two-point. The probabilities 2^-44, 2^-46, 2^-47 here, the 2^44 of the introduction, the 2^-44 of the conclusion and the 2^40 of the restart paragraph are computed from that 1.6.

source: audit. Marker for the authors: the Gram-matrix figures were not re-measured. Comment line only.

## 11. original line 1380 (dated copy 1382 to 1383)

old:

    and computed its successive minima for 4\,180 commitments at level I.

new:

    and computed its successive minima for 4\,180 commitments at level I
    (with an earlier version of the crate; not re-measured for this one).

source: 2026-11-14, at the author's word: the figure is marked as not re-measured. No instrumentation at head.

## 12. inserted before original line 1497 (dated copy 1500)

old:

    (nothing)

new:

    % TODO(author): "over two hundred hint pairs per signature" comes from an earlier version of the crate and was not re-measured; no test at the head of agent/dim4-two-point counts re-encodings with the hint check off.

source: audit. Marker for the authors: the hint-pair count was not re-measured. Comment line only.

## 13. original line 1498 (dated copy 1502 to 1503)

old:

    over two hundred hint pairs per signature at level I. With it, each

new:

    over two hundred hint pairs per signature at level I (with an earlier
    version of the crate; not re-measured for this one). With it, each

source: 2026-11-14, at the author's word: the figure is marked as not re-measured. No such test at head.

## 14. original lines 1536 to 1538 (dated copy 1541 to 1544)

old:

    with one vCPU of an Intel Xeon at 2.30~GHz and 3~GiB of memory, on
    2026-11-04. One binary measures the Rust code, compiled by \texttt{rustc}
    1.89 in release mode with fat link-time optimization and no target-CPU

new:

    with 32 vCPUs of an Intel Xeon at 2.80~GHz and 126~GiB of memory (Google
    Compute Engine, \texttt{n2-standard-32}), with the measuring process
    pinned to one core whose sibling was kept idle, on 2026-11-15. One binary measures the Rust code, compiled by \texttt{rustc}
    1.100 in release mode with fat link-time optimization and no target-CPU

source: session/machine.txt: 32 CPUs, Intel Xeon 2.80 GHz, MemTotal 132002464 kB, n2-standard-32, pinned core 13 with sibling 29, rustc 1.100.0; passes on 2026-11-15.

## 15. original line 1550 (dated copy 1556)

old:

    (\texttt{-{}-iterations=300} at level I, 100 at III and V), medians.

new:

    (\texttt{-{}-iterations=300} at level I, 100 at III and V), medians. Every row was measured in three passes of the session: key generation and signing of the C reference at levels III and V, whose benchmark draws a fresh seed on every run, differed by more than 1.5\,\% between the first two passes, so a third was run, and the tables report the median of the three passes' values; no row of the Rust code moved by more than 1.2\,\%.

source: numbers.json: three passes; the four C reference rows that differed by more than 1.5 % between passes 1 and 2. audit: 1.1 % became 1.2 % (dim4-III-sign moved by 1.13 % over the three passes).

## 16. original lines 1552 to 1557 (dated copy 1558 to 1562)

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

source: raw headers say "packed layout"; profile columns checks_hints_pairing and checks_goodness: 1.24, 3.18, 5.51 megacycles. audit: the shares, 1.6 / 1.6 / 1.1 % of a four-point verification and 2.6 / 2.6 / 1.9 % of a two-point one.

## 17. original lines 1589 to 1593 (dated copy 1594 to 1598)

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

source: Table 3's ratios from the unrounded medians: 6.84 and 3.60 at level I; 6.53 and 6.03; 4.03 and 3.64; 2.43 and 2.23.

## 18. original line 1595 (dated copy 1600)

old:

    the signatures took 300 to 1289 megacycles at level III and 567 to 1596 at

new:

    the signatures took 398 to 1471 megacycles at level III and 869 to 3392 at

source: quartiles of the 52 signatures of dim4-III-sign.csv and dim4-V-sign.csv, median of the three passes: 398.3 and 1470.9; 868.7 and 3391.95. audit: 3391 (pass 1 alone) became 3392.

## 19. original lines 1598 to 1601 (dated copy 1603 to 1606)

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

source: dim4-I-keygen 35.269, dim2-I-keygen 32.938; Rust over C at level I: 1.019, 1.056, 1.088; compressed over standard verify 1.186. audit: 1.05 became 1.06 (unrounded medians).

## 20. original lines 1618 to 1630 (dated copy 1623 to 1635)

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

source: every cell: medians recomputed from pass{1,2,3}/raw/*.csv and c-ref-*.txt, median of the three passes (audit/recompute.out).

## 21. original line 1640 (dated copy 1645)

old:

    xmin=120, xmax=240, ymin=0, ymax=60,

new:

    xmin=120, xmax=240, ymin=0, ymax=90,

source: figure layout: the y axis has to reach 76.5.

## 22. original line 1644 (dated copy 1649)

old:

    xtick={120,140,...,240}, ytick={0,20,40,60},

new:

    xtick={120,140,...,240}, ytick={0,30,60,90},

source: figure layout: ticks for the new axis.

## 23. original lines 1648 to 1654 (dated copy 1653 to 1659)

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

source: the four points are Table 2's bytes against Table 3's level-I cells (76.5, 13.3, 11.2, 48.0); labels moved next to their points.

## 24. original lines 1679 to 1682 (dated copy 1684 to 1687)

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

source: ns column of dim2-I-*.csv, median of three passes: 11.764, 35.724, 3.995, 4.738 ms.

## 25. original lines 1684 to 1687 (dated copy 1689 to 1692)

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

## 26. original lines 1697 to 1698 (dated copy 1702 to 1703)

old:

    the spread of signing times: the median signature needs 80 samples, the
    mean 107, and the range is 1 to 566.

new:

    the spread of signing times: the median signature needs 95 samples, the
    mean 109, and the range is 1 to 446.

source: extra column of dim4-I-sign.csv: median 95, mean 108.94, min 1, max 446 samples.

## 27. original line 1712 (dated copy 1717)

old:

    response samples per signature, mean & 107 & 121 \\

new:

    response samples per signature, mean & 109 & 121 \\

source: as entry 26.

## 28. original lines 1714 to 1717 (dated copy 1719 to 1722)

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

source: dim4-I-sign.csv: median 95 (1, 446); 54.87 primality tests; signing 360.54 megacycles, min 45.51, max 1528.9 (medians of the three passes); dim2-I-sign 100.02. audit: the max 1546 (pass 1 alone) became 1529.

## 29. original lines 1723 to 1724 (dated copy 1728 to 1729)

old:

    1673~KiB (1.63~MiB) of heap with four points and 1637~KiB with two,
    measured with a counting allocator around the call, and 3474 and 6271~KiB

new:

    1655~KiB (1.62~MiB) of heap with four points and 1631~KiB with two,
    measured with a counting allocator around the call, and 3443 and 6160~KiB

source: session/heap.csv, means of five signatures: 1654.7, 1630.7, 3442.8, 6159.9 KiB.

## 30. original lines 1774 to 1777 (dated copy 1779 to 1782)

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

source: Table 3's ratios, as entry 17, with 4.30 and 6.84 at level I.

## 31. original line 1779 (dated copy 1784)

old:

    the two-point verifier takes 0.63, 0.61 and 0.60 times the cycles of

new:

    the two-point verifier takes 0.63, 0.62 and 0.60 times the cycles of

source: as entry 3.

## 32. inserted before original line 1797 (dated copy 1802)

old:

    (nothing)

new:

    % TODO(author): the next two sentences set our megacycles per step, measured on this session's Xeon, against MIKE's 0.07, measured on another processor (and over F_p, in constant-time C). Megacycles from two machines do not compare. The audit of 2026-11-15 left the wording alone; an earlier pass had already replaced the original's "about 0.16 for ours (28.5 megacycles" by this session's 0.28 and 48.0.

source: audit. Marker for the authors: a comparison of megacycles across machines. Comment line only.

## 33. original line 1798 (dated copy 1804)

old:

    at level I, against about 0.16 for ours (28.5 megacycles over $2 \cdot

new:

    at level I, against about 0.28 for ours (48.0 megacycles over $2 \cdot

source: 2026-11-15 (before the audit): 48.048 / 174 = 0.276. The audit left the sentence as it found it and marked it (entry 32).

## 34. original line 2442 (dated copy 2448)

old:

    equations~\cite[Eq.~(21)]{Dar24} it solves along a single path, which the

new:

    equations~\cite[Eq.~(21)]{Dar24} it solves along a tree, which the

source: audit. Code: the equations are solved along a tree (TREE_EDGES, along_tree, two_point.rs).

## 35. original line 2454 (dated copy 2460)

old:

    1\,282 and 1\,708 plain steps of five vectors each, and the same verdict on 75 and 70 tampered variants.

new:

    1\,274 and 1\,706 plain steps of five vectors each, and the same verdict on 75 and 70 tampered variants.

source: tests at head: 1274 and 1706 plain steps at levels III and V (audit/tests-head-8747fa9.log).
