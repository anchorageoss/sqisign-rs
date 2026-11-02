# CHANGES: the paper update of SECENG-2442, written before any edit (P1, SECENG-2444)

File: paper-dim4/main_standalone.tex at commit 8552b9d (1895 lines, sha256 76dc480b1e96d02482976be6119772a7816559990bc5b2ab237cc7965c9d2138), the draft as the human committed it. Line numbers are those of 8552b9d. paper-dim4/appendix_symmetry.tex (441 lines, f5076fda) is not read by main_standalone.tex (no input or include of it; the appendix of the draft is inline from line 790) and carries no figure of a session, so it is not changed. main.tex (682 lines), refs.bib, tables/ and figures/ are not changed.

Sources. Measured: the session of record, paper-dim4/bench/sessions/2026-11-02-99c13b0/ (table.md, profile.md, costs.md, numbers.json, machine.txt, session.log, raw/), 2026-11-02 07:24:15Z to 07:28:39Z, build 99c13b0, core 3; every ratio is formed from the printed one-decimal cells of table.md. Counted: the SECENG-2433 hand-over (comment 8bfcb27f on SECENG-2433, "Where the floor is") and SECENG-2439's floor-after.md Table C. Tested: the committed tests at d45823d (crates/sqisign-rs/tests/two_point.rs) and the dumps of SECENG-2439. Text: ~/work/results/20261101T023738Z-SECENG-2425/paper-input.md (Part 1 variant 2 under option (i); Part 2's eleven ratio places; Part 3's proposed paragraph as item (f); Part 4's profile sentence) and the SECENG-2442 description. Code reading: the file and line at d45823d (the verifier code of 99c13b0 and 8552b9d is d45823d's).

Registers in the after text: a sentence is measured (it names the session of Section sec:bench), counted (it names the counters at the session's commit), tested (it counts the committed tests), or the draft's own words kept. No sentence of ours claims a proof or a priority. ASCII only; no double hyphen in our text (the draft's own `-{}-` and `Miller--Rabin` are in lines we do not write).

## A. The edits (applied in P2; the scratch diff logs/scratch-diff.txt is this list applied to a copy)

### E01 abstract (Part 2 place 1 of 11; the C reference clause out)

Where: line 48 (sub). Why: the headline ratio from the new cells; the comparison against the reference removed (Rust against Rust only).

Before:
```latex
        The price is verification: measured in one session next to the reference's assembly build, the dimension-4 verifier costs 9.9 times the dimension-2 verifier of the same code (120.3 against 12.2 megacycles).
```
After:
```latex
        The price is verification: measured in one session, the dimension-4 verifier costs 6.4 times the dimension-2 verifier of the same code (85.4 against 13.3 megacycles), or 4.6 times with two kernel points pushed instead of four (60.6 megacycles).
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): dim4 verify four points 85.4, two points with the check 60.6, dim2-I verify 13.3; 85.4/13.3 = 6.42 -> 6.4, 60.6/13.3 = 4.56 -> 4.6 (numbers.json ratios_from_cells: 6.4, 4.6); the two-point clause from paper-input.md Part 2 row 1.

### E02 Contributions bullet (the C reference clause out)

Where: lines 100 to 101 (range). Why: the bullet names a comparison against the reference that the paper no longer makes.

Before:
```latex
        \item measure everything in one session against the reference's
        assembly build, and compare our verifier's cost with MIKE's
```
After:
```latex
        \item measure everything in one session, and compare our verifier's cost with MIKE's
```
Source: SECENG-2442: Rust against Rust only; no number.

### E03 Section 4, The cost (Part 2 places 1 and 2 of 11)

Where: lines 266 to 269 (range). Why: the two ratios from the new cells; the 10.6 against the reference replaced by the two-point dimension-2 ratio (SECENG-2442: reworded to the dimension-2 Rust ratio).

Before:
```latex
    step, 87 steps from each end at level I. It costs 9.9 times the
    dimension-2 verifier of the same crate at the same prime (120.3 against
    12.2 megacycles), and 10.6 times the reference's assembly build
    (Section~\ref{sec:bench}).
```
After:
```latex
    step, 87 steps from each end at level I. It costs 6.4 times the
    dimension-2 verifier of the same crate at the same prime (85.4 against
    13.3 megacycles), or 4.6 times with the two-point verifier of
    Section~\ref{sec:conc} (60.6 megacycles; Section~\ref{sec:bench}).
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 85.4, 60.6, 13.3; ratios 6.4 and 4.6 as in E01.

### E04 Method paragraph (host, date, counter, run counts; the C reference sentence out)

Where: lines 487 to 500 (range). Why: the session of record replaces the 2026-09-26 session; the C reference build and iteration sentence removed with the C rows.

Before:
```latex
    All measurements come from one session on one machine, an Intel Xeon at
    2.80~GHz (8 vCPUs, 31~GiB), on 2026-09-26. One binary measures the Rust
    code, compiled by \texttt{rustc} in release mode with fat link-time
    optimization and no target-CPU flag; the code chooses its assembly field
    kernels at run time from CPUID. Each operation is timed with
    \texttt{rdtsc} and the wall clock around the call, and the tables report
    medians over 50 runs for key generation and signing, 200 for
    verification, and 100 for the compact signer and verifier. The C reference
    is the-sqisign at tag \texttt{nist-v3}, built with CMake in
    \texttt{Release} mode and \texttt{-DSQISIGN\_BUILD\_TYPE=broadwell}, the
    assembly-backed field arithmetic, and measured by its own benchmark
    (\texttt{-{}-iterations=300} at level I, 100 at III and V), medians.
    Verification is measured from bytes: decode the key, decode the signature,
    verify.
```
After:
```latex
    All measurements come from one session on one machine, an AMD EPYC 7B13
    guest (16 vCPUs, 63~GiB; the cycle counter at 3049.994~MHz), on
    2026-11-02. One binary measures the Rust code, compiled by \texttt{rustc}
    1.99.1 in release mode with fat link-time optimization and no target-CPU
    flag; the code chooses its assembly field kernels at run time from CPUID.
    Each operation is timed with \texttt{rdtsc} and the wall clock around the
    call, on one pinned core, in four rounds that visit every block in turn;
    the tables report medians over 52 runs for key generation and signing in
    dimension 2, 200 for verification, and, for the compact format, 52 for key
    generation, 100 for signing and 200 for each verifier. Verification is
    measured from bytes: decode the key, decode the signature, verify.
```
Source: machine.txt lines 7 (16 CPUs), 10 (AMD EPYC 7B13), 47 (tsc: Detected 3049.994 MHz processor), 54 (rustc 1.99.1), 60 (MemTotal 65854756 kB = 62.8 GiB); session.log line 1 (2026-11-02, core 3, rounds 4); table.md runs column (52, 200, 52, 100, 200); fat LTO and no target-CPU flag: paper-dim4/bench/Cargo.toml lines 26 to 28 and .cargo/config.toml line 1 at d45823d; the CPUID dispatch: crates/verify/src/fp/mod.rs line 28 at d45823d; the sentences kept are the draft's own.

### E05 Table 2 note (a) (a figure of the 2026-09-26 session)

Where: line 528 (sub). Why: the date replaced by the session of record, which measured the same 176-byte layout.

Before:
```latex
                \item[a] The compressed format's current layout. The session of 2026-09-26 observed an earlier layout of 176 / 269 / 353 bytes.
```
After:
```latex
                \item[a] The compressed format's current layout. The session of Section~\ref{sec:bench} measured the earlier layout of 176 / 269 / 353 bytes.
```
Source: crates/verify/src/compressed.rs line 39 and crates/sqisign-rs/tests/compressed.rs line 49 at d45823d (176 / 269 / 353 bytes; roundtrip::<Level1>(176)); the 174-byte layout is the draft's own claim and stays.

### E06 Cycles paragraph (Part 2 places 1, 3, 4, 5 and 6 of 11)

Where: lines 535 to 543 (range). Why: the four ratios from the new cells; the three ratios against the reference (1.08, 1.10, 1.12) dropped.

Before:
```latex
    Table~\ref{tab:cycles} gives the session in megacycles. At level I the
    dimension-4 verifier costs 9.9 times the dimension-2 verifier of the same
    crate, and dimension-4 signing costs 3.0 times dimension-2 signing, because
    the good degree is found by rejection (Table~\ref{tab:loop}). Key
    generation is the same in both formats up to the basis convention (41.7
    against 39.0 megacycles). The crate's dimension-2 code runs at 1.08 times
    the reference's assembly build for verification, 1.10 times for key
    generation and 1.12 times for signing; the compressed format's recovery
    adds 18\,\% to a verification.
```
After:
```latex
    Table~\ref{tab:cycles} gives the session in megacycles. At level I the
    dimension-4 verifier costs 6.4 times the dimension-2 verifier of the same
    crate with four kernel points and 4.6 times with two, and dimension-4
    signing costs 2.9 times dimension-2 signing, because the good degree is
    found by rejection (Table~\ref{tab:loop}). Key generation is the same in
    both formats up to the basis convention (40.8 against 38.5 megacycles);
    the compressed format's recovery adds 16\,\% to a verification.
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 85.4/13.3 = 6.4, 60.6/13.3 = 4.6, 329.0/112.5 = 2.92 -> 2.9, 40.8 and 38.5, 15.4/13.3 = 1.158 -> 16 percent.

### E07 Table 3 caption (the C reference sentence out)

Where: lines 547 to 549 (range). Why: no C rows in the table.

Before:
```latex
        \caption{Megacycles at the round-3 primes, medians of the session. The
            last column is the dimension-4 verifier over the dimension-2 verifier of
            the same code. The C reference is the \texttt{broadwell} build.}
```
After:
```latex
        \caption{Megacycles at the round-3 primes, medians of the session. The
            last column is the dimension-4 verifier over the dimension-2 verifier of
            the same code.}
```
Source: SECENG-2442: the C reference rows come out of Table 3.

### E08 Table 3 rows, wholesale (the C rows out; the two-point rows in)

Where: lines 555 to 562 (range). Why: every cell from the session of record; the x dim. 2 column is the dimension-4 verifier over the dimension-2 standard verify of the same session at level I.

Before:
```latex
            C reference, dimension 2 & I   & 35.5  & 101.4 & 11.3 & \\
            C reference, dimension 2 & III & 111.4 & 371.0 & 28.7 & \\
            C reference, dimension 2 & V   & 241.5 & 780.6 & 79.8 & \\
            Rust, dimension 2 & I   & 39.0  & 113.5 & 12.2 & \\
            Rust, dimension 2 & III & 137.8 & 448.1 & 32.7 & \\
            Rust, dimension 2 & V   & 280.1 & 922.6 & 83.0 & \\
            Rust, dimension 2, compressed verify & I & & & 14.4 & \\
            Rust, dimension 4 (compact) & I & 41.7 & 336.6 & 120.3 & 9.9 \\
```
After:
```latex
            Rust, dimension 2 & I   & 38.5  & 112.5 & 13.3 & \\
            Rust, dimension 2 & III & 142.2 & 392.9 & 33.7 & \\
            Rust, dimension 2 & V   & 257.0 & 787.5 & 80.7 & \\
            Rust, dimension 2, compressed verify & I & & & 15.4 & \\
            Rust, dimension 4 (compact), four points & I & 40.8 & 329.0 & 85.4 & 6.4 \\
            Rust, dimension 4, two points, with the image check & I & & & 60.6 & 4.6 \\
            Rust, dimension 4, two points, without the check & I & & & 59.6 & 4.5 \\
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): dim2 I 38.5 / 112.5 / 13.3, III 142.2 / 392.9 / 33.7, V 257.0 / 787.5 / 80.7, compressed I 15.4; dim4 40.8 / 329.0 / 85.4; two-check 60.6; two-nocheck 59.6; ratios 85.4/13.3 = 6.4, 60.6/13.3 = 4.6, 59.6/13.3 = 4.48 -> 4.5.

### E09 Figure trade, the points (Part 2 place 9 of 11)

Where: lines 580 to 584 (range). Why: the verification costs from the new cells; the two-point verifier added as a second mark at 142 bytes; the sizes 142, 174 and 200 are the draft's (the session did not record sizes).

Before:
```latex
                \addplot[only marks, mark=*, mark size=2.5pt, color=teal] coordinates {(142,120.3)};
                \addplot[only marks, mark=*, mark size=2.5pt, color=blue!45] coordinates {(174,14.4) (200,12.2)};
                \node[teal, anchor=west, font=\small] at (axis cs:145,128) {compact (dimension 4)};
                \node[blue!55, anchor=south, font=\small] at (axis cs:168,20) {compressed};
                \node[blue!55, anchor=south west, font=\small] at (axis cs:203,16) {standard};
```
After:
```latex
                \addplot[only marks, mark=*, mark size=2.5pt, color=teal] coordinates {(142,85.4)};
                \addplot[only marks, mark=square*, mark size=2.5pt, color=teal!60] coordinates {(142,60.6)};
                \addplot[only marks, mark=*, mark size=2.5pt, color=blue!45] coordinates {(174,15.4) (200,13.3)};
                \node[teal, anchor=west, font=\small] at (axis cs:145,92) {compact (dimension 4), four points};
                \node[teal!70!black, anchor=west, font=\small] at (axis cs:145,62) {two points};
                \node[blue!55, anchor=south, font=\small] at (axis cs:168,20) {compressed};
                \node[blue!55, anchor=south west, font=\small] at (axis cs:203,16) {standard};
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 85.4, 60.6, 15.4, 13.3; the label positions are layout only.

### E10 Figure trade caption

Where: line 588 (sub). Why: the caption names the added mark.

Before:
```latex
            against verification cost for the three formats of one crate, from the
```
After:
```latex
            against verification cost for the three formats of one crate, the compact one with four and with two kernel points, from the
```
Source: no number; the 176-byte sentence of the caption stays true for the session of record (E05's source).

### E11 Table 4 caption (Part 2 place 10 of 11)

Where: line 602 (sub). Why: the host of the session of record.

Before:
```latex
        \caption{Wall-clock milliseconds at level I on the round-3 prime, medians, on the Xeon at 2.80\,GHz of Section~\ref{sec:bench} and on an M4 Pro. The M4 Pro column is
```
After:
```latex
        \caption{Wall-clock milliseconds at level I on the round-3 prime, medians, on the EPYC 7B13 guest of Section~\ref{sec:bench}, as megacycles over its counter rate of 3049.994\,MHz, and on an M4 Pro. The M4 Pro column is
```
Source: machine.txt lines 10 and 47.

### E12 Table 4 header

Where: line 609 (sub). Why: the host column.

Before:
```latex
            format & operation & Xeon (ms) & M4 Pro (ms) \\
```
After:
```latex
            format & operation & EPYC (ms) & M4 Pro (ms) \\
```
Source: machine.txt line 10.

### E13 Table 4, dimension 2 key generation

Where: line 611 (sub). Why: recomputed.

Before:
```latex
            dimension 2 (standard) & key generation & 13.9 & {\color{red}pending} \\
```
After:
```latex
            dimension 2 (standard) & key generation & 12.6 & {\color{red}pending} \\
```
Source: 38.5 / 3.049994 = 12.6; the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt).

### E14 Table 4, dimension 2 signing

Where: line 612 (sub). Why: recomputed.

Before:
```latex
            & signing & 40.5 & {\color{red}pending} \\
```
After:
```latex
            & signing & 36.9 & {\color{red}pending} \\
```
Source: 112.5 / 3.049994 = 36.9; the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt).

### E15 Table 4, dimension 2 compressed verification

Where: line 614 (sub). Why: recomputed.

Before:
```latex
            & verification, compressed & 5.1 & {\color{red}pending} \\
```
After:
```latex
            & verification, compressed & 5.0 & {\color{red}pending} \\
```
Source: 15.4 / 3.049994 = 5.05 -> 5.0 (median_ms 5.04); the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt).

### E16 Table 4, dimension 4 key generation

Where: line 616 (sub). Why: recomputed.

Before:
```latex
            dimension 4 (compact) & key generation & 14.9 & {\color{red}pending} \\
```
After:
```latex
            dimension 4 (compact) & key generation & 13.4 & {\color{red}pending} \\
```
Source: 40.8 / 3.049994 = 13.4; the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt).

### E17 Table 4, dimension 4 signing

Where: line 617 (sub). Why: recomputed.

Before:
```latex
            & signing & 120.2 & {\color{red}pending} \\
```
After:
```latex
            & signing & 107.9 & {\color{red}pending} \\
```
Source: 329.0 / 3.049994 = 107.9; the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt).

### E18 Table 4, dimension 4 verification

Where: line 618 (sub). Why: recomputed (the four-point verifier, the table's row).

Before:
```latex
            & verification & 43.0 & {\color{red}pending} \\
```
After:
```latex
            & verification & 28.0 & {\color{red}pending} \\
```
Source: 85.4 / 3.049994 = 28.0; the session medians (table.md) over the counter rate 3049.994 MHz (machine.txt line 47), equal to table.md's median (ms) column and to the wall-clock medians of the raw CSVs within 0.1 ms (logs/session-numbers.txt); the two-point rows (19.9 with the check, 19.5 without) are not added, Table 4 keeps its rows.

### E19 The rejection loop paragraph (Part 2 place 11 of 11)

Where: lines 629 to 630 (range). Why: the loop statistics of the session of record.

Before:
```latex
    median signature needs 68 samples, the mean 114, and the range is 1 to
    835.
```
After:
```latex
    median signature needs 79.5 samples, the mean 111.5, and the range is 2 to
    523.
```
Source: dim4-sign extra column over 100 signatures: samples mean 111.5, median 79.5, min 2, max 523; primality tests mean 54.8; signing Mcycles median 329.0, min 48.8, max 1900.9 (raw/dim4-sign.csv, the extra column "response samples;primality tests", 25 signatures per round over 4 rounds).

### E20 Table 5, samples mean

Where: line 644 (sub). Why: the session of record; the expected column unchanged.

Before:
```latex
            response samples per signature, mean & 114 & 121 \\
```
After:
```latex
            response samples per signature, mean & 111.5 & 121 \\
```
Source: raw/dim4-sign.csv, mean of the samples 111.5.

### E21 Table 5, samples median (min, max)

Where: line 646 (sub). Why: the session of record.

Before:
```latex
            response samples per signature, median (min, max) & 68 (1, 835) & \\
```
After:
```latex
            response samples per signature, median (min, max) & 79.5 (2, 523) & \\
```
Source: raw/dim4-sign.csv: median 79.5, min 2, max 523.

### E22 Table 5, primality tests mean

Where: line 647 (sub). Why: the session of record; the expected column unchanged.

Before:
```latex
            primality tests per signature, mean & 57 & 60 \\
```
After:
```latex
            primality tests per signature, mean & 54.8 & 60 \\
```
Source: raw/dim4-sign.csv, mean of the tests 54.8.

### E23 Table 5, signing median (min, max)

Where: line 648 (sub). Why: the session of record.

Before:
```latex
            signing, Mcycles, median (min, max) & 337 (50, 3712) & \\
```
After:
```latex
            signing, Mcycles, median (min, max) & 329.0 (48.8, 1900.9) & \\
```
Source: table.md dim4 sign median 329.0, min 48.8, max 1900.9.

### E24 Table 5, dimension-2 signing median

Where: line 649 (sub). Why: the session of record.

Before:
```latex
            signing, dimension 2, Mcycles, median & 114 & \\
```
After:
```latex
            signing, dimension 2, Mcycles, median & 112.5 & \\
```
Source: table.md dim2-I sign median 112.5.

### E25 Conclusion (Part 2 place 1 of 11, Section 8)

Where: lines 710 to 712 (range). Why: the ratios from the new cells; the heap figure stays (see F1).

Before:
```latex
    costs 9.9 times a dimension-2 one at level I (120.3 against 12.2
    megacycles) and needs about 1.6~MiB of heap where the dimension-2 verifier
    needs none; signing costs three times as much, because of the search for a
```
After:
```latex
    costs 6.4 times a dimension-2 one at level I (85.4 against 13.3
    megacycles), 4.6 times with the two-point verifier below, and needs about
    1.6~MiB of heap where the dimension-2 verifier needs none; signing costs
    2.9 times as much, because of the search for a
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 6.4, 4.6 as in E01; 329.0/112.5 = 2.9.

### E26 The [TBD] passage (paper-input.md Part 1, variant 2, option (i))

Where: lines 735 to 739 (range). Why: variant 2 (the co-author's 24-signature sentence dropped) with the numbers of the session of record under option (i) (one session, so "in the same session" and the host clause collapse into Section sec:bench).

Before:
```latex
    the other two by permuting coordinates. We have implemented this variant;
    it checks the condition at the gluing on each signature and falls back to
    four points otherwise. On 24 signatures it reproduces every codomain of
    both half-chains and every verdict of the four-point verifier, on honest
    and tampered signatures alike; its cost is \textbf{[TBD]}.
```
After:
```latex
    the other two by permuting coordinates. We have implemented this variant
    as an option of our Rust port; it checks the condition at the gluing on
    each signature and falls back to four points otherwise. In the session of
    Section~\ref{sec:bench}, two-point verification with the check costs 60.6
    megacycles at level I against 85.4 for the four-point verifier, 0.71 of
    it, and 59.6 without the check, 0.70. In the port's tests, the two-point
    and the four-point verifier return the same verdict and the same
    rejection reason on the five level-I vectors and on 70 distinct tampered
    inputs, and the same codomain null point at each of the 844 plain steps
    of both half-chains over the five vectors, in both two-point modes.
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 60.6, 85.4, 59.6; 60.6/85.4 = 0.71, 59.6/85.4 = 0.698 -> 0.70 (numbers.json ratios_from_cells 0.71 and 0.7); tested: crates/sqisign-rs/tests/two_point.rs at d45823d, COUNT = 5 vectors (line 27), 15 tampered variants per vector (lines 88 to 122 and 144: 6 + 9), 75 variants of which 70 distinct (the flip of byte 82 repeats the q-bit flip); the plain steps from the codomain dump of SECENG-2439 (runs/codomains-r12-s0.txt, identical to the baseline): 844 two-check and 844 four-point plain-step lines over vectors 0 to 4, compared in both two-point modes by codomains_agree_at_every_step (lines 172 to 211).

### E27 Outlook, the per-step figure and the penalty sentence (Part 2 places 7 and 8 of 11)

Where: lines 766 to 771 (range). Why: the per-step figure restated from the new cells; the penalty sentence gets the measured number beside it; no MIKE sentence added or removed.

Before:
```latex
    reason to expect room. MIKE's chain costs about 0.07 megacycles per step
    at level I, against about 0.69 for ours (120.3 megacycles over $2 \cdot
    87$ steps). The comparison is loose, since MIKE works over $\mathbb{F}_p$
    on a different processor and our verifier also evaluates a point for its
    image check, but it suggests that much of our 9.9-fold penalty lies in the
    implementation rather than in the dimension. MIKE's computation of each
```
After:
```latex
    reason to expect room. MIKE's chain costs about 0.07 megacycles per step
    at level I, against about 0.49 for our four-point verifier and 0.35 for
    the two-point one (85.4 and 60.6 megacycles over $2 \cdot 87$ steps). The
    comparison is loose, since MIKE works over $\mathbb{F}_p$ on a different
    processor and our verifier also evaluates a point for its image check,
    but it suggests that much of our 6.4-fold penalty lies in the
    implementation rather than in the dimension (measured below: 1.18 of the
    counted floor at our own operation count). MIKE's computation of each
```
Source: the session table.md (paper-dim4/bench/sessions/2026-11-02-99c13b0/table.md), Rust rows, median (Mcycles): 85.4/174 = 0.491 -> 0.49, 60.6/174 = 0.348 -> 0.35, 6.4 as in E01; 1.18 from the SECENG-2433 hand-over (comment 8bfcb27f, "Where the floor is") and SECENG-2439 floor-after.md Table C, whole call 60.37 / 50.96 = 1.18.

### E28 Outlook, the profile sentence and the measured floor paragraph (paper-input.md Part 4 with the 99c13b0 profile; the SECENG-2433 hand-over)

Where: after line 777 (insert). Why: inserted after the MIKE paragraph, before "Three years separate".

Before:
```latex
    and vectorization tuned to the Hadamard transform.
```
After (the before lines kept, these lines added below them):
```latex

    In our two-point verifier the plain dimension-4 steps are 67.1 percent of
    the cost (point pushes 44.9 percent, codomain computations 22.2 percent),
    the setup before the chain 19.6 percent, the dimension-2 gluing chain 0.9
    percent, the image check 6.5 percent and the gluing 6.0 percent, measured
    in the session of Section~\ref{sec:bench} (basis generation 7.48
    megacycles, Cornacchia and kernel setup 4.23, dimension-2 steps 0.52,
    gluing 3.60, point pushes 26.88, codomain computations 13.28, image check
    3.87; the whole call 60.37); the codomain computations, 33.1 percent of
    the plain steps, are the part a codomain-from-two-points step would
    address. Its operation count, taken with counters in the field arithmetic
    at the session's commit, is 132\,921 multiplications, 91\,426 squarings,
    435\,363 additions, 54 inversions and 11 square roots in
    $\mathbb{F}_{p^2}$ per verification, 51.0 megacycles at the session's
    per-operation costs; the SQIsignHD library's verifier, counted from its
    code, would take 50.8. Measured in the same session, the verifier runs at
    1.18 of the former floor and 1.19 of the latter, 9.6 megacycles above the
    count, 5.8 of them in the two plain-step phases (1.15 for the point
    pushes, 1.21 for the codomains; 1.12 for the basis, 1.16 for the setup,
    1.35 for the gluing, 1.16 for the image check). What remains is the
    sixteen-coordinate arithmetic of the plain steps on general points: their
    operation count is the library's, and the gap sits inside the Hadamard
    transform and the squaring kernels.
```
Source: profile.md (two points with the image check): shares 44.9, 22.2, 12.5, 7.1, 0.9, 6.0, 6.5 percent (plain steps 44.9 + 22.2 = 67.1; setup before the chain 12.5 + 7.1 = 19.6), medians 7.48, 4.23, 0.52, 3.60, 26.88, 13.28, 3.87, whole call 60.37, "Inside the plain dimension-4 steps: codomain computation 33.1%"; the floor: SECENG-2433 hand-over comment 8bfcb27f ("Where the floor is": 1.18, 51.0, 132921 M, 91426 S, 435363 A, 54 I, 11 sqrt, 1.19, 50.8, 9.6, 5.8, per slot 1.15, 1.21, 1.12, 1.16, 1.35, 1.16, "the count is the library's; the gap is the Hadamard's stage structure and the squaring, inside the kernels") and SECENG-2439 floor-after.md Table C (counts from the counters build at 99c13b0, costs from the session's fp2-ops block; the reference count is a hand count of the SQIsignHD library's code, SECENG-2434 floor.md line 141); 1.18 = 60.37 / 50.96, 1.19 = 60.37 / 50.81, 9.57 -> 9.6, 3.45 + 2.33 = 5.78 -> 5.8.

### E29 Closing paragraph (a date of the 2026-09-26 session)

Where: line 783 (sub). Why: the code measured is the session of record's.

Before:
```latex
    September 2026, not a settled design; that is the reason for caution in
```
After:
```latex
    November 2026, not a settled design; that is the reason for caution in
```
Source: session.log line 1 (2026-11-02).

### E30 Remark conformance, item (f) (paper-input.md Part 3, the proposed paragraph, with the two sentences SECENG-2442 asks for)

Where: after line 1820 (insert). Why: inserted as the last item of the enumerate, after (e); the co-author's (a) to (e) untouched.

Before:
```latex
            Lemma~\ref{lem:basepoints}.
```
After (the before lines kept, these lines added below them):
```latex
            \item \emph{Our Rust port.} Of our Rust port we established the same
            items by reading its code at the commit of the session of
            Section~\ref{sec:bench}: the index convention of (a) line by line
            against the appendix, the index $\mathbf i_0$ of (b), the projective
            use of the kernel points of (d), the base points and the rule of (e),
            and, for (c), the odd $a_1$, the parity and residue tests on $2^e - q$,
            and the bound on $m$. The port's verifier does not repeat the goodness
            test before Cornacchia's algorithm: the only tests on $2^e - q$ are the
            parity, the residue modulo 4 and the success of the descent. It rejects
            a transmitted $a$ and $b$ both even before any inversion. The
            dimension-2 gluing is one chain applied to both pairs of copies;
            whether its structure, the dimension-4 gluing's change of basis and the
            plain steps are those of Definition~\ref{def:structures} is not
            established for the port either, and what the port's tests show is
            only that its two-point and four-point chains agree at every step and
            that the SQIsignHD library's verifier accepts its signatures.
```
Source: paper-input.md Part 3 (the proposed paragraph) with the no-goodness-test sentence expanded and the both-even sentence as SECENG-2442 asks; code reading at d45823d: crates/verify/src/hd/kani.rs lines 53 to 55 (n odd and n = 1 mod 4 or None) and 85 (the odd member returned first), self_contained.rs line 64 (MAX_GLUING_M = 40), 372 (norm_equation_2f_minus_q called from build_setup with no primality test before it; no primality function under crates/verify/src/hd, grep), 376 to 382 (m = trailing zeros of a2, m > 40 rejected as GluingBound), response.rs lines 155 to 159 (a and b both even rejected before any inversion, since 783f783), two_point.rs line 34 (PI = 0 4 8 12 1 5 9 13 2 6 10 14 3 7 11 15), dim4.rs lines 24 and 41 (the base change) and 66 to 81 (the first i0 with a non-zero row), gluing_chain.rs lines 164 and 238 (the kernel generators from columns 4 to 7), isogeny.rs lines 13 to 27 and 249 to 250 (the codomain from each kernel point's own ratios), hd_verify.rs lines 103 to 106 (the dropped entry); items (a), (b), (d) and (e) as read in paper-input.md Part 3 at 2541c3c with SECENG-2420's convention.md, the points and codomain dumps identical from 2541c3c to d45823d (SECENG-2439 REPORT item 6); the acceptance by the SQIsignHD library's verifier is paper-input.md Part 3's statement (the SECENG-2420 harness), not re-run here.

## B. Lines read and left as they are (context, not edits)

- Line 46 and line 82: "September 2026" is the date the round-3 primes were adopted, not the session; unchanged.
- Lines 185 to 186 (Terminology, megacycle): "as in the reference's own benchmark" names the unit's origin, no number; unchanged.
- Lines 207 to 210 (Related work, The submission): the reference's portable and assembly-backed arithmetic, no number; unchanged.
- Line 253: the GitHub search of 2026-09-26 is a search date, not a session figure; unchanged.
- Lines 109, 733 and 759: the red cross-references to Remark rem:conformance ("checked only in part"); the co-author's sentences, still true after item (f); unchanged.
- Lines 413 and 627 (the Miller-Rabin sentences of Section 4 and the loop paragraph), and the co-author's items (a) to (e) of the Remark (lines 1777 to 1820): unchanged, as SECENG-2442 asks.
- Line 645 (Table 5, "the same if q had random parity": the red ???): the co-author's open cell, not a session figure; unchanged.
- Lines 1843 to 1877 (the red Evidence and What the implementation does meanwhile paragraphs: 20 random points, 32 honest and 72 of 96 forged signatures, 5 472, 12 318 and 17 332 steps, 102 signatures): the co-author's own tests, not the session; unchanged.
- Line 613 (Table 4, dimension-2 verification 4.4 ms): recomputed from the session as 13.3 / 3.049994 = 4.36 -> 4.4, the same printed value; no edit.
- Lines 502 to 507 (Sizes) and Table 2's byte counts: sizes of the formats, not measured in the session; unchanged except note (a), E05.
- Lines 587 and 589 to 591 (the figure caption's 174 and 176 bytes): the draft's claim about the compressed layout; the session of record measured the 176-byte layout (E05's source), so the sentence stays true; unchanged.

## C. Open for the human before P2 (not applied; not in the decision list)

- F1, lines 655 to 657 (Memory) and line 711 (the Conclusion's "about 1.6~MiB of heap"): the 1674 KiB peak heap is a figure of the 2026-09-26 session; the session of record has no memory block, so it cannot be replaced from the records. Options: (a) keep the paragraph and name its session ("measured with a counting allocator around the call in a session of 2026-09-26"), which keeps one dated figure in the paper and needs an exception in P2's grep; (b) drop the paragraph and the clause at line 711. Recommendation: (a), the smaller change; the figure is the draft's and the verifier's heap use has not been re-measured.
- F2, lines 755 to 756 (Limitations): "tested at level I only, on 24 signatures, and is not yet the default" refers to the co-author's implementation, whose 24-signature sentence E26 drops. For the port the sentence would read "tested at level I only, on five vectors and their tampered variants, and is not the default" (paper-input.md Part 1, the note on the Limitations paragraph). Not applied: the human says whether the count stays the co-author's or becomes the port's.
- F3, Table 4: the two-point verifier's milliseconds (19.9 with the check, 19.5 without, from 60.6 and 59.6 over 3.049994) are not added as rows, since the decision names the column, not new rows; one line to add if wanted.

## D. Checks

- The session numbers recomputed from the raw CSVs (medians over all rounds; the counter-rate milliseconds against the wall-clock medians), logs/session-numbers.txt:
  - dim2-I-keygen: runs 52, median 38.5 Mcycles, 12.6 ms at the counter rate, wall-clock median 12.6 ms
  - dim2-I-sign: runs 52, median 112.5 Mcycles, 36.9 ms at the counter rate, wall-clock median 36.9 ms
  - dim2-I-verify: runs 200, median 13.3 Mcycles, 4.4 ms at the counter rate, wall-clock median 4.4 ms
  - dim2-I-verify-compressed: runs 200, median 15.4 Mcycles, 5.0 ms at the counter rate, wall-clock median 5.0 ms
  - dim2-III-keygen: runs 52, median 142.2 Mcycles, 46.6 ms at the counter rate, wall-clock median 46.6 ms
  - dim2-III-sign: runs 52, median 392.9 Mcycles, 128.8 ms at the counter rate, wall-clock median 128.8 ms
  - dim2-III-verify: runs 200, median 33.7 Mcycles, 11.0 ms at the counter rate, wall-clock median 11.0 ms
  - dim2-V-keygen: runs 52, median 257.0 Mcycles, 84.2 ms at the counter rate, wall-clock median 84.2 ms
  - dim2-V-sign: runs 52, median 787.5 Mcycles, 258.2 ms at the counter rate, wall-clock median 258.2 ms
  - dim2-V-verify: runs 200, median 80.7 Mcycles, 26.4 ms at the counter rate, wall-clock median 26.4 ms
  - dim4-keygen: runs 52, median 40.8 Mcycles, 13.4 ms at the counter rate, wall-clock median 13.4 ms
  - dim4-sign: runs 100, median 329.0 Mcycles, 107.9 ms at the counter rate, wall-clock median 107.9 ms
  - dim4-verify-four: runs 200, median 85.4 Mcycles, 28.0 ms at the counter rate, wall-clock median 28.0 ms
  - dim4-verify-two-check: runs 200, median 60.6 Mcycles, 19.9 ms at the counter rate, wall-clock median 19.9 ms
  - dim4-verify-two-nocheck: runs 200, median 59.6 Mcycles, 19.5 ms at the counter rate, wall-clock median 19.5 ms
  - dim4-sign extra column over 100 signatures: samples mean 111.5, median 79.5, min 2, max 523; primality tests mean 54.8; signing Mcycles median 329.0, min 48.8, max 1900.9
  - codomain dump: plain-step lines over the five vectors, two-check 844, four 844
- The scratch apply: scripts/gen_changes.py applied the edits of section A to a copy under the results folder (scratch/main_standalone.applied.tex); the unified diff against 8552b9d is logs/scratch-diff.txt; every changed original line is inside a listed range, and lines inside a listed range whose text is the same before and after (the diff skips them): 535 547 548 583 584 766. The branch's file is untouched (git status clean).
- Part 2's eleven ratio places, matched: 1 abstract line 48, Section 4 lines 266 to 268, Conclusion lines 710 to 711 (E01, E03, E25); 2 Section 4 line 268 (E03); 3 Cycles line 537 (E06); 4 Cycles lines 539 to 540 (E06); 5 Cycles lines 540 to 542 (E06, dropped); 6 Cycles line 543 (E06); 7 Outlook line 767 (E27); 8 Outlook line 770 (E27); 9 Figure trade lines 580 to 581 (E09); 10 Table 4 lines 602 to 618 (E11 to E18); 11 Table 5 lines 644 to 649 and the loop paragraph lines 629 to 630 (E19 to E24). The twelfth row of Part 2's table, the memory figure, is F1.
- The C reference in the draft: Table 3 rows 555 to 557 (E08, out), the caption's sentence at 549 (E07, out), the Method sentence at 494 to 498 with its iteration counts (E04, out), the ratios at 540 to 542 (E06, out), "10.6 times the reference's assembly build" at 268 (E03, reworded), "next to the reference's assembly build" at 48 (E01, out), "against the reference's assembly build" at 100 to 101 (E02, out). Every figure of the 2026-09-26 session or the Xeon: lines 48, 266 to 268, 487 to 498, 528, 536 to 543, 555 to 562, 580 to 581, 602, 609, 611 to 618, 629 to 630, 644 to 649, 655 to 657 (F1), 710 to 712, 767, 770, 783; the comment lines of the draft carry none.
