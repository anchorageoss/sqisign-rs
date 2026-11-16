# paper-dim4 on the branch agent/dim4-two-point

## What is on the branch

The code. 52 commits on top of main at 7a4a01b, from 950731d (2026-11-01)
to d45823d (2026-11-02), the code head: the dimension-4 two-point verifier
with the image check (crates/verify/src/hd/, four files of
crates/verify/src/fp/, the challenge ladder of crates/verify/src/sqisign.rs)
and the work on its cost down to the counted floor. The Linear goals, one
line each: SECENG-2419, the two-point verifier implemented from the paper's
appendix behind a verifier option, its tests against the four-point
verifier, one benchmark session and the phase profile; SECENG-2426, the
setup phase (the bases from the hints, Cornacchia and the kernel points) and
the field arithmetic path, one session before and one after; SECENG-2433,
the plain steps, the gluing, the image check and the dimension-2 ladder
worked to the counted floor slot by slot, the floor table, one session
after; SECENG-2442, this paper update. The tests of the two-point verifier
are crates/sqisign-rs/tests/two_point.rs and base_change.rs;
crates/kat/tests/tampered.rs tests the dimension-2 verifier's rejection of a
KAT signature with any one byte flipped (SECENG-2438). crates/verify/src/ec/
and the KAT files under crates/kat/kat/ are unchanged from main.

The paper. paper-dim4/main_standalone.tex is the paper's source (the
appendix is inline from its \appendix line). Its base on the branch is the
author's file of 2026-11-10, committed unedited at e29dc36 from
gs://eclipse-revision-ops/refs/dim4/ (UPDATE-2026-11-10.md; sha256
8d84d50ba716ed0a3df71bfa3f31be48dbebceec74ad9eac8fcae1edfacdfa40, 2227
lines) together with its bibliography, paper-dim4/refs.bib and
paper-dim4/Dim4SQI.bib, the same 26 entries under two names (sha256
05765b2710da7394f39387ee81f11fcb310df901badc46fc51b290e9441a0edf; the
\bibliography line names Dim4SQI). The author's file is their own edit of
the paper and is not on top of the branch's earlier edit 78c2106; that edit
and its base 8552b9d (the human's draft of Oct 31) are history.
paper-dim4/appendix_symmetry.tex is the appendix as a separate file, the
human's file of Oct 31 (f5076fda), which main_standalone.tex does not read.
paper-dim4/iacrtrans.cls, CC-by.pdf and COPYING-iacrtrans are the IACR
Transactions class and its licence, copied unmodified from
Cryptosaurus/iacrtrans at tag v0.94 (commit cbe3b01, the hashes in its
message), build support only. paper-dim4/CHANGES.md is version 2 of the
edit log: the disposition of the 33 edits of the earlier log against the
author's text (covered, superseded by the author's rewrite, superseded by
the coming session of record, or PENDING and not applied), the author's
text flagged where a build or a record of this branch disagrees with it,
and one PENDING item; it applies no edit, and no byte of any .tex or .bib
has changed since e29dc36.

Earlier material, not the paper. paper-dim4/main.tex (682 lines),
paper-dim4/tables/, paper-dim4/figures/ with figures/data/, BENCH.md,
DRAFT_NOTES.md and FIGURES.md are the files of an earlier generated
pipeline: figures/bundle.py wrote the earlier main_standalone.tex from
main.tex and the tables, and they carry the session of 2026-09-26 on an
Intel Xeon, which is not the session the paper now reports. They are not
the paper's figures. They stay as they are under the residue decision (a)
of SECENG-2442 (the DECIDED comment headed 11:13Z), which also removed the
generated main_standalone.pdf of 2026-10-05 that came with them from the
branch. paper-dim4/bench/ is the measurement harness
(its own README.md) and bench/sessions/ holds the session records.

## The head commit

The paper's base is commit e29dc36 (main_standalone.tex, refs.bib and
Dim4SQI.bib, the author's bytes of 2026-11-10); the edit log's version 2 is
d068284; the committed PDF and this README follow in one commit on top of
d068284, which cannot carry its own hash. Earlier, now history: the edit
78c2106 and its log fd5b937 (with a01a199, 515ce42, fec8401 and 739a400
before it), against the draft 8552b9d; the class files are cbe3b01. The
code head is d45823d; no paper commit changes code.

## The session of record

paper-dim4/bench/sessions/2026-11-02-99c13b0/ is the session of record on
the branch: the measurement of the branch's code (d45823d) on this host,
2026-11-02 07:24:15Z to 07:28:39Z, build 99c13b0, core 3, and the source of
every figure the earlier edit 78c2106 put in the paper. The author's text at
e29dc36 quotes a session of 2026-11-04 on another machine (its Setup
paragraph) of which no record is on the branch; the next session of record,
one session on one build of the branch head, is the goal SECENG-2617, and
which figures the paper quotes is decided on SECENG-2425. Its files: table.md (the rows in megacycles with the
run counts, medians, means, minima and maxima, and the ratios from the
printed cells), numbers.json (the same, machine-readable, with the TSC rate
3049.994 MHz), profile.md (the phase rows, the harness's decode row and the
code's eight slots, and the sixteen sub-slots of both dimension-4
verifiers), costs.md (the cycles per field operation of the
fp2-ops block), machine.txt (the host, the TSC flags, the kernel, the
toolchain, the memory), session.log (one start, one end, every block with
its exit code), raw/ (26 CSVs: one per table row with every iteration's cycles and
nanoseconds, one for the fp2-ops block, eight for the profiles), binary-sha256.txt (the session binary's hash and commit),
blocks.stdout and the three C reference outputs (tracking only, not cited
in the paper). The method, measured in that session: one build of the
harness from a wiped target, one pinned core, four rounds that visit every
block in turn, medians over the rounds' iterations; nothing but the session
ran on the host inside its window. The rule for every local measurement
between sessions is the same: no chain, build or dump runs while a
measurement runs. The two earlier sessions, 2026-11-01-327c5d2 (the state
before the setup work) and 2026-11-01-d91ab15 (before the floor work), sit
beside it with the same layout, plus session.stdout and without costs.md and
binary-sha256.txt.

## The vectors

The five level-I vectors the committed tests use are not stored files.
crates/sqisign-rs/tests/compact_vectors.rs regenerates them at test time
(lines 17 to 18 and 51 to 60 at d45823d): for i from 0 to 4, a ShakeRng
seeded with the label "compact level I vector i" and the domain "kat", a key
pair from compact_keygen, the message "message i" and a signature from
compact_sign, the port's own signer; the test vectors_level1 verifies each
and, when COMPACT_VECTORS_OUT is set, writes them in the SQIsignHD library's
text format for the Sage oracle. crates/sqisign-rs/tests/two_point.rs
rebuilds the same five from the same seeds (lines 51 to 74) and runs them,
with their tampered variants, through its four tests
(verdicts_agree_on_vectors_and_tampered_variants, codomains_agree_at_every_step,
fall_back_on_unpaired_gluing_images, fall_back_flags_through_rejections).
These vectors are not the reference implementation's known-answer tests;
those are the three files under crates/kat/kat/ for the dimension-2 formats,
unchanged on this branch.

The checked-in files under sqisignhd-harness/round3/Data (Public_keys_r3lvl1.txt,
Signatures_r3lvl1.txt, Signatures_r3lvl1_bad.txt and the two non-residue
tables) are the inputs of the Sage oracle sqisignhd-harness/round3/verify_r3.py,
in the library's text format, committed at 08e9f65 (2026-09-25). The record
(SECENG-2422's REPORT, 2026-11-01): their signatures are not what the signer
produces on this machine at 08e9f65, 5abfdff, 8718c97, 6e1c134 or c6a0c38,
nor with SQISIGN_FORCE_PORTABLE=1; the public keys, commitments, hints and
challenges are identical and the scalars a, b, c_or_d and q differ in every
vector; both sets are accepted by the library's verifier under passagemath
10.8.13. sqisignhd-harness/round3/Data/README.md
(line 6, unchanged since main) records their provenance as
compact_vectors.rs with the seeds it names; the runs of SECENG-2422 did not
reproduce the checked-in signatures from it, and where they came from is the
human's question. The human keeps them as they are (SECENG-2442, record
item 1).

## How to rebuild the PDF

In paper-dim4/, with the class files of this directory, four passes by
hand:

    pdflatex -interaction=nonstopmode -output-directory=build main_standalone.tex
    (cd build && BIBINPUTS=..: bibtex main_standalone)
    pdflatex -interaction=nonstopmode -output-directory=build main_standalone.tex
    pdflatex -interaction=nonstopmode -output-directory=build main_standalone.tex

bibtex runs inside build/ and needs the bib path pointed at paper-dim4/
(BIBINPUTS=..:); run without it, it opens no database and every citation
stays undefined. latexmk (latexmk -pdf -interaction=nonstopmode -f
-outdir=build main_standalone.tex) runs the same passes but stops after
bibtex on the db.png error below, so the four passes are the recipe. The
output is build/main_standalone.pdf (build/ is ignored by git).
paper-dim4/build/main.pdf is the committed PDF of the base at e29dc36 with
the log at d068284, built by these passes, sha256
a6c362a6703c5facfa71f250a560708a0cada5819b82f55954eadbc158cb4991, 29 pages:
0 undefined references, 0 undefined citations (21 entries from Dim4SQI.bib,
bibtex 0 warnings), one overfull box (4.2 pt, the author's paragraph at
lines 733 to 739), and one error: the author's appendix heading at line
1641 includes db.png, a file neither on the branch nor in the bucket, so
pdftex prints the file name in its place (the human has been asked for the
file; nothing of ours changes the .tex). Do not run make pdf, make figures,
make check or figures/bundle.py: the pdf target depends on figures, and
they regenerate main_standalone.tex from the earlier pipeline's main.tex
and tables and would overwrite the paper.
