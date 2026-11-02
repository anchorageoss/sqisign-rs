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
appendix is inline from its \appendix line) and paper-dim4/appendix_symmetry.tex
is the appendix as a separate file; both are the human's draft committed
unedited at 8552b9d and then edited by the change log. paper-dim4/iacrtrans.cls,
CC-by.pdf and COPYING-iacrtrans are the IACR Transactions class and its
licence, copied unmodified from Cryptosaurus/iacrtrans at tag v0.94 (commit
cbe3b01, the hashes in its message), build support only.
paper-dim4/CHANGES.md is the full list of edits to main_standalone.tex since
the draft at 8552b9d: every changed line with its text before, its text after
and the source of every number and sentence.

Earlier material, not the paper. paper-dim4/main.tex (682 lines),
paper-dim4/tables/, paper-dim4/figures/ with figures/data/, BENCH.md,
DRAFT_NOTES.md, FIGURES.md and the committed paper-dim4/main_standalone.pdf
of 2026-10-05 (12 pages) are the files of an earlier generated pipeline:
figures/bundle.py wrote the earlier main_standalone.tex from main.tex and
the tables, and they carry the session of 2026-09-26 on an Intel Xeon,
which is not the session the paper now reports. They are not the paper's
figures. They stay as they are until the orchestrator's pending decision on
them (SECENG-2442) is taken. paper-dim4/bench/ is the measurement harness
(its own README.md) and bench/sessions/ holds the session records.

## The head commit

The paper edit is commit 78c2106 (main_standalone.tex alone, applied from
CHANGES.md); the change log's commits are fd5b937 before it and a01a199,
515ce42, fec8401 and 739a400 earlier; the class files are cbe3b01; the
draft is 8552b9d; the code head is d45823d. This README's commit sits on top
of 78c2106 and cannot carry its own hash.

## The session of record

paper-dim4/bench/sessions/2026-11-02-99c13b0/ is the session every number
of the paper's Section 5 comes from: 2026-11-02 07:24:15Z to 07:28:39Z,
build 99c13b0, core 3. Its files: table.md (the rows in megacycles with the
run counts, medians, means, minima and maxima, and the ratios from the
printed cells), numbers.json (the same, machine-readable, with the TSC rate
3049.994 MHz), profile.md (the eight phases and sixteen sub-slots of both
dimension-4 verifiers), costs.md (the cycles per field operation of the
fp2-ops block), machine.txt (the host, the TSC flags, the kernel, the
toolchain, the memory), session.log (one start, one end, every block with
its exit code), raw/ (one CSV per block with every iteration's cycles and
nanoseconds), binary-sha256.txt (the session binary's hash and commit),
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

In paper-dim4/, with the class files of this directory:

    latexmk -pdf -interaction=nonstopmode -f -outdir=build main_standalone.tex

The output is build/main_standalone.pdf (build/ is ignored by git). The
draft's bibliography file Dim4SQI.bib is required and is not on the branch
yet; until it is committed beside the draft, the build ends with undefined
citations and no bibliography: at 78c2106, 0 errors, 0 undefined
references, 98 undefined citations over the draft's 23 cited keys, 23
pages. Do not run make figures, make check or figures/bundle.py: they
regenerate main_standalone.tex from the earlier pipeline's main.tex and
tables and would overwrite the paper.
