# paper-dim4 on the branch agent/dim4-two-point

What is on the branch. 52 commits on top of main at 7a4a01b: the
dimension-4 two-point verifier with the image check (crates/verify/src/hd/,
with four files of crates/verify/src/fp/ and the challenge ladder of
crates/verify/src/sqisign.rs; the Linear goals SECENG-2419, SECENG-2426 and
SECENG-2433), its tests (crates/sqisign-rs/tests/two_point.rs and
base_change.rs, crates/kat/tests/tampered.rs), the session harness under
paper-dim4/bench/ and the records of three sessions under
paper-dim4/bench/sessions/. crates/verify/src/ec/ and the KAT files are
unchanged from main. The paper files (main.tex, main_standalone.tex,
refs.bib, tables/) are unchanged from main; the paper update is the Linear
goal SECENG-2442.

The head commit. d45823d is the last code and records commit; this README
is the one commit above it.

Where the session of record lives. paper-dim4/bench/sessions/2026-11-02-99c13b0/
(2026-11-02 07:24:15Z to 07:28:39Z, build 99c13b0, the binary's sha256 in
binary-sha256.txt): table.md, profile.md, costs.md and numbers.json are
printed from the raw CSVs under raw/ by paper-dim4/bench/session_report.py;
machine.txt and session.log describe the host and the run. The two earlier
sessions, 2026-11-01-327c5d2 and 2026-11-01-d91ab15, sit beside it.
