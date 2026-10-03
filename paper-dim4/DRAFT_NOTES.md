# Draft notes (W05, SECENG-1997)

Every judgment call the author should look at, and every claim the
drafter was unsure was supportable. The section list is the ticket's; no
section was added.

## Judgment calls

1. **The abstract's size claim.** "The smallest level-I combination we are
   aware of at parameters that stand after the 2026 attacks": 84 + 142 =
   226 bytes (Table 2). The round-2 compact format was 64 + 108 = 172 at
   `5·2^248 − 1`, which the round-3 submission replaced; the abstract's
   qualifier is there for it. The author should decide whether the claim
   stands against schemes this paper does not cite.
2. **The session, not the crate's numbers.** Section 5 reads paper-dim4/BENCH.md's
   session of 2026-09-26, run by `bench/session.sh`: one binary, both
   crate versions, `rdtsc` medians (50 / 200 / 100 runs), the C reference
   at `nist-v3` and `nist-v2` built `broadwell`. The README and the crate's
   own BENCH table keep the criterion session of 2026-09-25; the two agree
   within their spread (verify 12.2 vs 12.3 Mcycles at level I). The
   ticket's "no historical figures" is met by quoting only this session.
3. **Round-2 rows.** The round-2 crate (`v0.4.31`, from GitHub as a
   dependency) is measured as it is: no assembly kernels, `num-bigint`
   signer, no `target-cpu=native` (its own `.cargo/config.toml` does not
   apply to a dependency). Section 5 states this in one sentence and
   draws no comparison between the two crate versions.
4. **The rejection loop.** Measured 114 samples per signature (median
   68, 1 to 835) against the doc's 241. The sampler's norms are odd, so
   the `mod 4` filter passes one in two, not one in four: `1/(e ln 2) =
   121` at level I. Table 5 shows both expectations; `docs/COMPACT_R3.md`
   Section 2 got a measured note saying so. Primality tests per signature,
   57 measured against 60, are unaffected.
5. **Signing median vs mean.** Compact signing is 337 Mcycles median and
   546 mean (50 to 3712): the loop's tail. Table 3 reports medians
   throughout, as the method paragraph says; Table 5 carries the spread.
6. **M4 Pro column.** Pending the author's run (`bench/README.md`). Table 4
   prints "pending" and its caption says so; `\MFourAvailable` is 0. The
   ticket's two-machine table exists in shape; the author fills it.
7. **Memory.** 1 674 KiB peak heap for one compact verification, from a
   counting allocator around the call (heap only; the stack is not
   counted). The doc's estimate for both half-chains kept was about
   1.53 MiB.
8. **Section 6.** ECLIPSE appears once in prose plus the footnote, whose
   text is `docs/NAMING.md` of prism-rs with the citation key changed; the
   section title carries the name because the ticket's section list does.
   The feasibility study is cited in prose ("a feasibility study carried
   out for prism-rs") with `\cite{PRISMrs}`; the repository path was
   removed from the text on 2026-09-26 at the author's request, as was
   the SQIsignHD library's script path in Section 3.
9. **The "only implementation" claim.** Sourced to a GitHub search on
   2026-09-26 (`SQIsignHD`, `SQISignHD`, `Theta_dim4`), which returned
   the SQIsignHD authors' two repositories; no W01 lit-review record on
   dimension-4 code was found in prism-rs. The paragraph says what was
   searched and what would change if something was missed.
10. **SQIsign v1 "signed in seconds".** Stated with the citation to
    DKLPW20 and no number; the author may want the paper's own figure.
11. **Bibliography.** DOIs of the nine venue entries were checked against
    CrossRef on 2026-09-26 (titles, page ranges, authors); the ePrint
    entries against their landing pages (titles, authors, revision dates).
    LNCS volume numbers were checked on 2026-09-26 against the round-3
    specification's own bibliography (12491, 14651, 15486, 14008 all
    agree) and, for SQIsignHD, against the chapter's running header
    (LNCS 14651, pp. 3–32).

12. **Figure 1, panel (b).** Redrawn as a Kani square on the doubled
    curves, the same shape as panel (a). Its degree labels are the
    polarized degrees of Kani's lemma (`sigma + sigma` has degree `q` on
    each factor, `alpha` has degree `2^e - q`, the two sum to `2^e` as
    in panel (a)); as maps, the degrees are `q^2` and `(2^e - q)^2`. The
    matrix of `alpha` and its sum of squares are written inside the
    square beside the left arrow, since outside it the figure exceeds
    the text width. The author may prefer `deg q^2` on the horizontal
    arrows.

13. **Reference audit of 2026-09-26.** Every DOI resolves (CrossRef
    metadata and registered resource URLs); every ePrint page and both
    sqisign.org PDFs resolve with matching titles, author lists and dates;
    the-sqisign tags, the SQIsignHD library commit and its parameter rules,
    and the round-2 crate tag all match. Corrected: the sqisign-rs entry
    paired tag `v0.6.34` with commit `14361c0`, but the tag is at
    `e6ccbc5` (version 0.6.34); paper-dim4/BENCH.md's session measured `14361c0`,
    and the note now says so. The round-3 specification calls itself
    "Version 3.0", so the entry's title now does too. The verifier's two
    acceptance conditions are Section 4.5 of the SQIsignHD ePrint
    (Algorithm 5, IsValid), not 4.3; fixed in two places. The `4λ/3`
    figure is Nakagawa and Onuki's conservative countermeasure for
    PRISM-id (their Section 5.3, `a = 4λ/3`; the aggressive size is
    `7λ/6`), not a general recommendation; the sentence now says so and
    marks the transfer to `e` as ours. The SQIsignHD pointers (Lemma 4,
    Lemma 12, Sections 4, 4.2, 4.4, 4.5, 6.1) were checked against both
    the ePrint of 2024-09-20 and the EUROCRYPT 2024 chapter (the author's
    copy); the numbering agrees in both. `anchorlabsinc/prism-rs` is an
    INTERNAL repository at the time of writing, so the URL in
    `\cite{PRISMrs}` does not yet resolve for readers; the author's
    decision is to keep it, since the repository is released shortly.
    Two entries that were never cited, `DMPR24` and `PRISM25`, were
    removed from `refs.bib` at the author's request.

## Claims the author should read twice

- Section 4, "Why pay it": the fourth scalar from `ad − bc ≡ kq mod 2^r`
  is cited to SQIsignHD Section 6.1; the doc says "equation (2)".
- Section 4, "Moving the parameters": the `2^(2·margin)` candidate count
  and the Nakagawa–Onuki bound are the doc's arguments, quoted not
  re-derived.
- Section 6, paragraph (c): "within a few bytes" and "far sparser" are the
  feasibility study's findings (239 bytes either way; 13 300 against 67
  candidates per hit at level I), stated without numbers as the ticket
  asks.
- Section 7, the SIDH timeline: "three years" from the 2022 attacks to the
  2025 round-2 submission, "one more" to September 2026.
