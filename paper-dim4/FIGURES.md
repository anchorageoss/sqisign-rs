# Figures and tables (W05, SECENG-1997)

Every number in a figure or table is read from a file of this repository
by a script in `figures/`; nothing is typed into a script or into
`main.tex`. `make figures` writes:

- `figures/data/values.tex`: the macros the prose and captions quote
  (sizes, parameters, medians, multiples), from paper-dim4/BENCH.md's session
  "2026-09-26 (the dimension-4 paper)" and `docs/COMPACT_R3.md`;
- `figures/data/fig02_tradeoff.csv`: the series pgfplots reads for
  Figure 2, and `figures/out/fig02_tradeoff.{pdf,svg,png}`, a matplotlib
  rendering of the same data;
- `tables/tab*.tex`: the five tables as booktabs `tabular`s;
- `main_standalone.tex`: `main.tex` with the tables, the macro file and
  the CSV inlined (`figures/bundle.py`), the one file to paste into
  Overleaf next to `refs.bib`.

`make check` runs the generation twice and verifies byte-identical output
(`SOURCE_DATE_EPOCH=0`, no timestamps in PDF or SVG metadata); `make pdf`
renders locally with TeX Live (pdflatex + bibtex, three passes) and
prints the error count and the page count; `make words` greps `main.tex`
for the banned words and for em dashes, comment lines excluded.
Requires `python3` with `matplotlib`.

Style: one accent colour per implementation, identical in
`figures/common.py` and the `\definecolor`s of `main.tex` (C reference
`#4C6A92`, this crate in dimension 2 `#8FB3DA`, this crate in dimension 4
`#2A9D8F`, the round-2 crate `#D98E2B`); sans-serif labels; no titles
inside figures. Figure 1 is TikZ in `main.tex`, with its degree labels in
the neutral colour (`neutralc`, `#6B6B6B`) and its braces in the light one
(`lightc`, `#BFBFBF`), the `neutral` and `light` values of `common.py`;
Figure 2 is pgfplots reading the CSV.

| # | label | kind | script | data source | judgment calls |
|---|---|---|---|---|---|
| 1 | `fig:idea` | TikZ | none (no numbers) | panel (a): the left panel of Figure 1 of the ECLIPSE paper (prism-rs `paper/prism-suf/main.tex`); panel (b): SQIsignHD's Kani square on the doubled curves (DLRW24, Section 4.1) | Padlocks removed. Both panels are Kani squares of one geometry (same node size, arrow style and degree-label style); in (b) the auxiliary side is the endomorphism `alpha = (a1 a2; -a2 a1)`, written inside the square beside the left arrow with `a1^2 + a2^2 = 2^e - q` under it, and light braces mark each column as one curve. Degree labels are the polarized degrees of Kani's lemma (`sigma + sigma`: deg `q` on each factor; `alpha`: deg `2^e - q`), so the two sides sum to `2^e` as in (a). The panel captions say what is sent, not what is computed: the verifier computes the whole dimension-4 chain. |
| 2 | `fig:tradeoff` | pgfplots | `fig02_tradeoff.py` | paper-dim4/BENCH.md session (sizes table; `verify from bytes` and `verify compressed from bytes` rows, level I) | Three points, one crate, medians; labels from the CSV's `format` column through `point meta=explicit symbolic`. |
| T1 | `tab:parameters` | booktabs | `tab1_parameters.py` | `docs/COMPACT_R3.md` Sections 1 and 3 (tables; the public-key sentence); paper-dim4/BENCH.md session sizes table (measured level I) | "not built" for the measured cells at III and V. The `lambda` row is typed from the specification's levels (128/192/256) since no doc carries it as a table. |
| T2 | `tab:sizes` | booktabs | `tab2_sizes.py` | paper-dim4/BENCH.md session, sizes table | Cells are `pk + sig = total` so the level-I combination the abstract claims is visible. |
| T3 | `tab:bench` | booktabs | `tab3_bench.py` | paper-dim4/BENCH.md session, Rust rows and C rows | Medians only. The multiple column is filled for dimension-4 rows only, against the same crate's dimension-2 verify at the same round; the compressed verify is its own row. |
| T4 | `tab:ms` | booktabs | `tab4_ms.py` | paper-dim4/BENCH.md session, `median (ms)` column; the "M4 Pro" table of the same section | The M4 Pro cells read "pending" until the author runs `bench/session.sh` on that machine and pastes its `median (ms)` values into the M4 Pro table (`bench/README.md`). |
| T5 | `tab:loop` | booktabs | `tab5_loop.py` | paper-dim4/BENCH.md session, rejection-loop row and the `sign` rows; `docs/COMPACT_R3.md` Section 2 table | The expected count for odd `q` is `1/(e ln 2)`, computed in the script from `e`; the doc's random-parity figure is shown on its own row. |

Macros (`values.tex`): names carry no digits (TeX), levels are spelled
`I`, `III`, `V`, round 2 is `Rtwo`. `\MFourAvailable` is 0 until the M4
Pro table has numbers.
