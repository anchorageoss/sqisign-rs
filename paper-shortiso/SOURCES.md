# Item 6: the claims about other schemes, with sources

Looked up on 2026-11-27. The MIKE paper itself (ePrint 2026/2377, "MIKE: a
fast and compact post-quantum NIKE") could not be read from this machine:
eprint.iacr.org answered HTTP 403 to every fetch. Every MIKE finding below
rests on the authors' public repositories (github.com/tensor-MIKE), not on
the paper's text. Lines marked "fetched" were read here from the raw files;
the others were reported by a search agent and not re-read.

## 1. SQIsign starts from y^2 = x^3 + x. Confirmed.

- Reference code, SQISign/the-sqisign, main,
  `scripts/precomp/precompute_E0_basis.sage` line 13 (fetched):
  `E0 = EllipticCurve(Fp2, [1, 0])`.
- Specification 2.0.1 (sqisign.org/spec/sqisign-20250707.pdf), Section 3.2,
  as reported: "In SQIsign, we will always use a prime p = 3 (mod 4). In
  this case, the curve E0 : y^2 = x^3 + x is supersingular, with
  endomorphism ring isomorphic to O0".
- This repository: `crates/verify/src/precomp/mod.rs` line 2, and
  `crates/sqisign-rs/src/precomp/p324_3.rs` line 11 ("Action of
  `i: (x, y) -> (-x, i y)`") and line 41 ("`j: (x, y) -> (x^p, y^p)`"),
  which are the paper's iota and pi.

## 2. MIKE starts from y^2 = x^3 - x. Confirmed from the code only.

- tensor-MIKE/mike_py, main, `mike/mike.py` line 117 (fetched):
  `E0 = EllipticCurve(self.Fp2, [-1, 0])`.
- tensor-MIKE/mike_c, main, `scripts/precompute_e0_basis.py` lines 38 and
  41 (fetched): `E0 = EllipticCurve(Fp2, [1, 0])` for the even torsion
  basis and `E0p = EllipticCurve(Fp2, [-1, 0])` for the secret-key basis.
  Both have j = 1728; they are twists over F_p and isomorphic over F_{p^2}.
- The sentence of the paper that states the curve: not seen.

## 3. "The three SQIsign primes are also proposed for MIKE." Not confirmed.

- None of 3*2^324-1, 27*2^500-1, 17*2^664-1 is a parameter set in mike_c,
  mike_rs or mike_py (`mike_c/scripts/mike_parameters.py`, fetched, lists
  six primes and no others).
- The paper's text was not readable, so the claim is neither confirmed nor
  refuted for the paper. If the sentence stays, it needs a section
  reference into ePrint 2026/2377.
- For what it is worth, the dimension-4 draft in this repository (branch
  agent/dim4-two-point, paper-dim4/main_standalone.tex lines 341 to 343)
  says the opposite direction: "MIKE also argues that the round-3 primes
  are conservative and proposes log2 p = 317, 484 and 651 instead".

## 4. The six MIKE primes. Confirmed from the code.

`mike_c/scripts/mike_parameters.py` lines 3 to 8 (fetched) give `(f, c)`
for `p = c*2^f - 1`, and `mike_c/README.md` assigns the families and
levels ("fast: for primes p308_633, p474_593 and p628_317 for NIST levels
I, III and V respectively", and the same for "conservative").

| prime | family | level |
|---|---|---|
| 633*2^308 - 1 | fast | I |
| 593*2^474 - 1 | fast | III |
| 317*2^628 - 1 | fast | V |
| 117*2^374 - 1 | conservative | I |
| 77*2^566 - 1 | conservative | III |
| 41*2^758 - 1 | conservative | V |

No other prime appears in the three repositories.

## The (D, m, h(-D)) table, recomputed

From `crates/sqisign-rs/tests/shortiso_table.rs` (results/item6.log): D
and m follow from each prime by Corollary B and agree with the nine pairs
given for the paper; `p = D m^2 - 1` holds exactly; h(-D) is the number of
reduced primitive forms. The paper's own h column was not available for
comparison.

| prime | D | m | h(-D) | order |
|---|---|---|---|---|
| 3*2^324 - 1 | 3 | 2^162 | 1 | maximal |
| 27*2^500 - 1 | 27 | 2^250 | 1 | conductor 3 in Q(sqrt(-3)) |
| 17*2^664 - 1 | 68 | 2^331 | 4 | maximal |
| 117*2^374 - 1 | 468 | 2^186 | 8 | conductor 3 in Q(sqrt(-13)) |
| 77*2^566 - 1 | 308 | 2^282 | 8 | maximal |
| 41*2^758 - 1 | 164 | 2^378 | 8 | maximal |
| 633*2^308 - 1 | 2532 | 2^153 | 20 | maximal |
| 593*2^474 - 1 | 2372 | 2^236 | 24 | maximal |
| 317*2^628 - 1 | 1268 | 2^313 | 10 | maximal |
