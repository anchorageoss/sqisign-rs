# The SQIsignHD reference at the round-3 level I prime (P29 step 2)

The oracle for the compact format is the SQIsignHD library's Sage verifier,
unmodified, run at `p = 3·2^324 − 1` with the parameters of
`docs/COMPACT_R3.md` (`e = 174`, `r = 89`, `λ = 128`). Nothing here edits the
library; these files supply parameters, tables and drivers around it.

## Setup

```sh
git clone https://github.com/Pierrick-Dartois/SQISignHD-lib   # master 0287566 was used
cd SQISignHD-lib
git config submodule.Verification/Theta_dim4.url https://github.com/Pierrick-Dartois/Theta_dim4
git submodule update --init                                     # .gitmodules pins an ssh URL
```

Sage 10.8 was used (`sage -python`). The library's own self-test must pass
first: `cd Verification && sage -python Verify.py -lvl 1 -test -i 0`
(level 1 of the library is the round-2 SQIsign prime `5·2^248 − 1`).

## Files

- `gen_tables_r3.sage`: the two non-residue tables the library's
  `torsion_basis_2f_from_hint` indexes, regenerated for the round-3 prime with
  the library's rule (`Signature/scripts/precompute_gf_constants.sage`: seed 0,
  twenty entries each, `F_p^2 = F_p[i]/(i^2 + 1)`). Output in `Data/`, the
  library's text format.
- `verify_r3.py`: the oracle driver. A subclass of the library's parameter
  object with `p, c = 3, e = 324` (its `e` is the 2-adic exponent), `f = 174`
  (its `f` is the embedding exponent, our `e`), `r = 89`, the tables above;
  reads public keys and signatures in the library's `Data/` text format
  (three lines per key: `A_pk`, two hints; eight per signature: `A_com`,
  `a`, `b`, `c_or_d`, `q`, two hints, `chal`) and runs `verify()`, which is
  the middle-codomain match and the image check. `--expect-reject` for the
  negatives.
- `kani_r3_smoke.py`: the real-prime check of the dimension-4 half-chain at
  `e = 174`, `r = 89`, built as the library's own `Tests.py` builds its toy
  cases: a walk off `E_0` to a generic domain, `σ` a walk of `k = 109`
  rational 3-isogenies (`q = 3^109`, 173 bits, `≡ 3 mod 4`), `2^e − q =
  a_1^2 + a_2^2` by Cornacchia on the factored value, a basis of the
  `2^89`-torsion pushed through `σ`, then `KaniEndoHalf` and the two checks
  of `Verify.py`. Run from `SQISignHD-lib/Verification`:
  `sage -python kani_r3_smoke.py <seed> 174`.

## Results, 2026-09-25

- Library self-test at its level 1: passes (1.4 s per signature).
- Library toy tests of the half-chain (`Tests.py --KaniEndoHalf -l_B 3 -i N`):
  index 0 (a 26-bit prime) fails with the library's "product of abelian
  varieties" exception, as its own message anticipates for tiny fields;
  indices 2 and 5 (54- and 208-bit primes) pass.
- `kani_r3_smoke.py` at the round-3 prime, `e = 174`, `r = 89`, seeds 1 and
  2: half-chains of 87 + 87 steps in 1.5 s, middle codomain match `True`,
  image check `True`. The same script at `e = 40` (seed 7) also passes. So
  the library's dimension-4 arithmetic runs unmodified at our field size and
  chain length, and the parameters `e = 174`, `r = 89` are consistent with
  it (`2f ≥ e + 4` in its terms).

## Two facts learned on the way

- The goodness rule (`2^e − q` a prime `≡ 1 mod 4`) forces `q ≡ 3 (mod 4)`.
  An odd sum of two squares is `≡ 1 (mod 4)`, so no endomorphism `a + b·ι`
  of `E_0` can be a response degree, and neither can `q` with `2^e − q` a
  sum of two squares when `q` itself is one (the two would sum to
  `2 mod 4`). A test isogeny must have degree `≡ 3 mod 4`: an odd power of 3
  does.
- A chain whose `σ` passes through `E_0` (an endomorphism conjugated by
  short walks) hit products along the chain at every size tried; a plain
  3-power walk from a generic curve does not. The library's own test walks
  off `E_0` first for the same reason.
