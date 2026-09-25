# Experiment C — wrong cross-isogeny rejection test

**Claim under test (ePrint 2026/1305).** The correct cross isogeny
`c: E_aux -> E_com` is not merely sufficient but **necessary**: a wrong `c''`
(same degree `2^a - q`, different kernel) should give a wrong `psi` via eq. (5),
hence a wrong `M`, hence verification failure (codomain fails to split with
`E_com`). This experiment tests that converse.

## Headline finding

**At toy parameters, a genuine wrong cross isogeny EXISTS and PASSES the
split-with-`E_com` check.** The condition "the `(2^a,2^a)`-gluing splits with
`E_com` as a factor" does **not**, by itself, single out the correct cross
isogeny. This is not a new break — it is the same Kani / CCAI-vacuity mechanism
seen in `CCAI_EXHAUSTIVE.md`, and the wrong `c''` yields the **same** `E_com`
(not a different commitment). But it means the necessity of `c` cannot rest on
the split condition alone (see *Verdict*).

## Why (analytic)

With `a = alpha.dual` and `hat(c) : E_com -> E_aux` the degree-`(2^a-q)` cross
isogeny, eq. (5) gives, for `R in E_com[2^a]`:

```
psi(alpha(R)) = - hat(c)( [(2^a-q)^{-1} * q] R ) = hat(c)(R)      (since -(2^a-q)^{-1}*q = 1 mod 2^a)
```

so the Kani kernel is `graph{(alpha(R), hat(c)(R))}`. A wrong cross isogeny just
replaces `hat(c)` by another degree-`(2^a-q)` isogeny `hat(c'') : E_com -> E_aux`.
Then `(alpha, hat(c''))` is *still* a pair of isogenies from the common source
`E_com` with coprime degrees `q + (2^a-q) = 2^a`, so by Kani (1997, Thm 2.3) the
gluing on `E_A x E_aux` **still splits with `E_com` as the source factor**. The
split check therefore cannot distinguish `c` from `c''`.

## Results (real theta (2,2)-chain, ThetaIsogenies/two-isogenies, Sage 10.9)

All instances non-degenerate (`E_A != E_com`). "# cross isogenies" counts
**distinct-kernel** degree-`(2^a-q)` isogenies `E_com -> E_aux` (dedup by kernel;
verified Sage-native for prime degree).

| instance | p | gluing | comp. `q` | cross deg | `E_A != E_com` | # distinct-kernel cross isog. | split with `E_com` | rejected |
|----------|--:|--------|:--:|:--:|:--:|:--:|:--:|:--:|
| small     | 479   | `(2^3,2^3)` | 3  | 5  | yes | **2** | **2 / 2** | 0 |
| larger    | 51839 | `(2^5,2^5)` | 5  | 27 | yes | 2 (image-dedup) | 2 / 2 | 0 |
| prime-cross | 51839 | `(2^5,2^5)` | 27 | 5  | yes | **2** | **2 / 2** | 0 |

For the prime cross-degree cases (`5`), a Sage-native count deduping by
`kernel_polynomial` confirms **exactly 2 distinct-kernel** degree-5 isogenies
`E_com -> E_aux` for *every* reachable target — so the second one is a genuine
different-kernel `c''`, not `c` composed with an automorphism. Both pass the
theta-chain split-with-`E_com` check. **Zero** wrong `c''` were rejected.

So: a wrong `c''` produces a **different** `M` (different anti-isometry — "wrong
`c` -> wrong `M`" holds), but that different `M` **still splits with `E_com`**
("wrong `M` -> failure" does **not** hold at the split level).

## Interpretation / verdict

- **The split-with-`E_com` check does NOT enforce cross-isogeny necessity.** Any
  degree-`(2^a-q)` cross isogeny `E_com -> E_aux` gives an `M` that passes it
  (Kani). At the split level, `c` is *not* necessary.
- **This is not a forgery of a new commitment.** Every alternative `c''` yields
  the same `E_com` (the source factor). It is the cross-isogeny-side view of the
  CCAI malleability already documented — it does not let an adversary reach a
  *different* `E_com`.
- **Where necessity must actually come from.** Two candidates, neither of which
  is the split condition:
  1. **Uniqueness at cryptographic scale.** The multiplicity `2` here is a
     toy-parameter artifact: `E_com = E0` has `j = 1728` (extra automorphisms)
     and the supersingular graph is tiny, so short cycles create multiple
     connecting isogenies of a given degree. In the cryptographic regime the
     graph is a Ramanujan expander and two *generic* curves have `<= 1`
     connecting isogeny of a fixed degree — so no wrong `c''` exists and `c` is
     effectively unique. The reduction should state this uniqueness explicitly.
  2. **Torsion-image binding in full verification.** Real verification also fixes
     the response's action on chosen torsion points; a wrong `c''` gives a
     different `M` that need not be consistent with those images. If necessity is
     meant to hold at *fixed* parameters, it must be argued from this binding,
     **not** from the split condition. (This experiment models only the split
     check — the CCAI/Kani condition — so it cannot see a torsion-image
     rejection; that is a limitation to close.)

**Bottom line for the paper:** the statement "any `M` passing verification must
encode the `psi'` of eq. (5)" is **too strong if "passing verification" means
only the split-with-`E_com` condition** — wrong-kernel cross isogenies also pass
(demonstrated). The necessity claim needs to invoke either uniqueness of the
connecting isogeny at scale, or the torsion-image binding, and the write-up
should say which. The isogeny-path intuition ("computing `c` is hard") is sound
only once uniqueness is established; where `c` is non-unique (as at these toy
primes), multiple valid cross isogenies coexist for the same `E_com`.

## Caveats

- Split oracle = real theta `(2,2)`-chain; "verification" here means the
  split-with-`E_com` condition only (as in Experiments 1/A/B). Torsion-image
  consistency is not modeled.
- Distinct-kernel counts verified Sage-native for **prime** cross degree (5);
  the degree-27 count uses image-level dedup from the library run (naive
  Sage-native chaining over-counts via backtracking).
- Toy multiplicity `2` reflects `j(E_com)=1728` and the small graph; not
  representative of cryptographic parameters.

## Reproduction

```
# SageMath >= 10 + ThetaIsogenies/two-isogenies
cd two-isogenies/Theta-SageMath
sage verify_wrong_cross.py     # scripts/experiments/verify_wrong_cross.py
```
