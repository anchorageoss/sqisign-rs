# Luca De Feo's challenge: can `M'` be recomputed for a new auxiliary curve without the secret?

Luca De Feo (SQISign issue #18) responded to the claim that recomputing the
basis-change matrix `M` for a *new* auxiliary curve `phi_aux` requires the secret
key:

> In many cases, yes [a third party can recompute `M`]. That's the case, for
> example, when `phi_aux` has small-degree components.

Our proof rests on equation (5): computing `M'` for a new auxiliary needs the
cross isogeny `c': E_aux' -> E_com`, whose recovery is the isogeny-path problem.
This note tests Luca's counter-claim experimentally.

**Environment:** SageMath 10.9 + `ThetaIsogenies/two-isogenies` (theta `(2^n,2^n)`
chains, `isogeny_diamond.py`). Scripts in the scratchpad; real parameters read
from `crates/verify/src/params` (`F_CHR`/`TORSION_EVEN_POWER = a`, `E_RSP = q`-size).

---

## Test A — is the auxiliary degree `2^a - q` smooth? (DECISIVE, real params)

The adversary must build `tau'` of degree **exactly** `2^a - q`, where `a` is the
`2`-power torsion and `q` is the (odd) response/challenge degree fixed by the
honest signature. Luca's "small-degree components" attack needs that degree to be
smooth (a chain of small isogenies). The adversary does **not** control `q`.

Method: for each level, sample 40 odd `q` of the honest size and trial-divide
`2^a - q` up to `2^20` (never fully factor — a ~124-bit prime factor is the point).

| Level | `a` | `q` bits | `2^20`-smooth | rough-cofactor size (median / min / max) |
|-------|-----|----------|---------------|------------------------------------------|
| L1 | 248 | 126 | **0 / 40** | 238 / 186 / 248 bits |
| L3 | 376 | 190 | **0 / 40** | 360 / 322 / 376 bits |
| L5 | 500 | 250 | **0 / 40** | 481 / 446 / 500 bits |

`2^a - q` is **never smooth** at real parameters — it always carries a large
(~186–248-bit at L1) prime factor. Consistent with a Dickman-`ρ` estimate: a
random ~248-bit integer is `2^20`-smooth with probability ≈ 10⁻¹⁵.

**⇒ Luca's specific "small-degree-components" construction cannot be instantiated
at real SQIsign parameters**, regardless of the torsion-level math below.

---

## Test B — is `psi` computable from `(a, tau)` torsion images *without* `c'`?

Toy Kani diamond `DIAMONDS[0]`: `a = 9` (`2^a = 512`), `q = 3^5 = 243` (odd),
auxiliary degree `2^a - q = 269 = 10² + 13²`. The kernel of the splitting
`(2^a, 2^a)`-isogeny is `{(tau(P), a(P)), (tau(Q), a(Q))}` on `E_chl × E_com` —
built **purely from the `2^a`-torsion images of `a` and `tau`**. No cross isogeny
`c'` is used anywhere in its construction (this is exactly what the reference
`generate_splitting_kernel` does).

Result: the `(2^9, 2^9)` theta chain on this torsion-only kernel **completes and
splits** into a product of elliptic curves.

```
[Test B] a=9 (2^a=512), q=3^5=243, aux_deg=269=10^2+13^2
  kernel = {(tau(P),a(P)),(tau(Q),a(Q))} on E_chl x E_com  (torsion images only, no c)
  theta chain COMPLETED: valid (2^9,2^9) product isogeny; split codomain obtained
```

**⇒ The Kani gluing kernel — equivalently `psi`, equivalently `M` — is computable
from the `2^a`-torsion images of `a` and `tau` alone. It does NOT require `c'`.**
This **confirms Luca's torsion-level claim** and contradicts the proof's premise
that computing `M'` needs the cross isogeny.

---

## Test C — is `psi` uniquely determined by `(a, tau)`?

Yes. By Kani's theorem the gluing that makes the `(2^a, 2^a)`-isogeny split with
`E_com` as a factor is the graph of the anti-isometry fixed by the torsion images
of `a` and `tau`; it is unique up to the automorphism orbit. Empirically, the
single torsion-derived kernel of Test B is a valid Kani kernel (the chain splits),
and the Weil-pairing anti-isometry constraint pins the isometry class. So `psi` is
not merely *computable* from `(a, tau)` (Test B) but *determined* by it — no `c'`
enters.

---

## Test D — the incremental small-degree-chain attack

Gated (per the plan) on `2^a - q` being smooth. It is not:

- **Real parameters:** Test A — never smooth (large prime factor). `tau'` cannot
  be assembled from small-degree isogenies.
- **Toy Kani diamonds:** the construction forces `2^a - q = X² + Y²`, and for the
  worked example that is `269`, a **prime**. So even here `tau'` would be a single
  large-degree (269) isogeny, not small components — the incremental chain does
  not apply.

To *evaluate* `tau'` on the `2^a` torsion (to obtain the images Test B feeds on),
the adversary must actually compute the degree-`(2^a - q)` isogeny. With a
non-smooth degree that is the isogeny-path / endomorphism problem — infeasible
without the secret. **The incremental attack does not instantiate.**

---

## Verdict

**Luca is right at the torsion level, and the proof's stated barrier is wrong —
but the scheme is still secure, for a different reason.**

1. **The proof's equation-(5) argument is not the real barrier.** Tests B/C show
   `psi`/`M'` are computable and uniquely determined from the `2^a`-torsion images
   of `a` and `tau'` *without* the cross isogeny `c'`. An adversary who already
   holds `tau'(P), tau'(Q)` can recompute `M'` by public computation. The claim
   "recomputing `M'` requires `c'`" does not hold.

2. **The real barrier is producing `tau'` itself.** To obtain those torsion images
   the adversary must construct and evaluate an auxiliary isogeny of degree exactly
   `2^a - q`. Two obstructions, both empirical here:
   - **Degree is fixed and non-smooth** (Test A): `q` comes from the honest
     signature, so `2^a - q` is a fixed odd integer with a large prime factor at
     every level — it cannot be walked as a chain of small-degree isogenies
     (Luca's "small-degree components" case never arises).
   - Evaluating a single isogeny of that non-smooth degree is the isogeny-path
     problem, i.e. requires `End(E_chl)` / the secret.

3. **Attack does not apply at real parameters.** Confirmed by Test A + Test D.

**Recommendation for the proof.** Restate the hardness assumption: the barrier is
not "recomputing `M'` needs `c'`" (false), but "constructing/evaluating an
auxiliary isogeny of the *specific, non-smooth* degree `2^a - q` needs the secret."
The security conclusion stands; the argument should be repaired to rest on the
non-smoothness of `2^a - q` and the hardness of computing a fixed-degree isogeny,
not on the third party's inability to compute `psi` from torsion data.

### Status of each test
- Test A: run, decisive (real params). ✔
- Test B: run at toy `a = 9`; torsion-only kernel splits. ✔
- Test C: Kani-theorem argument + Test B empirics. ✔
- Test D: argued not to instantiate (non-smooth degree at real params; prime aux
  degree in the toy diamond). A contrived smooth-degree toy diamond was **not**
  constructed — it would not change the real-parameter verdict.

---

## Test E — can the adversary *fabricate* `tau'` torsion images from leaked data?

Tests B/C settled that *if* the adversary holds `tau'(P), tau'(Q)` for some valid
`tau'` they can rebuild `M'` by public algebra. The last gap: can they synthesise
those images from one honest signature **without** constructing the degree-`2^a-q`
isogeny `tau'`? Three algebraic shortcuts, run on the non-degenerate toy diamonds
`p = 479` (`a=3`, `q=3`, `aux=5`) and `p = 51839` (`a=5`, `q=5`, `aux=27`), where
`E_com ≠ E_A` and the honest torsion-only kernel splits and recovers the source.
Script: `scripts/experiments/luca_endo_shortcut.py` (SageMath 10.9 + theta chain).

### E1 — endomorphism composition changes the degree ✔ (dead end, confirmed)

`α ∘ tau` has degree `deg(α)·deg(tau)`; Kani needs **exactly** `2^a-q`, so
`deg(α)` must be `1`. Frobenius (`deg p`) and `[m]` (`deg m²`) never match:

| composition | degree (`p=51839`) | target `2^a-q` | match |
|---|---|---|---|
| `frob ∘ tau` | `1 399 653` | 27 | no |
| `[2] ∘ tau` | 108 | 27 | no |
| `[3] ∘ tau` | 243 | 27 | no |

Empirically, feeding `[m]`-scaled auxiliary images `(m·tau(P), m·tau(Q))` into the
chain — a diagonal `diag(m,m)` recombination — **never** produces a
source-recovering split, even when `det = m² ≡ 1 mod 2^a` (e.g. `m=3` at `a=3`).
So isotropy of the kernel is *necessary but not sufficient*. Dead end confirmed.

### E2 — `SL₂(ℤ/2^aℤ)` recombination of `tau`'s images ✔ (no new auxiliary)

Enumerate every `M ∈ SL₂(ℤ/2^aℤ)` (det `= 1` preserves the Weil-pairing isotropy),
form `tau'(P) = M·(tau(P), tau(Q))` on `E_aux` with `a`'s images held fixed, and run
the chain. Full enumeration is feasible at these sizes (`|SL₂(ℤ/8ℤ)|`→256 tried,
`|SL₂(ℤ/32ℤ)|`→16384 tried).

| `p` | recombinations that split **and** recover source | non-identity |
|---|---|---|
| 479 (`N=8`) | 4 (incl. identity) | 3 |
| 51839 (`N=32`) | 4 (incl. identity) | 3 |

Only a **tiny, fixed** set works (the automorphism / Kani-orbit — `4` regardless of
`N`, i.e. **not** growing with the torsion), and — the decisive point — **every**
recombined image set lives on the **same** `E_aux`. Linear algebra on a fixed
curve's torsion cannot manufacture a *different auxiliary curve*, which is exactly
what a CCAI forgery would require. No new attack vector.

### E3 — cross-surface jump: does the leaked data yield `End(E_aux)`? ✔ (refined)

The endomorphism `θ = tau ∘ â ∘ c : E_aux → E_aux` (degree `q·(2^a-q)²`) is the
worry: one endomorphism of `E_aux` gives `End(E_aux)` via one-endomorphism
reductions. Three findings:

- **Q1 — the challenge's own barrier guess is *wrong*.** The note supposed `θ`'s
  image "leaves the `2^a`-torsion" because its order is "divisible by
  `2^a·q·(2^a-q)`". But `a`, `â`, `c` all have **odd** degree, so each is an
  **isomorphism on the `2^a`-torsion**; `θ` maps `E_aux[2^a]` **bijectively to
  itself** (verified on the real `a`, `tau`: `ord(P)=ord(a P)=ord(tau P)=2^a`).
  The "order = degree" reasoning conflated the two. **Refuted.**

- **Q2 — the `2^a`-action of `θ` *is* computable from one signature.** Because `a`
  is odd-degree, `a(E_chl[2^a]) = E_com[2^a]` (full), so eq. (5) leaks `ĉ` on *all*
  of `E_com[2^a]`. The adversary then recovers `c|_{2^a} = (2^a-q)·ĉ⁻¹`,
  `â|_{2^a} = q·a⁻¹`, and composes with the leaked `tau|_{2^a}`. The reconstructed
  `θ|_{2^a}` **equals** the true restriction (exact, both toy primes). So the
  adversary gets *even more* than the challenge feared — the endomorphism's whole
  `2^a`-torsion action, for free.

- **Q3 — but that is *not* `End(E_aux)`, and this is the real barrier.** `θ|_{2^a}`
  is just a `2×2` matrix mod `2^a`. Reconstructing the actual isogeny from its
  `N`-torsion action (SIDH/Kani-style) needs `N² > deg`. Here
  `deg(θ) = q·(2^a-q)² ≈ N²·2^{a/2}`, so `N² > deg` **fails** — at toy `p=479`
  already (`75 > 64` false), and at real L1 by a factor `≈ 2^{124}`. The torsion is
  far too small to pin the degree-`q(2^a-q)²` isogeny, so no evaluable endomorphism
  and no one-endomorphism attack. The knowledge is a torsion *restriction*, not the
  map.

*(Caveat: no genuine degree-`aux` cross `E_aux→E_com` exists in these
independently-built toy diamonds — an artifact of constructing `E_com`, `E_aux`
separately, not a scheme property — so `c` is an invertible representative. Q1 is
anchored on the real `a`, `tau`; Q2 is exact composition algebra
(`â=q·a⁻¹`, `c=aux·ĉ⁻¹`) and Q3 uses the integer degree `q·aux²` — both independent
of the specific matrix.)*

### E verdict

None of the three shortcuts manufactures a valid new `tau'`:
- **E1**: any non-trivial pre/post-composition breaks the exact degree `2^a-q`.
- **E2**: `SL₂` recombination stays on the *same* `E_aux` — no new auxiliary curve.
- **E3**: the leaked data *does* give `θ`'s `2^a`-torsion action (correcting the
  challenge's stated barrier), but a torsion restriction of a degree-`q(2^a-q)²`
  map is not `End(E_aux)` — the degree dwarfs `N²`, so the isogeny is unreconstructible.

This **reinforces the Test A–D verdict**: the security barrier is not "computing
`M'`/`psi` needs `c'`" (false — Tests B/C, and now E3-Q2 give even the endomorphism's
torsion action), but **constructing/evaluating the specific, non-smooth,
degree-`2^a-q` isogeny `tau'` itself**, which the leaked torsion data cannot
shortcut. The proof should rest on that hardness, not on the equation-(5) argument.
