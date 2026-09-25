# Exhaustive CCAI search at toy parameters

**Question (ePrint 2026/1305, Assumption 5.7 / CCAI).** For a fixed
`(E_A, E_com, q)` and a fixed component isogeny `a: E_A -> E_com` of degree
`q`, how many distinct auxiliary isogenies `tau: E_A -> E_aux` of degree
`2^a - q` (distinct `j(E_aux)`) produce a Kani diamond whose `(2^a, 2^a)`-isogeny
codomain **splits with `E_com` as a factor**?

- If always exactly 1 (only the honest `tau`) → CCAI is a vacuous assumption.
- If it can be `> 1` **because splitting is automatic from Kani's theorem** →
  the `E_com` factor is forced by `a` alone, so CCAI still needs no bespoke
  hardness (see *Verdict*).

## TL;DR verdict

**CCAI is vacuous.** For a fixed component `a` (hence a fixed `E_com`), the
`(2^a, 2^a)` gluing splits with `E_com` as a codomain factor for **every**
auxiliary `tau`, because `E_com` is the common *source* of the two isogenies
being glued and Kani's theorem (Kani 1997, Thm 2.3) forces the source to
reappear as a codomain factor. The recovered `E_com` is therefore **invariant
under the choice of `tau`** — it is determined by `a`, not by the auxiliary.
Producing a different commitment-compatible `tau'` is trivially always possible
but yields the *same* `E_com`, so it confers no forgery advantage and the
SUF-CMA proof needs no dedicated CCAI hardness assumption.

This is the outcome anticipated in the experiment brief's "critical note", now
confirmed by direct computation of the `(2^a, 2^a)`-isogeny codomains.

## Method

The direct abelian-surface splitting check is essential — the shortcut of
returning `codomain(a).j == j(E_com)` is a **tautology** and proves nothing.
Stock SageMath 9.5 (Debian bookworm) has no `(2,2)`-isogeny / Richelot / product
gluing machinery (verified: no `sage.rings.generic.ProductTree`, no Jacobian
Richelot, no `torsion_basis`), so the codomain of the gluing cannot be computed
there.

We therefore used the SQIsign2D-West theta `(2,2)`-isogeny library
(`ThetaIsogenies/two-isogenies`, Dartois–Maino–Pope–Robert, "An Algorithmic
Approach to (2,2)-isogenies in the Theta Model"), which computes the full
`(2^n, 2^n)`-isogeny chain between elliptic products and returns the split
codomain factors. It requires SageMath ≥ 10; we installed **SageMath 10.9** via
conda-forge for this.

`compute_richelot_chain(kernel, n, N_constant, strategy)` takes a kernel on
`E1 × E2 = codomain(a) × codomain(tau)` given as the graph of the Kani
anti-isometry, computes the `(2^n, 2^n)`-chain, and returns the two elliptic
codomain factors `(F1, F2)` — or reports non-splitting.

**Diamond construction (reference-faithful).** Both constituent isogenies
emanate from a common source `S`; the kernel is the graph
`{(a-image(R), tau-image(R)) : R in S[2^(a+2)]}`. In the constructible toy
setup the source `S` is the commitment curve, realized as `E0` (`j = 1728`):
the component `a = X + Y·iota` is an endomorphism of `E0` of degree
`C = X^2 + Y^2`, and the auxiliary is a variable isogeny `phi_B` of degree
`B = 3^{eb}` to `E_aux`, with `C + B = 2^{ea}`. Varying the auxiliary secret
gives many distinct `E_aux`; the component (hence `E_com = E0`) is held fixed.

**Toy primes** are of the form `p = f · 4 · 2^{ea} · 3^{eb} − 1` (so the
`2^{ea+2}`- and `3^{eb}`-torsion are `F_{p^2}`-rational and the component degree
`C = 2^{ea} − 3^{eb}` is a sum of two squares, realizable by `X + Y·iota`).

## Results

For each parameter set: a fixed component/`E_com`, many distinct auxiliaries
`tau` (→ distinct `E_aux`), the real `(2^a,2^a)`-chain, and a **control** using
a non-diamond (random) anti-isometry kernel.

| param | prime size | gluing        | component deg `C` | aux deg `B` | distinct `E_aux` tested | split with `E_com` factor | control (non-diamond) splits? |
|------:|-----------:|---------------|-------------------|-------------|-------------------------|---------------------------|-------------------------------|
| 0     | 19-bit     | `(2^9, 2^9)`   | 269               | 243 = 3^5   | 16                      | **16 / 16 (100%)**        | no                            |
| 1     | 139-bit    | `(2^72, 2^72)` | 4.68e21           | 3^41        | 16                      | **16 / 16 (100%)**        | no                            |
| 2     | 254-bit    | `(2^126,2^126)`| 8.45e37           | 3^75        | 4                       | **4 / 4 (100%)**          | no                            |

- **Every** distinct auxiliary `tau` tested produces a split codomain with
  `E_com` as a factor. The compatible-`tau` count equals the number of distinct
  `E_aux` reached — i.e. **all** of them, not 1.
- The recovered `E_com` (`j = 1728` here) is **identical across all auxiliaries**
  — the auxiliary only changes the *other* codomain factor.
- The **control** confirms splitting is not trivially automatic for *arbitrary*
  kernels: a random (non-diamond) anti-isometry kernel of the same shape fails
  to split (the chain aborts / returns no product). Splitting is a genuine
  property of *genuine `(component, auxiliary)` diamonds* — which is exactly
  what Kani's theorem characterizes.

## Interpretation

Kani's theorem (1997, Thm 2.3): given two isogenies `alpha, beta` from a common
source `S` with `deg(alpha) + deg(beta) = 2^a` and `gcd(deg alpha, deg beta) = 1`
(here `q` odd, `2^a − q` odd, `gcd = gcd(q, 2^a) = 1` — **always**), the induced
`(2^a, 2^a)`-gluing on `codomain(alpha) × codomain(beta)` is reducible and its
codomain splits as `S × E'`. The source `S = E_com` is therefore *always* a
codomain factor, for *every* auxiliary.

Consequently:

- The number of commitment-compatible `tau` is **not 1** — it is essentially
  every auxiliary of the correct degree.
- But the `E_com` factor is **determined by the component `a`** (as the source),
  **independent of `tau`**. An adversary who swaps `tau` for a different `tau'`
  recovers the *same* `E_com`.

So the CCAI "problem" (produce a second compatible `tau'`) is trivially always
solvable, yet gives no control over `E_com` and hence no forgery leverage: the
commitment is pinned by `a`, and this pinning is a *theorem* (Kani), not a
hardness assumption. **Assumption 5.7 (CCAI) is vacuous**, and the SUF-CMA proof
needs no bespoke CCAI hardness assumption for recovery of the `E_com` factor.

## Non-degenerate verification (E_A != E_com)

The result above uses a component that is an endomorphism (`E_A = E_com = E0`).
Experiment B removes that caveat by choosing primes whose odd torsion has **two
coprime smooth factors**, so both the component and the auxiliary are genuine
smooth-degree isogenies from a common source `E_com`:

```
    alpha : E_com -> E_A     (degree q,       FIXED)   => E_A = codomain(alpha) != E_com
    beta  : E_com -> E_aux   (degree 2^a - q, VARIED)  => many distinct E_aux
```

so the component `a = alpha.dual()` is a genuine isogeny between *distinct*
curves, not an endomorphism. The gluing on `E_A x E_aux` has kernel
`graph{(alpha(R), beta(R))}`; by Kani the source `E_com` is a codomain factor.
We confirm this on the real theta `(2,2)`-chain, varying the auxiliary, with
`E_A != E_com` throughout.

| `p` | `p+1` | gluing | comp. deg `q` | aux deg | `j(E_A)` | `j(E_com)` | `E_A != E_com` | distinct `E_aux` | all split w/ `E_com` | control splits |
|----:|-------|--------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 479   | `2^5·3·5`   | `(2^3,2^3)` | 3 | 5  | 365   | 291 (`=1728 mod 479`) | yes | 3  | **yes** | no |
| 51839 | `2^7·3^4·5` | `(2^5,2^5)` | 5 | 27 | 41633 | 1728 | yes | 18 | **yes** | no |

Across all **21 non-degenerate auxiliaries** (`E_A != E_com` and `E_aux != E_A`),
the `(2^a,2^a)`-gluing splits with `E_com` as a codomain factor; every
**non-diamond control** kernel fails to split. This matches the degenerate case
and confirms the Kani source-factor mechanism does **not** depend on
`E_A = E_com`. The E_A = E_com caveat is therefore removed.

## Caveats and scope

1. **`E_A = E_com` — RESOLVED (Experiment B, above).** The main result's
   constructible primes force the component to be an *endomorphism* of `E0`
   (degree `C = X^2 + Y^2`, via `X + Y·iota`), giving `E_A = E_com = E0`. The
   *Non-degenerate verification* section removes this: with primes carrying two
   coprime smooth odd-torsion factors, the component is a genuine isogeny and
   `E_A != E_com`, yet the theta chain still splits with `E_com` as a factor for
   every auxiliary. (At toy scale `End(E_A)` is always computable; the structural
   point that matters here — genuine-isogeny component, `E_A != E_com` — holds.)
2. **Special toy primes.** We used `p = f·4·2^{ea}·3^{eb} − 1` (required for
   rational torsion + sum-of-two-squares component in the theta library), not the
   plain `p = f·2^a − 1` supersingular-graph enumeration originally sketched.
   This is the price of using validated `(2,2)` machinery rather than a
   from-scratch (error-prone) implementation.
3. **"Exhaustive" = over sampled auxiliaries** (16/16/4 distinct `E_aux`), not
   over all supersingular curves. The Kani argument covers all auxiliaries; the
   computation is a representative, uniformly-passing sample plus a
   negative control.
4. **Tooling.** Stock Sage 9.5 could not compute the gluing; SageMath 10.9 +
   `ThetaIsogenies/two-isogenies` was used.

## Reproduction

```
# SageMath >= 10 (we used 10.9 via conda-forge)
git clone https://github.com/ThetaIsogenies/two-isogenies.git
cd two-isogenies/Theta-SageMath
# place ccai_experiment.py here (scripts/experiments/ccai_experiment.py)
sage ccai_experiment.py 16 0,1     # 16 auxiliaries each, params 0 and 1
sage ccai_experiment.py 4 2        # larger (2^126,2^126) gluing, 4 auxiliaries
```

The experiment script (`ccai_experiment.py`) fixes the component, varies the
auxiliary over many secrets, builds the reference Kani kernel, runs
`compute_richelot_chain`, and checks whether `E_com` is among the split codomain
factors; it also runs the non-diamond control.
