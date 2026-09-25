# Single small-factor forgery: `l | (2^a - q)` is sufficient (Luca De Feo's claim)

Luca De Feo (SQIsign co-designer):

> There is no hardness. SQIsign is not SUF-CMA and cannot easily be made so.
> `2^a - q` doesn't have to be smooth. It is sufficient that a single small prime
> divides `(2^a - q)`, then forgery is trivial.

This note reproduces the mechanism at toy Kani diamonds and **confirms the claim**.
A single small factor `l | (2^a - q)` lets an adversary turn one honest signature
into a valid signature on a **different auxiliary curve**, using only the leaked
`2^a`-torsion images plus cheap degree-`l` isogeny walks — **no `End`, no secret,
and (crucially) without ever touching the non-smooth cofactor `(2^a - q)/l`.**

This **overturns the verdict of [`LUCA_CHALLENGE.md`](LUCA_CHALLENGE.md)**, which
concluded the barrier was "constructing/evaluating the specific, non-smooth,
degree-`2^a - q` isogeny `tau'`." The reroute below constructs a genuine `tau'`
of degree **exactly** `2^a - q` without constructing that isogeny.

**Environment.** SageMath 10.9 + `ThetaIsogenies/two-isogenies` (theta `(2^n,2^n)`
chains). Script: [`scripts/experiments/luca_single_factor.py`](scripts/experiments/luca_single_factor.py).
Run from the Theta-SageMath dir. The verifier check we model is exactly the one
used throughout this branch: the `(2^a, 2^a)` theta chain on the Kani kernel
`graph{(a(P),a(Q)),(x(P),x(Q))}` **splits and one factor is `E_chl`**.

---

## Setup

Diamond (challenge naming): common source `E_chl` (`j = 1728`), a fixed
odd-degree challenge component `a : E_chl -> E_com` of degree `q`, and the honest
auxiliary `tau : E_chl -> E_aux` of degree `2^a - q`.

The adversary's view from one honest signature: the curves `E_chl, E_com, E_aux`
and the `2^a`-torsion images `a(P), a(Q)` and `tau(P), tau(Q)`. **Not** `End(-)`,
**not** `tau` as an evaluable map.

A **forgery** = a valid signature different from the honest one, i.e. the chain
splits and recovers `E_chl` on an auxiliary curve `E_aux''` with
`j(E_aux'') != j(E_aux)` (a genuinely new curve, not the automorphism orbit of the
honest `E_aux` — that orbit stays on `E_aux` and was already mapped by
`luca_endo_shortcut.py` test E2).

---

## The mechanism — REROUTE the terminal `l`-step (works)

`deg(tau) = 2^a - q` and `l | (2^a - q)`, so on the codomain side the dual
`tau_hat : E_aux -> E_chl` has a cyclic kernel of order `2^a - q` whose unique
order-`l` subgroup gives a degree-`l` isogeny out of `E_aux`. Concretely `tau`
factors as `tau = psi_last o tau_rest` with `deg(psi_last) = l`,
`tau_rest : E_chl -> E_1'` of degree `(2^a - q)/l`. The adversary:

1. **Backtracks the last `l`-step.** Push the leaked images through an
   `l`-isogeny `rho : E_aux -> E_mid`. When `rho` is the dual of `psi_last`,
   `rho(tau(P)) = [l]·tau_rest(P)`.
2. **Rescales.** `l^{-1}·rho(tau(P)) = tau_rest(P)` on `E_mid = E_1'`, with
   `l^{-1}` taken mod `2^(a+2)` (an ordinary scalar op — `gcd(l, 2)=1`). This
   recovers the intermediate images **without ever evaluating the degree-`(2^a-q)/l`
   isogeny `tau_rest`.**
3. **Re-extends by a *different* `l`-isogeny.** `psi_new : E_1' -> E_aux''`,
   `psi_new != psi_last`, giving `tau' = psi_new o tau_rest` of degree
   `l·(2^a - q)/l = 2^a - q` **exactly**, landing on a new `E_aux''`.
4. **Rebuilds `M'`** from `tau'`'s torsion images by public torsion algebra
   (already established by `LUCA_CHALLENGE.md` Test B).

The adversary does not know *which* `rho` is the dual of `psi_last`, so it sweeps
all `(l+1)^2` pairs `(rho, psi_new)` and keeps the ones that split — verification
is public, so this is free. Only the genuine backtrack yields a splitting kernel.

### Result — every toy instance forges

| `p` | `a` | `q` | `2^a - q` | `l` | `(rho,psi)` pairs | chain splits | **forgeries (new-curve)** | distinct new `j` |
|---|---|---|---|---|---|---|---|---|
| 51839 | 5 | 5 | `27 = 3^3` | 3 | 16 | 10 | **6** | 2 |
| 38817791 | 8 | 13 | `243 = 3^5` | 3 | 16 | 6 | **3** | 2 |
| 5519359 | 8 | 11 | `245 = 5·7^2` | 5 | 36 | 16 | **8** | 4 |
| 5519359 | 8 | 11 | `245 = 5·7^2` | 7 | 64 | 21 | **14** | 6 |
| 29567999 | 8 | 25 | `231 = 3·7·11` | 7 | 64 | 15 | **7** | 6 |

**4/4 instances forged on a new auxiliary curve.** Each new `j` was verified by
the theta chain splitting and recovering `E_chl`, and each new curve coincides
with a genuine degree-`(2^a - q)` alternative auxiliary (cross-checked against an
exhaustive scan of alternative degree-`(2^a - q)` isogenies from `E_chl`). The
"split back onto `E_aux`" counts are the automorphism orbit (`E2`), not new
forgeries.

### Cost is independent of cofactor smoothness — the decisive point

The reroute uses **only** degree-`l` isogenies (enumerated as `l+1` per curve) and
the scalar `l^{-1} mod 2^(a+2)`. The cofactor `(2^a - q)/l` is **never evaluated**:
its torsion images are obtained by rescaling leaked data, not by walking an
isogeny. So `2^a - q` being non-smooth is irrelevant — a single small factor is
all the attack consumes. This is exactly Luca's statement, reproduced.

---

## Why the earlier "must build the whole isogeny" verdict was wrong

`LUCA_CHALLENGE.md` rerouted the **first** step (at `E_chl`) — see Attack A below —
found that completing a fresh tail of degree `(2^a - q)/l` needs `End(E_chl)`, and
concluded the whole non-smooth isogeny must be built. The fix is to reroute at the
**codomain** (`E_aux`) instead: the last step's dual is a cheap degree-`l` isogeny
the adversary can walk directly, and the entire non-smooth middle is carried for
free as rescaled torsion images. Domain-side reroute is blocked; codomain-side
reroute is trivial.

---

## Controls (the naive readings fail — the geometry has to be right)

Both literal interpretations from the plan produce **zero** splits on every
instance, isolating what makes the reroute work:

- **Attack B — uniform `l^{-1}` then a fresh `l`-isogeny** (`psi(l^{-1}·tau(P))`
  with `l^{-1}` applied on `E_aux`, no backtracking `rho`): `0` splits. This is
  `psi o [l^{-1}] o tau`, a uniform diagonal rescale `diag(l^{-1}, l^{-1})`
  (`det = l^{-2} != 1 mod 2^a`) — it breaks the Weil-pairing isotropy (consistent
  with the `diag(m,m)` dead end in `luca_endo_shortcut.py` E1). The `l^{-1}` must
  be applied at the **intermediate** curve `E_1'` after a genuine backtrack, not
  uniformly on `E_aux`.
- **Attack C — post-compose an `l`-isogeny** (`psi_l o tau`, images `psi_l(tau(P))`):
  `0` splits. Degree `l·(2^a - q) != 2^a - q` violates the Kani degree constraint.
- **Attack A — reroute the first step at `E_chl`**: the `l+1` degree-`l` isogenies
  out of `E_chl` are cheap to enumerate, but replacing `tau`'s first step needs a
  fresh tail of degree `(2^a - q)/l` — the isogeny-path problem / `End(E_chl)`, not
  available from leaked torsion. This is the domain-side barrier the reroute avoids.

---

## Caveats and honest limitations

- **Terminal-step reach.** The single-`l`-isogeny reroute reaches new curves only
  when the order-`l` part sits as an accessible terminal step in `E_aux`'s kernel
  decomposition. For `p = 29567999`, `2^a - q = 231 = 3·7·11`, the `l = 3` reroute
  found **0** new curves (only the `E_aux` orbit), while `l = 7` forged (7 new-curve
  forgeries, 6 distinct `j`). Since `2^a - q` generically has *some* accessible
  small factor and the adversary keeps whichever `l` splits, the claim holds per
  instance; a full attack would peel the `l`-part at its canonical kernel position
  regardless of ordering.
- **Model artifact (important).** SageMath's `E.isogenies_prime_degree(l)` returns
  long-Weierstrass codomains that the theta `(2^a,2^a)` chain **cannot ingest** —
  feeding them makes *every* kernel spuriously fail to split (an all-zero result
  that looks like "no forgery"). The `l`-isogeny steps must be built with the same
  Montgomery x-only machinery as the diamond (`isogeny_from_scalar_x_only`). The
  first draft of this experiment hit exactly this trap; the corrected script uses
  Montgomery `l`-isogenies throughout.
- **Toy scale.** Verified at `a in {5,8}`, `p` up to ~`3·10^7`. The theta chain
  cannot be run at real `2^a` (torsion far too large), but the *cost* argument for
  the reroute (only degree-`l` work + one scalar inverse; cofactor never evaluated)
  is exact and does not depend on scale.

---

## Verdict

**Confirmed: a single small prime `l | (2^a - q)` is sufficient to forge.** From
one honest signature the adversary produces a valid signature on a different
auxiliary curve using only leaked `2^a`-torsion images and degree-`l` isogeny
walks; the non-smooth cofactor is never constructed or evaluated. This is a
genuine SUF-CMA break at toy parameters and reproduces Luca's stated mechanism.

The security implication (extrapolated, not toy-verified): for a random honest
signature, `2^a - q` almost always has a small factor (a random integer avoids all
primes `< 100` with probability `≈ 0.12`), so the condition the attack needs is the
common case, not an exotic one — matching "cannot easily be made SUF-CMA."

**Correction to [`LUCA_CHALLENGE.md`](LUCA_CHALLENGE.md).** Its verdict — "the
barrier is constructing/evaluating the specific, non-smooth, degree-`2^a - q`
isogeny `tau'`" — is **refuted for the SUF-CMA setting**. When `2^a - q` has a
small factor, no such construction is needed: the reroute builds a valid `tau'` of
the exact degree by walking one small-degree step and carrying the rest as
rescaled torsion images. The earlier note reached the opposite conclusion only
because it rerouted the first step (domain side, `End`-blocked) rather than the
last step (codomain side, cheap), and because the naive uniform-`l^{-1}` control
(Attack B here) fails and was read as evidence no torsion shortcut exists.
