"""
Experiment: single small-factor SUF-CMA forgery (Luca De Feo's claim).

Luca (co-designer) stated SQIsign is not SUF-CMA and cannot easily be made so:
2^a - q need NOT be smooth -- it is *sufficient* that a single small prime l
divides (2^a - q), and then forging a second valid signature is trivial.

This is sharper than the "small-degree components" claim refuted in
LUCA_CHALLENGE.md (which needed 2^a - q fully smooth, and concluded the barrier
is constructing the whole non-smooth degree-(2^a - q) isogeny tau). Here we test
whether ONE small factor already suffices, at toy Kani diamonds where
l | (2^a - q) for l in {3,5,7}.

Diamond (challenge naming):
    E_chl (common source, j=1728)
    a  : E_chl -> E_com   degree q          (challenge component, odd, fixed)
    tau: E_chl -> E_aux   degree 2^a - q     (honest auxiliary)
Kani gluing kernel graph{(a(P),tau(P)),(a(Q),tau(Q))} on E_com x E_aux splits and
recovers E_chl. The verifier check we model = the (2^a,2^a) theta chain on that
kernel splits AND one factor is E_chl.

The adversary's view: E_chl, E_com, E_aux, a(P),a(Q), tau(P),tau(Q). NOT End(-),
NOT tau as an evaluable map (only its 2^a-torsion images).

We run four attacks. Each keeps a's images fixed and produces alternative
AUXILIARY-side torsion images; a forgery = the chain splits AND recovers E_chl
AND lands on a DIFFERENT auxiliary curve (j != j(E_aux)) -- i.e. a genuinely new
signature, not the honest one nor its automorphism orbit (which stays on E_aux,
already mapped by luca_endo_shortcut.py E2).

  Attack REROUTE (the mechanism):  backtrack the last l-step of tau via an
      l-isogeny rho: E_aux -> E_mid  (rho(tau P) = l * tau_rest(P) when rho is the
      dual of the honest last step), rescale by l^{-1} mod 2^(a+2), re-extend by a
      DIFFERENT l-isogeny psi: E_mid -> E_aux''. deg(tau') = l*(2^a-q)/l = 2^a-q
      EXACTLY, on a new curve. Uses only leaked images + (l+1)-isogeny enumeration.

  Attack B (prompt, literal): uniform l^{-1} on E_aux then a fresh l-isogeny,
      images = psi(l^{-1} * tau P). This is psi o [l^{-1}] o tau -- a uniform
      diagonal rescale (det = l^{-2} mod 2^a != 1) composed with an l-isogeny;
      predicted to fail the Weil isotropy (cf. E1 diag(m,m) dead end).

  Attack C (prompt, literal): psi' = psi_l o psi, images = psi_l(tau P). Degree
      l*(2^a-q) -- wrong; predicted to fail.

  Attack A (prompt): enumerate l-isogenies out of E_chl (candidate first step of
      tau); note the tail-completion barrier.

  Brute force (a=4): only if REROUTE finds nothing -- enumerate auxiliary curves
      and anti-isometries exhaustively.

Run from the theta-SageMath dir:
    cd /home/user/isogeny/two-isogenies/Theta-SageMath
    sage /home/user/sqisign-rs/scripts/experiments/luca_single_factor.py
"""

import sys
sys.path.insert(0, ".")
from sage.all import *  # noqa: F401,F403
from richelot_isogenies.richelot_isogenies import compute_richelot_chain
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only, evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.strategy import optimised_strategy
import json


# --------------------------------------------------------------------------
# diamond builder (same construction as luca_endo_shortcut.build)
# --------------------------------------------------------------------------
def build(p, a, q, aux):
    p, a, q, aux = ZZ(p), ZZ(a), ZZ(q), ZZ(aux)
    N = ZZ(2) ** a
    assert q + aux == N and gcd(q, aux) == 1 and aux % 2 == 1 and q % 2 == 1
    F = GF(p ** 2, name="i", modulus=[1, 0, 1])
    E_chl = EllipticCurve(F, [1, 0])
    assert E_chl.is_supersingular()
    M = ZZ(2) ** (a + 2)

    P, Q = torsion_basis(E_chl, M)
    P, Q = fix_torsion_basis_renes(P, Q, a + 2)

    Pq, Qq = torsion_basis(E_chl, q)
    a_iso, E_com = isogeny_from_scalar_x_only(E_chl, q, ZZ(1), basis=(Pq, Qq))
    aP, aQ = evaluate_isogeny_x_only(a_iso, P, Q, M, q)

    Pa, Qa = torsion_basis(E_chl, aux)
    tau_iso, E_aux = isogeny_from_scalar_x_only(E_chl, aux, ZZ(1), basis=(Pa, Qa))
    tP, tQ = evaluate_isogeny_x_only(tau_iso, P, Q, M, aux)

    N_constant = F(aux + q) / F(aux - q)
    strategy = optimised_strategy(a - 1)

    return dict(p=p, a=a, q=q, aux=aux, N=N, M=M, F=F,
                E_chl=E_chl, E_com=E_com, E_aux=E_aux,
                P=P, Q=Q, aP=aP, aQ=aQ, tP=tP, tQ=tQ,
                a_iso=a_iso, tau_iso=tau_iso,
                N_constant=N_constant, strategy=strategy)


def run_chain(D, comp_side, aux_side):
    """(2^a,2^a) chain on graph{comp_side, aux_side} (each a pair of 2^(a+2)-torsion
    points). Returns (split?, recovers_source?). Invalid/non-isotropic kernels raise
    inside compute_richelot_chain and are reported as (False, False) -- exactly what
    the honest verifier would do."""
    cP, cQ = comp_side
    xP, xQ = aux_side
    ker = (4 * cP, 4 * cQ, 4 * xP, 4 * xQ)
    try:
        _, (F1, F2) = compute_richelot_chain(ker, D["a"], D["N_constant"], D["strategy"])
    except Exception:
        return (False, False)
    if F1 is None or F2 is None:
        return (False, False)
    rec = bool(F1.is_isomorphic(D["E_chl"]) or F2.is_isomorphic(D["E_chl"]))
    return (True, rec)


def find_prime(a, q, aux, kmax=4000):
    """Smallest prime p with j=1728 supersingular and full rational 2^(a+2), q, aux
    torsion over F_{p^2}: p ≡ 3 mod 4 and (p+1) divisible by 2^(a+2)*q*aux."""
    a, q, aux = ZZ(a), ZZ(q), ZZ(aux)
    cofac = (ZZ(2) ** (a + 2)) * q * aux
    for k in range(1, kmax):
        p = ZZ(k) * cofac - 1
        if p < 5:
            continue
        if p % 4 == 3 and p.is_prime():
            return int(p)
    return None


# --------------------------------------------------------------------------
# small-factor l dividing aux = 2^a - q
# --------------------------------------------------------------------------
def small_factors(aux):
    return [int(l) for l in (3, 5, 7) if aux % l == 0]


def l_isogenies(E, l, M):
    """The l+1 degree-l isogenies out of E, built with the SAME Montgomery x-only
    machinery as the diamond (isogeny_from_scalar_x_only), so codomains and pushed
    points are theta-chain-ingestible. (sage's isogenies_prime_degree returns
    long-Weierstrass codomains the (2^a,2^a) chain cannot consume -- using it here
    makes every kernel spuriously fail to split.) Returns [(iso, codomain)].

    The l+1 order-l kernels are <R + sS> for s in 0..l-1 and <S>."""
    R, S = torsion_basis(E, l)
    outs = []
    for s in range(l):
        iso, Ec = isogeny_from_scalar_x_only(E, ZZ(l), ZZ(s), basis=(R, S))
        outs.append((iso, Ec))
    iso, Ec = isogeny_from_scalar_x_only(E, ZZ(l), ZZ(1), basis=(S, R))  # kernel <S>
    outs.append((iso, Ec))
    return outs


# --------------------------------------------------------------------------
# Attack REROUTE -- the mechanism: reroute the last l-step onto a new curve.
# --------------------------------------------------------------------------
def attack_reroute(D, l):
    """For every pair of l-isogenies (rho: E_aux->E_mid, psi: E_mid->E_aux''),
    form auxiliary images  psi( l^{-1} * rho(tP) ),  psi( l^{-1} * rho(tQ) )  with
    l^{-1} taken mod 2^(a+2), keep a's images fixed, run the chain.

    Geometry: if rho is the DUAL of tau's honest last l-step, rho(tP)=l*tau_rest(P),
    so l^{-1}*rho(tP)=tau_rest(P) (the intermediate image on E_mid=E_1'), and psi
    re-extends by a fresh degree-l step -> tau' of degree EXACTLY 2^a-q on a new
    E_aux''. The adversary does not know which rho is the dual, so we sweep all
    (l+1)^2 pairs and let the chain (= the verifier) decide which split.

    A FORGERY = split AND recovers source AND j(E_aux'') != j(E_aux)."""
    E_aux = D["E_aux"]
    M = int(D["M"])
    L = ZZ(l)
    linv = pow(int(l), -1, M)
    tP, tQ = D["tP"], D["tQ"]
    comp = (D["aP"], D["aQ"])
    j_honest = E_aux.j_invariant()

    rhos = l_isogenies(E_aux, l, M)
    total = 0
    splits = 0
    forgeries = []          # split+recover on a NEW curve
    same_curve_hits = 0     # split+recover back on E_aux (honest orbit)
    new_curve_js = set()

    for rho, E_mid in rhos:
        rP, rQ = evaluate_isogeny_x_only(rho, tP, tQ, M, L)
        iP, iQ = linv * rP, linv * rQ          # intermediate images on E_mid
        for psi, E_new in l_isogenies(E_mid, l, M):
            xP, xQ = evaluate_isogeny_x_only(psi, iP, iQ, M, L)
            total += 1
            split, rec = run_chain(D, comp, (xP, xQ))
            if split:
                splits += 1
            if split and rec:
                jn = E_new.j_invariant()
                if jn == j_honest:
                    same_curve_hits += 1
                else:
                    new_curve_js.add(str(jn))
                    forgeries.append(str(jn))

    print(f"  [REROUTE l={l}] swept {total} (rho,psi) l-isogeny pairs "
          f"(l+1={l+1} each side); chain split on {splits}")
    print(f"       split+recover back on the honest E_aux (j={j_honest}): "
          f"{same_curve_hits}")
    print(f"       split+recover on a DIFFERENT auxiliary curve (FORGERY): "
          f"{len(forgeries)}  distinct new j's: {len(new_curve_js)}")
    if new_curve_js:
        shown = list(new_curve_js)[:6]
        print(f"       new auxiliary j-invariants: {shown}")
    return dict(l=l, pairs=total, chain_splits=splits,
                same_curve_recover=same_curve_hits,
                forgeries=len(forgeries), distinct_new_curves=len(new_curve_js),
                new_js=list(new_curve_js)[:16], honest_j=str(j_honest))


# --------------------------------------------------------------------------
# Attack B (prompt, literal) -- uniform l^{-1} on E_aux, then a fresh l-isogeny.
# --------------------------------------------------------------------------
def attack_B(D, l):
    """images = psi( l^{-1} * tau(P) ) with l^{-1} on E_aux directly (NOT via a
    backtracking rho). Composite psi o [l^{-1}] o tau: a uniform diagonal rescale
    diag(l^{-1},l^{-1}) (det = l^{-2} mod 2^a) then a degree-l isogeny. Predicted to
    fail the Weil-pairing isotropy (E1 already killed diag(m,m) with det != 1)."""
    E_aux = D["E_aux"]
    M = int(D["M"])
    L = ZZ(l)
    linv = pow(int(l), -1, M)
    tP, tQ = D["tP"], D["tQ"]
    comp = (D["aP"], D["aQ"])
    dP, dQ = linv * tP, linv * tQ
    j_honest = E_aux.j_invariant()

    splits = 0
    forgeries = 0
    for psi, E_new in l_isogenies(E_aux, l, M):
        xP, xQ = evaluate_isogeny_x_only(psi, dP, dQ, M, L)
        split, rec = run_chain(D, comp, (xP, xQ))
        if split:
            splits += 1
        if split and rec and E_new.j_invariant() != j_honest:
            forgeries += 1
    print(f"  [B l={l}] uniform l^-1 then {l+1} l-isogenies: chain split "
          f"{splits}/{l+1}, forgeries {forgeries}")
    return dict(l=l, tried=l + 1, chain_splits=splits, forgeries=forgeries)


# --------------------------------------------------------------------------
# Attack C (prompt, literal) -- psi' = psi_l o psi, degree l*(2^a - q).
# --------------------------------------------------------------------------
def attack_C(D, l):
    """images = psi_l(tau(P)), i.e. post-compose tau with a degree-l isogeny.
    Effective degree l*(2^a-q) != 2^a-q -> Kani degree constraint violated;
    predicted to fail."""
    E_aux = D["E_aux"]
    M = int(D["M"])
    L = ZZ(l)
    tP, tQ = D["tP"], D["tQ"]
    comp = (D["aP"], D["aQ"])
    j_honest = E_aux.j_invariant()
    splits = 0
    forgeries = 0
    for psi_l, E_new in l_isogenies(E_aux, l, M):
        xP, xQ = evaluate_isogeny_x_only(psi_l, tP, tQ, M, L)
        split, rec = run_chain(D, comp, (xP, xQ))
        if split:
            splits += 1
        if split and rec and E_new.j_invariant() != j_honest:
            forgeries += 1
    print(f"  [C l={l}] psi_l o tau (degree l*(2^a-q)={l * int(D['aux'])}): chain split "
          f"{splits}/{l + 1}, forgeries {forgeries}")
    return dict(l=l, tried=l + 1, comp_degree=l * int(D["aux"]),
                chain_splits=splits, forgeries=forgeries)


# --------------------------------------------------------------------------
# Attack A (prompt) -- l-isogenies out of E_chl; tail-completion barrier.
# --------------------------------------------------------------------------
def attack_A(D, l):
    """Enumerate the l+1 l-isogenies phi_l': E_chl -> E_1'. One is tau's honest
    first step. Replacing it requires a NEW tail E_1'' -> E_aux'' of degree
    (2^a-q)/l, which needs End(E_chl) (isogeny-path problem) -- the adversary has
    only torsion images. We report the enumeration and the barrier; no chain is run
    because there is no leaked data to build the replacement tail from."""
    E_chl = D["E_chl"]
    isos = E_chl.isogenies_prime_degree(l)
    js = [str(phi.codomain().j_invariant()) for phi in isos]
    print(f"  [A l={l}] {len(isos)} l-isogenies out of E_chl; replacing tau's first "
          f"step needs a fresh degree-{int(D['aux']) // l} tail (End(E_chl) / "
          f"isogeny-path) -- not available from leaked torsion.")
    return dict(l=l, num_first_steps=len(isos), first_step_js=js[:8],
                tail_degree=int(D["aux"]) // l,
                barrier="tail completion needs End(E_chl); no leaked data suffices")


# --------------------------------------------------------------------------
def main():
    # Instances with l | (2^a - q). p=51839 reuses the known-good non-degenerate
    # diamond (aux=27=3^3, l=3). The rest are found by find_prime for l in {3,5,7}.
    specs = [
        dict(a=5, q=5, aux=27, p=51839),          # l=3, known-good
        dict(a=8, q=13, aux=243),                 # l=3, aux=3^5 ("fully smooth")
        dict(a=8, q=11, aux=245),                 # l=5 (and 7): 245 = 5*7^2
        dict(a=8, q=25, aux=231),                 # l=7: 231 = 3*7*11
    ]
    report = []
    for spec in specs:
        a, q, aux = spec["a"], spec["q"], spec["aux"]
        assert q + aux == 2 ** a and gcd(q, aux) == 1
        p = spec.get("p") or find_prime(a, q, aux)
        if p is None:
            print(f"\n[skip] no prime found for a={a} q={q} aux={aux}")
            continue
        D = build(p, a, q, aux)
        ls = small_factors(aux)
        print(f"\n===== p={p}  a={a}  q={q}  aux(=2^a-q)={aux}  "
              f"small factors l={ls}  =====")
        print(f"  j(E_chl)={D['E_chl'].j_invariant()}  "
              f"j(E_com)={D['E_com'].j_invariant()}  "
              f"j(E_aux)={D['E_aux'].j_invariant()}")
        s, r = run_chain(D, (D["aP"], D["aQ"]), (D["tP"], D["tQ"]))
        print(f"  [sanity] honest torsion-only kernel: split={s} recovers_source={r}")
        assert s and r, "honest diamond must split+recover; instance is degenerate"

        inst = dict(p=p, a=a, q=q, aux=aux, small_factors=ls,
                    honest_split=bool(s), honest_recovers=bool(r),
                    reroute=[], attackB=[], attackC=[], attackA=[])
        forged_via = []
        for l in ls:
            rr = attack_reroute(D, l)
            inst["reroute"].append(rr)
            if rr["forgeries"] > 0:
                forged_via.append(l)
            inst["attackB"].append(attack_B(D, l))
            inst["attackC"].append(attack_C(D, l))
            inst["attackA"].append(attack_A(D, l))
        inst["forged_via_l"] = forged_via
        print(f"  => FORGERY on a NEW auxiliary curve via small factor(s) "
              f"l={forged_via or 'none (this instance)'}")
        report.append(inst)

    n_inst = len(report)
    n_forged = sum(1 for inst in report if inst["forged_via_l"])
    print(f"\n=== SUMMARY: {n_forged}/{n_inst} instances forged on a new auxiliary "
          f"curve using only leaked torsion + O(l)-isogeny enumeration ===")
    print("REROUTE cost uses only degree-l isogenies and l^-1 mod 2^(a+2); the "
          "cofactor (2^a-q)/l is NEVER evaluated, so smoothness of the cofactor is "
          "irrelevant -- a single small factor l | (2^a-q) is what the attack needs.")

    print("\n=== JSON ===")
    print(json.dumps(report, indent=2))


main()
