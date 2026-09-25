"""
Experiment C: wrong cross-isogeny rejection test.

Claim under test (ePrint 2026/1305): the correct cross isogeny c: E_aux -> E_com
is NECESSARY -- a wrong c'' should give a wrong psi (eq. 5), hence a wrong kernel,
hence verification failure (codomain does not split with E_com).

Analytic prediction (from eq. 5, see EQ5_VERIFICATION.md): with a = alpha.dual
and hat(c) = the degree-(N-q) isogeny E_com -> E_aux, eq. (5) gives
    psi(alpha(R)) = hat(c)(R)             for R in E_com[N].
So the kernel is graph{(alpha(R), hat(c)(R))}. A wrong cross isogeny c'' just
replaces hat(c) by another degree-(N-q) isogeny hat(c'') : E_com -> E_aux, and
(alpha, hat(c'')) is STILL a common-source diamond -> by Kani it still splits
with E_com. Prediction: wrong c'' also passes the split check. This experiment
checks that prediction on the real theta (2,2)-chain, and counts how many
distinct cross isogenies E_com -> E_aux exist per diamond (the uniqueness
question: if 1, c is unique and necessity is trivial).

Setup (non-degenerate, E_A != E_com): p, a, component deg q, cross deg N-q.
"""

import sys
sys.path.insert(0, ".")
from sage.all import *
from richelot_isogenies.richelot_isogenies import compute_richelot_chain
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only, evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.strategy import optimised_strategy


def enumerate_cross_isogenies(E_com, cross_deg, j_target, M, T1, T2):
    """All degree-cross_deg isogenies E_com -> (curve with j = j_target),
    returned as (codomain_j, images (cT1,cT2), kernel_id). Enumerate the two
    scalar families <P+sQ> and <Q+sP> and dedup by (codomain_j, image data)."""
    P, Q = torsion_basis(E_com, cross_deg)
    seen = {}
    out = []
    families = [((P, Q), range(cross_deg)), ((Q, P), range(cross_deg))]
    for (basis, srange) in families:
        for s in srange:
            try:
                phi, Ec = isogeny_from_scalar_x_only(E_com, cross_deg, ZZ(s), basis=basis)
            except Exception:
                continue
            if Ec.j_invariant() != j_target:
                continue
            cT1, cT2 = evaluate_isogeny_x_only(phi, T1, T2, M, cross_deg)
            key = (str(cT1[0]), str(cT1[1]), str(cT2[0]), str(cT2[1]))
            if key in seen:
                continue
            seen[key] = True
            out.append((cT1, cT2))
    return out


def run(p, a, q, cross, label):
    p, a, q, cross = ZZ(p), ZZ(a), ZZ(q), ZZ(cross)
    N = 2**a
    assert q + cross == N and gcd(q, cross) == 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E_com = EllipticCurve(F, [1, 0])
    M = 2**(a + 2)
    T1, T2 = torsion_basis(E_com, M)
    T1, T2 = fix_torsion_basis_renes(T1, T2, a + 2)

    # fixed component alpha: E_com -> E_A (deg q); a = alpha.dual (genuine, E_A != E_com)
    Pq, Qq = torsion_basis(E_com, q)
    alpha, E_A = isogeny_from_scalar_x_only(E_com, q, ZZ(1), basis=(Pq, Qq))
    aT1, aT2 = evaluate_isogeny_x_only(alpha, T1, T2, M, q)
    jA, jcom = E_A.j_invariant(), E_com.j_invariant()

    # pick the "correct" cross isogeny hat(c): E_com -> E_aux (deg cross), fixing E_aux
    Pc, Qc = torsion_basis(E_com, cross)
    beta, E_aux = isogeny_from_scalar_x_only(E_com, cross, ZZ(1), basis=(Pc, Qc))
    j_aux = E_aux.j_invariant()

    # enumerate ALL degree-cross isogenies E_com -> (j = j_aux): these are the
    # hat(c'') candidates (correct + wrong), all landing on the same auxiliary curve.
    candidates = enumerate_cross_isogenies(E_com, cross, j_aux, M, T1, T2)
    k = len(candidates)

    N_constant = F(cross + q) / F(cross - q)
    strategy = optimised_strategy(a - 1)

    n_split_ecom = 0
    n_split_notecom = 0
    n_nosplit = 0
    for (cT1, cT2) in candidates:
        ker = (4 * aT1, 4 * aT2, 4 * cT1, 4 * cT2)   # graph{(alpha(R), hat(c'')(R))}
        try:
            Phi, (F1, F2) = compute_richelot_chain(ker, a, N_constant, strategy)
            if F1 is None or F2 is None:
                n_nosplit += 1
            elif F1.is_isomorphic(E_com) or F2.is_isomorphic(E_com):
                n_split_ecom += 1
            else:
                n_split_notecom += 1
        except Exception:
            n_nosplit += 1

    print(f"\n===== {label}: p={p}, a={a}, q={q}, cross deg N-q={cross}; "
          f"j(E_A)={jA} j(E_com)={jcom} [E_A!=E_com:{jA!=jcom}] j(E_aux)={j_aux} =====")
    print(f"  # distinct cross isogenies E_com->E_aux (deg {cross}) found: {k}")
    print(f"  of these: split w/ E_com = {n_split_ecom}, "
          f"split w/o E_com = {n_split_notecom}, no-split = {n_nosplit}")
    if k <= 1:
        print("  => cross isogeny is UNIQUE: no wrong c'' exists; necessity is trivial/tight.")
    else:
        n_wrong = k - 1
        wrong_pass = n_split_ecom - 1  # correct one also splits w/ E_com
        print(f"  => {n_wrong} WRONG c'' exist; of these, split-w/-E_com (i.e. 'pass'): "
              f"{max(wrong_pass,0)}; genuinely rejected: {n_wrong - max(wrong_pass,0)}")
    return dict(label=label, p=int(p), a=int(a), q=int(q), cross=int(cross),
                EA_ne_Ecom=bool(jA != jcom), n_cross=int(k),
                split_ecom=int(n_split_ecom), split_not_ecom=int(n_split_notecom),
                nosplit=int(n_nosplit))


import json
res = []
# non-degenerate (E_A != E_com): p=51839 (a=5, q=5, cross=27) and p=479 (a=3,q=3,cross=5)
res.append(run(479, 3, 3, 5, "non-degenerate small"))
res.append(run(51839, 5, 5, 27, "non-degenerate larger"))
# also the other split of a: component=cross swap gives cross deg = 5 (prime): more chance of uniqueness
res.append(run(51839, 5, 27, 5, "non-degenerate (cross deg 5, prime)"))
print("\n=== JSON ===")
print(json.dumps(res, indent=2))
