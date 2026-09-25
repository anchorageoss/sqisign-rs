"""
Experiment A (Check 1): verify eq. (5) of ePrint 2026/1305 pointwise, broadly.

    psi'(P) = - hat(c') * [(2^a - q)]^{-1} * a(P)          for P in E_A[2^a]

Diamond from a common source E_com:
    alpha: E_com -> E_A   (deg q)     => a = alpha.dual(),  a o alpha = [q]
    beta : E_com -> E_aux (deg N-q)   => c' = beta.dual(),  hat(c') = beta
Actual anti-isometry (graph of ker Phi):  psi'_true(alpha(S)) = beta(S).

This check is model-independent (pure isogeny/point arithmetic in Sage), so it
runs across many primes, powers a, degrees q, commitment curves, and isogeny
choices. A single disagreement would mean the derivation of eq. (5) is wrong.
"""

import sys
from sage.all import *


def isogeny_of_degree(E, d, pick=0):
    phi, cur = None, E
    for (ell, e) in factor(d):
        for _ in range(e):
            opts = cur.isogenies_prime_degree(ell)
            step = opts[pick % len(opts)]
            phi = step if phi is None else step * phi
            cur = step.codomain()
    return phi


def torsion_basis_2a(E, a, p):
    """Basis (P, Q) of E[2^a] (rational since 2^a | p+1)."""
    N = 2**a
    order = (p + 1)**2
    assert E.order() == order and order % (N * N) == 0
    cof = (p + 1) // N

    def rand_2a():
        while True:
            R = cof * E.random_point()
            if R.order() == N:
                return R
    P = rand_2a()
    for _ in range(200):
        Q = rand_2a()
        # independent iff Weil pairing has full order N
        if P.weil_pairing(Q, N).multiplicative_order() == N:
            return P, Q
    raise RuntimeError("no basis")


def check_instance(p, a, q, com_pick, alpha_pick, beta_pick):
    N = 2**a
    Fp2 = GF(p**2, name="i", modulus=[1, 0, 1])
    E_com = EllipticCurve(Fp2, [1, 0])
    for _ in range(com_pick):
        E_com = E_com.isogenies_prime_degree(2)[0].codomain()

    alpha = isogeny_of_degree(E_com, q, alpha_pick)
    beta = isogeny_of_degree(E_com, N - q, beta_pick)
    a_iso = alpha.dual()          # E_A -> E_com
    chat = beta                   # hat(c')
    inv = int(inverse_mod(int(N - q) % N, N))

    R1, R2 = torsion_basis_2a(E_com, a, p)
    # test on a spread of torsion points S; P = alpha(S) ranges over E_A[2^a]
    S_list = [R1, R2, R1 + R2, R1 - R2, 3 * R1 + R2, R1 + 5 * R2]
    all_ok = True
    for S in S_list:
        P = alpha(S)                          # in E_A[2^a]
        psi_true = beta(S)                    # geometric anti-isometry
        psi_formula = -chat(inv * a_iso(P))   # eq. (5)
        if psi_formula != psi_true:
            all_ok = False
            print(f"    DISAGREE at S: p={p} a={a} q={q} "
                  f"com={com_pick} ap={alpha_pick} bp={beta_pick}")
    return all_ok, str(E_com.j_invariant()), str(alpha.codomain().j_invariant()), \
        str(beta.codomain().j_invariant())


def main():
    configs = [
        # (p, a) with 2^a | p+1 ; then sweep q, com, isogeny picks
        (191, 4), (191, 5), (191, 6),
        (1279, 6), (1279, 8),
    ]
    n_ok = n_total = 0
    for (p, a) in configs:
        N = 2**a
        # odd q with 1 < q < N, coprime-to-2 (all odd), a few values
        qs = [x for x in (3,5,7,9,15,27,63,127) if 1 < x < N-1]
        def maxpf(n):
            return max([pf for pf,_ in factor(n)]) if n>1 else 1
        for q in qs:
            if max(maxpf(q), maxpf(N-q)) > 19:
                continue
            for com_pick in [0, 1, 2]:
                for beta_pick in [0, 1]:
                    ok, jc, jA, jaux = check_instance(p, a, q, com_pick, 0, beta_pick)
                    n_total += 1
                    n_ok += int(ok)
                    tag = "OK " if ok else "FAIL"
                    print(f"[{tag}] p={p:>5} a={a} q={q:>3} (N-q={N-q:>3}) "
                          f"com={com_pick} bp={beta_pick}  jE_com={jc[:16]}")
    print(f"\nSUMMARY: eq.(5) pointwise identity holds on {n_ok}/{n_total} instances "
          f"(each tested at 6 torsion points)")
    print("VERDICT:", "VALIDATED" if n_ok == n_total else "FAILED — DISAGREEMENT")


main()
