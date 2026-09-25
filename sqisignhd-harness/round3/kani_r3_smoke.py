# Real-prime check of the dimension-4 half-chain at the round-3 level I prime, with the
# SQIsignHD library's Theta_dim4 package unmodified (as Verify.py uses it), built like the
# library's own test (Tests.py, test_kani_endomorphism_half): walk off E0 to a generic domain,
# then sigma is a walk of k rational 3-isogenies, q = 3^k (odd k gives q = 3 mod 4), and
# 2^e - q = a1^2 + a2^2 (a prime = 1 mod 4 when one exists among the usable k, else any
# two-squares value found by factoring).
from sage.all import *
proof.all(False)
import sys, os
from time import time
sys.path.insert(0, os.getcwd())
from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_endomorphism import KaniEndoHalf
from Theta_dim4.Theta_dim4_sage.pkg.utilities.strategy import precompute_strategy_with_first_eval

p = 3*2**324 - 1
e = int(sys.argv[2]) if len(sys.argv) > 2 else 174
r = ceil(e/2) + 2
print(f"p = 3*2^324 - 1, e = {e}, r = {r}", flush=True)
Fp2 = GF(p**2, 'i', modulus=[1,0,1], proof=False)
E0 = EllipticCurve(Fp2, [0,0,0,1,0])
Q0 = BinaryQF([1,0,1])
set_random_seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

# the degree: q = 3^k < 2^e, k odd, with 2^e - q a sum of two squares
t0 = time()
kmax = floor(e * log(2) / log(3))
kmin = max(1, ceil((e - 8) * log(2) / log(3)))   # q within 8 bits of 2^e, as a response would be
ks = [k for k in range(kmax, kmin - 1, -1) if k % 2 == 1]
choice = None
for k in ks:
    N = 2**e - 3**k
    if N % 4 == 1 and ZZ(N).is_pseudoprime():
        choice = (k, N, "prime"); break
if choice is None:
    for k in ks:
        N = ZZ(2**e - 3**k)
        if N % 4 != 1: continue
        t_f = time()
        try:
            sol = Q0.solve_integer(N)   # factors N; None if not a sum of two squares
        except Exception:
            sol = None
        print(f"  k = {k}: factoring 2^e - 3^k took {time()-t_f:.1f} s, two squares: {sol is not None}", flush=True)
        if sol is not None:
            choice = (k, N, "composite two-squares"); break
if choice is None:
    print("no usable k", flush=True); sys.exit(2)
k, N, kind = choice
q = 3**k
a1, a2 = Q0.solve_integer(ZZ(N))
assert a1**2 + a2**2 + q == 2**e
print(f"q = 3^{k} ({ZZ(q).nbits()} bits, q mod 4 = {q % 4}); 2^e - q {kind}, decomposed in {time()-t0:.2f} s", flush=True)

def rational_3_point(E, avoid=None):
    cof3 = (p + 1) // 3
    while True:
        T = cof3 * E.random_point()
        if T != E(0) and (avoid is None or avoid(T) != avoid.codomain()(0)):
            return T
def montgomery(E):
    M, iso = E.montgomery_model(morphism=True)
    ai = M.a_invariants()
    assert ai[0] == 0 and ai[2] == 0 and ai[3] == 1 and ai[4] == 0
    return M, iso

# a generic domain: 2-power walk of length 8 from E0 (its degree does not enter sigma)
t0 = time()
E = E0
for _ in range(8):
    cof2 = (p + 1) // 2
    while True:
        T = cof2 * E.random_point()
        if T != E(0): break
    E = E.isogeny(T).codomain()
Edom, iso_dom = montgomery(E)
# sigma: k rational 3-isogenies without backtracking, evaluated step by step
steps = []
E = Edom
prev = None
for _ in range(k):
    T = rational_3_point(E, avoid=prev.dual() if prev is not None else None)
    phi = E.isogeny(T)
    steps.append(phi)
    prev = phi
    E = phi.codomain()
E1, iso_1 = montgomery(E)
def sigma(P):
    for phi in steps:
        P = phi(P)
    return iso_1(P)
print(f"domain walk, {k} 3-isogeny steps and Montgomery models: {time()-t0:.2f} s", flush=True)

cof = (p+1) // 2**r
def random_point_order_2r():
    while True:
        P = cof * Edom.random_point()
        if (2**(r-1))*P != Edom(0):
            return P
t0 = time()
P = random_point_order_2r()
while True:
    Qp = random_point_order_2r()
    if (2**(r-1))*P != (2**(r-1))*Qp:
        break
R2, S2 = sigma(P), sigma(Qp)
assert (2**r)*R2 == E1(0) and (2**(r-1))*R2 != E1(0)
print(f"basis of Edom[2^{r}] and sigma images on E1: {time()-t0:.2f} s", flush=True)

m = 0
ai = a2 if a2 % 2 == 0 else a1
while ai % 2 == 0:
    m += 1; ai //= 2
e1 = ceil(e/2); e2 = e - e1
s1 = precompute_strategy_with_first_eval(e1, m, M=1, S=0.8, I=100)
s2 = s1 if e2 == e1 else precompute_strategy_with_first_eval(e2, m, M=1, S=0.8, I=100)
print(f"strategies ready, m = {m}", flush=True)
t0 = time()
F = KaniEndoHalf(P, Qp, R2, S2, q, a1, a2, e, r, s1, s2)
t1 = time()
C1 = F.F1._isogenies[-1]._codomain
C2 = F.F2_dual._isogenies[-1]._codomain
match = C1.zero() == C2.hadamard().zero()
t2 = time()
T4 = TuplePoint(P, Edom(0), E1(0), E1(0))
FT = F(T4)
img = ((FT[0] == a1*P) or (FT[0] == -a1*P)) and ((FT[1] == a2*P) or (FT[1] == -a2*P)) and (FT[3] == E1(0))
t3 = time()
print(f"half-chains ({e1} + {e2} steps, m = {m}): {t1-t0:.2f} s; middle codomain match: {match} ({t2-t1:.3f} s); image check: {img} ({t3-t2:.2f} s)", flush=True)
