#!/usr/bin/env python3
"""The constant of "at least about 0.07 X^2 / p vectors of norm at most X
outside L_2" (Section "Parameters for the round-3 primes"), and the counts
2^18.6, 2^19.4, 2^20.1 and 2^9.7 that follow from it.

Setting. q is a positive definite quaternary form of determinant p^2/16
(covolume p/4), N1 <= N2 <= N3 <= N4 its successive minima (norms), L_2 the
sublattice of the first two, L' the projection of the lattice orthogonal to
L_2, X = 2^e, and N3 <= X/2 is assumed.

Step 1, vectors per translate. A translate v + L_2 whose projection has
norm at most X/2 meets the ball of norm X in a disc of squared radius at
least X/2. The covering radius of L_2 is at most sqrt((N1 + N2)/4) <=
sqrt(X/4), so the disc contains a disc of radius sqrt(X/2) - sqrt(X/4)
around a lattice point of the translate, hence at least
    c1 * X / det(L_2),   c1 = (pi/2) * (1 - 1/sqrt 2)^2
points of it (area over covolume), and det(L_2) <= sqrt(N1 N2).

Step 2a, N4 <= X/2. The minima of L' are at most N3, N4 <= X/2, its
covering radius at most sqrt(X/4), and the same count gives at least
c1 * X / det(L') classes of projection norm at most X/2 (the zero class,
which is L_2 itself, is one of them and is negligible against the rest).
With det(L_2) det(L') = p/4:
    at least c1^2 * X^2 / (p/4) = 4 c1^2 * X^2 / p   vectors.

Step 2b, N4 > X/2. The classes k * b3 with 0 < k^2 N3 <= X/2 number at
least sqrt(X / (2 N3)), each with c1 * X / sqrt(N1 N2) vectors. Minkowski's
second theorem, N1 N2 N3 N4 <= 64 p^2 / pi^4, with N4 > X/2 gives
N1 N2 N3 < 128 p^2 / (pi^4 X), hence
    at least c1 * (pi^2 / 16) * X^2 / p   vectors.

The constant is the smaller of the two.
"""
import math

c1 = (math.pi / 2) * (1 - 1 / math.sqrt(2)) ** 2
case_a = 4 * c1 * c1
case_b = c1 * math.pi**2 / 16
c = min(case_a, case_b)
print(f"c1 = (pi/2)(1 - 1/sqrt 2)^2 = {c1:.4f}   (the paper's 0.13)")
print(f"N4 <= X/2: 4 c1^2      = {case_a:.4f}")
print(f"N4 >  X/2: c1 pi^2/16  = {case_b:.4f}")
print(f"constant: {c:.4f}   (the paper's 'about 0.07')")
levels = {"I": (math.log2(3) + 324, 174), "III": (math.log2(27) + 500, 264), "V": (math.log2(17) + 664, 346)}
for name, (log2p, e) in levels.items():
    n = math.log2(0.07) + 2 * e - log2p
    print(f"level {name}: 0.07 X^2/p = 2^{n:.2f}   (with {c:.4f}: 2^{math.log2(c) + 2*e - log2p:.2f})")
log2p, e = levels["I"]
vec = math.log2(0.07) + 2 * e - log2p
good = vec - 1 - 1 - math.log2(e * math.log(2))
print(f"level I, good norms: half of the norms odd (-1), v and -v counted once (-1), one odd norm in e ln 2 = {e*math.log(2):.1f} good: 2^{good:.2f}")
