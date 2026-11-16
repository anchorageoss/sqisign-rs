#!/usr/bin/env python3
"""Successive minima of the response lattices recorded by paper_instrument.

    python3 gram_minima.py <dir with gram.csv> [--first N]

gram.csv has one Gram matrix G per commitment, with x^T G x equal to twice
the degree of the response x. The norm form is q = G / 2, of determinant
p^2 / 16 (checked on every row). For each lattice the script reduces the
basis (LLL on the Gram matrix, exact rational arithmetic), enumerates the
coefficient vectors in a box around the reduced basis and reads off the
four successive minima N1 <= N2 <= N3 <= N4 of q. It prints what the
paper's Section "Parameters" quotes: the third minimum normalised by
sqrt(p)/2 (median, maximum, the maximum as a power of X = 2^e), the tail
constant y^4 Pr[N3 > y sqrt(p)/2] over the range where the sample has
between 38 percent of its points and 6 points, the first minimum's
Pr[nu1 <= x] / x^2 over 2^-5 <= x <= 2^-1, and the probabilities that
follow from the tail constant. Pure Python, no dependency.
"""
import csv, itertools, math, sys
from fractions import Fraction

P = 3 * 2**324 - 1
E = {"I": 174, "III": 264, "V": 346}
LOG2P = {"I": math.log2(3) + 324, "III": math.log2(27) + 500, "V": math.log2(17) + 664}
BOX = 3  # coefficients of the reduced basis searched: -BOX..BOX


def lll_gram(g):
    """LLL (delta 3/4) on a Gram matrix; returns the reduced Gram matrix."""
    n = len(g)
    b = [[int(i == j) for j in range(n)] for i in range(n)]  # rows: new basis in the old

    def gram(u, v):
        return sum(u[i] * g[i][j] * v[j] for i in range(n) for j in range(n))

    def gso():
        mu = [[Fraction(0)] * n for _ in range(n)]
        bb = [Fraction(0)] * n
        gg = [[gram(b[i], b[j]) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                s = Fraction(gg[i][j]) - sum(mu[i][k] * mu[j][k] * bb[k] for k in range(j))
                mu[i][j] = s / bb[j]
            bb[i] = Fraction(gg[i][i]) - sum(mu[i][k] ** 2 * bb[k] for k in range(i))
        return mu, bb

    k = 1
    mu, bb = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            r = round(mu[k][j])
            if r:
                b[k] = [x - r * y for x, y in zip(b[k], b[j])]
                mu, bb = gso()
        if bb[k] >= (Fraction(3, 4) - mu[k][k - 1] ** 2) * bb[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            mu, bb = gso()
            k = max(k - 1, 1)
    return [[gram(b[i], b[j]) for j in range(n)] for i in range(n)]


def rank_increases(vs, v):
    """Whether v is independent of the integer vectors vs (fraction elimination)."""
    rows = [[Fraction(x) for x in u] for u in vs + [v]]
    rank = 0
    for col in range(4):
        piv = next((r for r in range(rank, len(rows)) if rows[r][col] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(rank + 1, len(rows)):
            f = rows[r][col] / rows[rank][col]
            rows[r] = [a - f * c for a, c in zip(rows[r], rows[rank])]
        rank += 1
    return rank == len(vs) + 1


COEFFS = [c for c in itertools.product(range(-BOX, BOX + 1), repeat=4) if c > (0, 0, 0, 0)]


def minima(g):
    """The successive minima of x^T g x (g reduced), as integers."""
    vals = []
    for c in COEFFS:
        v = sum(c[i] * g[i][j] * c[j] for i in range(4) for j in range(4))
        vals.append((v, c))
    vals.sort()
    out, vs = [], []
    for v, c in vals:
        if rank_increases(vs, list(c)):
            vs.append(list(c))
            out.append(v)
            if len(out) == 4:
                break
    # the box is wide enough when the fourth minimum's vectors are interior
    assert all(abs(x) < BOX for u in vs for x in u), "box too small"
    return out


def det4(m):
    m = [[Fraction(x) for x in r] for r in m]
    d = Fraction(1)
    for i in range(4):
        p = next(r for r in range(i, 4) if m[r][i] != 0)
        if p != i:
            m[i], m[p] = m[p], m[i]
            d = -d
        d *= m[i][i]
        for r in range(i + 1, 4):
            f = m[r][i] / m[i][i]
            m[r] = [a - f * c for a, c in zip(m[r], m[i])]
    return d


def main():
    d = sys.argv[1]
    first = int(sys.argv[sys.argv.index("--first") + 1]) if "--first" in sys.argv else None
    rows = [r for r in csv.reader(l for l in open(f"{d}/gram.csv") if not l.startswith("#"))][1:]
    if first is not None:
        rows = [r for r in rows if int(r[0]) < first]
    unit = math.sqrt(P) / 2  # the fourth root of the determinant
    nus = []
    for r in rows:
        g = [[int(x, 16) for x in r[4 + 4 * i : 8 + 4 * i]] for i in range(4)]
        assert det4(g) == P * P, "determinant of G is not p^2 (that of G/2 is p^2/16)"
        m = minima(lll_gram(g))
        nus.append([v / 2 / unit for v in m])  # q = G/2
    n = len(nus)
    print(f"commitments: {n} (signatures {rows[0][0]} to {rows[-1][0]})")
    for i in range(4):
        col = sorted(v[i] for v in nus)
        print(f"nu{i+1}: median 2^{math.log2(col[n//2]):.2f}, min 2^{math.log2(col[0]):.2f}, max 2^{math.log2(col[-1]):.2f}")
    nu3 = sorted(v[2] for v in nus)
    mx = nu3[-1]
    e = E["I"]
    y_x = e - LOG2P["I"] / 2  # log2 of X / sqrt(p)
    print(f"max nu3 = 2^{math.log2(mx):.2f}, that is N3 = 2^{math.log2(mx) - 1 - y_x:.1f} X")
    # tail: y^4 * Pr[nu3 > y] for y at each order statistic from 38 % of the sample down to 6 points
    lo_k, hi_k = 6, int(0.38 * n)
    prods = []
    for k in range(lo_k, hi_k + 1):
        y = nu3[n - k]  # k points are >= y... strictly above the one before
        y = nu3[n - k - 1]  # exactly k points lie above this value
        prods.append((y, k, y**4 * k / n))
    ys = [p[0] for p in prods]
    ps = sorted(p[2] for p in prods)
    c_med = ps[len(ps) // 2]
    # least squares of log Pr = log c - a log y, to show the exponent
    lx = [math.log(p[0]) for p in prods]
    ly = [math.log(p[1] / n) for p in prods]
    mxx, myy = sum(lx) / len(lx), sum(ly) / len(ly)
    slope = sum((a - mxx) * (b - myy) for a, b in zip(lx, ly)) / sum((a - mxx) ** 2 for a in lx)
    print(f"tail of nu3 over {min(ys):.2f} <= y <= {max(ys):.2f} (between {hi_k} and {lo_k} points above y):")
    print(f"  y^4 Pr[nu3 > y]: min {ps[0]:.2f}, median {c_med:.2f}, max {ps[-1]:.2f}; fitted exponent {slope:.2f}")
    for kmin in (20, 50, 100):
        sub = sorted(p[2] for p in prods if p[1] >= kmin)
        if sub:
            print(f"  with at least {kmin} points above y (y <= {max(p[0] for p in prods if p[1] >= kmin):.2f}): min {sub[0]:.2f}, max {sub[-1]:.2f}")
    # the first minimum
    nu1 = sorted(v[0] for v in nus)
    rat = []
    for k in range(-10, -1):
        x = 2.0 ** (k / 2)
        cnt = sum(1 for v in nu1 if v <= x)
        if cnt:
            rat.append((x, cnt, cnt / n / x**2))
    if rat:
        print("Pr[nu1 <= x] / x^2 over 2^-5 <= x <= 2^-1: " + ", ".join(f"{r[2]:.2f} ({r[1]} pts at 2^{math.log2(r[0]):.1f})" for r in rat))
        print(f"  Siegel: pi^2/4 = {math.pi**2/4:.3f}")
    print(f"probabilities from Pr[N3 > X/2] = c (X/sqrt p)^-4 with c = {c_med:.2f} (and with the extremes {ps[0]:.2f}, {ps[-1]:.2f}):")
    for lvl in ("I", "III", "V"):
        yx = E[lvl] - LOG2P[lvl] / 2
        vals = [math.log2(c) - 4 * yx for c in (c_med, ps[0], ps[-1])]
        print(f"  level {lvl}: y = 2^{yx:.2f}; 2^{vals[0]:.1f} (2^{vals[1]:.1f} to 2^{vals[2]:.1f}); random-lattice bound 690 y^-4 = 2^{math.log2(690) - 4*yx:.1f}")
    yx = E["I"] - LOG2P["I"] / 2
    p1 = c_med * 2 ** (-4 * yx)
    print(f"level I: a signature restarts with probability about 2^{math.log2(p1):.1f}; 64 commitments in a row fail with probability 2^{64*math.log2(p1):.0f}")


main()
