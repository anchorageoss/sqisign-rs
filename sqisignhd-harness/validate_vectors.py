#!/usr/bin/env python3
r"""
validate_vectors.py  --  internal-consistency checker for the SQIsignHD
dim-4 verification test vectors produced by extract_vectors.py.

Pure Python 3 (no SageMath required).  It re-derives, from the JSON alone:

  * structural invariants  : 16 coords per theta null point; every F_{p^2}
                             coordinate in [0,p); chain length == e - m;
  * arithmetic invariants  : a1^2 + a2^2 == N, q + a1^2 + a2^2 == 2^f,
                             N % 4 == 1, N (probably) prime;
  * the verifier's own checks, recomputed from the dumped coordinates:
        - codomain match : C1.zero() == HC2.zero()      (projective)
        - Hadamard link  : HC2.zero() == H(C2.null_point())   (projective)
        - C1.zero()      == last F1 chain codomain        (identity)
  * result == "accept" for every vector.

Exit code 0 iff all checks pass.

Usage:  python3 validate_vectors.py [test_vectors_l1.json]
"""
import json
import sys
import os


def load(path):
    with open(path) as fh:
        return json.load(fh)


# ----- F_{p^2} = F_p[i]/(i^2+1) arithmetic on (re, im) integer tuples --------
def parse_coord(pair, p):
    a = int(pair[0], 16) if isinstance(pair[0], str) else int(pair[0])
    b = int(pair[1], 16) if isinstance(pair[1], str) else int(pair[1])
    return (a % p, b % p)


def fp2_mul(x, y, p):
    (a, b), (c, d) = x, y
    return ((a * c - b * d) % p, (a * d + b * c) % p)


def fp2_sub(x, y, p):
    return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)


def fp2_add(x, y, p):
    return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)


def is_zero(x):
    return x[0] == 0 and x[1] == 0


def hadamard16(P, p):
    """Unnormalized tensor Hadamard: H[chi] = sum_j (-1)^popcount(chi&j) P[j]."""
    out = []
    for chi in range(16):
        acc = (0, 0)
        for j in range(16):
            term = P[j]
            if bin(chi & j).count("1") & 1:
                acc = fp2_sub(acc, term, p)
            else:
                acc = fp2_add(acc, term, p)
        out.append(acc)
    return out


def proj_equal(P, Q, p):
    """Projective equality of two 16-coordinate theta points (cross-multiply),
    matching ThetaPointDim4.__eq__ in the reference."""
    if len(P) != len(Q):
        return False
    k0 = 0
    while k0 < len(P) - 1 and is_zero(P[k0]):
        k0 += 1
    for l in range(len(P)):
        lhs = fp2_mul(P[l], Q[k0], p)
        rhs = fp2_mul(Q[l], P[k0], p)
        if lhs != rhs:
            return False
    return True


def coords_of(rec, p):
    cc = rec["coords"]
    assert len(cc) == 16, "theta point must have 16 coords, got %d" % len(cc)
    out = []
    for pair in cc:
        a, b = parse_coord(pair, p)
        assert 0 <= a < p and 0 <= b < p, "coordinate out of range [0,p)"
        out.append((a, b))
    return out


def miller_rabin(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
    if n < 2:
        return False
    for q in bases:
        if n % q == 0:
            return n == q
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


class Checker:
    def __init__(self):
        self.fail = 0
        self.ok = 0

    def check(self, cond, msg):
        if cond:
            self.ok += 1
        else:
            self.fail += 1
            print("  FAIL: " + msg)
        return cond


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_vectors_l1.json")
    doc = load(path)
    c = Checker()

    p = int(doc["prime_decimal"])
    lvl = doc["level"]
    params = doc["params"]
    f = params["f"]
    e1, e2 = params["e1"], params["e2"]
    print("File: %s" % path)
    print("level=%s  p=%s  f=%d  e1=%d e2=%d  sage=%s" % (
        lvl, doc.get("prime"), f, e1, e2, doc.get("sage_version")))
    print("reference_commit=%s" % doc.get("reference_commit"))
    print("theta_dim4_commit=%s" % doc.get("theta_dim4_commit"))

    # Top-level
    c.check(p == 5 * 2 ** 248 - 1, "prime is 5*2^248-1 (level 1)")
    c.check(e1 == (f + 1) // 2 and e2 == f - (f + 1) // 2, "e1,e2 split of f")
    vectors = doc["test_vectors"]
    c.check(len(vectors) == doc["n_vectors"], "n_vectors matches array length")
    c.check(len(vectors) >= 5, "at least 5 vectors present")

    for v in vectors:
        idx = v["index"]
        print("\n[vector %d]" % idx)

        # result
        c.check(v["result"] == "accept", "v%d result == accept" % idx)

        s4 = v["stage4_compute_hd"]
        q = int(s4["q"]); N = int(s4["N"])
        a1 = int(s4["a1"]); a2 = int(s4["a2"])
        m = s4["m"]; lF1 = s4["len_F1"]; lF2 = s4["len_F2_dual"]

        # arithmetic invariants
        c.check(a1 * a1 + a2 * a2 == N, "v%d a1^2+a2^2 == N" % idx)
        c.check(q + a1 * a1 + a2 * a2 == 2 ** f, "v%d q+a1^2+a2^2 == 2^f" % idx)
        c.check(N % 4 == 1, "v%d N %% 4 == 1" % idx)
        c.check(s4["N_mod_4"] == 1 and s4["a1_sq_plus_a2_sq_eq_N"], "v%d recorded N flags" % idx)
        c.check(miller_rabin(N), "v%d N is (probably) prime" % idx)

        # chain lengths == e - m
        c.check(lF1 == e1 - m, "v%d len_F1 == e1-m (%d == %d-%d)" % (idx, lF1, e1, m))
        c.check(lF2 == e2 - m, "v%d len_F2_dual == e2-m (%d == %d-%d)" % (idx, lF2, e2, m))

        F1 = s4["F1"]; F2 = s4["F2_dual"]
        c.check(len(F1["chain_codomains"]) == lF1, "v%d F1 chain array length" % idx)
        c.check(len(F2["chain_codomains"]) == lF2, "v%d F2_dual chain array length" % idx)

        # 16 coords + range for every theta point on both chains + starts
        all_ok = True
        for name, node in (("F1_start", F1["theta_null_start"]),
                           ("F1_prod", F1["theta_null_domain_product"]),
                           ("F2_start", F2["theta_null_start"]),
                           ("F2_prod", F2["theta_null_domain_product"])):
            try:
                coords_of(node, p)
            except AssertionError as e:
                all_ok = False
                print("    %s: %s" % (name, e))
        for chain_name, chain in (("F1", F1), ("F2_dual", F2)):
            for step, node in enumerate(chain["chain_codomains"]):
                try:
                    coords_of(node, p)
                except AssertionError as e:
                    all_ok = False
                    print("    %s step %d: %s" % (chain_name, step, e))
        c.check(all_ok, "v%d every theta null point has 16 in-range coords" % idx)

        # stage 5: recompute the verifier's codomain checks from coordinates
        s5 = v["stage5_codomain_check"]
        C1z = coords_of(s5["C1_zero"], p)
        HC2z = coords_of(s5["HC2_zero"], p)
        c.check(s5["match"] is True, "v%d recorded codomain match flag" % idx)
        c.check(proj_equal(C1z, HC2z, p),
                "v%d C1.zero() == HC2.zero() (recomputed, projective)" % idx)

        # C1.zero() must equal the last F1 chain codomain (same structure)
        lastF1 = coords_of(F1["chain_codomains"][-1], p)
        c.check(proj_equal(C1z, lastF1, p),
                "v%d C1.zero() == last F1 chain codomain" % idx)

        # Hadamard link: HC2.zero() == H(C2.null_point()) where C2.null_point()
        # is the last F2_dual chain codomain (projective; off by factor 16).
        lastF2 = coords_of(F2["chain_codomains"][-1], p)
        c.check(proj_equal(HC2z, hadamard16(lastF2, p), p),
                "v%d HC2.zero() == Hadamard(last F2_dual codomain)" % idx)

        # normalized middle codomains must be byte-identical (same projective
        # point normalized at the same pivot)
        n1 = s5["C1_zero"].get("normalized"); n2 = s5["HC2_zero"].get("normalized")
        p1 = s5["C1_zero"].get("pivot"); p2 = s5["HC2_zero"].get("pivot")
        if n1 is not None and n2 is not None:
            c.check(p1 == p2 and n1 == n2,
                    "v%d normalized C1.zero()==HC2.zero() and same pivot" % idx)

        # stage 6
        s6 = v["stage6_image_check"]
        c.check(s6["correct"] is True, "v%d recorded image check flag" % idx)
        c.check(s6["FT"][3].get("inf") is True, "v%d FT[3] is identity on E_chal" % idx)

    print("\n=== %d checks passed, %d failed ===" % (c.ok, c.fail))
    sys.exit(1 if c.fail else 0)


if __name__ == "__main__":
    main()
