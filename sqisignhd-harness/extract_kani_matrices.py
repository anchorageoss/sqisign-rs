#!/usr/bin/env sage
r"""
extract_kani_matrices.py  --  additive oracle for Phase 5b.4 (Kani embedding
integer matrices + sum-of-two-squares).

For each Level-1 signature this dumps the integer linear-algebra intermediates
of the Kani embedding (`basis_change/kani_base_change.py`), driven by the
reference's own functions on the recorded `(a1, a2, q)`:

  * the post-swap `(a1, a2)` (a1 odd, a2 even) and `m = v2(a2)` actually used by
    `KaniEndoHalf`, plus `q`, `f = r = 70` (the matrix modulus exponent), `N=2^f`;
  * the DETERMINISTIC matrices (closed-form, no linear solve):
      - `matrix_F`, `matrix_F_dual`              (8x8 mod 2^f)
      - kernel blocks `C,D` of `complete_kernel_matrix_F1` / `_F2_dual`
                                                  (4x4 mod 2^f -- the columns
                                                   `kernel_basis` consumes)
      - `gluing_base_change_matrix_dim2_F1`/`_F2` (4x4 mod 4)
  * the COMPLETION-dependent matrices (route through sage `solve_right` =
    PARI `matsolvemod`):
      - `M1`, `M2` from `starting_two_symplectic_matrices`   (8x8 mod 2^f)
      - `M_gluing_1`, `M_gluing_2` from `gluing_base_change_matrix_dim2_dim4_F1/F2`
                                                              (8x8 mod 4)

The reference parameters at Level 1 are `e=248`, `lamb=128`, `r=70`, `f_sig=136`.
`compute_HD` sets `N_HD = 2^{f_sig} - q = a1^2 + a2^2`. `KaniEndoHalf` is then
called with its parameter `e := f_sig = 136` and `f := r = 70`, so every
symplectic matrix here is over Z/2^{70}.

Pure read of the reference; nothing in it (or any Phase 0/3/4/5 deliverable) is
modified. The recorded `(a1, a2, q, m)` are cross-checked against
`test_vectors_l1.json`'s `stage4_compute_hd`.

Usage:  sage extract_kani_matrices.py [--out kani_matrices_l1.json] [--n N] [--lib PATH]
Env  :  SQISIGNHD_LIB
"""

import os
import sys
import json
import argparse


def find_lib(cli_path):
    cands = []
    if cli_path:
        cands.append(cli_path)
    if os.environ.get("SQISIGNHD_LIB"):
        cands.append(os.environ["SQISIGNHD_LIB"])
    cands += [
        os.path.join(os.getcwd(), "SQISignHD-lib"),
        os.path.expanduser("~/SQISignHD-lib"),
        "/home/user/SQISignHD-lib",
    ]
    for c in cands:
        if c and os.path.isdir(os.path.join(c, "Verification", "Theta_dim4",
                                            "Theta_dim4_sage")):
            return c
    raise SystemExit("Could not locate SQISignHD-lib; set --lib or $SQISIGNHD_LIB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--lib", default=None)
    ap.add_argument("--vectors", default=None,
                    help="path to test_vectors_l1.json (for a1,a2,q)")
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "kani_matrices_l1.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    here = os.path.dirname(os.path.abspath(__file__))
    vec_path = args.vectors or os.path.join(here, "test_vectors_l1.json")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    from sage.all import Integers, matrix
    import Theta_dim4.Theta_dim4_sage.pkg.basis_change.kani_base_change as kbc

    doc = json.load(open(vec_path))

    def mat_strs(M, nrows, ncols):
        return [[str(int(M[i, j])) for j in range(ncols)] for i in range(nrows)]

    out = {"level": 1, "f_matrix": 70, "note": "all symplectic matrices over Z/2^70; gluing dim2/dim4 over Z/4",
           "vectors": []}

    for v in doc["test_vectors"][: args.n]:
        s4 = v["stage4_compute_hd"]
        a1r = int(s4["a1"])
        a2r = int(s4["a2"])
        q = int(s4["q"])
        # post-swap canonicalisation used by KaniEndoHalf: a1 odd, a2 even.
        if a1r % 2 == 1:
            a1, a2 = a1r, a2r
        else:
            a1, a2 = a2r, a1r
        m = 0
        t = a2
        while t % 2 == 0:
            m += 1
            t //= 2
        f = 70  # = r ; the matrix modulus exponent inside KaniEndoHalf
        N = 2 ** f

        assert m == int(s4["m"]), "m mismatch vs oracle"

        # Deterministic matrices.
        MatF = kbc.matrix_F(a1, a2, q, f)
        MatFd = kbc.matrix_F_dual(a1, a2, q, f)
        ZN = Integers(N)
        C1 = matrix(ZN, [[a1, 0, -a2, 0], [a2, 0, a1, 0], [1, 0, 0, 0], [0, 0, 1, 0]])
        D1 = matrix(ZN, [[0, a1, 0, -a2], [0, a2, 0, a1], [0, q, 0, 0], [0, 0, 0, q]])
        C2 = matrix(ZN, [[a1, 0, a2, 0], [-a2, 0, a1, 0], [-1, 0, 0, 0], [0, 0, -1, 0]])
        D2 = matrix(ZN, [[0, a1, 0, a2], [0, -a2, 0, a1], [0, -q, 0, 0], [0, 0, 0, -q]])
        G2F1 = kbc.gluing_base_change_matrix_dim2_F1(a1, a2, q)
        G2F2 = kbc.gluing_base_change_matrix_dim2_F2(a1, a2, q)

        # Completion-dependent matrices (PARI matsolvemod inside).
        M1, M2 = kbc.starting_two_symplectic_matrices(a1, a2, q, f)
        MG1 = kbc.gluing_base_change_matrix_dim2_dim4_F1(a1, a2, q, m, M1)
        MG2 = kbc.gluing_base_change_matrix_dim2_dim4_F2(a1, a2, q, m, M2)

        out["vectors"].append({
            "index": v["index"],
            "a1": str(a1), "a2": str(a2), "m": m, "q": str(q), "f": f,
            "matrix_F": mat_strs(MatF, 8, 8),
            "matrix_F_dual": mat_strs(MatFd, 8, 8),
            "ckm_F1_C": mat_strs(C1, 4, 4),
            "ckm_F1_D": mat_strs(D1, 4, 4),
            "ckm_F2dual_C": mat_strs(C2, 4, 4),
            "ckm_F2dual_D": mat_strs(D2, 4, 4),
            "gluing_dim2_F1": mat_strs(G2F1, 4, 4),
            "gluing_dim2_F2": mat_strs(G2F2, 4, 4),
            "M1": mat_strs(M1, 8, 8),
            "M2": mat_strs(M2, 8, 8),
            "M_gluing_1": mat_strs(MG1, 8, 8),
            "M_gluing_2": mat_strs(MG2, 8, 8),
        })

    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote %s (%d vectors)" % (out_path, len(out["vectors"])))


if __name__ == "__main__":
    main()
