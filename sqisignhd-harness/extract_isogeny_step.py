#!/usr/bin/env sage
r"""
extract_isogeny_step.py  --  per-step oracle for Phase 2 (dim-4 isogeny step).

The Phase 0 vectors record per-step codomain theta null points but NOT the
kernel that drives each step, so a single `(2,2,2,2)`-isogeny step cannot be
validated in isolation from them. This script supplies that missing ground
truth: for a real verification it captures, for individual *plain* `IsogenyDim4`
steps, the kernel 8-torsion points `K_8`, the resulting codomain theta null
point, and a few (input, image) pairs.

It does this by wrapping `IsogenyDim4.__init__` at runtime (a new script; no
reference file and no Phase 0 deliverable is edited). The gluing/splitting
isogenies (`GluingIsogenyDim4`) override `__init__` and are therefore not
captured - exactly the steps that are out of scope for this phase.

Usage:  sage extract_isogeny_step.py [--out isogeny_step_vectors.json]
                                     [--n NVEC] [--lib PATH]
Env  :  SQISIGNHD_LIB
"""
import sys
import os
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
    ap.add_argument("--n", type=int, default=5, help="number of signatures")
    ap.add_argument("--lib", default=None)
    args = ap.parse_args()
    out_path = os.path.abspath(args.out or "isogeny_step_vectors.json")

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)
    sys.path.insert(0, VERIF)

    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    import Theta_dim4.Theta_dim4_sage.pkg.isogenies.isogeny_dim4 as isogmod

    pp = SQIsignHD(1)
    p = int(pp.p)
    Fp2 = pp.Fp2

    def fp2_pair(x):
        c = list(Fp2(x).list())
        a = int(c[0]) % p if len(c) > 0 else 0
        b = int(c[1]) % p if len(c) > 1 else 0
        return ["0x%x" % a, "0x%x" % b]

    def theta(coords):
        co = list(coords)
        assert len(co) == 16
        return [fp2_pair(c) for c in co]

    # Capture references to every plain IsogenyDim4 as it is built.
    captures = []  # list of (domain, list(K_8), iso)
    orig_init = isogmod.IsogenyDim4.__init__

    def patched_init(self, domain, K_8, codomain=None, precomputation=None):
        orig_init(self, domain, K_8, codomain, precomputation)
        if K_8 is not None:
            captures.append((domain, list(K_8), self))

    isogmod.IsogenyDim4.__init__ = patched_init

    def serialize(domain, K_8, iso):
        prec = iso._precomputation
        has_zero = any(x is None for x in prec)
        rec = {
            "domain_null": theta(domain.null_point().coords()),
            "K_8": [theta(k.coords()) for k in K_8],
            "codomain_null": theta(iso._codomain.null_point().coords()),
            "has_zero_dual": bool(has_zero),
        }
        if not has_zero:
            rec["image_pairs"] = [
                {"in": theta(k.coords()), "out": theta(iso.image(k).coords())}
                for k in K_8
            ]
        return rec

    all_cases = []
    for idx in range(args.n):
        captures.clear()
        v = SQIsignHD_verif(pp, idx)
        v.recover_pk_and_com()
        v.recover_chal()
        v.image_response()
        v.compute_HD()

        # Plain steps are _isogenies[1..]; F1 is built before F2_dual, so the
        # first (n1 - 1) captures belong to F1, the rest to F2_dual.
        n1 = len(v.F.F1._isogenies)
        n_f1_plain = n1 - 1
        total = len(captures)

        # Sample early F1 steps and the first F2_dual step(s).
        if idx == 0:
            picks = [(0, "F1", 1), (1, "F1", 2), (2, "F1", 3),
                     (n_f1_plain, "F2_dual", 1), (n_f1_plain + 1, "F2_dual", 2)]
        else:
            picks = [(0, "F1", 1), (n_f1_plain, "F2_dual", 1)]

        for (ci, chain, step) in picks:
            if 0 <= ci < total:
                domain, K_8, iso = captures[ci]
                rec = serialize(domain, K_8, iso)
                rec["vector"] = idx
                rec["chain"] = chain
                rec["plain_step"] = step
                all_cases.append(rec)
        print("vector %d: captured %d plain steps (n1=%d), selected %d"
              % (idx, total, n1, len([1 for q in picks if 0 <= q[0] < total])))

    import sage.version
    doc = {
        "purpose": "Phase 2 oracle: kernel (K_8), codomain theta null point, and "
                   "(input,image) pairs for individual plain dim-4 (2,2,2,2)-"
                   "isogeny steps, captured from real Level-1 verifications.",
        "level": 1,
        "prime_decimal": str(p),
        "field": "F_{p^2}=F_p[i]/(i^2+1); a+b*i serialized as [hex(a),hex(b)]",
        "sage_version": str(sage.version.version),
        "note": "K_8 are four 8-torsion theta points (4*K_8 is a (2,2,2,2) kernel "
                "basis). codomain_null is the standard theta null point H(O). "
                "image_pairs use IsogenyDim4.image; absent when the codomain has "
                "a zero dual theta-null coordinate (has_zero_dual).",
        "n_cases": len(all_cases),
        "cases": all_cases,
    }
    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=1)
    print("wrote %d isogeny-step cases to %s" % (len(all_cases), out_path))
    nz = sum(1 for c in all_cases if c["has_zero_dual"])
    print("cases with zero dual theta-null coordinate: %d / %d" % (nz, len(all_cases)))


if __name__ == "__main__":
    main()
