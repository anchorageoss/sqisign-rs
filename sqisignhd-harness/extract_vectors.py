#!/usr/bin/env sage
r"""
extract_vectors.py  --  Phase 0 validation-harness extractor for the dim-4
                        SQIsignHD verifier.

This script imports the *unmodified* SQIsignHD sage verifier
(Pierrick-Dartois/SQISignHD-lib, file Verification/Verify.py) together with its
Theta_dim4 submodule, re-runs the verification of N level-`lvl` signatures
stage by stage, and dumps every stage-boundary intermediate value to JSON.

It does NOT modify the reference: it imports `Verify.SQIsignHD` /
`Verify.SQIsignHD_verif`, calls the public stage methods in the same order as
`SQIsignHD_verif.verify()`, and introspects the resulting objects (the
KaniEndoHalf chain `F`, its sub-chains `F1` / `F2_dual`, the per-step dim-4
codomain theta-structures, etc.).

Run with SageMath 10.0+ :

    sage extract_vectors.py --lvl 1 --n 5 --out test_vectors_l1.json \
         --lib /path/to/SQISignHD-lib

The reference checkout is located via, in order: --lib, $SQISIGNHD_LIB,
./SQISignHD-lib, ~/SQISignHD-lib, /home/user/SQISignHD-lib.  See
setup_reference.sh for how to obtain it (the Theta_dim4 submodule must be
cloned via https because .gitmodules pins an ssh URL).

Serialization conventions are documented at length in HARNESS_NOTES.md.  In
short: an element of F_{p^2} = F_p[i]/(i^2+1), written a + b*i with
0 <= a,b < p, is serialized as the 2-tuple ["0x<a>", "0x<b>"].  A theta null
point is 16 such coordinates and is *projective*; see HARNESS_NOTES.md.
"""

import sys
import os
import json
import argparse
from time import time


# --------------------------------------------------------------------------
# Locate the reference checkout and put it on the path.  The SQIsignHD verifier
# reads its data tables (Data/*.txt) via *relative* paths, so we must chdir
# into the Verification directory before importing it.
# --------------------------------------------------------------------------
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
    raise SystemExit(
        "Could not locate a SQISignHD-lib checkout with its Theta_dim4 "
        "submodule.\nTried:\n  " + "\n  ".join(cands) +
        "\nPass --lib PATH or set $SQISIGNHD_LIB, and run setup_reference.sh.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lvl", type=int, default=1, choices=[1, 3, 5])
    ap.add_argument("--n", type=int, default=5, help="number of signatures")
    ap.add_argument("--out", default=None, help="output JSON path")
    ap.add_argument("--lib", default=None, help="path to SQISignHD-lib checkout")
    args = ap.parse_args()

    lvl = args.lvl
    n = args.n
    out_path = args.out or "test_vectors_l{}.json".format(lvl)
    out_path = os.path.abspath(out_path)

    LIB = find_lib(args.lib)
    VERIF = os.path.join(LIB, "Verification")
    os.chdir(VERIF)                 # required for the relative Data/*.txt reads
    sys.path.insert(0, VERIF)

    # ---- imports that need the path / cwd set up above ----
    from sage.all import EllipticCurve, ZZ, GF, is_prime
    import cypari2
    import Verify
    from Verify import SQIsignHD, SQIsignHD_verif
    from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
    from Theta_dim4.Theta_dim4_sage.pkg.utilities.discrete_log import (
        weil_pairing_pari, discrete_log_pari)

    # Resolve the reference git commits for provenance (best effort).
    def git_head(path):
        try:
            import subprocess
            return subprocess.check_output(
                ["git", "-C", path, "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return None

    pp = SQIsignHD(lvl)
    p = int(pp.p)
    Fp2 = pp.Fp2
    iU = pp.i
    order_r = ZZ(2 ** pp.r)

    # ----------------------------------------------------------------------
    # Serialization helpers.
    # ----------------------------------------------------------------------
    def to_fp2(x):
        """Coerce x (sage F_{p^2} elt, pari Gen weil-pairing value, F_p elt, or
        python int) into a sage F_{p^2} element of pp.Fp2."""
        if isinstance(x, cypari2.gen.Gen):
            return Fp2(x)
        if not hasattr(x, "list"):
            return Fp2(x)
        # already an Fp2/Fp/sage element with .list()
        return Fp2(x)

    def fp2_pair(x):
        """F_{p^2} element a+b*i  ->  ["0x<a>", "0x<b>"], 0<=a,b<p."""
        e = to_fp2(x)
        c = list(e.list())
        a = int(c[0]) % p if len(c) > 0 else 0
        b = int(c[1]) % p if len(c) > 1 else 0
        return ["0x%x" % a, "0x%x" % b]

    def int_str(x):
        return str(int(x))

    def theta_raw(coords):
        """coords: iterable of 16 F_{p^2} elts -> list of 16 [re,im] pairs."""
        L = list(coords)
        assert len(L) == 16, "expected 16 theta coordinates, got %d" % len(L)
        return [fp2_pair(c) for c in L]

    def first_nonzero(coords):
        for k in range(len(coords)):
            if coords[k] != 0:
                return k
        return -1

    def theta_record(theta_point, normalized=False):
        """Build a JSON record for a theta null point (ThetaPointDim4 or the
        output of .zero()/.null_point()).  Always stores the raw 16 projective
        coordinates exactly as sage computed them.  If normalized=True, also
        stores a canonical representative (divide by first non-zero coord, so
        that coordinate becomes 1) plus the pivot index used."""
        # Coerce to Fp2 first: dual null points mix sage Integer 0/1 with Fp2
        # elements, which breaks bare arithmetic.
        coords = [to_fp2(c) for c in theta_point.coords()]
        rec = {"coords": theta_raw(coords)}
        if normalized:
            piv = first_nonzero(coords)
            inv = ~coords[piv]
            norm = [c * inv for c in coords]
            rec["pivot"] = piv
            rec["normalized"] = theta_raw(norm)
        return rec

    def ec_point(P):
        """Affine (x,y) of an elliptic-curve point, or {'inf':True}."""
        if P.is_zero():
            return {"inf": True}
        return {"x": fp2_pair(P[0]), "y": fp2_pair(P[1])}

    def chain_codomains(chain):
        """For an IsogenyChainDim4, the codomain theta null point after each of
        its len = (e - m) steps.  chain._isogenies[0] is the gluing chain
        (KaniGluingIsogenyChainDim4Half); the rest are IsogenyDim4."""
        return [theta_record(iso._codomain.null_point()) for iso in chain._isogenies]

    def chain_start(chain):
        """Theta null points at the entry of a half-chain's dim-4 part.
        domain_product : product theta structure A_m^2  (before base change)
        domain_start   : dim-4 theta structure that is the domain of the
                         gluing isogeny (after the N_dim4 base change)."""
        gl = chain._isogenies[0]
        return {
            "theta_null_domain_product": theta_record(gl.domain_product.null_point()),
            "theta_null_start": theta_record(gl.domain_base_change.null_point(),
                                             normalized=True),
        }

    # ----------------------------------------------------------------------
    # Raw verbatim lines from the data files (faithful record of the inputs).
    # ----------------------------------------------------------------------
    def raw_lines(path, start, count):
        out = []
        with open(path) as fh:
            lines = fh.readlines()
        for j in range(start, start + count):
            out.append(lines[j].rstrip("\n"))
        return out

    pk_file = "Data/Public_keys_lvl{}.txt".format(lvl)
    sig_file = "Data/Signatures_lvl{}.txt".format(lvl)

    # ----------------------------------------------------------------------
    # Per-signature extraction.
    # ----------------------------------------------------------------------
    vectors = []
    t_all = time()
    for idx in range(n):
        print("Extracting vector {}/{} ...".format(idx, n - 1), flush=True)
        t0 = time()
        v = SQIsignHD_verif(pp, idx)

        rec = {"index": idx}

        # --- inputs ---
        rec["message_hex"] = ""   # the sage verifier consumes no message; chal is appended
        rec["public_key"] = {
            "A_pk": fp2_pair(v.A_pk),
            "hint_pk_P": int(v.h_pk_P),
            "hint_pk_Q": int(v.h_pk_Q),
            "raw_lines": raw_lines(pk_file, idx * 3, 3),
        }
        rec["signature"] = {
            "A_com": fp2_pair(v.A_com),
            "a": int_str(v.a), "b": int_str(v.b),
            "c_or_d": int_str(v.c_or_d), "q": int_str(v.q),
            "hint_com_P": int(v.h_com_P), "hint_com_Q": int(v.h_com_Q),
            "chal": int_str(v.chal),
            "raw_lines": raw_lines(sig_file, idx * 8, 8),
        }

        # --- stage 1: recover_pk_and_com ---
        v.recover_pk_and_com()
        rec["stage1_recover_basis"] = {
            "E_pk_A": fp2_pair(v.E_pk.a2()),
            "P_pk": ec_point(v.P_pk), "Q_pk": ec_point(v.Q_pk),
            "E_com_A": fp2_pair(v.E_com.a2()),
            "P_com": ec_point(v.P_com), "Q_com": ec_point(v.Q_com),
        }

        # --- stage 2: recover_chal ---
        v.recover_chal()
        rec["stage2_recover_chal"] = {
            "E_chal_A": fp2_pair(v.E_chal.a2()),
            "w_chal": fp2_pair(v.w_chal),
            "P_chal_resc": ec_point(v.P_chal_resc),
            "Q_chal_resc": ec_point(v.Q_chal_resc),
        }

        # --- stage 3: image_response ---
        v.image_response()
        # w_com and k are local to image_response(); recompute them with the
        # exact same inputs/functions for the record (R_com, S_com, w_chal are
        # all now stored on v).
        w_com = weil_pairing_pari(v.R_com, v.S_com, order_r)
        k = discrete_log_pari(w_com, v.w_chal, order_r)
        rec["stage3_image_response"] = {
            "w_com": fp2_pair(w_com),
            "k": int_str(k),
            "c": int_str(v.c), "d": int_str(v.d),
            "R_com": ec_point(v.R_com), "S_com": ec_point(v.S_com),
            "phi_rsp_R_com": ec_point(v.phi_rsp_R_com),
            "phi_rsp_S_com": ec_point(v.phi_rsp_S_com),
        }

        # --- stage 4: compute_HD  (the dim-4 isogeny core) ---
        v.compute_HD()
        N = ZZ(2) ** pp.f - ZZ(v.q)
        a1, a2 = ZZ(v.a1), ZZ(v.a2)
        # m as computed inside KaniEndoHalf (after putting the odd one first)
        aa1, aa2 = a1, a2
        if aa1 % 2 == 0:
            aa1, aa2 = aa2, aa1
        m = 0
        t = aa2
        while t % 2 == 0:
            m += 1
            t //= 2
        e1 = (pp.f + 1) // 2
        e2 = pp.f - e1
        F = v.F
        F1, F2d = F.F1, F.F2_dual

        s4 = {
            "q": int_str(v.q), "N": int_str(N),
            "N_is_prime": bool(is_prime(N)), "N_mod_4": int(N % 4),
            "a1": int_str(a1), "a2": int_str(a2),
            "a1_sq_plus_a2_sq_eq_N": bool(a1 * a1 + a2 * a2 == N),
            "f": int(pp.f), "m": int(m), "e1": int(e1), "e2": int(e2),
            "len_F1": len(F1._isogenies), "len_F2_dual": len(F2d._isogenies),
            "swap": bool(F.swap),
        }
        s4["F1"] = chain_start(F1)
        s4["F1"]["chain_codomains"] = chain_codomains(F1)
        s4["F2_dual"] = chain_start(F2d)
        s4["F2_dual"]["chain_codomains"] = chain_codomains(F2d)
        rec["stage4_compute_hd"] = s4

        # --- stage 5: codomain matching check ---
        C1 = F1._isogenies[-1]._codomain
        C2 = F2d._isogenies[-1]._codomain
        HC2 = C2.hadamard()
        match = v.verify_middle_codomain()
        rec["stage5_codomain_check"] = {
            "C1_zero": theta_record(C1.zero(), normalized=True),
            "HC2_zero": theta_record(HC2.zero(), normalized=True),
            "match": bool(match),
        }

        # --- stage 6: HD image check ---
        T = TuplePoint(v.P_com, v.E_com(0), v.E_chal(0), v.E_chal(0))
        FT = F(T)
        a1P = v.a1 * v.P_com
        a2P = v.a2 * v.P_com
        correct = v.verify_HD_image()
        rec["stage6_image_check"] = {
            "input_T": [ec_point(T[0]), ec_point(T[1]), ec_point(T[2]), ec_point(T[3])],
            "FT": [ec_point(FT[0]), ec_point(FT[1]), ec_point(FT[2]), ec_point(FT[3])],
            "a1_P_com": ec_point(a1P), "a2_P_com": ec_point(a2P),
            "correct": bool(correct),
        }

        rec["result"] = "accept" if (match and correct) else "reject"
        rec["extract_seconds"] = time() - t0
        vectors.append(rec)
        print("  -> result={} ({:.2f}s, chain len F1={} F2_dual={})".format(
            rec["result"], rec["extract_seconds"],
            s4["len_F1"], s4["len_F2_dual"]), flush=True)

    import sage.version
    doc = {
        "scheme": "F-SQIsignHD (dimension-4 2-isogeny verification, theta level 2)",
        "reference": "https://github.com/Pierrick-Dartois/SQISignHD-lib (Verification/Verify.py)",
        "theta_dim4": "https://github.com/Pierrick-Dartois/Theta_dim4",
        "reference_commit": git_head(LIB),
        "theta_dim4_commit": git_head(os.path.join(VERIF, "Theta_dim4")),
        "sage_version": str(sage.version.version),
        "level": lvl,
        "prime": "{}*2^{} - 1".format(pp.c, pp.e),
        "prime_decimal": str(p),
        "field": ("F_{p^2} = F_p[i]/(i^2+1); an element a+b*i with 0<=a,b<p is "
                  "serialized as the pair [hex(a), hex(b)]"),
        "theta_convention": (
            "A theta null point has 16 coordinates in F_{p^2} and is PROJECTIVE "
            "(defined up to a common non-zero scalar). Compare two null points "
            "projectively (cross-multiply) or via the 'normalized' field, which "
            "divides every coordinate by the coordinate at index 'pivot' (the "
            "first non-zero one) so that coordinate becomes 1. The per-step "
            "'chain_codomains' store the raw representative sage computes "
            "(codomain.null_point() = Hadamard of a dual null point whose pivot "
            "coordinate is fixed to 1, hence deterministic)."),
        "params": {"lamb": int(pp.lamb), "f": int(pp.f), "r": int(pp.r),
                   "e1": int((pp.f + 1) // 2), "e2": int(pp.f - (pp.f + 1) // 2),
                   "n_bytes": int(pp.n_bytes)},
        "stages": [
            "stage1_recover_basis: build E_pk,E_com (y^2=x^3+A x^2+x) and their "
            "2^f torsion bases from hints",
            "stage2_recover_chal: challenge isogeny phi_chal:E_pk->E_chal of "
            "degree 2^lamb; rescaled E_chal basis + Weil pairing w_chal",
            "stage3_image_response: w_com, k=dlog(w_com,w_chal); recover c,d; "
            "response images phi_rsp_R_com, phi_rsp_S_com",
            "stage4_compute_hd: N=2^f-q=a1^2+a2^2; KaniEndoHalf builds two dim-4 "
            "(2,2,2,2)-isogeny half-chains F1 and F2_dual",
            "stage5_codomain_check: F1 and F2_dual must reach the same codomain: "
            "C1.zero() == C2.hadamard().zero() (projective)",
            "stage6_image_check: F(P_com,0,0,0) must equal (+-a1 P_com, +-a2 "
            "P_com, *, 0_Echal)",
        ],
        "note_message": ("The sage reference does NOT recompute the challenge as "
                         "a hash of (pk, commitment, message): the SHAKE256 in "
                         "common python libs does not match Signature/src/common/"
                         "fips202.c, so the challenge scalar is appended to the "
                         "signature ('chal'). message_hex is therefore empty."),
        "n_vectors": len(vectors),
        "test_vectors": vectors,
    }

    def _json_default(o):
        # Safety net: any stray sage Integer-like value -> int, else str.
        try:
            return int(o)
        except Exception:
            return str(o)

    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=1, default=_json_default)
    print("\nWrote {} vectors to {} in {:.1f}s".format(
        len(vectors), out_path, time() - t_all))
    print("results: " + ", ".join(v["result"] for v in vectors))


if __name__ == "__main__":
    main()
