#!/usr/bin/env python3
"""The statistics of a paper session (paper_session.sh), from its raw files.

    paper_report.py <session-dir> [--check-noise]

Per pass and per row: the median of the cycles and of the nanoseconds over
the runs the paper prescribes (for the compact verifier at levels III and
V, the mean over the six signatures of the median of the 15 runs of each).
The session's value of a row is the value of pass 1 when two passes agree
within 3 percent on every row, and otherwise the median of the three
passes' values. Writes numbers.json and table.md into the session
directory. With --check-noise: exit 1 when two passes disagree by more
than 3 percent on some row (the session script then runs a third pass). The
threshold is the environment variable PAPER_NOISE (default 0.03).
"""
import glob
import json
import os
import re
import statistics
import sys

NOISE = float(os.environ.get("PAPER_NOISE", "0.03"))


def read(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line[0].isdigit():
                continue
            it, cyc, ns, extra = line.rstrip("\n").split(",", 3)
            rows.append((None if cyc == "NA" else int(cyc), int(ns), extra))
    return rows


def med(xs):
    return statistics.median(xs)


def med_cyc(xs):
    """Median of a cycle column; None when the timer had no cycles (NA)."""
    return None if any(x is None for x in xs) else statistics.median(xs)


def row_stats(name, rows):
    """(megacycles, milliseconds) of a row by the paper's rule."""
    m = re.match(r"dim4-(III|V)-verify-", name)
    if m:
        by_sig = {}
        for cyc, ns, extra in rows:
            by_sig.setdefault(extra.split(";")[0], []).append((cyc, ns))
        mcs = [med_cyc([c for c, _ in v]) for v in by_sig.values()]
        mc = None if any(x is None for x in mcs) else statistics.mean(mcs)
        ms = statistics.mean(med([n for _, n in v]) for v in by_sig.values())
    else:
        mc = med_cyc([c for c, _, _ in rows])
        ms = med([n for _, n, _ in rows])
    return (None if mc is None else mc / 1e6), ms / 1e6


def c_ref(path):
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*(keypair|sign|verify)\s*\|.*median\s+([0-9.]+)", line)
            if m:
                out[m.group(1)] = float(m.group(2))
    return out


def load_pass(pdir):
    vals = {}
    for path in sorted(glob.glob(os.path.join(pdir, "raw", "*.csv"))):
        name = os.path.basename(path)[:-4]
        if "profile" in name:
            continue
        vals[name] = row_stats(name, read(path))
    for lvl in ("I", "III", "V"):
        p = os.path.join(pdir, f"c-ref-{lvl}.txt")
        if os.path.exists(p):
            for op, v in c_ref(p).items():
                key = f"c-{lvl}-{'keygen' if op == 'keypair' else op}"
                # built with -DNO_CYCLE_COUNTER the reference's "megacycles"
                # are nanoseconds / 1e6, that is milliseconds
                vals[key] = (None, v) if os.environ.get("C_REF_UNIT") == "ms" else (v, None)
    return vals


def loop_stats(pdir):
    rows = read(os.path.join(pdir, "raw", "dim4-I-sign.csv"))
    samples = [int(e.split(";")[0]) for _, _, e in rows]
    tests = [int(e.split(";")[1]) for _, _, e in rows]
    restarts = sum(int(e.split(";")[2]) for _, _, e in rows)
    has_cyc = all(c is not None for c, _, _ in rows)
    cyc = [c / 1e6 for c, _, _ in rows] if has_cyc else [float("nan")]
    msv = [n / 1e6 for _, n, _ in rows]
    return {
        "sign_ms_median": med(msv),
        "sign_ms_min": min(msv),
        "sign_ms_max": max(msv),
        "signatures": len(rows),
        "samples_mean": statistics.mean(samples),
        "samples_median": med(samples),
        "samples_min": min(samples),
        "samples_max": max(samples),
        "primality_tests_mean": statistics.mean(tests),
        "restarts": restarts,
        "sign_mc_median": med(cyc),
        "sign_mc_min": min(cyc),
        "sign_mc_max": max(cyc),
    }


def quartiles(pdir, lvl):
    rows = read(os.path.join(pdir, "raw", f"dim4-{lvl}-sign.csv"))
    unit = "mc" if all(c is not None for c, _, _ in rows) else "ms"
    cyc = sorted((c if unit == "mc" else n) / 1e6 for c, n, _ in rows)
    q = statistics.quantiles(cyc, n=4)
    return {"unit": unit, "q1": q[0], "q3": q[2], "min": cyc[0], "max": cyc[-1]}


def profile(pdir, name):
    path = os.path.join(pdir, "raw", name + ".csv")
    names, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("# slot"):
                names.append(line.split("=", 1)[1].strip())
            elif line[0].isdigit():
                p = line.strip().split(",")
                rows.append([int(x) for x in p[2:]])
    total = med([r[0] for r in rows])
    # the unit is the timer's: megacycles with rdtsc, milliseconds with the wall clock
    out = {"total_mc": total / 1e6, "decode_share": med([r[1] for r in rows]) / total}
    for k, n in enumerate(names):
        out[n] = med([r[2 + k] for r in rows]) / total
    if len(rows[0]) >= 2 + len(names) + 2:
        # the pre-chain checks, megacycles per verification
        out["decode_mc"] = med([r[1] for r in rows]) / 1e6
        out["checks_hints_pairing_mc"] = med([r[-2] for r in rows]) / 1e6
        out["checks_goodness_mc"] = med([r[-1] for r in rows]) / 1e6
    return out


def heap(sdir):
    path = os.path.join(sdir, "heap.csv")
    if not os.path.exists(path):
        return {}
    by = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("row"):
                continue
            row, _, peak, n = line.strip().split(",")
            by.setdefault(row, []).append(int(peak))
    return {k: {"peak_kib_mean": statistics.mean(v) / 1024, "peak_kib_max": max(v) / 1024, "signatures": len(v)} for k, v in by.items()}


def main():
    sdir = sys.argv[1]
    check = "--check-noise" in sys.argv
    passes = sorted(glob.glob(os.path.join(sdir, "pass[0-9]")))
    vals = [load_pass(p) for p in passes]
    names = sorted(vals[0])
    noisy = []
    if len(vals) >= 2:
        for n in names:
            k = 0 if vals[0][n][0] is not None else 1
            a, b = vals[0][n][k], vals[1][n][k]
            if abs(a - b) / min(a, b) > NOISE:
                noisy.append((n, a, b))
    if check:
        for n, a, b in noisy:
            print(f"noise: {n}: pass 1 {a:.3f}, pass 2 {b:.3f} ({abs(a - b) / min(a, b) * 100:.1f} percent)")
        sys.exit(1 if noisy else 0)

    use_three = len(vals) >= 3
    session = {}
    for n in names:
        mcs = [v[n][0] for v in vals if v[n][0] is not None]
        mss = [v[n][1] for v in vals if v[n][1] is not None]
        if use_three:
            session[n] = {"mc": med(mcs) if mcs else None, "ms": med(mss) if mss else None}
        else:
            session[n] = {"mc": mcs[0] if mcs else None, "ms": mss[0] if mss else None}
        session[n]["passes_mc"] = mcs
        session[n]["passes_ms"] = mss
    out = {
        "passes": len(vals),
        "rule": "median of the three passes' values" if use_three else "pass 1 (two passes within 3 percent on every row)",
        "noisy_rows_between_pass_1_and_2": [{"row": n, "pass1": a, "pass2": b} for n, a, b in noisy],
        "rows": session,
        "loop": loop_stats(passes[0]),
        "sign_quartiles": {lvl: quartiles(passes[0], lvl) for lvl in ("I", "III", "V")},
        "profile": {lvl + "-" + m: profile(passes[0], f"dim4-{lvl}-profile-{m}") for lvl in ("I", "III", "V") for m in ("four", "two")},
        "profile_per_pass": [{lvl + "-" + m: profile(p, f"dim4-{lvl}-profile-{m}") for lvl in ("I", "III", "V") for m in ("four", "two")} for p in passes],
        "heap": heap(sdir),
    }
    with open(os.path.join(sdir, "numbers.json"), "w") as f:
        json.dump(out, f, indent=1)
    lines = [f"passes: {len(vals)}; session value: {out['rule']}", "", "| row | megacycles | ms | per pass (megacycles, or ms where there are no cycles) |", "|---|---:|---:|---|"]
    for n in names:
        r = session[n]
        ms = "" if r["ms"] is None else f"{r['ms']:.2f}"
        mcv = "" if r["mc"] is None else f"{r['mc']:.2f}"
        per = r["passes_mc"] if r["passes_mc"] else r["passes_ms"]
        lines.append(f"| {n} | {mcv} | {ms} | {', '.join(f'{x:.2f}' for x in per)} |")
    with open(os.path.join(sdir, "table.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
