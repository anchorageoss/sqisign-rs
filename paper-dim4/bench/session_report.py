#!/usr/bin/env python3
"""Turn one session of paper-dim4/bench/session.sh into its table, its
profile and numbers.json, every statistic recomputed from the raw files
(SECENG-2423, SECENG-2424).

    session_report.py <session dir>

Reads, under the session directory:
  raw/<block>.csv            one line per iteration (round, iteration, cycles, ns, extra)
  raw/profile-<mode>.csv     one line per iteration with the phase cycles
  sqisign-c-<prime>.txt      the C reference's own benchmark output
  session.log                the commit and the TSC line

Writes, under the session directory:
  table.md                   the rows (medians; signing with mean and maximum)
  profile.md                 the phase profile of both verifiers
  numbers.json               every value, and the ratios from the printed cells

Medians are the statistic; the signing rows carry the mean and the maximum
too; a megacycle is a million ticks of the TSC.
"""
import csv
import glob
import json
import os
import re
import statistics
import sys

# The phases, by the functions the timers enclose (both half-chains added).
# The timers sit on the shipped code's phase boundaries; what no slot
# encloses (optimised_strategy, the first precompute of the gluing codomain,
# gluing_images_paired in the checked two-point mode, the timing calls) is
# the remainder, the whole call's median minus the sum of the phase medians.
PROFILE_PHASES = [
    ("decode", "decode, in the harness (parse_signature, parse_public_key, challenge_limbs)"),
    ("slot0", "basis generation (recover_challenge, recover_response, hd_torsion_basis)"),
    ("slot1", "Cornacchia and kernel setup (norm_equation_2f_minus_q, make_canonical, the product null, starting_two_symplectic_matrices, the gluing base changes, b_kpp)"),
    ("slot2", "dimension-2 steps (IsogenyChainDim2::new)"),
    ("slot3", "gluing (KaniGluingChainHalf::new except the dimension-2 chain: the base changes, the kernel points, the five evaluations, the gluing step, the translates; and the gluing images of the generators, four in the four-point and the checked two-point modes, two without the check)"),
    ("slot4", "plain dimension-4 steps, point pushes (double_iter and image of the saved bases or pairs)"),
    ("slot5", "plain dimension-4 steps, codomain (kernel_from_pair in the two-point mode, IsogenyDim4::from_kernel, the next structure's precompute)"),
    ("slot6", "final comparisons (middle_codomain_matches)"),
    ("slot7", "the stage-6 image check (hd_image_check), run in every mode; the two-point option's own comparison gluing_images_paired is in the remainder"),
]

ROWS = [
    # (raw file, implementation, dimension, level, operation, signing-style stats)
    ("dim2-I-keygen", "Rust, dimension 2", 2, "I", "keygen", False),
    ("dim2-I-sign", "Rust, dimension 2", 2, "I", "sign", True),
    ("dim2-I-verify", "Rust, dimension 2", 2, "I", "verify from bytes", False),
    ("dim2-I-verify-compressed", "Rust, dimension 2", 2, "I", "verify compressed from bytes", False),
    ("dim2-III-keygen", "Rust, dimension 2", 2, "III", "keygen", False),
    ("dim2-III-sign", "Rust, dimension 2", 2, "III", "sign", True),
    ("dim2-III-verify", "Rust, dimension 2", 2, "III", "verify from bytes", False),
    ("dim2-III-verify-compressed", "Rust, dimension 2", 2, "III", "verify compressed from bytes", False),
    ("dim2-V-keygen", "Rust, dimension 2", 2, "V", "keygen", False),
    ("dim2-V-sign", "Rust, dimension 2", 2, "V", "sign", True),
    ("dim2-V-verify", "Rust, dimension 2", 2, "V", "verify from bytes", False),
    ("dim2-V-verify-compressed", "Rust, dimension 2", 2, "V", "verify compressed from bytes", False),
    ("dim4-keygen", "Rust, dimension 4 (compact)", 4, "I", "keygen", False),
    ("dim4-sign", "Rust, dimension 4 (compact)", 4, "I", "sign", True),
    ("dim4-verify-four", "Rust, dimension 4 (compact)", 4, "I", "verify from bytes, four points", False),
    ("dim4-verify-two-check", "Rust, dimension 4 (compact)", 4, "I", "verify from bytes, two points with the image check", False),
    ("dim4-verify-two-nocheck", "Rust, dimension 4 (compact)", 4, "I", "verify from bytes, two points without the image check", False),
]

C_LEVELS = [("p324_3", "I"), ("p500_27", "III"), ("p664_17", "V")]


def read_raw(path):
    """[(round, iteration, cycles, ns, extra)]."""
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("round,"):
                continue
            r, it, cyc, ns, extra = line.rstrip("\n").split(",", 4)
            out.append((int(r), int(it), int(cyc), int(ns), extra))
    return out


def stats(cycles, ns):
    mc = [c / 1e6 for c in cycles]
    return {
        "n": len(mc),
        "median": statistics.median(mc),
        "mean": statistics.fmean(mc),
        "min": min(mc),
        "max": max(mc),
        "median_ms": statistics.median(ns) / 1e6,
    }


def read_c(path):
    """{op: {runs, median, mean, min, max}} from the C benchmark's summary lines."""
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*(keypair|sign|verify)\s*\|(.*)", line)
            if not m:
                continue
            op = {"keypair": "keygen", "sign": "sign", "verify": "verify from bytes"}[m.group(1)]
            rest = m.group(2)
            vals = {}
            for key in ("average", "median", "min", "max"):
                mm = re.search(key + r"\s+([0-9.]+)", rest)
                if mm:
                    vals[key] = float(mm.group(1))
            out[op] = vals
    return out


def fmt1(x):
    return f"{x:.1f}"


def main():
    sdir = sys.argv[1]
    numbers = {"rows": {}, "c_reference": {}, "profile": {}, "session": {}}
    with open(os.path.join(sdir, "session.log")) as f:
        log = f.read()
    m = re.search(r"commit ([0-9a-f]{40})", log)
    numbers["session"]["commit"] = m.group(1) if m else None
    m = re.search(r"Detected ([0-9.]+) MHz processor", log)
    numbers["session"]["tsc_mhz"] = float(m.group(1)) if m else None
    m = re.search(r"(\S+) session start", log)
    numbers["session"]["start"] = m.group(1) if m else None
    m = re.search(r"(\S+) session end", log)
    numbers["session"]["end"] = m.group(1) if m else None

    lines = []
    lines.append("#### Rust rows (this crate; megacycles over `runs`, medians; the signing rows with mean and maximum; one pinned core, round-robin over the session's rounds)\n")
    lines.append("| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max | median (ms) |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for raw, impl_, dim, level, op, _signing in ROWS:
        path = os.path.join(sdir, "raw", raw + ".csv")
        if not os.path.exists(path):
            continue
        rows = read_raw(path)
        st = stats([r[2] for r in rows], [r[3] for r in rows])
        rounds = sorted({r[0] for r in rows})
        st["rounds"] = rounds
        st["per_round_median"] = {
            str(rd): statistics.median([r[2] / 1e6 for r in rows if r[0] == rd]) for rd in rounds
        }
        if raw == "dim4-sign":
            samples = [int(r[4].split(";")[0]) for r in rows]
            tests = [int(r[4].split(";")[1]) for r in rows]
            st["loop"] = {
                "signatures": len(samples),
                "samples_mean": statistics.fmean(samples),
                "samples_median": statistics.median(samples),
                "samples_min": min(samples),
                "samples_max": max(samples),
                "tests_mean": statistics.fmean(tests),
            }
        if raw.startswith("dim4-verify"):
            fb = [r[4].split(";") for r in rows]
            st["fell_back"] = sum(int(x[1]) + int(x[2]) for x in fb)
            per_sig = {}
            for r in rows:
                j = r[4].split(";")[0]
                per_sig.setdefault(j, []).append(r[2] / 1e6)
            st["per_signature_median"] = {j: statistics.median(v) for j, v in sorted(per_sig.items())}
        numbers["rows"][raw] = st
        lines.append(
            f"| {impl_} | {dim} | {level} | {op} | {st['n']} | {fmt1(st['median'])} | {fmt1(st['mean'])} | {fmt1(st['min'])} | {fmt1(st['max'])} | {st['median_ms']:.2f} |"
        )

    lines.append("\n#### The C reference (the-sqisign nist-v3, Release, broadwell; its own benchmark and cycle counter, medians over iterations)\n")
    lines.append("| implementation | dimension | level | operation | runs | median (Mcycles) | mean | min | max |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for prime, level in C_LEVELS:
        path = os.path.join(sdir, f"sqisign-c-{prime}.txt")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            txt = f.read()
        m = re.search(r"--iterations=(\d+)", txt)
        runs = None
        mm = re.search(r"(\d+)\s+iterations", txt)
        if mm:
            runs = int(mm.group(1))
        ops = read_c(path)
        numbers["c_reference"][level] = {"runs": runs, "ops": ops}
        for op in ("keygen", "sign", "verify from bytes"):
            if op in ops:
                v = ops[op]
                lines.append(
                    f"| C reference, dimension 2 | 2 | {level} | {op} | {runs if runs else ''} | {v.get('median', float('nan')):.3f} | {v.get('average', float('nan')):.3f} | {v.get('min', float('nan')):.3f} | {v.get('max', float('nan')):.3f} |"
                )

    # Ratios from the printed cells (one decimal), as the REPORT quotes them.
    def cell(raw):
        return float(fmt1(numbers["rows"][raw]["median"])) if raw in numbers["rows"] else None

    def c_cell(level, op):
        v = numbers["c_reference"].get(level, {}).get("ops", {}).get(op, {}).get("median")
        return round(v, 1) if v is not None else None

    ratios = {}
    four, two_c, two_n = cell("dim4-verify-four"), cell("dim4-verify-two-check"), cell("dim4-verify-two-nocheck")
    d2v = cell("dim2-I-verify")
    cv = c_cell("I", "verify from bytes")
    if four and two_c:
        ratios["two_point_with_check_over_four_point"] = round(two_c / four, 2)
    if four and two_n:
        ratios["two_point_without_check_over_four_point"] = round(two_n / four, 2)
    if four and cv:
        ratios["four_point_over_c_reference_verify_I"] = round(four / cv, 1)
    if two_c and cv:
        ratios["two_point_with_check_over_c_reference_verify_I"] = round(two_c / cv, 1)
    if four and d2v:
        ratios["four_point_over_dim2_standard_verify_I"] = round(four / d2v, 1)
    if two_c and d2v:
        ratios["two_point_with_check_over_dim2_standard_verify_I"] = round(two_c / d2v, 1)
    numbers["ratios_from_cells"] = ratios
    lines.append("\n#### Ratios from the printed cells\n")
    for k, v in ratios.items():
        lines.append(f"- {k}: {v}")
    with open(os.path.join(sdir, "table.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # The profile.
    plines = []
    for mode, label in (("two-check", "two points with the image check"), ("four", "four points")):
        path = os.path.join(sdir, "raw", f"profile-{mode}.csv")
        if not os.path.exists(path):
            continue
        rows = []
        with open(path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("round,"):
                    continue
                rows.append([int(x) for x in line.rstrip("\n").split(",")])
        # columns: round, iteration, signature, total, decode, slot0..slot7
        cols = {"total": 3, "decode": 4}
        for k in range(8):
            cols[f"slot{k}"] = 5 + k
        med = {name: statistics.median([r[i] / 1e6 for r in rows]) for name, i in cols.items()}
        phase_sum = sum(med[name] for name, _ in PROFILE_PHASES)
        unslotted = med["total"] - phase_sum
        prof = {
            "n": len(rows),
            "median_total_mcycles": med["total"],
            "phases": {},
            "sum_of_phase_medians": phase_sum,
            "unslotted_median_total_minus_sum": unslotted,
        }
        plines.append(f"#### Phase profile, {label} (n = {len(rows)}; medians in megacycles, both half-chains added; share = median of the phase over the sum of the phase medians)\n")
        plines.append("| phase | median (Mcycles) | share |")
        plines.append("|---|---|---|")
        for name, desc in PROFILE_PHASES:
            share = med[name] / phase_sum if phase_sum else 0.0
            prof["phases"][name] = {"description": desc, "median_mcycles": med[name], "share": share}
            plines.append(f"| {desc} | {med[name]:.2f} | {100 * share:.1f}% |")
        plain = med["slot4"] + med["slot5"]
        pushes_frac = med["slot4"] / plain if plain else 0.0
        prof["plain_steps"] = {
            "median_mcycles": plain,
            "pushes_fraction": pushes_frac,
            "codomain_fraction": 1 - pushes_frac,
        }
        plines.append(f"| sum of the phase medians | {phase_sum:.2f} | 100% |")
        plines.append(f"| median of the whole call (decode and verify) | {med['total']:.2f} | |")
        plines.append(f"| the remainder (the whole minus the sum: optimised_strategy, the first precompute, gluing_images_paired in the checked mode, the timing calls) | {unslotted:.2f} | |")
        plines.append(f"\nInside the plain dimension-4 steps ({plain:.2f} Mcycles): point pushes {100 * pushes_frac:.1f}%, codomain computation {100 * (1 - pushes_frac):.1f}%.\n")
        numbers["profile"][mode] = prof
    with open(os.path.join(sdir, "profile.md"), "w") as f:
        f.write("\n".join(plines) + "\n")
    with open(os.path.join(sdir, "numbers.json"), "w") as f:
        json.dump(numbers, f, indent=1, sort_keys=True)
        f.write("\n")
    print(open(os.path.join(sdir, "table.md")).read())
    print(open(os.path.join(sdir, "profile.md")).read())


if __name__ == "__main__":
    main()
