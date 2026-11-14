#!/usr/bin/env python3
"""TABLES.md in milliseconds from a session's numbers.json, in the row
layout of the paper's Table 3: for machines without a cycle column (Apple
Silicon). The C reference rows read "not run" when the session has none.

    paper_tables_ms.py <session-dir> > TABLES.md
"""
import json
import os
import sys

N = json.load(open(os.path.join(sys.argv[1], "numbers.json")))
R = N["rows"]


def ms(name):
    r = R.get(name)
    return None if r is None or r["ms"] is None else r["ms"]


def cell(x):
    return "not run" if x is None else f"{x:.2f}"


print("# Wall-clock milliseconds, medians (the paper's Table 3 layout)")
print()
print(f"Passes: {N['passes']}; session value of a row: {N['rule']}.")
print("Compact verification at levels III and V: the mean over six signatures of the median of 15 runs.")
print()
print("| implementation | level | keygen | sign | verify | x dim. 2 |")
print("|---|---|---:|---:|---:|---:|")
for lvl in ("I", "III", "V"):
    print(f"| C reference, dimension 2 (portable build) | {lvl} | {cell(ms(f'c-{lvl}-keygen'))} | {cell(ms(f'c-{lvl}-sign'))} | {cell(ms(f'c-{lvl}-verify'))} | |")
for lvl in ("I", "III", "V"):
    print(f"| Rust, dimension 2 | {lvl} | {cell(ms(f'dim2-{lvl}-keygen'))} | {cell(ms(f'dim2-{lvl}-sign'))} | {cell(ms(f'dim2-{lvl}-verify'))} | |")
for lvl in ("I", "III", "V"):
    print(f"| Rust, dimension 2, compressed verify | {lvl} | | | {cell(ms(f'dim2-{lvl}-verify-compressed'))} | |")
for lvl in ("I", "III", "V"):
    d2 = ms(f"dim2-{lvl}-verify")
    four, two, nc = ms(f"dim4-{lvl}-verify-four"), ms(f"dim4-{lvl}-verify-two"), ms(f"dim4-{lvl}-verify-two-nocheck")
    print(f"| Rust, dimension 4 (compact) | {lvl} | {cell(ms(f'dim4-{lvl}-keygen'))} | {cell(ms(f'dim4-{lvl}-sign'))} | {cell(four)} | {four / d2:.1f} |")
    print(f"| two points | {lvl} | | | {cell(two)} | {two / d2:.1f} |")
    print(f"| two points, no check at the gluing | {lvl} | | | {cell(nc)} | {nc / d2:.1f} |")
print()
print("## Two-point over four-point verification")
print()
print("| level | two / four | ratio |")
print("|---|---|---:|")
for lvl in ("I", "III", "V"):
    four, two = ms(f"dim4-{lvl}-verify-four"), ms(f"dim4-{lvl}-verify-two")
    print(f"| {lvl} | {two:.2f} / {four:.2f} | {two / four:.2f} |")
L = N["loop"]
print()
print("## The good-degree loop at level I")
print()
print("| | value |")
print("|---|---|")
print(f"| signatures measured | {L['signatures']} |")
print(f"| response samples per signature, mean | {L['samples_mean']:.0f} |")
print(f"| response samples per signature, median (min, max) | {L['samples_median']:.0f} ({L['samples_min']}, {L['samples_max']}) |")
print(f"| primality tests per signature, mean | {L['primality_tests_mean']:.0f} |")
print(f"| signing, ms, median (min, max) | {L['sign_ms_median']:.1f} ({L['sign_ms_min']:.1f}, {L['sign_ms_max']:.1f}) |")
print(f"| commitments abandoned for a fresh one | {L['restarts']} |")
if N.get("heap"):
    print()
    print("## Heap peak of one verification, KiB (counting allocator around the call)")
    print()
    print("| row | mean over 5 signatures (max) |")
    print("|---|---|")
    for k, v in sorted(N["heap"].items()):
        print(f"| {k} | {v['peak_kib_mean']:.0f} ({v['peak_kib_max']:.0f}) |")
