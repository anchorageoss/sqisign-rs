"""The paper's session in paper-dim4/BENCH.md, parsed once for every script.

Section "## Session 2026-09-26 (the dimension-4 paper)" of paper-dim4/BENCH.md: the
sizes table, the Rust rows, the C reference rows, the rejection-loop row,
the memory row and, when the author has run it, the M4 Pro table.
"""
from common import *

SESSION_HEADING = r"^## Session 2026-09-26"
bench = read("paper-dim4/BENCH.md")
sess = section(bench, SESSION_HEADING)

R3 = "SQIsign (Rust, this crate)"
R3C = "SQIsign compressed (Rust, this crate)"
R2 = "SQIsign (Rust, v0.4.31)"
C3 = "SQIsign (C reference, nist-v3, broadwell)"
C2 = "SQIsign (C reference, nist-v2, broadwell)"

_, _size_rows = md_table(section(sess, r"^#### Sizes observed"), r"^implementation$")
SIZES = {(r[0], int(r[1]), int(r[2]), r[3]): (int(r[4]), int(r[5])) for r in _size_rows}


def _rows(sec_regex):
    _, rows = md_table(section(sess, sec_regex), r"^implementation$")
    out = {}
    for r in rows:
        key = (r[0], int(r[1]), int(r[2]), r[3], r[4])
        ms = num(r[10]) if r[10] else None
        out[key] = dict(runs=int(r[5]), median=num(r[6]), mean=num(r[7]), min=num(r[8]), max=num(r[9]), ms=ms)
    return out


BENCH = _rows(r"^#### Key generation, signing, verification")
BENCH.update(_rows(r"^#### The C reference"))

_, _loop = md_table(section(sess, r"^#### The compact signer's rejection loop"), r"^signatures$")
LOOP = dict(signatures=int(_loop[0][0]), samples_mean=num(_loop[0][1]), samples_median=num(_loop[0][2]),
            samples_min=int(_loop[0][3]), samples_max=int(_loop[0][4]), tests_mean=num(_loop[0][5]))

_, _mem = md_table(section(sess, r"^#### Memory of one compact verification"), r"^peak heap")
MEM = dict(peak_heap_bytes=int(_mem[0][0]), heap_before=int(_mem[0][1]), vm_hwm_kib=_mem[0][2])

try:
    _, _m4 = md_table(section(sess, r"^#### M4 Pro"), r"^implementation$")
    M4 = {(r[0], int(r[1]), int(r[2]), r[3], r[4]): num(r[6]) for r in _m4 if r[6] and r[6] != "pending"}
except KeyError:
    M4 = {}

MACHINE = grab(sess, r"Machine: ([^;\n]+)")
XEON_GHZ = grab(sess, r"@ ([0-9.]+)GHz")


def b(impl, rnd, dim, level, op):
    return BENCH[(impl, rnd, dim, level, op)]
