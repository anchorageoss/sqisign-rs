"""The macros main.tex quotes in prose and captions: sizes, parameters,
costs and multiples.  Sources: paper-dim4/BENCH.md (the paper's session),
docs/COMPACT_R3.md Sections 1 to 3 (parameters, predicted sizes, expected
rejection-loop counts)."""
from session import *

c3 = read("docs/COMPACT_R3.md")
_, prow = md_table(section(c3, r"^## 1\."), r"^level$")
P = {r[0]: r for r in prow}
_, srow = md_table(section(c3, r"^## 3\."), r"^level$")
SZ = {r[0]: r for r in srow}
_, lrow = md_table(section(c3, r"^## 2\."), r"^level$")
LP = {r[0]: r for r in lrow}
lam = {"I": 128, "III": 192, "V": 256}
pk_pred = grab(c3, r"convention \(\*\*theirs\*\*\), (\d+) / (\d+) / (\d+) bytes")
pk_pred = re.search(r"(\d+) / (\d+) / (\d+) bytes", grab(c3, r"(hints in the library's\s+convention \(\*\*theirs\*\*\), \d+ / \d+ / \d+ bytes)"))

pairs = []
add = lambda k, v: pairs.append((k, v))
lv = {"I": "I", "III": "III", "V": "V"}
for L in lv:
    add(f"Embed{L}", num(P[L][5].replace("**", "")))
    add(f"TorsionR{L}", num(P[L][6].replace("**", "")))
    add(f"HalfSteps{L}", P[L][7])
    add(f"LogTwoP{L}", P[L][2])
    add(f"PredSig{L}", SZ[L][5].replace("**", ""))
    add(f"PredSigAll{L}", SZ[L][7])
    add(f"SpecSig{L}", SZ[L][8])
    add(f"Lambda{L}", lam[L])
    add(f"ExpectedSamples{L}", LP[L][3])
    add(f"ExpectedTests{L}", LP[L][4])
    import math
    add(f"ExpectedSamplesOdd{L}", f"{int(LP[L][1]) * math.log(2):.0f}")
for L, g in zip(lv, (1, 2, 3)):
    add(f"PredPk{L}", pk_pred.group(g))
pairs = [(k, (str(int(v)) if isinstance(v, float) and v == int(v) else str(v))) for k, v in pairs]
add = lambda k, v: pairs.append((k, str(v)))

# sizes observed in the session
for L in lv:
    pk, sig = SIZES[(R3, 3, 2, L)]
    add(f"PkStd{L}", pk); add(f"SigStd{L}", sig); add(f"ComboStd{L}", pk + sig)
    add(f"SigComp{L}", SIZES[(R3C, 3, 2, L)][1])
    pk2, sig2 = SIZES[(R2, 2, 2, L)]
    add(f"PkRtwoStd{L}", pk2); add(f"SigRtwoStd{L}", sig2)
pk, sig = SIZES[(R3, 3, 4, "I")]
add("PkCompactI", pk); add("SigCompactI", sig); add("ComboCompactI", pk + sig)
pk, sig = SIZES[(R2, 2, 4, "I")]
add("PkRtwoCompactI", pk); add("SigRtwoCompactI", sig); add("ComboRtwoCompactI", pk + sig)
add("SigSavingPct", round(100 * (1 - SIZES[(R3, 3, 4, "I")][1] / SIZES[(R3, 3, 2, "I")][1])))

# costs, megacycles (medians) and the multiples
def M(impl, rnd, dim, L, op, key="median"):
    return b(impl, rnd, dim, L, op)[key]
f1 = lambda x: f"{x:.1f}"
for L in lv:
    add(f"VerStd{L}", f1(M(R3, 3, 2, L, "verify from bytes")))
    add(f"VerComp{L}", f1(M(R3, 3, 2, L, "verify compressed from bytes")))
    add(f"SignStd{L}", f1(M(R3, 3, 2, L, "sign")))
    add(f"KeygenStd{L}", f1(M(R3, 3, 2, L, "keygen")))
    add(f"VerCref{L}", f1(M(C3, 3, 2, L, "verify from bytes")))
    add(f"SignCref{L}", f1(M(C3, 3, 2, L, "sign")))
    add(f"KeygenCref{L}", f1(M(C3, 3, 2, L, "keygen")))
    add(f"VerRtwoStd{L}", f1(M(R2, 2, 2, L, "verify from bytes")))
    add(f"SignRtwoStd{L}", f1(M(R2, 2, 2, L, "sign")))
    add(f"KeygenRtwoStd{L}", f1(M(R2, 2, 2, L, "keygen")))
    if (C2, 2, 2, L, "verify from bytes") in BENCH:
        add(f"VerCrefRtwo{L}", f1(M(C2, 2, 2, L, "verify from bytes")))
        add(f"SignCrefRtwo{L}", f1(M(C2, 2, 2, L, "sign")))
        add(f"KeygenCrefRtwo{L}", f1(M(C2, 2, 2, L, "keygen")))
    add(f"VerStdOverCref{L}", f"{M(R3, 3, 2, L, 'verify from bytes') / M(C3, 3, 2, L, 'verify from bytes'):.2f}")
    add(f"KeygenStdOverCref{L}", f"{M(R3, 3, 2, L, 'keygen') / M(C3, 3, 2, L, 'keygen'):.2f}")
    add(f"SignStdOverCref{L}", f"{M(R3, 3, 2, L, 'sign') / M(C3, 3, 2, L, 'sign'):.2f}")
v4 = M(R3, 3, 4, "I", "verify from bytes"); v2 = M(R3, 3, 2, "I", "verify from bytes")
add("VerCompactI", f1(v4)); add("VerMultI", f1(v4 / v2))
add("VerMultCrefI", f1(v4 / M(C3, 3, 2, "I", "verify from bytes")))
s4 = M(R3, 3, 4, "I", "sign"); s2 = M(R3, 3, 2, "I", "sign")
add("SignCompactI", f1(s4)); add("SignMultI", f1(s4 / s2))
add("SignCompactMeanI", f1(M(R3, 3, 4, "I", "sign", "mean")))
add("SignCompactMinI", f1(M(R3, 3, 4, "I", "sign", "min"))); add("SignCompactMaxI", f1(M(R3, 3, 4, "I", "sign", "max")))
k4 = M(R3, 3, 4, "I", "keygen"); k2 = M(R3, 3, 2, "I", "keygen")
add("KeygenCompactI", f1(k4)); add("KeygenMultI", f"{k4 / k2:.2f}")
rv4 = M(R2, 2, 4, "I", "verify from bytes"); rv2 = M(R2, 2, 2, "I", "verify from bytes")
add("VerRtwoCompactI", f1(rv4)); add("VerRtwoMultI", f1(rv4 / rv2))
add("SignRtwoCompactI", f1(M(R2, 2, 4, "I", "sign"))); add("KeygenRtwoCompactI", f1(M(R2, 2, 4, "I", "keygen")))
add("VerCompOverStdPct", round(100 * (M(R3, 3, 2, "I", "verify compressed from bytes") / v2 - 1)))
# milliseconds on the bench box
for (name, impl, rnd, dim, op) in [("MsVerStdI", R3, 3, 2, "verify from bytes"), ("MsVerCompactI", R3, 3, 4, "verify from bytes"),
                                   ("MsSignStdI", R3, 3, 2, "sign"), ("MsSignCompactI", R3, 3, 4, "sign")]:
    add(name, f"{M(impl, rnd, dim, 'I', op, 'ms'):.1f}")
# the rejection loop and memory
add("LoopSignatures", LOOP["signatures"]); add("LoopSamplesMean", f"{LOOP['samples_mean']:.0f}")
add("LoopSamplesMedian", f"{LOOP['samples_median']:.0f}"); add("LoopSamplesMin", LOOP["samples_min"]); add("LoopSamplesMax", LOOP["samples_max"])
add("LoopTestsMean", f"{LOOP['tests_mean']:.0f}")
add("PeakHeapKiB", round(MEM["peak_heap_bytes"] / 1024)); add("PeakHeapMiB", f"{MEM['peak_heap_bytes'] / 2**20:.2f}")
add("MachineXeon", MACHINE.replace("(R)", "")); add("XeonGHz", XEON_GHZ)
add("MFourAvailable", "1" if M4 else "0")
write_macros("values.tex", pairs, "values.py")
