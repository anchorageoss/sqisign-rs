"""Fig 2: signature bytes against verification megacycles for the three
round-3 formats of one crate at level I (standard, compressed, compact).
Source: paper-dim4/BENCH.md, the paper's session (sizes table; verify rows)."""
from session import *

pts = [
    ("standard", SIZES[(R3, 3, 2, "I")][1], b(R3, 3, 2, "I", "verify from bytes")["median"], "rs_dim2"),
    ("compressed", SIZES[(R3C, 3, 2, "I")][1], b(R3, 3, 2, "I", "verify compressed from bytes")["median"], "rs_dim2"),
    ("compact (dimension 4)", SIZES[(R3, 3, 4, "I")][1], b(R3, 3, 4, "I", "verify from bytes")["median"], "rs_dim4"),
]
write_csv("fig02_tradeoff.csv", ["format", "bytes", "mcyc"], [[n, by, f"{m:.1f}"] for n, by, m, _ in pts])
# one file per point so main.tex can anchor each label on its own side
for fname, (n, by, m, _) in zip(("fig02_standard.csv", "fig02_compressed.csv", "fig02_compact.csv"), pts):
    write_csv(fname, ["format", "bytes", "mcyc"], [[n, by, f"{m:.1f}"]])
fig, ax = plt.subplots(figsize=(4.4, 3.0))
for n, by, m, c in pts:
    ax.plot([by], [m], "o", color=COLORS[c], ms=7)
    ax.annotate(n, (by, m), textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel("signature bytes")
ax.set_ylabel("verification, Mcycles (median)")
ax.set_xlim(120, 220)
ax.set_ylim(0, max(p[2] for p in pts) * 1.15)
fig.tight_layout()
save(fig, "fig02_tradeoff")
