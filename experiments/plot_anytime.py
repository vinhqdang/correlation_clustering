"""Figure: certified gap of the lower bound over time (cold vs warm start)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BEST_UB = {"ca-HepTh": 15814, "com-DBLP": 606752}
COLORS = {"block LP (cold start)": "#eb6834", "packing + block LP + gap blocks": "#2a78d6"}
LABELS = {"block LP (cold start)": "block LP, cold start",
          "packing + block LP + gap blocks": "packing warm start + block LP"}

files = sys.argv[1:] or [os.path.join(ROOT, "results", f) for f in
                         ("anytime_hepth.csv", "anytime_dblp.csv")]
d = pd.concat([pd.read_csv(f) for f in files])
graphs = list(dict.fromkeys(d.graph))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, len(graphs), figsize=(3.3 * len(graphs), 2.9))
for ax, gname in zip(axes, graphs):
    ub = BEST_UB[gname]
    for var, s in d[d.graph == gname].groupby("variant", sort=False):
        s = s.sort_values("time")
        s = s[s.lb > 0]
        gap = 100 * (ub - s.lb) / s.lb
        ax.step(s.time, gap, where="post", color=COLORS[var], lw=2, label=LABELS[var])
        ax.annotate(f"{gap.iloc[-1]:.1f}%", (s.time.iloc[-1], gap.iloc[-1]),
                    textcoords="offset points", xytext=(3, 0), va="center", fontsize=8,
                    color="#333333")
    ax.set_yscale("log")
    ax.set_title(gname, fontsize=9)
    ax.set_xlabel("time (s)")
    ax.grid(True, which="major", color="#e5e5e5", lw=0.6)
    ax.set_xlim(0, 700)
axes[0].set_ylabel("certified gap (%)")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=2)
fig.tight_layout(rect=(0, 0.08, 1, 1))
out = os.path.join(ROOT, "paper", "figures", "anytime.pdf")
fig.savefig(out)
fig.savefig(out.replace(".pdf", ".png"), dpi=150)
print(out)
