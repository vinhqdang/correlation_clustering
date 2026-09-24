"""LaTeX tables for the PXMem experiments (head-to-head against KaPoCE, budget
dependence, parallel runs combined by partition crossover)."""
import csv, json, glob, os
from collections import defaultdict

R = os.path.join(os.path.dirname(__file__), "..", "results")
T = os.path.join(os.path.dirname(__file__), "..", "paper", "tables")
NAMES = {"BitcoinAlpha+": "Bitcoin-Alpha$^+$", "BitcoinOTC+": "Bitcoin-OTC$^+$"}
ORDER = ["ca-GrQc", "BitcoinAlpha+", "BitcoinOTC+", "ca-HepTh", "ca-HepPh", "ca-AstroPh",
         "ca-CondMat", "email-Enron", "loc-Brightkite", "com-DBLP", "com-Amazon"]


def tt(g):
    return NAMES.get(g, r"\texttt{%s}" % g)


def read(path):
    return {r["graph"]: r for r in csv.DictReader(open(path))}


# per-graph mean costs from the raw runs
runs = defaultdict(list)
for f in glob.glob(os.path.join(R, "colab", "m3_*.json")):
    r = json.load(open(f))
    if r["tag"] == "m3":
        runs[r["graph"]].append(r)

h = read(os.path.join(R, "headtohead_m3.csv"))
with open(os.path.join(T, "headtohead.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrcr}\n\\toprule\n")
    fh.write("graph & $n$ & PXMem & KaPoCE & $\\Delta$ (\\%) & W/T/L & unfin. & $p$ \\\\\n\\midrule\n")
    for g in ORDER:
        r = h[g]
        rs = runs[g]
        ok = [x for x in rs if x["kapoce"] <= 1.5 * x["ours"]]
        ours = sum(x["ours"] for x in rs) / len(rs)
        kap = sum(x["kapoce"] for x in ok) / len(ok)
        p = float(r["sign_p"])
        ps = ("\\textbf{%.3f}" % p) if p < 0.05 else "%.2f" % p
        fh.write(f"{tt(g)} & \\num{{{rs[0]['n']}}} & \\num{{{ours:.1f}}} & \\num{{{kap:.1f}}} & "
                 f"${float(r['mean_rel_diff_pct_valid_runs']):+.4f}$ & {r['wins']}/{r['ties']}/{r['losses']} & "
                 f"{r['kapoce_unfinished'] if r['kapoce_unfinished'] != '0' else '--'} & {ps} \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# budget dependence
b60 = read(os.path.join(R, "headtohead_m3t60.csv"))
b150 = read(os.path.join(R, "headtohead_m3t150.csv"))
with open(os.path.join(T, "budget.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrr}\n\\toprule\n")
    fh.write("graph & \\SI{60}{s} & \\SI{150}{s} & \\SI{600}{s} \\\\\n\\midrule\n")
    for g in ORDER:
        cells = []
        for d in (b60, b150, h):
            r = d[g]
            cells.append(f"${float(r['mean_rel_diff_pct_valid_runs']):+.3f}$")
        fh.write(f"{tt(g)} & " + " & ".join(cells) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# parallel runs and PX post-processing
par = defaultdict(dict)
for r in csv.DictReader(open(os.path.join(R, "parallel_px.csv"))):
    par[r["graph"]][(r["k"], r["method"])] = float(r["cost"])
cols = [("1", "best-of-k"), ("2", "best-of-k"), ("2", "PX-of-k"), ("4", "best-of-k"),
        ("4", "PX-of-k"), ("1", "KaPoCE"), ("1", "PX(KaPoCE, ours)")]
with open(os.path.join(T, "parallel.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrrr}\n\\toprule\n")
    fh.write(" & 1 run & \\multicolumn{2}{c}{2 runs} & \\multicolumn{2}{c}{4 runs} & & PX(KaPoCE, \\\\\n")
    fh.write("\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n")
    fh.write("graph & & best & PX & best & PX & KaPoCE & PXMem) \\\\\n\\midrule\n")
    for g in [x for x in ORDER if x in par]:
        vals = [par[g][c] for c in cols]
        best = min(vals)
        cells = []
        for v in vals:
            s = f"\\num{{{v:.1f}}}" if v != int(v) else f"\\num{{{int(v)}}}"
            cells.append(f"\\textbf{{{s}}}" if v == best else s)
        fh.write(f"{tt(g)} & " + " & ".join(cells) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
print("tables written")
