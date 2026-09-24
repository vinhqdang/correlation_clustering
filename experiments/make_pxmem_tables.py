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


ORDER15 = ["ca-GrQc", "BitcoinAlpha+", "BitcoinOTC+", "ca-HepTh", "ca-HepPh", "ca-AstroPh",
           "ca-CondMat", "email-Enron", "loc-Brightkite", "Slashdot+", "soc-Epinions",
           "Epinions+", "com-DBLP", "com-Amazon", "com-Youtube"]


def pval(p):
    p = float(p)
    return ("\\textbf{%.3f}" % p) if p < 0.05 else ("%.2f" % p if p >= 0.01 else "%.3f" % p)


# head-to-head, sequential protocol (experiments/seq_stats.py s1)
h = read(os.path.join(R, "headtohead_s1.csv"))
with open(os.path.join(T, "headtohead.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrrrr}\n\\toprule\n")
    fh.write("graph & $n$ & PXMem & KaPoCE & $\\bar\\Delta$ [95\\% CI] & W/T/L & $p$ & "
             "$p_{\\mathrm{Holm}}$ & $t_{=}$ (s) \\\\\n\\midrule\n")
    for g in ORDER15:
        r = h[g]
        d = float(r["diff_mean"])
        ci = f"$[{float(r['ci_lo']):+.1f}, {float(r['ci_hi']):+.1f}]$" if d != 0 else ""
        ttt = float(r["ttt"])
        ts = "--" if ttt == float("inf") else f"{ttt:.0f}"
        fh.write(f"{tt(g)} & \\num{{{r['n']}}} & \\num{{{float(r['ours_mean']):.1f}}} & "
                 f"\\num{{{float(r['kapoce_mean']):.1f}}} & ${d:+.1f}$ {ci} & "
                 f"{r['wins']}/{r['ties']}/{r['losses']} & {pval(r['sign_p'])} & "
                 f"{pval(r['holm_p'])} & {ts} \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# budget dependence
b60 = read(os.path.join(R, "headtohead_m3t60.csv"))
b150 = read(os.path.join(R, "headtohead_m3t150.csv"))
b600 = read(os.path.join(R, "headtohead_m3.csv"))
with open(os.path.join(T, "budget.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrr}\n\\toprule\n")
    fh.write("graph & \\SI{60}{s} & \\SI{150}{s} & \\SI{600}{s} \\\\\n\\midrule\n")
    for g in ORDER:
        cells = []
        for d in (b60, b150, b600):
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
