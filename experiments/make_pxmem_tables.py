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

# budget dependence, sequential protocol (tags s1t60, s1t150, s1): per graph the
# median paired relative difference (%) and wins/ties/losses of PXMem
def runs(tag):
    by = defaultdict(list)
    for f in glob.glob(os.path.join(R, "colab", f"{tag}_*.json")):
        r = json.load(open(f))
        if r.get("tag") == tag:
            by[r["graph"]].append(r)
    return by


def median(xs):
    xs = sorted(xs)
    k = len(xs)
    return (xs[k // 2] + xs[(k - 1) // 2]) / 2


budgets = [runs(t) for t in ("s1t60", "s1t150", "s1")]
with open(os.path.join(T, "budget.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
    fh.write(" & \\multicolumn{2}{c}{\\SI{60}{s}} & \\multicolumn{2}{c}{\\SI{150}{s}} & "
             "\\multicolumn{2}{c}{\\SI{600}{s}} \\\\\n")
    fh.write("\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n")
    fh.write("graph & med.\\ (\\%) & W/T/L & med.\\ (\\%) & W/T/L & med.\\ (\\%) & W/T/L "
             "\\\\\n\\midrule\n")
    for g in ORDER15:
        cells = []
        for by in budgets:
            rs = by.get(g, [])
            if not rs:
                cells += ["--", "--"]
                continue
            # a run without a valid KaPoCE solution counts as a win (shown in W)
            valid = [r for r in rs if r["kapoce_valid"]]
            d = [100 * (r["ours"] - r["kapoce"]) / r["kapoce"] for r in valid]
            w = sum(x < 0 for x in d) + len(rs) - len(valid)
            l = sum(x > 0 for x in d)
            m = round(median(d), 3) + 0.0 if d else None   # + 0.0 turns -0.0 into 0.0
            cells += [f"${m:+.3f}$" if d else "--", f"{w}/{len(rs) - w - l}/{l}"]
        fh.write(f"{tt(g)} & " + " & ".join(cells) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# parallel runs and PX post-processing
par = defaultdict(dict)
for r in csv.DictReader(open(os.path.join(R, "parallel_px.csv"))):
    par[r["graph"]][(r["k"], r["method"])] = float(r["cost"])
cols = [("1", "best-of-k"), ("2", "best-of-k"), ("2", "PX-of-k"), ("4", "best-of-k"),
        ("4", "PX-of-k"), ("1", "KaPoCE"), ("1", "PX(KaPoCE, ours)"),
        ("4", "PX(KaPoCE, 4 x ours)")]
with open(os.path.join(T, "parallel.tex"), "w") as fh:
    fh.write("\\begin{tabular}{lrrrrrrrr}\n\\toprule\n")
    fh.write(" & 1 run & \\multicolumn{2}{c}{2 runs} & \\multicolumn{2}{c}{4 runs} & & "
             "\\multicolumn{2}{c}{PX with KaPoCE} \\\\\n")
    fh.write("\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){8-9}\n")
    fh.write("graph & & best & PX & best & PX & KaPoCE & 1 run & 4 runs \\\\\n\\midrule\n")
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
