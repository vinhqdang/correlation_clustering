"""Turn benchmark CSVs into LaTeX tables and summary statistics."""
from __future__ import annotations

import argparse
import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
TAB = os.path.join(ROOT, "paper", "tables")

NAMES = {"pivot": "Pivot", "pivot50": "Pivot-50", "vote": "Vote", "louvain": "Louvain-CC",
         "leiden": "Leiden-CPM", "mfp": "MatchFlipPivot", "ins": "Insertion LS",
         "flip": "IteratedFlip", "kapoce": "KaPoCE", "certiflip": "CertiFlip"}
ORDER = ["pivot", "pivot50", "mfp", "vote", "louvain", "leiden", "ins", "flip", "kapoce",
         "certiflip"]


def load(path):
    d = pd.read_csv(path)
    if "error" in d:
        d = d[d["error"].isna()]
    d["cost"] = pd.to_numeric(d["cost"], errors="coerce")
    d["lb"] = pd.to_numeric(d["lb"], errors="coerce")
    d["time"] = pd.to_numeric(d["time"], errors="coerce")
    return d


def pace_exact(path):
    d = load(path)
    d["inst"] = d["instance"].str.split("/").str[-1]
    ref = pd.read_csv(os.path.join(ROOT, "data", "pace2021_exact_kapoce_bounds.csv"))
    ref = ref.rename(columns={"instance": "inst"})
    ref["opt"] = np.where(ref["solved_1h"] == 1, ref["opt_or_empty"], np.nan)
    ref["pub_lb"] = np.where(ref["solved_1h"] == 1, ref["opt_or_empty"], ref["low_star"])
    # published per-instance lower bounds from the branch-and-bound log:
    # low_star / low_p3 are its root packing bounds (see data README)
    prim = d[~d["algo"].str.startswith("lb:")]
    best_ours = prim.groupby("inst")["cost"].min()
    ref = ref.set_index("inst")
    ref["best_found"] = best_ours
    ref["best_known"] = np.fmin(ref["upper"], ref["best_found"])
    rows = []
    for a in ORDER:
        s = prim[prim.algo == a]
        if s.empty:
            continue
        g = s.groupby("inst").agg(cost=("cost", "mean"), best=("cost", "min"), time=("time", "mean"))
        g = g.join(ref[["best_known", "opt"]])
        gap = 100 * (g["cost"] - g["best_known"]) / g["best_known"].clip(lower=1)
        rows.append(dict(algo=NAMES[a], n_inst=len(g),
                         hit=int((g["best"] <= g["best_known"]).sum()),
                         mean_gap=gap.mean(), max_gap=gap.max(), time=g["time"].mean()))
    prim_tab = pd.DataFrame(rows)
    # lower bounds
    lbs = {}
    cf = d[d.algo == "certiflip"].set_index("inst")
    lbs["CertiFlip bound"] = np.ceil(cf["lb"] - 1e-6)
    for a, nm in [("lb:tri-mwu", "Triangle packing (MWU)"), ("lb:greedy", "Greedy packing")]:
        s = d[d.algo == a].set_index("inst")
        lbs[nm] = np.ceil(s["lb"] - 1e-6)
    lbs["B\\&B star packing"] = ref["low_star"]
    lbs["B\\&B $P_3$ packing"] = ref["low_p3"]
    solved = ref[ref["solved_1h"] == 1]
    lrows = []
    for nm, v in lbs.items():
        v = v.reindex(solved.index)
        ok = v.notna()
        ratio = (v[ok] / solved["opt"][ok].clip(lower=1))
        lrows.append(dict(bound=nm, n=int(ok.sum()), mean_ratio=ratio.mean(),
                          min_ratio=ratio.min(), tight=int((v[ok] >= solved["opt"][ok]).sum())))
    lb_tab = pd.DataFrame(lrows)
    # open instances
    openi = ref[ref["solved_1h"] == 0].copy()
    openi["our_ub"] = best_ours.reindex(openi.index)
    openi["our_lb"] = lbs["CertiFlip bound"].reindex(openi.index)
    return prim_tab, lb_tab, openi, ref, d


def fmt(x, nd=2):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    return f"{x:.{nd}f}"


def write_pace_tables(path):
    prim, lbt, openi, ref, d = pace_exact(path)
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, "pace_primal.tex"), "w") as fh:
        fh.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
        fh.write("Algorithm & best known hit & mean gap (\\%) & max gap (\\%) & time (s)\\\\\n\\midrule\n")
        for _, r in prim.iterrows():
            fh.write(f"{r.algo} & {r.hit}/{r.n_inst} & {fmt(r.mean_gap, 3)} & {fmt(r.max_gap, 2)} & {fmt(r.time, 2)}\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    with open(os.path.join(TAB, "pace_bounds.tex"), "w") as fh:
        fh.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
        fh.write("Lower bound & instances & mean LB/OPT & min LB/OPT & LB $=$ OPT\\\\\n\\midrule\n")
        for _, r in lbt.iterrows():
            fh.write(f"{r.bound} & {r.n} & {fmt(r.mean_ratio, 4)} & {fmt(r.min_ratio, 3)} & {r.tight}\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    with open(os.path.join(TAB, "pace_open.tex"), "w") as fh:
        fh.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        fh.write("Instance & $n$ & $m$ & published UB & published LB & our UB & our LB\\\\\n\\midrule\n")
        for inst, r in openi.iterrows():
            nm = inst.replace(".gr", "")
            ub = int(r.our_ub) if not math.isnan(r.our_ub) else "--"
            lb = int(r.our_lb) if not math.isnan(r.our_lb) else "--"
            ubs = f"\\textbf{{{ub}}}" if ub != "--" and ub < r.upper else f"{ub}"
            lbs = f"\\textbf{{{lb}}}" if lb != "--" and lb > r.low_star else f"{lb}"
            fh.write(f"{nm} & {r.n} & {r.m} & {int(r.upper)} & {int(r.low_star)} & {ubs} & {lbs}\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(prim.to_string())
    print(lbt.to_string())
    print(openi[["n", "m", "upper", "low_star", "our_ub", "our_lb"]].to_string())


def write_snap_tables(path):
    d = load(path)
    os.makedirs(TAB, exist_ok=True)
    prim = d[~d["algo"].str.startswith("lb:")]
    algos = [a for a in ORDER if a in set(prim.algo)]
    insts = list(dict.fromkeys(d["instance"]))
    with open(os.path.join(TAB, "snap_primal.tex"), "w") as fh:
        fh.write("\\begin{tabular}{l" + "r" * (len(algos) + 2) + "}\n\\toprule\n")
        fh.write("Graph & $n$ & $m$ & " + " & ".join(NAMES[a] for a in algos) + "\\\\\n\\midrule\n")
        for inst in insts:
            s = prim[prim.instance == inst]
            if s.empty:
                continue
            n, m = int(s.n.iloc[0]), int(s.m.iloc[0])
            vals = {a: s[s.algo == a]["cost"].min() for a in algos}
            best = np.nanmin(list(vals.values()))
            cells = []
            for a in algos:
                v = vals[a]
                if v is None or math.isnan(v):
                    cells.append("--")
                elif v == best:
                    cells.append(f"\\textbf{{{int(v)}}}")
                else:
                    cells.append(f"{int(v)}")
            fh.write(f"\\texttt{{{inst}}} & {n} & {m} & " + " & ".join(cells) + "\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    with open(os.path.join(TAB, "snap_bounds.tex"), "w") as fh:
        fh.write("\\begin{tabular}{lrrrrr}\n\\toprule\n")
        fh.write("Graph & best UB & greedy packing & triangle packing & CertiFlip LB & certified gap (\\%)\\\\\n\\midrule\n")
        for inst in insts:
            s = d[d.instance == inst]
            ub = s[~s.algo.str.startswith("lb:")]["cost"].min()
            g = s[s.algo == "lb:greedy"]["lb"].max()
            t = s[s.algo == "lb:tri-mwu"]["lb"].max()
            c = s[s.algo == "certiflip"]["lb"].max()
            gap = 100 * (ub - math.ceil(c - 1e-6)) / max(1, math.ceil(c - 1e-6)) if not math.isnan(c) else float("nan")
            fh.write(f"\\texttt{{{inst}}} & {int(ub) if not math.isnan(ub) else '--'} & "
                     f"{int(g) if not math.isnan(g) else '--'} & {int(math.ceil(t-1e-6)) if not math.isnan(t) else '--'} & "
                     f"{int(math.ceil(c-1e-6)) if not math.isnan(c) else '--'} & {fmt(gap, 2)}\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    with open(os.path.join(TAB, "snap_time.tex"), "w") as fh:
        fh.write("\\begin{tabular}{l" + "r" * len(algos) + "}\n\\toprule\n")
        fh.write("Graph & " + " & ".join(NAMES[a] for a in algos) + "\\\\\n\\midrule\n")
        for inst in insts:
            s = prim[prim.instance == inst]
            if s.empty:
                continue
            cells = [fmt(s[s.algo == a]["time"].mean(), 1) for a in algos]
            fh.write(f"\\texttt{{{inst}}} & " + " & ".join(cells) + "\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")


def write_ablation(path):
    d = pd.read_csv(path)
    d = d[d.method != "support"]
    graphs = list(dict.fromkeys(d.graph))
    methods = list(dict.fromkeys(d.method))
    best_ub = {}
    snap = os.path.join(ROOT, "results", "snap.csv")
    if os.path.exists(snap):
        s2 = load(snap)
        s2 = s2[~s2.algo.str.startswith("lb:")]
        for gname in graphs:
            v = s2[s2.instance == gname]["cost"].min()
            if not math.isnan(v):
                best_ub[gname] = v
    with open(os.path.join(TAB, "ablation_lb.tex"), "w") as fh:
        fh.write("\\begin{tabular}{l" + "rr" * len(graphs) + "}\n\\toprule\n")
        fh.write(" & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{\\texttt{{{gname}}}}}" for gname in graphs) + "\\\\\n")
        fh.write("Bound & " + " & ".join("gap (\\%) & s" for _ in graphs) + "\\\\\n\\midrule\n")
        for mth in methods:
            cells = []
            for gname in graphs:
                r = d[(d.graph == gname) & (d.method == mth)]
                if r.empty:
                    cells += ["--", "--"]
                    continue
                ub = best_ub.get(gname, r.ub.iloc[0])
                lb = math.ceil(r.lb.iloc[0] - 1e-6)
                cells += [fmt(100 * (ub - lb) / max(lb, 1), 2), fmt(r.time.iloc[0], 0)]
            fh.write(mth + " & " + " & ".join(cells) + "\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(d.pivot_table(index="method", columns="graph", values="lb").to_string())


def write_scaling(path):
    d = pd.read_csv(path)
    tab = d.pivot_table(index="step", columns="n", values="time")
    with open(os.path.join(TAB, "scaling.tex"), "w") as fh:
        cols = list(tab.columns)
        fh.write("\\begin{tabular}{l" + "r" * len(cols) + "}\n\\toprule\n")
        fh.write("Step & " + " & ".join(f"$n={c:,}$".replace(",", "\\,") for c in cols) + "\\\\\n\\midrule\n")
        for step, row in tab.iterrows():
            fh.write(step + " & " + " & ".join(fmt(v, 1) for v in row.values) + "\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(tab.to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pace-exact")
    ap.add_argument("--snap")
    ap.add_argument("--ablation")
    ap.add_argument("--scaling")
    a = ap.parse_args()
    if a.ablation:
        write_ablation(a.ablation)
    if a.scaling:
        write_scaling(a.scaling)
    if a.pace_exact:
        write_pace_tables(a.pace_exact)
    if a.snap:
        write_snap_tables(a.snap)


if __name__ == "__main__":
    main()
