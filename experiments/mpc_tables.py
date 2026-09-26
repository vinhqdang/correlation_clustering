"""Tables and numbers of the MPC manuscript (paper_mpc/) from the raw results.

Every table is written to paper_mpc/tables/, every number used in the text to
paper_mpc/numbers.tex.  Results that are not yet available are printed as
``pending'' so that the manuscript always compiles.

usage: python experiments/mpc_tables.py [instances|all]"""
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments")]
RES = os.path.join(ROOT, "results")
COLAB = os.path.join(RES, "colab")
TAB = os.path.join(ROOT, "paper_mpc", "tables")
NUMBERS = os.path.join(ROOT, "paper_mpc", "numbers.tex")

DEV = ["ca-GrQc", "ca-HepTh", "ca-HepPh", "ca-AstroPh", "ca-CondMat", "email-Enron",
       "loc-Brightkite", "soc-Epinions", "Slashdot+", "Epinions+", "BitcoinOTC+",
       "BitcoinAlpha+", "com-Amazon", "com-DBLP", "com-Youtube"]
HELD = ["email-Eu-core", "facebook", "wiki-Vote", "cit-HepTh", "cit-HepPh", "p2p-Gnutella31",
        "soc-Slashdot0902", "loc-Gowalla", "email-EuAll", "amazon0302", "web-NotreDame",
        "roadNet-PA"]
# PACE 2021 heuristic-track instances with the same graph (experiments/heldout.md)
PACE_TWIN = {"ca-GrQc": 167, "ca-HepTh": 173, "ca-HepPh": 175, "ca-AstroPh": 176,
             "ca-CondMat": 178, "email-Enron": 179, "soc-Epinions": 180, "com-DBLP": 191,
             "com-Amazon": 193, "com-Youtube": 198, "email-Eu-core": 94, "facebook": 166,
             "soc-Slashdot0902": 186, "loc-Gowalla": 188, "amazon0302": 189,
             "web-NotreDame": 192, "roadNet-PA": 197}
MAX_PAIRS = 3e7

NUM = {}
TIMES = {}
# every number used in the text; results not yet available stay \pending
DEFAULTS = {k: r"\pending{}" for k in (
    "numPaceOptOurs", "numPaceMeanStar", "numPaceMeanOurs", "numSnapGapLo", "numSnapGapHi",
    "numBaseFacLo", "numBaseFacHi", "numCheckMaxGraph", "numCheckMaxRows", "numCheckMaxTime",
    "numCpuModels")}
DEFAULTS["numKapoceCommit"] = "64e2101"
DEFAULTS.update({k: r"\pending{}" for k in (
    "numPHn", "numPHBetter", "numPHEqual", "numPHWorse", "numPHWilcoxon", "numPHSign",
    "numPHKapoceInvalid", "numPHNoSolution", "numPDn", "numPDBetter", "numPDEqual",
    "numPDWorse")})
DEFAULTS.update({f"numPace{k}{f}": r"\pending{}" for k in ("Pack", "Warm", "Tri")
                 for f in ("N", "Mean", "Opt", "OursMean", "StarMean", "OursOpt")})


def tt(name):
    return r"\texttt{" + name.replace("_", r"\_") + "}"


def write(name, body):
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, name + ".tex"), "w") as fh:
        fh.write(body)


def instances():
    """n, m and |P| of every SNAP graph (|P| computed when the estimate
    m + sum_v C(deg v, 2) is at most MAX_PAIRS), cached in results/mpc."""
    cache_path = os.path.join(RES, "mpc", "instances.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    import datasets as D
    for name in DEV + HELD:
        if name in cache:
            continue
        g = D.load(name)
        deg = np.diff(g.indptr).astype(np.float64)
        est = g.m + float((deg * (deg - 1) / 2).sum())
        rec = {"n": int(g.n), "m": int(g.m), "est_pairs": est, "pairs": None}
        if est <= MAX_PAIRS:
            from ccbench.support import build_support
            rec["pairs"] = int(build_support(g).npairs)
        cache[name] = rec
        print(name, rec, flush=True)
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        json.dump(cache, open(cache_path, "w"), indent=1)
    rows = []
    for role, names in (("dev", DEV), ("held-out", HELD)):
        for name in names:
            r = cache[name]
            pairs = r"\num{%d}" % r["pairs"] if r["pairs"] is not None else "--"
            twin = "heur%03d" % PACE_TWIN[name] if name in PACE_TWIN else ""
            rows.append(f"{tt(name)} & \\num{{{r['n']}}} & \\num{{{r['m']}}} & {pairs} & "
                        f"{role} & {twin} \\\\")
        rows.append(r"\midrule")
    body = ("\\begin{tabular}{lrrrll}\n\\toprule\ngraph & $n$ & $m$ & $|P|$ & role & PACE \\\\\n"
            "\\midrule\n" + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
    write("instances", body)
    cert = [n for n in DEV + HELD if cache[n]["pairs"] is not None]
    NUM["numSnapCert"] = str(len(cert))
    NUM["numSnapCertHeld"] = str(len([n for n in cert if n in HELD]))
    big = max(cert, key=lambda n: cache[n]["n"])
    NUM["numMaxCertN"] = r"\num{%d}" % cache[big]["n"]
    NUM["numMaxCertM"] = r"\num{%d}" % max(cache[n]["m"] for n in cert)
    return cache


def runs(tag, directory=COLAB):
    """All result records of one tag, keyed by (graph, seed)."""
    out = {}
    for f in glob.glob(os.path.join(directory, f"{tag}_*.json")):
        d = json.load(open(f))
        if d.get("tag", tag) == tag or directory != COLAB:
            out[(d["graph"], d["seed"])] = d
    return out


def checked(rec):
    """Checked bound of a record, or None."""
    if rec and rec.get("check") == "ok":
        return rec["certified"]
    return None


def fmt_pct(x, nd=2):
    return "--" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def pace_exact():
    import pandas as pd
    ref = pd.read_csv(os.path.join(ROOT, "data", "pace2021_exact_kapoce_bounds.csv"))
    ref["inst"] = ref["instance"]
    ref = ref.set_index("inst")
    ref["density"] = ref["m"] / (ref["n"] * (ref["n"] - 1) / 2)
    ref["opt"] = np.where(ref["solved_1h"] == 1, ref["opt_or_empty"], np.nan)
    key = lambda i: f"pace-exact/{i}"
    # bounds per instance
    B = {}
    old = pd.read_csv(os.path.join(RES, "pace_exact.csv"))
    old["inst"] = old["instance"].str.split("/").str[-1]
    for a, nm in [("lb:greedy", "greedy triangle packing"), ("lb:tri-mwu", "fractional triangle packing")]:
        B[nm] = np.ceil(old[old.algo == a].set_index("inst")["lb"] - 1e-6)
    for tag, seed, nm in [("lpack", 1, "sparse star packing"), ("c2x", 0, "CertiFlip bound"),
                          ("ltri", 0, "metric LP on $P$, triangle rows"),
                          ("lwarm", 0, "CertiFlip bound, 1800 s")]:
        R = runs(tag)
        vals = {}
        for i in ref.index:
            r = R.get((key(i), seed))
            if r is None:
                continue
            if tag == "ltri":
                vals[i] = np.ceil(r["value"] - 1e-6) if r.get("converged") else np.nan
            else:
                v = checked(r)
                vals[i] = np.nan if v is None else v
        B[nm] = pd.Series(vals, dtype=float)
    B["B\\&B star packing~\\cite{BlasiusEtAl22sea}"] = ref["low_star"].astype(float)
    B["B\\&B $P_3$ packing~\\cite{BlasiusEtAl22sea}"] = ref["low_p3"].astype(float)
    solved = ref[ref["solved_1h"] == 1]
    rows = []
    for nm, v in B.items():
        v = v.reindex(solved.index)
        ok = v.notna()
        if not ok.any():
            rows.append(f"{nm} & \\pending{{}} & & & \\\\")
            continue
        r = v[ok] / solved["opt"][ok].clip(lower=1)
        rows.append(f"{nm} & {int(ok.sum())} & {r.mean():.4f} & {r.min():.3f} & "
                    f"{int((v[ok] >= solved['opt'][ok]).sum())} \\\\")
    write("pace_bounds", "\\begin{tabular}{lrrrr}\n\\toprule\nbound & instances & mean LB/OPT & "
          "min LB/OPT & LB $=$ OPT \\\\\n\\midrule\n" + "\n".join(rows) +
          "\n\\bottomrule\n\\end{tabular}\n")
    ours = B["CertiFlip bound"].reindex(solved.index)
    star = B["B\\&B star packing~\\cite{BlasiusEtAl22sea}"].reindex(solved.index)
    if ours.notna().sum() == len(solved):
        NUM["numPaceOptOurs"] = str(int((ours >= solved["opt"]).sum()))
        NUM["numPaceMeanOurs"] = f"{(ours / solved['opt'].clip(lower=1)).mean():.4f}"
        NUM["numPaceMeanStar"] = f"{(star / solved['opt'].clip(lower=1)).mean():.4f}"
        NUM["numPaceOursAbove"] = str(int((ours > star).sum()))
        NUM["numPaceOursBelow"] = str(int((ours < star).sum()))
        NUM["numPaceOursEqual"] = str(int((ours == star).sum()))
    # density bands
    bands = [(0, 0.1), (0.1, 0.3), (0.3, 0.5), (0.5, 1.01)]
    brow = []
    pk = B["sparse star packing"].reindex(solved.index)
    for lo, hi in bands:
        sel = solved[(solved.density >= lo) & (solved.density < hi)]
        cells = [f"$[{lo:g},{min(hi, 1):g})$", str(len(sel))]
        for v in (ours, pk, star):
            v = v.reindex(sel.index)
            cells.append("--" if v.isna().all() else f"{(v / sel['opt'].clip(lower=1)).mean():.4f}")
        brow.append(" & ".join(cells) + " \\\\")
    write("pace_density", "\\begin{tabular}{lrrrr}\n\\toprule\nedge density & instances & CertiFlip "
          "& sparse star packing & B\\&B star packing \\\\\n\\midrule\n" + "\n".join(brow) +
          "\n\\bottomrule\n\\end{tabular}\n")
    # further bounds on the solved instances, each on the instances where it ran
    for nm, kk in (("sparse star packing", "Pack"), ("CertiFlip bound, 1800 s", "Warm"),
                    ("metric LP on $P$, triangle rows", "Tri")):
        v = B[nm].reindex(solved.index)
        ok = v.notna()
        if ok.any():
            o = solved["opt"][ok].clip(lower=1)
            NUM[f"numPace{kk}N"] = str(int(ok.sum()))
            NUM[f"numPace{kk}Mean"] = f"{(v[ok] / o).mean():.4f}"
            NUM[f"numPace{kk}Opt"] = str(int((v[ok] >= solved["opt"][ok]).sum()))
            NUM[f"numPace{kk}OursMean"] = f"{(ours[ok] / o).mean():.4f}"
            NUM[f"numPace{kk}StarMean"] = f"{(star[ok] / o).mean():.4f}"
            NUM[f"numPace{kk}OursOpt"] = str(int((ours[ok] >= solved["opt"][ok]).sum()))
    # open instances: best checked bound against the published bounds
    unsolved = ref[ref["solved_1h"] != 1]
    C = runs("c2x")
    orows, improved = [], 0
    for i in unsolved.index:
        vals = [B[nm].get(i) for nm in ("CertiFlip bound", "sparse star packing",
                                        "CertiFlip bound, 1800 s")]
        vals = [x for x in vals if x is not None and np.isfinite(x)]
        if not vals:
            continue
        lb = int(max(vals))
        pub = int(max(ref.loc[i, "low_star"], ref.loc[i, "low_p3"]))
        ub = C.get((key(i), 0), {}).get("cost")
        if lb > pub:
            improved += 1
            ubs = "--" if ub is None else f"\\num{{{int(ub)}}}"
            orows.append(f"{tt(i.replace('.gr', ''))} & {ref.loc[i, 'n']} & {ref.loc[i, 'm']} & "
                         f"\\num{{{int(ref.loc[i, 'upper'])}}} & \\num{{{pub}}} & "
                         f"{ubs} & \\num{{{lb}}} \\\\")
    write("pace_open", "\\begin{tabular}{lrrrrrr}\n\\toprule\n"
          "instance & $n$ & $m$ & published UB & published LB & our UB & checked LB \\\\\n"
          "\\midrule\n" + ("\n".join(orows) if orows else "\\multicolumn{7}{c}{none} \\\\") +
          "\n\\bottomrule\n\\end{tabular}\n")
    NUM["numPaceOpenImproved"] = str(improved)
    NUM["numPaceOpen"] = str(len(unsolved))
    return ref, B


def best_costs():
    """Best cost per SNAP graph: (archived, i.e. with a re-checkable clustering;
    any run)."""
    import pandas as pd
    arch, anyrun = {}, {}

    def upd(dct, g, c):
        if c is not None and np.isfinite(c):
            dct[g] = min(dct.get(g, np.inf), int(c))
    for (g, sd), r in runs("c2").items():
        upd(arch, g, r.get("cost"))
    for tag in ("s2h", "s2h150", "s2h60"):
        for (g, sd), r in runs(tag).items():
            if tag == "s2h" and sd == 0:
                upd(arch, g, r.get("ours_final"))
                if r.get("kapoce_valid"):
                    upd(arch, g, r.get("kapoce"))
            upd(anyrun, g, r.get("ours_final"))
            if r.get("kapoce_valid"):
                upd(anyrun, g, r.get("kapoce"))
    for tag in ("s1", "m3", "c1"):
        for (g, sd), r in runs(tag).items():
            for k in ("ours", "kapoce", "cost"):
                if k == "kapoce" and r.get("kapoce_valid") is False:
                    continue
                upd(anyrun, g, r.get(k))
    old = pd.read_csv(os.path.join(RES, "snap.csv"))
    for _, r in old[old.cost.notna()].iterrows():
        upd(anyrun, r["instance"], r["cost"])
    for g, c in arch.items():
        upd(anyrun, g, c)
    return arch, anyrun


def snap():
    import pandas as pd
    cache = json.load(open(os.path.join(RES, "mpc", "instances.json")))
    arch, anyrun = best_costs()
    C2 = runs("c2")
    PK = runs("lpack")
    LT = runs("ltri")
    LS = runs("lstar")
    tri = pd.concat([pd.read_csv(os.path.join(RES, f)) for f in
                     ("snap.csv", "snap_heldout_bounds.csv")
                     if os.path.exists(os.path.join(RES, f))])
    tri = tri[tri.algo.isin(["lb:greedy", "lb:tri-mwu"]) & tri.lb.notna()]
    tri = tri.groupby("instance")["lb"].max().apply(lambda x: np.ceil(x - 1e-6))
    rows, gaps, facs, spreads, cf_better, pk_better, lp_rows = [], [], [], [], [], [], []
    role_gaps = {}

    def lp_gap(rec, ub):
        # an LP value is reported only when cutting planes converged: then it
        # is the optimum of the relaxation (in floating point, not certified)
        if rec is None or not rec.get("converged"):
            return None
        return 100 * (ub - np.ceil(rec["value"] - 1e-6)) / max(1.0, np.ceil(rec["value"] - 1e-6))

    for role, names in (("dev", DEV), ("held-out", HELD)):
        for name in names:
            if cache[name]["pairs"] is None:
                continue
            lbs = [checked(C2.get((name, sd))) for sd in range(3)]
            lbs = [x for x in lbs if x is not None]
            ub_a, ub = arch.get(name), anyrun.get(name)
            if not lbs or ub_a is None:
                rows.append(f"{tt(name)} & \\multicolumn{{7}}{{c}}{{\\pending{{}}}} \\\\")
                continue
            cf = max(lbs)
            spreads.append(100 * (max(lbs) - min(lbs)) / cf)
            pk = checked(PK.get((name, 0)))
            lb = cf if pk is None else max(cf, pk)
            if pk is not None:
                (cf_better if cf > pk else pk_better).append(name)
            gap_a = 100 * (ub_a - lb) / lb
            gaps.append(gap_a)
            role_gaps.setdefault(role, []).append(gap_a)
            tb = tri.get(name)
            fac = None if tb is None or not np.isfinite(tb) else ub / tb
            if fac is not None:
                facs.append(fac)
            t_lb = np.median([C2[(name, sd)]["lb_time"] for sd in range(3) if (name, sd) in C2])
            TIMES[name] = t_lb
            gap_b = 100 * (ub - lb) / lb
            g_tri = lp_gap(LT.get((name, 0)), ub)
            g_star = lp_gap(LS.get((name, 0)), ub)
            if (name, 0) in LT or (name, 0) in LS:
                lp_rows.append((name, LT.get((name, 0)), LS.get((name, 0)), ub, pk, cf))
            b = lambda v, other: (f"\\textbf{{\\num{{{v}}}}}" if other is not None and v > other
                                  else f"\\num{{{v}}}")
            rows.append(f"{tt(name)} & \\num{{{ub_a}}} & {b(cf, pk)} & "
                        f"{'--' if pk is None else b(pk, cf)} & {gap_a:.2f} & {gap_b:.2f} & "
                        f"{fmt_pct(g_tri, 2) if g_tri is not None else '--'} & "
                        f"{fmt_pct(fac, 2) if fac else '--'} \\\\")
        rows.append("\\midrule")
    write("snap_bounds", "\\begin{tabular}{lrrrrrrr}\n\\toprule\n"
          "graph & UB & \\multicolumn{2}{c}{checked LB} & gap & gap$^*$ & tri.\\ LP & tri.\\\\\n"
          "\\cmidrule(lr){3-4}\n"
          " & & CertiFlip & packing & (\\%) & (\\%) & gap (\\%) & factor \\\\\n\\midrule\n"
          + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
    # metric LP on P against the checked bounds, where it was run
    lrows = []
    for name, rt, rs, ub, pk, cf in lp_rows:
        def cell(r):
            if r is None:
                return "--", "--"
            v = f"\\num{{{np.ceil(r['value'] - 1e-6):.0f}}}"
            if not r.get("converged"):
                v = f"({v})"
            return v, f"{r.get('time', float('nan')):.0f}"
        (vt, tt_), (vs, ts) = cell(rt), cell(rs)
        best = max(x for x in (pk, cf) if x is not None)
        lrows.append(f"{tt(name)} & \\num{{{ub}}} & \\num{{{best}}} & {vt} & {tt_} & {vs} & {ts} \\\\")
    write("snap_lp", "\\begin{tabular}{lrrrrrr}\n\\toprule\n"
          "graph & best & best checked & \\multicolumn{2}{c}{LP$_P$, triangle rows} & "
          "\\multicolumn{2}{c}{LP$_P$, + star rows} \\\\\n"
          "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n"
          " & known & LB & value & time (s) & value & time (s) \\\\\n\\midrule\n"
          + ("\n".join(lrows) if lrows else "\\multicolumn{7}{c}{\\pending{}} \\\\") +
          "\n\\bottomrule\n\\end{tabular}\n")
    if gaps:
        NUM["numSnapGapLo"] = f"{min(gaps):.2f}\\%"
        NUM["numSnapGapHi"] = f"{max(gaps):.1f}\\%"
        NUM["numSnapGapMedian"] = f"{np.median(gaps):.1f}\\%"
        NUM["numSnapGapMedDev"] = f"{np.median(role_gaps['dev']):.1f}\\%"
        NUM["numSnapGapMedHeld"] = f"{np.median(role_gaps['held-out']):.1f}\\%"
        NUM["numSnapGapMinHeld"] = f"{min(role_gaps['held-out']):.2f}\\%"
        NUM["numSnapSpreadMax"] = f"{max(spreads):.2f}\\%"
    if facs:
        NUM["numBaseFacLo"] = f"{min(facs):.2f}"
        NUM["numBaseFacHi"] = f"{max(facs):.2f}"
    pt = [PK[(n, 0)]["time"] for n in DEV + HELD if (n, 0) in PK]
    ct = [TIMES[n] for n in TIMES]
    if pt:
        NUM["numPackTimeLo"] = f"{min(pt):.0f}"
        NUM["numPackTimeHi"] = f"{max(pt):.0f}"
        NUM["numPackTimeMedian"] = f"{np.median(pt):.0f}"
    if ct:
        NUM["numCfLbTimeMedian"] = f"{np.median(ct):.0f}"
    for name, rt, rs, ub, pk, cf in lp_rows:
        if name == "ca-GrQc":
            if rt is not None and rt.get("converged"):
                v = np.ceil(rt["value"] - 1e-6)
                NUM["numGrqcTri"] = f"\\num{{{v:.0f}}}"
                NUM["numGrqcTriGap"] = f"{100 * (ub - v) / v:.1f}\\%"
                NUM["numGrqcTriTime"] = f"{rt['time']:.0f}"
            if rs is not None:
                v = np.ceil(rs["value"] - 1e-6)
                NUM["numGrqcStar"] = f"\\num{{{v:.0f}}}"
                NUM["numGrqcStarGap"] = f"{100 * (ub - v) / v:.2f}\\%"
                NUM["numGrqcStarTime"] = f"{rs['time']:.0f}"
                NUM["numGrqcBest"] = f"\\num{{{max(pk, cf)}}}"
    NUM["numPackBetter"] = str(len(pk_better))
    NUM["numCfBetter"] = str(len(cf_better))
    NUM["numCfBetterList"] = ", ".join(tt(x) for x in cf_better)


def signed_rank_cdf(n):
    """Exact null distribution of the Wilcoxon signed-rank statistic W+ for n
    untied, non-zero differences: P(W+ <= w) for w = 0..n(n+1)/2."""
    counts = np.zeros(n * (n + 1) // 2 + 1)
    counts[0] = 1
    for r in range(1, n + 1):
        counts[r:] = counts[r:] + counts[:-r].copy()
    return np.cumsum(counts) / 2.0 ** n


def hodges_lehmann(d, alpha=0.05):
    """Hodges-Lehmann estimate of the location of the paired differences d and
    its exact (1 - alpha) confidence interval from the Walsh averages."""
    d = np.asarray(d, dtype=float)
    n = len(d)
    w = np.sort([(d[i] + d[j]) / 2 for i in range(n) for j in range(i, n)])
    cdf = signed_rank_cdf(n)
    k = int(np.searchsorted(cdf, alpha / 2, side="right"))  # P(W+ <= k-1) <= alpha/2
    lo, hi = (w[k - 1], w[len(w) - k]) if k >= 1 else (-np.inf, np.inf)
    return float(np.median(w)), float(lo), float(hi)


def holm(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for i, j in enumerate(order):
        run = max(run, min(1.0, (len(p) - i) * p[j]))
        adj[j] = run
    return adj


BIG = 1.0  # relative difference used for a run that one solver lost by forfeit


def paired(tag, field="ours"):
    """Per graph: list of relative differences (PXMem - KaPoCE)/KaPoCE in %,
    with forfeits (invalid KaPoCE output: PXMem wins; no PXMem solution at the
    cut: PXMem loses) coded as -+100%."""
    out, forfeits = {}, {}
    for (g, sd), r in runs(tag).items():
        k, o = r.get("kapoce"), r.get(field)
        if not r.get("kapoce_valid", False):
            d = -100.0 * BIG
            forfeits.setdefault(g, []).append("kapoce")
        elif o is None:
            d = 100.0 * BIG
            forfeits.setdefault(g, []).append("pxmem")
        else:
            d = 100.0 * (o - k) / max(1, k)
        out.setdefault(g, {})[sd] = d
    return out, forfeits


def h2h(tag, label):
    from scipy.stats import binomtest, wilcoxon
    P, forfeits = paired(tag)
    graphs = [g for g in HELD if g in P]
    rows, meds, ps = [], [], []
    for g in graphs:
        d = np.array([P[g][s] for s in sorted(P[g])])
        pos, neg = int((d > 0).sum()), int((d < 0).sum())
        p = binomtest(neg, pos + neg).pvalue if pos + neg else 1.0
        ps.append(p)
        meds.append(float(np.median(d)))
        hl = hodges_lehmann(d) if len(d) >= 6 else (np.nan, np.nan, np.nan)
        rows.append([g, len(d), neg, int((d == 0).sum()), pos, float(np.median(d)), p, hl,
                     len(forfeits.get(g, []))])
    if not rows:
        write(f"h2h_{label}", "\\pending{}\n")
        return None
    adj = holm(ps)
    lines = []
    for r, a in zip(rows, adj):
        g, n, w, t, l, med, p, hl, ff = r
        ci = "--" if not np.isfinite(hl[0]) else f"${hl[0]:+.3f}$ [${hl[1]:+.3f}$, ${hl[2]:+.3f}$]"
        mark = r"$^\dagger$" if ff else ""
        lines.append(f"{tt(g)} & {n} & {w}/{t}/{l} & ${med:+.3f}$ & {ci} & {p:.3f} & {a:.3f}"
                     f"{mark} \\\\")
    body = ("\\begin{tabular}{lrcrlrr}\n\\toprule\ngraph & runs & W/T/L & median (\\%) & "
            "HL (\\%) [95\\% CI] & $p$ & $p_{\\mathrm{Holm}}$ \\\\\n\\midrule\n" + "\n".join(lines))
    m = np.array(meds)
    nz = m[m != 0]
    wp = wilcoxon(nz).pvalue if len(nz) >= 1 else 1.0
    sp = binomtest(int((m < 0).sum()), int((m != 0).sum())).pvalue if (m != 0).any() else 1.0
    body += (f"\n\\midrule\n\\multicolumn{{7}}{{l}}{{graph medians: {int((m < 0).sum())} better, "
             f"{int((m == 0).sum())} equal, {int((m > 0).sum())} worse; Wilcoxon $p={wp:.3f}$, "
             f"sign test $p={sp:.3f}$}} \\\\\n\\bottomrule\n\\end{{tabular}}\n")
    write(f"h2h_{label}", body)
    key = {"600": "Six", "150": "OneFifty", "60": "Sixty"}[label]
    NUM[f"numHH{key}Better"] = str(int((m < 0).sum()))
    NUM[f"numHH{key}Equal"] = str(int((m == 0).sum()))
    NUM[f"numHH{key}Worse"] = str(int((m > 0).sum()))
    NUM[f"numHH{key}Wilcoxon"] = f"{wp:.3f}"
    NUM[f"numHH{key}Sign"] = f"{sp:.3f}"
    NUM[f"numHH{key}Graphs"] = str(len(m))
    NUM[f"numHH{key}Runs"] = str(sum(r[1] for r in rows))
    NUM[f"numHH{key}ForfeitK"] = str(sum(v.count("kapoce") for v in forfeits.values()))
    NUM[f"numHH{key}ForfeitP"] = str(sum(v.count("pxmem") for v in forfeits.values()))
    NUM[f"numHH{key}MaxAbs"] = f"{max(abs(x) for x in meds):.2f}\\%"
    # sensitivity: PXMem scored at T instead of at KaPoCE's elapsed time
    PT, _ = paired(tag, "ours_at_T")
    mt = np.array([float(np.median([PT[g][s] for s in PT[g]])) for g in graphs])
    NUM[f"numHH{key}AtTBetter"] = str(int((mt < 0).sum()))
    NUM[f"numHH{key}AtTWorse"] = str(int((mt > 0).sum()))
    nzt = mt[mt != 0]
    NUM[f"numHH{key}AtTWilcoxon"] = f"{wilcoxon(nzt).pvalue:.3f}" if len(nzt) else "1.000"
    return rows


PACE_DEV = {167, 173, 175, 176, 178, 179, 180, 191, 193, 198}


def pace_heur():
    """Protocol v2 on the PACE heuristic track, one run per instance at 600 s;
    the instance is the unit.  Development twins are reported separately."""
    from scipy.stats import binomtest, wilcoxon
    R = runs("s2p")
    held, dev = [], []
    for (g, sd), r in R.items():
        i = int(g.split("heur")[-1].split(".")[0])
        k, o = r.get("kapoce"), r.get("ours")
        if not r.get("kapoce_valid", False):
            d, ratio = -100.0, (1.0, np.inf)
        elif o is None:
            d, ratio = 100.0, (np.inf, 1.0)
        else:
            d = 100.0 * (o - k) / max(1, k)
            best = max(1, min(o, k))
            ratio = (o / best, k / best)
        (dev if i in PACE_DEV else held).append((i, d, ratio, r))
    if not held:
        write("pace_heur2", "\\begin{tabular}{l}\n\\pending{}\n\\end{tabular}\n")
        return
    d = np.array([x[1] for x in held])
    nz = d[d != 0]
    wp = wilcoxon(nz).pvalue if len(nz) else 1.0
    sp = binomtest(int((d < 0).sum()), int((d != 0).sum())).pvalue if (d != 0).any() else 1.0
    NUM.update({"numPHn": str(len(d)), "numPHBetter": str(int((d < 0).sum())),
                "numPHEqual": str(int((d == 0).sum())), "numPHWorse": str(int((d > 0).sum())),
                "numPHWilcoxon": f"{wp:.3g}", "numPHSign": f"{sp:.3g}",
                "numPHKapoceInvalid": str(sum(1 for x in held if not x[3].get("kapoce_valid", False))),
                "numPHNoSolution": str(sum(1 for x in held if x[3].get("ours") is None
                                           and x[3].get("kapoce_valid", False)))})
    # performance profile: fraction of instances within factor tau of the better solver
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    taus = np.logspace(0, np.log10(1.02), 200)
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for j, (nm, col) in enumerate((("PXMem", "#2A78D6"), ("KaPoCE", "#EB6834"))):
        r = np.array([x[2][j] for x in held])
        ax.step(taus, [(r <= t).mean() for t in taus], where="post", label=nm, color=col)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau$ (cost / best of the two)")
    ax.set_ylabel("fraction of instances")
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.join(ROOT, "paper_mpc", "figures"), exist_ok=True)
    fig.savefig(os.path.join(ROOT, "paper_mpc", "figures", "profile_pace_heur.pdf"))
    plt.close(fig)
    # W/T/L by instance size; development twins in a separate row
    def row(label, xs):
        dx = np.array([x[1] for x in xs])
        diff = sum(x[3]["ours"] - x[3]["kapoce"] for x in xs
                   if x[3].get("ours") is not None and x[3].get("kapoce_valid", False))
        return (f"{label} & {len(xs)} & {int((dx < 0).sum())}/{int((dx == 0).sum())}/"
                f"{int((dx > 0).sum())} & ${diff:+d}$ & ${np.median(dx):+.4f}$ \\\\")
    bins = [(0, 1e3, "$<10^3$"), (1e3, 1e4, "$10^3$--$10^4$"), (1e4, 1e5, "$10^4$--$10^5$"),
            (1e5, np.inf, "$\\ge 10^5$")]
    lines = [row(lab, [x for x in held if lo <= x[3]["n"] < hi]) for lo, hi, lab in bins
             if any(lo <= x[3]["n"] < hi for x in held)]
    lines.append("\\midrule")
    lines.append(row("without dev.\\ twins", held))
    if dev:
        lines.append(row("development twins", dev))
    write("pace_heur2", "\\begin{tabular}{lrcrr}\n\\toprule\n"
          "vertices & instances & W/T/L & $\\sum$ PXMem$-$KaPoCE & median (\\%) \\\\\n"
          "\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n")
    dd = np.array([x[1] for x in dev])
    if len(dd):
        NUM.update({"numPDn": str(len(dd)), "numPDBetter": str(int((dd < 0).sum())),
                    "numPDEqual": str(int((dd == 0).sum())), "numPDWorse": str(int((dd > 0).sum()))})


def timing():
    """Budget overruns, CPU models and the cost of checking."""
    models = set()
    for tag, key in (("c2", "Snap"), ("c2x", "Pace")):
        R = runs(tag)
        if not R:
            continue
        T = next(iter(R.values()))["T"]
        t = np.array([r["time"] for r in R.values()])
        NUM[f"numOver{key}Count"] = str(int((t > 1.01 * T).sum()))
        NUM[f"numOver{key}Runs"] = str(len(t))
        NUM[f"numOver{key}Max"] = f"{t.max():.0f}"
        NUM[f"numOver{key}Median"] = f"{np.median(t):.0f}"
        NUM[f"numOver{key}Large"] = str(int((t > 1.2 * T).sum()))
        if key == "Snap":
            e = [r["lb_time"] for (g, sd), r in R.items() if g == "email-Enron"]
            if e:
                NUM["numEnronLbTime"] = f"{max(e):.0f}"
        models |= {r.get("cpu") for r in R.values() if r.get("cpu")}
        ok = [r for r in R.values() if r.get("check") == "ok"]
        rej = [r for r in R.values() if str(r.get("check", "")).startswith("rejected")]
        NUM[f"numCert{key}Ok"] = str(len(ok))
        NUM[f"numCert{key}Rejected"] = str(len(rej))
        if ok and key == "Snap":
            big = max(ok, key=lambda r: r["rows"])
            NUM["numCheckMaxGraph"] = big["graph"]
            NUM["numCheckMaxRows"] = r"\num{%d}" % big["rows"]
            NUM["numCheckMaxTime"] = f"{big['check_time']:.0f}"
            NUM["numCheckTimeMax"] = f"{max(r['check_time'] for r in ok):.0f}"
    for tag in ("s2h", "s2p"):
        for r in runs(tag).values():
            m = r.get("machine", {}).get("Model name")
            if m:
                models.add(m)
    if models:
        NUM["numCpuModels"] = "; ".join(sorted(models)).replace("(R)", r"\textsuperscript{\textregistered}")


def pace_primal(ref):
    """Primal quality on the PACE exact track: CertiFlip (checked runs), SCC on
    the complete encoding, and the earlier local runs of the other methods."""
    import pandas as pd
    costs = {}
    for (g, sd), r in runs("c2x").items():
        costs.setdefault("CertiFlip (60 s)", {})[g.split("/")[-1]] = r["cost"]
    for (g, sd), r in runs("sccevo", os.path.join(RES, "scc")).items():
        costs.setdefault("SCC~\\cite{hausberger2025scc} (60 s)", {})[g.split("/")[-1]] = r["cost"]
    old = pd.read_csv(os.path.join(RES, "pace_exact.csv"))
    old["inst"] = old["instance"].str.split("/").str[-1]
    for a, nm in (("kapoce", "KaPoCE (60 s)$^a$"), ("flip", "IteratedFlip$^a$"),
                  ("leiden", "Leiden-CPM$^a$"), ("pivot", "Pivot$^a$")):
        s_ = old[(old.algo == a) & old.cost.notna()].groupby("inst")["cost"].min()
        costs[nm] = s_.to_dict()
    best = ref["upper"].astype(float).copy()
    for nm, c in costs.items():
        for i, v in c.items():
            best[i] = min(best[i], v)
    rows = []
    for nm, c in costs.items():
        if not c:
            continue
        v = pd.Series(c)
        b = best.reindex(v.index)
        gap = 100 * (v - b) / b.clip(lower=1)
        rows.append(f"{nm} & {int((v <= b).sum())}/{len(v)} & {gap.mean():.3f} & {gap.max():.2f} \\\\")
    write("pace_primal", "\\begin{tabular}{lrrr}\n\\toprule\nmethod & best known & mean gap (\\%) & "
          "max gap (\\%) \\\\\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
    sc = costs.get("SCC~\\cite{hausberger2025scc} (60 s)", {})
    if sc:
        NUM["numSccHit"] = str(int(sum(sc[i] <= best[i] for i in sc)))
        NUM["numSccRuns"] = str(len(sc))
    cf = costs.get("CertiFlip (60 s)", {})
    if cf:
        NUM["numCfHit"] = str(int(sum(cf[i] <= best[i] for i in cf)))


def static_numbers():
    with open(os.path.join(ROOT, "experiments", "check_certificate.py")) as fh:
        NUM["numCheckerLines"] = str(sum(1 for _ in fh))


def write_numbers():
    old = dict(DEFAULTS)
    if os.path.exists(NUMBERS):
        for line in open(NUMBERS):
            if line.startswith(r"\newcommand{"):
                k = line[len(r"\newcommand{\\") - 1:].split("}", 1)[0]
                old[k] = line.split("}{", 1)[1].rstrip()[:-1]
    old.update(NUM)
    with open(NUMBERS, "w") as fh:
        fh.write("% generated by experiments/mpc_tables.py; do not edit\n")
        for k in sorted(old):
            fh.write(f"\\newcommand{{\\{k}}}{{{old[k]}}}\n")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    instances()
    ref_, _ = pace_exact()
    pace_primal(ref_)
    snap()
    for tag, label in (("s2h", "600"), ("s2h150", "150"), ("s2h60", "60")):
        h2h(tag, label)
    pace_heur()
    timing()
    static_numbers()
    write_numbers()
