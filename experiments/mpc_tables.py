"""Tables and numbers of the MPC manuscript (paper_mpc/) from the raw results.

Every table is written to paper_mpc/tables/, every number used in the text to
paper_mpc/numbers.tex.  Results that are not yet available are printed as
``pending'' so that the manuscript always compiles.

usage: python experiments/mpc_tables.py [instances|all]"""
import glob
import json
import subprocess
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
    "numCpuModels", "numRecheckN", "numRecheckOk", "numRecheckMax", "numRecheckTotal")}
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
    PK = runs("lpack")
    for role, names in (("dev", DEV), ("held-out", HELD)):
        for name in names:
            r = cache[name]
            pairs = r"\num{%d}" % r["pairs"] if r["pairs"] is not None else "--"
            if r["pairs"] is None and (name, 0) in PK:
                # measured by the later packing run, which ignored the threshold
                pairs = r"\num{%d}$^\ast$" % PK[(name, 0)]["pairs"]
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


LOCAL = os.path.join(RES, "local")


def runs(tag, directory=COLAB):
    """All result records of one tag, keyed by (graph, seed); for the Colab
    tags also the runs of the same jobs on a local machine (results/local)."""
    out = {}
    dirs = [directory] + ([LOCAL] if directory == COLAB else [])
    for f in (f for d in dirs for f in glob.glob(os.path.join(d, f"{tag}_*.json"))):
        d = json.load(open(f))
        if d.get("tag", tag) == tag or directory != COLAB:
            out[(d["graph"], d.get("seed", 0))] = d
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
    # the larger of two valid bounds is a valid bound
    B["max of CertiFlip and B\\&B star"] = pd.concat(
        [B["CertiFlip bound"], ref["low_star"].astype(float)], axis=1).max(axis=1, skipna=False)
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
        mx = np.maximum(ours, star)
        NUM["numPaceMaxOpt"] = str(int((mx >= solved["opt"]).sum()))
        NUM["numPaceMaxMean"] = f"{(mx / solved['opt'].clip(lower=1)).mean():.4f}"
    fr = B["fractional triangle packing"].reindex(solved.index)
    NUM["numPaceFracMean"] = f"{(fr / solved['opt'].clip(lower=1)).mean():.4f}"
    # density bands
    bands = [(0, 0.1), (0.1, 0.3), (0.3, 0.5), (0.5, 1.01)]
    brow = []
    pk = B["sparse star packing"].reindex(solved.index)
    for lo, hi in bands:
        sel = solved[(solved.density >= lo) & (solved.density < hi)]
        cells = [f"$[{lo:g},{min(hi, 1):g})$", str(len(sel))]
        for v in (ours, pk, fr, star):
            v = v.reindex(sel.index)
            cells.append("--" if v.isna().all() else f"{(v / sel['opt'].clip(lower=1)).mean():.4f}")
        brow.append(" & ".join(cells) + " \\\\")
    write("pace_density", "\\begin{tabular}{lrrrrr}\n\\toprule\nedge density & instances & CertiFlip "
          "& sparse star & fractional & B\\&B star \\\\\n"
          " & & bound & packing & triangle packing & packing \\\\\n\\midrule\n" + "\n".join(brow) +
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
    orows, improved, below = [], 0, []
    for i in unsolved.index:
        vals = [B[nm].get(i) for nm in ("CertiFlip bound", "sparse star packing",
                                        "CertiFlip bound, 1800 s")]
        vals = [x for x in vals if x is not None and np.isfinite(x)]
        if not vals:
            continue
        lb = int(max(vals))
        pub = int(max(ref.loc[i, "low_star"], ref.loc[i, "low_p3"]))
        ub = C.get((key(i), 0), {}).get("cost")
        if lb < pub:
            below.append(100 * (pub - lb) / pub)
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
    NUM["numPaceOpenBelow"] = str(len(below))
    NUM["numPaceOpenMaxBelow"] = f"{max(below):.1f}\\%" if below else "0"
    NUM["numPaceOpen"] = str(len(unsolved))
    # the root bound program of KaPoCE rerun on the input graphs, without its
    # data reductions (tag lkroot), against the published root bounds
    K = {i: runs("lkroot").get((key(i), 0)) for i in ref.index}
    K = {i: r for i, r in K.items() if r is not None}
    if K:
        p3 = [i for i, r in K.items() if r.get("p3") is not None]
        st = [i for i, r in K.items() if r.get("star") is not None]
        ds = [K[i]["star"] - int(ref.loc[i, "low_star"]) for i in st]
        NUM["numKrootPaceN"] = str(len(K))
        NUM["numKrootPaceStarN"] = str(len(st))
        NUM["numKrootPaceTimeout"] = str(len(K) - len(st))
        NUM["numKrootPaceTimeoutList"] = ", ".join(
            tt(i.replace(".gr", "")) for i in sorted(K) if K[i].get("star") is None) or "none"
        if os.path.exists(os.path.join(COLAB, "queue.json")):
            NUM["numKrootPaceColab"] = str(sum(
                1 for j in json.load(open(os.path.join(COLAB, "queue.json")))
                if j["id"].startswith("lkroot_pace") and j["status"] == "done"))
        NUM["numKrootPacePthreeEqual"] = str(sum(K[i]["p3"] == int(ref.loc[i, "low_p3"]) for i in p3))
        NUM["numKrootPacePthreeN"] = str(len(p3))
        NUM["numKrootPaceStarEqual"] = str(sum(d == 0 for d in ds))
        NUM["numKrootPaceStarAbove"] = str(sum(d > 0 for d in ds))
        NUM["numKrootPaceStarBelow"] = str(sum(d < 0 for d in ds))
        NUM["numKrootPaceStarDiffMax"] = str(max(map(abs, ds))) if ds else "0"
        NUM["numKrootPaceStarRelMax"] = f"{max(abs(d) / max(1, int(ref.loc[i, 'low_star'])) for d, i in zip(ds, st)) * 100:.2f}\\%"
        ts = [K[i]["star_time"] for i in st]
        NUM["numKrootPaceTimeMedian"] = f"{np.median(ts):.1f}"
        NUM["numKrootPaceTimeMax"] = f"{max(ts):.0f}"
        # the rerun star packing against our sparse star packing, same instances
        pk = B["sparse star packing"]
        both = [i for i in st if i in solved.index and np.isfinite(pk.get(i, np.nan))]
        if both:
            o = solved["opt"].reindex(both).clip(lower=1)
            NUM["numKrootPaceStarMean"] = f"{np.mean([K[i]['star'] for i in both] / o):.4f}"
            NUM["numKrootPacePackMean"] = f"{np.mean(pk.reindex(both) / o):.4f}"
            NUM["numKrootPaceBothN"] = str(len(both))
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
    PL = runs("lplp")
    LT = runs("ltri")
    LS = runs("lstar")
    EQ = runs("lpack+eq")   # the packing alone with the budget of pack.+LP
    LK = runs("lkstar")     # the root star packing of KaPoCE, exported and checked
    tri = pd.concat([pd.read_csv(os.path.join(RES, f)) for f in
                     ("snap.csv", "snap_heldout_bounds.csv")
                     if os.path.exists(os.path.join(RES, f))])
    tri = tri[tri.algo.isin(["lb:greedy", "lb:tri-mwu"]) & tri.lb.notna()]
    tri = tri.groupby("instance")["lb"].max().apply(lambda x: np.ceil(x - 1e-6))
    root, root_to = {}, set()
    qstat_root = {}
    if os.path.exists(os.path.join(COLAB, "queue.json")):
        qstat_root = {j["id"]: j["status"] for j in json.load(open(os.path.join(COLAB, "queue.json")))}
    for f in sorted(glob.glob(os.path.join(RES, "kapoce_root_snap_*.csv"))):
        d = pd.read_csv(f)
        for _, r in d[d.bound == "star"].iterrows():
            root[r.graph] = (int(r.lb), float(r.time))
    # the same program as a fleet job (colab_run_lp.py, mode kroot)
    for (g, sd), r in runs("lkroot").items():
        if g not in DEV + HELD:
            continue
        if r.get("star") is not None:
            root[g] = (int(r["star"]), float(r["star_time"]))
        elif r.get("status") == "timeout":
            root_to.add(g)
    rows, gaps, facs, spreads, cf_better, pk_better, lp_rows = [], [], [], [], [], [], []
    role_gaps, gaps_cf, gaps_pk, longer, seed_below, root_better, root_n = {}, [], [], [], [], [], 0
    plp_better, plp_below, plp_gain, pack_only, no_bound = [], [], [], [], []
    root_cmp = []  # (graph, B&B root star bound, our best checked bound, its time)
    root_cert, root_nocert, eq_rows = [], [], []

    def lp_gap(rec, ub):
        # an LP value is reported only when cutting planes converged: then it
        # is the optimum of the relaxation (in floating point, not certified)
        if rec is None or not rec.get("converged"):
            return None
        return 100 * (ub - np.ceil(rec["value"] - 1e-6)) / max(1.0, np.ceil(rec["value"] - 1e-6))

    for role, names in (("dev", DEV), ("held-out", HELD)):
        for name in names:
            if cache[name]["pairs"] is None:
                # above the support threshold: no CertiFlip run; the packing
                # alone was run without the threshold where memory allowed
                pk = checked(PK.get((name, 0)))
                ub_a, ub = arch.get(name), anyrun.get(name)
                ucell = "--" if ub_a is None else f"\\num{{{ub_a}}}"
                if pk is None or ub is None:
                    no_bound.append(name)
                    rows.append(f"{tt(name)} & {ucell} & \\multicolumn{{9}}{{c}}"
                                "{no bound: the packing job died three times (memory)} \\\\")
                    continue
                pack_only.append((name, role, pk, ub_a, ub))
                t_pk = PK[(name, 0)]["time"]
                tb = tri.get(name)
                fac = None if tb is None or not np.isfinite(tb) else ub / tb
                ga = "--" if ub_a is None else f"{100 * (ub_a - pk) / pk:.2f}"
                rows.append(f"{tt(name)} & {ucell} & -- & -- & "
                            f"\\num{{{pk}}} & {t_pk:.0f} & -- & -- & {ga} & "
                            f"{100 * (ub - pk) / pk:.2f} & {fmt_pct(fac, 2) if fac else '--'} \\\\")
                continue
            lbs = [checked(C2.get((name, sd))) for sd in range(3)]
            lbs = [x for x in lbs if x is not None]
            ub_a, ub = arch.get(name), anyrun.get(name)
            if not lbs or ub_a is None:
                rows.append(f"{tt(name)} & \\multicolumn{{10}}{{c}}{{\\pending{{}}}} \\\\")
                continue
            cf = max(lbs)
            spreads.append(100 * (max(lbs) - min(lbs)) / cf)
            pk = checked(PK.get((name, 0)))
            plp = checked(PL.get((name, 0)))
            prev = cf if pk is None else max(cf, pk)
            lb = prev if plp is None else max(prev, plp)
            eq = checked(EQ.get((name, 0)))
            lb_ours = lb  # the largest checked bound of the columns shown
            if eq is not None and eq > lb:
                lb = eq
            if eq is not None:
                eq_rows.append((name, role, eq, EQ[(name, 0)]["time"], plp, lb))
            if pk is not None:
                (cf_better if cf > pk else pk_better).append(name)
            if plp is not None:
                (plp_better if plp > prev else plp_below).append((name, plp - prev))
                # the block LPs against the packing they started from
                plp_gain.append((name, 100 * (plp - PL[(name, 0)]["pack_value"]) /
                                 PL[(name, 0)]["pack_value"]))
            gap_a = 100 * (ub_a - lb) / lb
            gaps.append(gap_a)
            # each procedure on its own: CertiFlip at its median seed, the packing alone
            gaps_cf.append(100 * (ub_a - np.median(lbs)) / np.median(lbs))
            if pk is not None:
                gaps_pk.append(100 * (ub_a - pk) / pk)
                if min(lbs) < pk < max(lbs):
                    seed_below.append(name)
            role_gaps.setdefault(role, []).append(gap_a)
            tb = tri.get(name)
            fac = None if tb is None or not np.isfinite(tb) else ub / tb
            if fac is not None:
                facs.append(fac)
            t_lbs = [C2[(name, sd)]["lb_time"] for sd in range(3) if (name, sd) in C2]
            t_lb = np.median(t_lbs)
            TIMES[name] = t_lb
            t_pk = PK[(name, 0)]["time"] if (name, 0) in PK else None
            if t_pk is not None and pk is not None and pk > cf and t_pk > max(t_lbs):
                longer.append(name)
            rt = root.get(name)
            kc = checked(LK.get((name, 0)))
            if rt is not None:
                root_n += 1
                root_cmp.append((name, rt[0], lb, rt[1]))
                if rt[0] > lb:
                    root_better.append(name)
                (root_cert if kc == rt[0] else root_nocert).append(name)
            gap_b = 100 * (ub - lb) / lb
            g_tri = lp_gap(LT.get((name, 0)), ub)
            g_star = lp_gap(LS.get((name, 0)), ub)
            lp_rows.append((name, LT.get((name, 0)), LS.get((name, 0)), ub, pk, lb))
            b = lambda v: (f"\\textbf{{\\num{{{v}}}}}" if v == lb else f"\\num{{{v}}}")
            running = qstat_root.get(f"lkroot_{name}_0") in ("pending", "running")
            # out of memory: the export job died with its machine three times
            oom = qstat_root.get(f"lkstar_{name}_0") == "failed"
            rcell = ("t.o." if name in root_to else "\\pending{}" if running else
                     "mem." if oom else "--") if rt is None else (
                f"\\num{{{rt[0]}}}" + ("$^\\dagger$" if rt[0] > lb else "") +
                ("" if kc == rt[0] else "$^\\circ$"))
            smark = "$^\\S$" if lb > lb_ours else ""
            rows.append(f"{tt(name)} & \\num{{{ub_a}}} & {b(cf)} & {t_lb:.0f} & "
                        f"{'--' if pk is None else b(pk)} & "
                        f"{'--' if t_pk is None else f'{t_pk:.0f}'} & "
                        f"{'--' if plp is None else b(plp)} & {rcell} & "
                        f"{gap_a:.2f}{smark} & {gap_b:.2f} & "
                        f"{fmt_pct(fac, 2) if fac else '--'} \\\\")
        rows.append("\\midrule")
    write("snap_bounds", "\\begin{tabular}{lrrrrrrrrrr}\n\\toprule\n"
          "graph & UB & \\multicolumn{5}{c}{checked LB} & B\\&B & gap & gap$^*$ & tri.\\\\\n"
          "\\cmidrule(lr){3-7}\n"
          " & & CertiFlip & $t$ (s) & packing & $t$ (s) & pack.+LP & root & (\\%) & (\\%) & factor \\\\\n"
          "\\midrule\n"
          + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
    # metric LP on P against the checked bounds, where it was run
    lrows = []
    # jobs that never returned: killed with their machine three times (the
    # manager then gives up on them), most likely for lack of memory
    qstat = {}
    qf = os.path.join(COLAB, "queue.json")
    if os.path.exists(qf):
        for j in json.load(open(qf)):
            qstat[j["id"]] = j["status"]
    rama = {}
    rf = os.path.join(RES, "rama_snap_user.csv")
    if os.path.exists(rf):
        for _, r in pd.read_csv(rf).iterrows():
            rama[r.graph] = r
    rama_ok, rama_fail = [], []
    for name, rt, rs, ub, pk, cf in lp_rows:
        def cell(r, tag=None):
            if r is None:
                st = qstat.get(f"{tag}_{name}_0")
                return ("killed", "") if st == "failed" else ("n.r.", "")
            v = f"\\num{{{np.ceil(r['value'] - 1e-6):.0f}}}"
            if not r.get("converged"):
                v = f"({v})"
            return v, f"{r.get('time', float('nan')):.0f}"
        (vt, tt_), (vs, ts) = cell(rt, "ltri"), cell(rs, "lstar")
        best = max(x for x in (pk, cf) if x is not None)
        ra = rama.get(name)
        if ra is None:
            vr, tr = "n.r.", ""
        elif ra.status == "ok":
            # RAMA prints six significant digits
            vr, tr = f"\\num{{{ra.cc_lb:.0f}}}", f"{ra.time:.0f}"
            rama_ok.append((name, ra.cc_lb, best))
        else:
            vr, tr = "failed", ""
            rama_fail.append(name)
        lrows.append(f"{tt(name)} & \\num{{{ub}}} & \\num{{{best}}} & {vt} & {tt_} & {vs} & {ts} & "
                     f"{vr} & {tr} \\\\")
    write("snap_lp", "\\begin{tabular}{lrrrrrrrr}\n\\toprule\n"
          "graph & best & best checked & \\multicolumn{2}{c}{LP$_P$, triangle rows} & "
          "\\multicolumn{2}{c}{LP$_P$, + star rows} & \\multicolumn{2}{c}{RAMA} \\\\\n"
          "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\n"
          " & known & LB & value & time (s) & value & time (s) & value & time (s) \\\\\n\\midrule\n"
          + ("\n".join(lrows) if lrows else "\\multicolumn{9}{c}{\\pending{}} \\\\") +
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
        NUM["numBaseGapLo"] = f"{100 * (min(facs) - 1):.0f}\\%"
        NUM["numBaseGapHi"] = f"{100 * (max(facs) - 1):.0f}\\%"
        NUM["numBaseGapMedian"] = f"{100 * (np.median(facs) - 1):.0f}\\%"
    # the smallest n above which the packing alone wins on every graph
    ns = sorted((cache[x]["n"], x) for x in cf_better + pk_better)
    thr = None
    for n_, x in reversed(ns):
        if x in cf_better:
            break
        thr = n_
    if thr is not None:
        NUM["numPackAllAbove"] = f"\\num{{{thr}}}"
        NUM["numPackAllAboveCount"] = str(sum(1 for n_, x in ns if n_ >= thr))
        NUM["numCfBetterMaxN"] = f"\\num{{{max(n_ for n_, x in ns if x in cf_better)}}}"
    NUM["numSnapTotal"] = str(len(DEV) + len(HELD))
    NUM["numSnapNoBound"] = str(len(DEV) + len(HELD) - len(gaps))
    NUM["numSnapHeldTotal"] = str(len(HELD))
    # over the graphs with a support, the ones compared with CertiFlip
    pt = [PK[(n, 0)]["time"] for n in DEV + HELD if (n, 0) in PK and cache[n]["pairs"] is not None]
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
                NUM["numGrqcTriBelow"] = f"{100 * (ub - v) / ub:.1f}\\%"
                NUM["numGrqcTriTime"] = f"{rt['time']:.0f}"
            if rs is not None:
                v = np.ceil(rs["value"] - 1e-6)
                NUM["numGrqcStar"] = f"\\num{{{v:.0f}}}"
                NUM["numGrqcStarGap"] = f"{100 * (ub - v) / ub:.2f}\\%"
                NUM["numGrqcStarTime"] = f"{rs['time']:.0f}"
                NUM["numGrqcBest"] = f"\\num{{{max(pk, cf)}}}"
    if gaps_cf:
        NUM["numSnapGapCfMedian"] = f"{np.median(gaps_cf):.1f}\\%"
        NUM["numSnapGapPkMedian"] = f"{np.median(gaps_pk):.1f}\\%"
    NUM["numPackLonger"] = str(len(longer))
    NUM["numPackLongerList"] = ", ".join(tt(x) for x in longer)
    NUM["numCfSeedBelowPack"] = ", ".join(tt(x) for x in seed_below) or "none"
    NUM["numRootN"] = str(root_n)
    NUM["numRootTimeout"] = str(len(root_to))
    NUM["numRootTimeoutList"] = ", ".join(tt(x) for x in sorted(root_to)) or "none"
    if root_cmp:
        up = [100 * (r - b) / b for _, r, b, _ in root_cmp if r > b]
        down = [100 * (b - r) / b for _, r, b, _ in root_cmp if r <= b]
        NUM["numRootBetterLo"] = f"{min(up):.1f}\\%" if up else "--"
        NUM["numRootBetterHi"] = f"{max(up):.1f}\\%" if up else "--"
        NUM["numRootBelowMax"] = f"{max(down):.2f}\\%" if down else "--"
        NUM["numRootBelowList"] = ", ".join(tt(x) for x, r, b, _ in root_cmp if r <= b) or "none"
        ts = [t for *_, t in root_cmp]
        for x, r, b_, _ in root_cmp:
            if arch.get(x):
                key_ = "".join(c for c in x if c.isalpha())
                NUM[f"numRootGap{key_}"] = f"{100 * (arch[x] - r) / r:.1f}\\%"
        NUM["numRootTimeLo"] = f"{min(ts):.0f}"
        NUM["numRootTimeHi"] = f"{max(ts) / 3600:.1f}"
        # the closest certified-looking gap the B&B root would give
        best = max(root_cmp, key=lambda x: x[1] / x[2])
        ubb = arch.get(best[0])
        if ubb:
            NUM["numRootBestGraph"] = tt(best[0])
            NUM["numRootBestGap"] = f"{100 * (ubb - best[1]) / best[1]:.2f}\\%"
            NUM["numRootBestOurGap"] = f"{100 * (ubb - best[2]) / best[2]:.2f}\\%"
    NUM["numRootBetter"] = str(len(root_better))
    NUM["numRootBetterList"] = ", ".join(tt(x) for x in root_better) or "none"
    NUM["numPackBetter"] = str(len(pk_better))
    NUM["numCfBetter"] = str(len(cf_better))
    NUM["numCfBetterList"] = ", ".join(tt(x) for x in cf_better)
    # a long packing followed by block LPs (tag lplp)
    if plp_gain:
        NUM["numPlpN"] = str(len(plp_gain))
        NUM["numPlpBetter"] = str(len(plp_better))
        # wins by less than 0.2 % of the previous best bound
        NUM["numPlpSmallWins"] = str(sum(1 for x, d in plp_better
                                         if d < 0.002 * (PL[(x, 0)]["certified"] - d)))
        NUM["numPlpBelow"] = str(len(plp_below))
        NUM["numPlpBelowList"] = ", ".join(tt(x) for x, _ in plp_below) or "none"
        NUM["numPlpBelowMax"] = f"\\num{{{-min(d for _, d in plp_below)}}}" if plp_below else "0"
        gains = [g for _, g in plp_gain]
        NUM["numPlpGainPos"] = str(sum(g > 0 for g in gains))
        NUM["numPlpGainMedian"] = f"{np.median(gains):.2f}\\%"
        top = max(plp_gain, key=lambda x: x[1])
        NUM["numPlpGainMax"] = f"{top[1]:.1f}\\%"
        NUM["numPlpGainMaxGraph"] = tt(top[0])
        pt_ = [PL[(x, 0)]["time"] for x, _ in plp_gain]
        NUM["numPlpTimeMedian"] = f"{np.median(pt_):.0f}"
        NUM["numPlpPackTimeMedian"] = f"{np.median([PL[(x, 0)]['pack_time'] for x, _ in plp_gain]):.0f}"
        NUM["numPlpTimeMax"] = f"{max(pt_):.0f}"
        NUM["numPlpTimeOver"] = str(sum(t > 1.05 * PL[(x, 0)]["T"] for (x, _), t in zip(plp_gain, pt_)))
    # the packing alone with the nominal budget of pack.+LP (tag lpack+eq)
    if eq_rows:
        erows = []
        for role in ("dev", "held-out"):
            for name, r_, eq, t_eq, plp, lb in eq_rows:
                if r_ != role:
                    continue
                pl = PL.get((name, 0))
                d = "--" if plp is None else f"${100 * (eq - plp) / plp:+.2f}$"
                vp = "--" if plp is None else "\\num{%d}" % plp
                tp = "--" if pl is None else "%.0f" % pl["time"]
                top = "yes" if eq == lb and (plp is None or eq > plp) else ""
                erows.append(f"{tt(name)} & {cache[name]['m'] / cache[name]['n']:.1f} & "
                             f"{vp} & {tp} & \\num{{{eq}}} & {t_eq:.0f} & {d} & {top} \\\\")
            erows.append("\\midrule")
        write("snap_eqtime", "\\begin{tabular}{lrrrrrrc}\n\\toprule\n"
              "graph & $m/n$ & \\multicolumn{2}{c}{pack.+LP} & \\multicolumn{2}{c}{packing, "
              "\\SI{1200}{s}} & difference & largest \\\\\n"
              "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n"
              " & & LB & $t$ (s) & LB & $t$ (s) & (\\%) & checked \\\\\n\\midrule\n"
              + "\n".join(erows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
        both = [(n_, eq, plp) for n_, _, eq, _, plp, _ in eq_rows if plp is not None]
        pw = [x for x in both if x[2] > x[1]]
        ew = [x for x in both if x[1] > x[2]]
        NUM["numEqN"] = str(len(both))
        NUM["numEqPlpBetter"] = str(len(pw))
        NUM["numEqPackBetter"] = str(len(ew))
        NUM["numEqPackBetterList"] = ", ".join(tt(x[0]) for x in ew) or "none"
        NUM["numEqPlpBetterMax"] = f"{max(100 * (p_ - e) / e for _, e, p_ in pw):.1f}\\%" if pw else "--"
        NUM["numEqPackBetterMax"] = f"{max(100 * (e - p_) / p_ for _, e, p_ in ew):.1f}\\%" if ew else "--"
        new_ = [n_ for n_, _, eq, _, _, lb in eq_rows if eq == lb and all(
            eq > (v or 0) for v in (checked(PL.get((n_, 0))), checked(PK.get((n_, 0)))))]
        NUM["numEqNew"] = str(len(new_))
        NUM["numEqNewList"] = ", ".join(tt(x) for x in new_) or "none"
        te = [t for *_, t, _, _ in eq_rows]
        NUM["numEqTimeMedian"] = f"{np.median(te):.0f}"
        NUM["numEqTimeMax"] = f"{max(te):.0f}"
        NUM["numEqTimeOver"] = str(sum(t > 1200 for t in te))
    NUM["numRootCert"] = str(len(root_cert))
    NUM["numRootNoCert"] = str(len(root_nocert))
    NUM["numRootNoCertList"] = ", ".join(tt(x) for x in root_nocert) or "none"
    NUM["numRootMemList"] = ", ".join(tt(x) for x in DEV + HELD
                                      if qstat_root.get(f"lkstar_{x}_0") == "failed"
                                      and x not in root_to) or "none"
    # graphs above the support threshold: the packing alone
    NUM["numPackOnlyN"] = str(len(pack_only))
    NUM["numPackOnlyHeld"] = str(sum(1 for x in pack_only if x[1] == "held-out"))
    NUM["numPackOnlyList"] = ", ".join(tt(x[0]) for x in pack_only) or "none"
    if pack_only:
        gs = [100 * (x[4] - x[2]) / x[2] for x in pack_only]
        NUM["numPackOnlyGapLo"] = f"{min(gs):.1f}\\%"
        NUM["numPackOnlyGapHi"] = f"{max(gs):.1f}\\%"
        NUM["numPackOnlyTimeMax"] = f"{max(PK[(x[0], 0)]['time'] for x in pack_only):.0f}"
        NUM["numPackOnlyMaxN"] = f"\\num{{{max(cache[x[0]]['n'] for x in pack_only)}}}"
    NUM["numNoBoundN"] = str(len(no_bound))
    NUM["numNoBoundList"] = ", ".join(tt(x) for x in no_bound) or "none"
    NUM["numSnapAnyBound"] = str(len(gaps) + len(pack_only))
    # RAMA (multicut dual on P, not certified)
    NUM["numRamaOk"] = str(len(rama_ok))
    NUM["numRamaFailed"] = str(len(rama_fail))
    NUM["numRamaFailedList"] = ", ".join(tt(x) for x in rama_fail) or "none"
    if rama_ok:
        fr = [(x, v / b_) for x, v, b_ in rama_ok]
        NUM["numRamaFracMax"] = f"{max(f for _, f in fr):.2f}"
        NUM["numRamaFracMaxGraph"] = tt(max(fr, key=lambda z: z[1])[0])
        # below the better triangle packing, a feasible dual of the same cycle LP
        NUM["numRamaBelowTri"] = str(sum(1 for x, v, _ in rama_ok
                                        if tri.get(x) is not None and v < tri.get(x)))
        NUM["numRamaNearZero"] = str(sum(1 for _, f in fr if f < 0.01))
        NUM["numRamaNearZeroList"] = ", ".join(tt(x) for x, f in fr if f < 0.01)


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
    if len(R) < 200:
        # partial results are not reported
        for k in DEFAULTS:
            if k.startswith(("numPH", "numPD")):
                NUM[k] = DEFAULTS[k]
        write("pace_heur2", "\\begin{tabular}{l}\n\\pending{} "
              f"({len(R)} of 200 runs done)\n\\end{{tabular}}\n")
        fig = os.path.join(ROOT, "paper_mpc", "figures", "profile_pace_heur.pdf")
        if os.path.exists(fig):
            os.remove(fig)
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
    ax.set_xlim(1.0, 1.02)
    ax.set_xticks([1.0, 1.005, 1.01, 1.015, 1.02])
    ax.set_xticklabels(["1.000", "1.005", "1.010", "1.015", "1.020"])
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
    # CPU model of every run reported in this section, by tag
    per = {}
    local = set()
    for tag in ("c2", "c2x", "lpack", "lplp", "ltri", "lstar", "lwarm", "s2h", "s2h150", "s2h60", "s2p"):
        for r in runs(tag).values():
            m = r.get("cpu") or (r.get("machine") or {}).get("Model name")
            if m and "Core(TM)" in m:  # the laptop of results/local, not a Colab machine
                local.add(m)
                continue
            if m:
                models.add(m)
                per.setdefault(tag, {}).setdefault(m, 0)
                per[tag][m] += 1
    if models:
        NUM["numCpuModels"] = "; ".join(sorted(models)).replace("(R)", r"\textsuperscript{\textregistered}")
    NUM["numCpuLocal"] = ("; ".join(sorted(local)).replace("(R)", r"\textsuperscript{\textregistered}")
                          or "none")
    # the 1800 s bound runs (lwarm): machine, overruns
    W = runs("lwarm")
    if W:
        wt = [r["time"] for r in W.values()]
        NUM["numWarmRuns"] = str(len(wt))
        NUM["numWarmLocal"] = str(sum("Core(TM)" in (r.get("cpu") or "") for r in W.values()))
        NUM["numWarmOver"] = str(sum(t > r["T"] for t, r in zip(wt, W.values())))
        NUM["numWarmMax"] = f"{max(wt):.0f}"
    amd = [(tag, c[m], sum(c.values())) for tag, c in per.items() for m in c if "AMD" in m]
    names = {"lpack": "star-packing runs", "s2h60": "head-to-head runs at \\SI{60}{s}",
             "s2h150": "head-to-head runs at \\SI{150}{s}", "s2h": "head-to-head runs at \\SI{600}{s}",
             "s2p": "PACE heuristic-track runs"}
    NUM["numCpuAmd"] = ", ".join(f"{k} of the {n} {names.get(t, t + ' runs')}" for t, k, n in amd) or "none"
    # time of the bound phase of CertiFlip (target: half of the remaining budget)
    # the thresholds are the ones named in the text ("more than 30 s", "300 s")
    for tag, key, thr in (("c2x", "Pace", 30), ("c2", "Snap", 300)):
        R = runs(tag)
        if not R:
            continue
        lt = np.array([r["lb_time"] for r in R.values()])
        NUM[f"numLb{key}Median"] = f"{np.median(lt):.0f}"
        NUM[f"numLb{key}Over"] = str(int((lt > thr).sum()))
        NUM[f"numLb{key}Max"] = f"{lt.max():.0f}"
        NUM[f"numLb{key}Runs"] = str(len(lt))
    # packing runs on PACE (T = 60 s): measured times
    pp = [r for (g, sd), r in runs("lpack").items() if g.startswith("pace")]
    if pp:
        t = np.array([r["time"] for r in pp])
        T = pp[0]["T"]
        NUM["numPackPaceMedian"] = f"{np.median(t):.0f}"
        NUM["numPackPaceOver"] = str(int((t > T).sum()))
        NUM["numPackPaceMax"] = f"{t.max():.0f}"
        NUM["numPackPaceRuns"] = str(len(t))
    pk = runs("lpack")
    NUM["numCertPackOk"] = str(sum(r.get("check") == "ok" for r in pk.values()))
    NUM["numCertPackRuns"] = str(len(pk))
    # line 7 of Algorithm 1 (the LP-rounding start)
    R = runs("c2x")
    used = imp = 0
    for r in R.values():
        h = r.get("history", [])
        for k, e in enumerate(h):
            if e[0] == "lp-seed":
                used += 1
                prev = [x for x in h[:k] if x[0] != "bound"]
                if prev and e[2] < prev[-1][2]:
                    imp += 1
    NUM["numLpSeedPace"] = str(used)
    NUM["numLpSeedPaceImproved"] = str(imp)
    NUM["numLpSeedSnap"] = str(sum(bool(r.get("lp_seed_used")) for r in runs("c2").values()))
    # KaPoCE stopped by SIGTERM in the head-to-head
    for tag, key in (("s2h", "Six"), ("s2h150", "OneFifty"), ("s2h60", "Sixty")):
        R = runs(tag)
        if R:
            NUM[f"numHH{key}Sigterm"] = str(sum(bool(r.get("kapoce_sigterm")) for r in R.values()))
    # the complete recheck against freshly downloaded files: the Colab and the
    # local certificates; every accepted value must equal the in-run check
    import pandas as pd
    parts = []
    for f, d_ in (("recheck.csv", COLAB), ("recheck_local.csv", LOCAL)):
        fp = os.path.join(RES, "mpc", f)
        if os.path.exists(fp):
            x = pd.read_csv(fp)
            x["resdir"] = d_
            parts.append(x)
    if parts:
        d = pd.concat(parts, ignore_index=True)
        lab = d.certificate.str.endswith(".labels.npz")
        NUM["numRecheckN"] = str(len(d))
        NUM["numRecheckCerts"] = str(int((~lab).sum()))
        NUM["numRecheckLabels"] = str(int(lab.sum()))
        NUM["numRecheckOk"] = str(int((d.status == "ok").sum()))
        NUM["numRecheckMax"] = f"{d.seconds.max():.0f}"
        NUM["numRecheckTotal"] = f"{d.seconds.sum() / 3600:.1f}"
        mism = []
        for _, r in d[~lab & (d.status == "ok")].iterrows():
            stem = os.path.basename(r.certificate)[:-4]
            jf = os.path.join(r.resdir, stem + ".json")
            if os.path.exists(jf):
                rec = json.load(open(jf))
                if rec.get("check") == "ok" and int(rec["certified"]) != int(r.certified):
                    mism.append(stem)
        assert not mism, f"recheck differs from the run records: {mism}"
        NUM["numRecheckMatched"] = "all"
        big = d[~lab].sort_values("rows").iloc[-1]
        NUM["numCheckMaxGraph"] = str(big.instance).split("/")[-1].replace(".txt.gz", "").replace(".gr", "")
        NUM["numCheckMaxRows"] = r"\num{%d}" % int(big.rows)
    r = subprocess.run(["git", "-C", ROOT, "log", "-1", "--format=%h", "--",
                        "experiments/check_certificate.py"], capture_output=True, text=True)
    NUM["numCheckerCommit"] = r.stdout.strip() or "unknown"


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
