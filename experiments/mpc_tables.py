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
    tri = pd.concat([pd.read_csv(os.path.join(RES, f)) for f in
                     ("snap.csv", "snap_heldout_bounds.csv")
                     if os.path.exists(os.path.join(RES, f))])
    tri = tri[tri.algo.isin(["lb:greedy", "lb:tri-mwu"]) & tri.lb.notna()]
    tri = tri.groupby("instance")["lb"].max().apply(lambda x: np.ceil(x - 1e-6))
    rows, gaps, facs, pk_share = [], [], [], []
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
            lb = max(lbs)
            spread = 100 * (max(lbs) - min(lbs)) / lb
            gap_a = 100 * (ub_a - lb) / lb
            gaps.append(gap_a)
            pk = checked(PK.get((name, 0)))
            pk_gap = "--" if pk is None else f"{100 * (ub - pk) / pk:.2f}"
            if pk is not None:
                pk_share.append((lb - pk) / max(1, ub - pk))
            tb = tri.get(name)
            fac = None if tb is None or not np.isfinite(tb) else ub / tb
            if fac is not None:
                facs.append(fac)
            t_lb = np.median([C2[(name, sd)]["lb_time"] for sd in range(3) if (name, sd) in C2])
            gap_b = 100 * (ub - lb) / lb
            TIMES[name] = t_lb
            rows.append(f"{tt(name)} & \\num{{{ub_a}}} & \\num{{{lb}}} & {spread:.2f} & "
                        f"{gap_a:.2f} & {gap_b:.2f} & {pk_gap} & {fmt_pct(fac, 2) if fac else '--'} \\\\")
        rows.append("\\midrule")
    write("snap_bounds", "\\begin{tabular}{lrrrrrrr}\n\\toprule\n"
          "graph & UB & LB & spread & gap & gap$^*$ & star & tri.\\\\\n"
          " & & & (\\%) & (\\%) & (\\%) & (\\%) & factor \\\\\n\\midrule\n" + "\n".join(rows[:-1]) +
          "\n\\bottomrule\n\\end{tabular}\n")
    if gaps:
        NUM["numSnapGapLo"] = f"{min(gaps):.2f}\\%"
        NUM["numSnapGapHi"] = f"{max(gaps):.1f}\\%"
    if facs:
        NUM["numBaseFacLo"] = f"{min(facs):.2f}"
        NUM["numBaseFacHi"] = f"{max(facs):.2f}"
    if pk_share:
        NUM["numLpShareLo"] = f"{100 * min(pk_share):.0f}\\%"
        NUM["numLpShareHi"] = f"{100 * max(pk_share):.0f}\\%"


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
    pace_exact()
    snap()
    static_numbers()
    write_numbers()
