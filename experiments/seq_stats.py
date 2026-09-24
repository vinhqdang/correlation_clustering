"""Statistics of a sequential head-to-head experiment (colab_run_seq.py).

Per graph: KaPoCE and PXMem runs with the same seeds, budget and machine.
A run in which KaPoCE returned no valid solution counts as a win for PXMem and
is shown as "-" (KaPoCE did not finish); it is left out of the differences.
Reported per graph: mean costs, the mean paired difference PXMem - KaPoCE with
a 95% bootstrap interval, the Hodges-Lehmann estimate, wins/ties/losses, the
two-sided sign test and its Holm adjustment over all graphs of the tag, and
the median time PXMem needed to reach the final KaPoCE cost.

usage: seq_stats.py TAG [OUT.csv] [RUNS.csv]
       seq_stats.py --pace TAG [OUT.csv]   (one run per instance, pooled test)"""
import csv, glob, json, os, sys
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "colab")


def load(tag):
    runs = []
    for f in sorted(glob.glob(os.path.join(ROOT, f"{tag}_*.json"))):
        r = json.load(open(f))
        if r["tag"] == tag:
            runs.append(r)
    return runs


def sign_p(w, l):
    m, x = w + l, min(w, l)
    return min(1.0, 2 * sum(comb(m, i) for i in range(x + 1)) / 2 ** m) if m else 1.0


def holm(ps):
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - k) * ps[i]))
        adj[i] = run
    return adj


def hodges_lehmann(d):
    d = np.asarray(d, dtype=float)
    w = (d[:, None] + d[None, :]) / 2
    return float(np.median(w[np.triu_indices(len(d))]))


def boot_ci(d, B=20000, seed=0):
    d = np.asarray(d, dtype=float)
    rng = np.random.default_rng(seed)
    means = d[rng.integers(0, len(d), size=(B, len(d)))].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def time_to(r, target):
    for t, c in r.get("hist", []):
        if c <= target:
            return t
    return float("inf")


def per_graph(tag):
    by = {}
    for r in load(tag):
        by.setdefault(r["graph"], []).append(r)
    out = []
    for g, rs in sorted(by.items(), key=lambda kv: kv[1][0]["n"]):
        valid = [r for r in rs if r["kapoce_valid"]]
        nf = len(rs) - len(valid)
        d = [r["ours"] - r["kapoce"] for r in valid]
        w = sum(x < 0 for x in d) + nf
        l = sum(x > 0 for x in d)
        t = len(rs) - w - l
        row = dict(graph=g, n=rs[0]["n"], runs=len(rs), kapoce_unfinished=nf,
                   kapoce_mean=np.mean([r["kapoce"] for r in valid]) if valid else float("nan"),
                   ours_mean=np.mean([r["ours"] for r in rs]),
                   diff_mean=np.mean(d) if d else float("nan"),
                   rel_pct=100 * np.mean([x / r["kapoce"] for x, r in zip(d, valid)]) if d
                   else float("nan"),
                   hl=hodges_lehmann(d) if d else float("nan"),
                   wins=w, ties=t, losses=l, sign_p=sign_p(w, l),
                   ttt=float(np.median([time_to(r, r["kapoce"]) for r in valid])) if valid
                   else float("nan"),
                   kapoce_time=np.mean([r["kapoce_time"] for r in rs]),
                   ours_time=np.mean([r["ours_time"] for r in rs]),
                   kapoce_rss=np.mean([r["kapoce_rss_mb"] for r in rs]),
                   ours_rss=np.mean([r["ours_rss_mb"] for r in rs]))
        row["ci_lo"], row["ci_hi"] = boot_ci(d) if len(d) > 1 else (float("nan"),) * 2
        out.append(row)
    if out:
        for row, a in zip(out, holm([r["sign_p"] for r in out])):
            row["holm_p"] = float(a)
    return out


def pace(tag):
    rs = load(tag)
    valid = [r for r in rs if r["kapoce_valid"]]
    nf = len(rs) - len(valid)
    d = [r["ours"] - r["kapoce"] for r in valid]
    w = sum(x < 0 for x in d) + nf
    l = sum(x > 0 for x in d)
    return dict(instances=len(rs), kapoce_unfinished=nf, wins=w, ties=len(rs) - w - l,
                losses=l, sign_p=sign_p(w, l), sum_ours=sum(r["ours"] for r in valid),
                sum_kapoce=sum(r["kapoce"] for r in valid))


def main(argv):
    if argv[0] == "--pace":
        res = pace(argv[1])
        print(res)
        if len(argv) > 2:
            with open(argv[2], "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(res))
                w.writeheader()
                w.writerow(res)
        return
    tag = argv[0]
    out = per_graph(tag)
    for r in out:
        print(f"{r['graph']:15s} runs={r['runs']:2d} unfin={r['kapoce_unfinished']} "
              f"diff {r['diff_mean']:+8.2f} [{r['ci_lo']:+.1f},{r['ci_hi']:+.1f}] "
              f"({r['rel_pct']:+.4f}%) W/T/L {r['wins']}/{r['ties']}/{r['losses']} "
              f"p={r['sign_p']:.4f} holm={r['holm_p']:.4f} ttt={r['ttt']:.0f}s")
    if len(argv) > 1 and out:
        with open(argv[1], "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
    if len(argv) > 2:
        with open(argv[2], "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["graph", "seed", "ours", "kapoce", "kapoce_time", "ours_time"])
            for r in sorted(load(tag), key=lambda r: (r["graph"], r["seed"])):
                w.writerow([r["graph"], r["seed"], r["ours"],
                            r["kapoce"] if r["kapoce_valid"] else "-",
                            f"{r['kapoce_time']:.1f}", f"{r['ours_time']:.1f}"])


if __name__ == "__main__":
    main(sys.argv[1:])
