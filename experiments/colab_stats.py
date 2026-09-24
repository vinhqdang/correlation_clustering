"""Per-graph statistics of a head-to-head experiment against KaPoCE on the same
machine and budget: mean relative difference, wins/ties/losses and the two-sided
sign test.

A run in which KaPoCE did not return a valid solution within the time limit
(it keeps a solution only after its first full repetition; when stopped
earlier it outputs an empty clustering) counts as a win for our solver, as in
the PACE rules.  Such runs are reported separately ("KaPoCE unfinished", shown
as "-" in the per-run tables) and are left out of the mean relative difference.

usage: colab_stats.py TAG [OUT.csv] [RUNS.csv]"""
import csv, glob, json, os, sys
from math import comb

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "colab")
tag = sys.argv[1]
rows, runs = {}, []
for f in sorted(glob.glob(os.path.join(ROOT, f"{tag}_*.json"))):
    r = json.load(open(f))
    if r["tag"] != tag:
        continue
    failed = r["kapoce"] > 1.5 * r["ours"]  # empty clustering: KaPoCE did not finish
    rows.setdefault(r["graph"], []).append((r["ours"], None if failed else r["kapoce"]))
    runs.append((r["graph"], r["seed"], r["ours"], "-" if failed else r["kapoce"]))
out = []
for g, v in sorted(rows.items()):
    valid = [(o, k) for o, k in v if k is not None]
    nf = len(v) - len(valid)
    rel = [100 * (o - k) / k for o, k in valid]
    w = sum(o < k for o, k in valid) + nf
    l = sum(o > k for o, k in valid)
    t = len(v) - w - l
    m, x = w + l, min(w, l)
    p = min(1.0, 2 * sum(comb(m, i) for i in range(x + 1)) / 2 ** m) if m else 1.0
    mean = sum(rel) / len(rel) if rel else float("nan")
    out.append((g, len(v), nf, mean, w, t, l, p))
    print(f"{g:15s} n={len(v):2d} KaPoCE unfinished={nf}  mean {mean:+.4f}%  "
          f"W/T/L {w}/{t}/{l}  sign-test p={p:.3f}")
if len(sys.argv) > 2:
    with open(sys.argv[2], "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["graph", "runs", "kapoce_unfinished", "mean_rel_diff_pct_valid_runs",
                     "wins", "ties", "losses", "sign_p"])
        for r in out:
            wr.writerow([r[0], r[1], r[2], f"{r[3]:.5f}", r[4], r[5], r[6], f"{r[7]:.4f}"])
if len(sys.argv) > 3:
    with open(sys.argv[3], "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["graph", "seed", "ours", "kapoce"])
        for r in sorted(runs):
            wr.writerow(r)
