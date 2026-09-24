"""Per-graph statistics of a head-to-head experiment against KaPoCE:
mean relative difference, wins/ties/losses and the two-sided sign test.
Runs where KaPoCE returned no valid solution are excluded.

usage: colab_stats.py TAG [OUT.csv]"""
import csv, glob, json, os, sys
from math import comb

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "colab")
tag = sys.argv[1]
rows = {}
for f in glob.glob(os.path.join(ROOT, f"{tag}_*.json")):
    r = json.load(open(f))
    if r["tag"] != tag or r["kapoce"] > 1.5 * r["ours"]:
        continue
    rows.setdefault(r["graph"], []).append((r["ours"], r["kapoce"]))
out = []
for g, v in sorted(rows.items()):
    rel = [100 * (o - k) / k for o, k in v]
    w = sum(o < k for o, k in v)
    l = sum(o > k for o, k in v)
    m, x = w + l, min(w, l)
    p = min(1.0, 2 * sum(comb(m, i) for i in range(x + 1)) / 2 ** m) if m else 1.0
    out.append((g, len(v), sum(rel) / len(rel), w, len(v) - w - l, l, p))
    print(f"{g:15s} n={len(v):2d} mean {sum(rel) / len(rel):+.4f}%  W/T/L {w}/{len(v) - w - l}/{l}"
          f"  sign-test p={p:.3f}")
if len(sys.argv) > 2:
    with open(sys.argv[2], "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["graph", "runs", "mean_rel_diff_pct", "wins", "ties", "losses", "sign_p"])
        for r in out:
            wr.writerow([r[0], r[1], f"{r[2]:.5f}", r[3], r[4], r[5], f"{r[6]:.4f}"])
