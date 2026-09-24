"""Summarise head-to-head results in results/colab: ours - KaPoCE per graph and tag."""
import glob, json, os, sys
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "colab")
rows = defaultdict(dict)
invalid = []
for f in sorted(glob.glob(os.path.join(ROOT, "*.json"))):
    if os.path.basename(f) in ("queue.json", "fleet_state.json"):
        continue
    r = json.load(open(f))
    if r["kapoce"] > 1.5 * r["ours"]:
        # KaPoCE returned no solution within the time limit: a win for us
        invalid.append(os.path.basename(f))
        rows[r["tag"]].setdefault(r["graph"], []).append(None)
        continue
    rows[r["tag"]].setdefault(r["graph"], []).append(r["ours"] - r["kapoce"])
tags = [t for t in rows if not sys.argv[1:] or any(t.startswith(a) for a in sys.argv[1:])]
graphs = sorted({g for t in tags for g in rows[t]})
for t in tags:
    d = rows[t]
    wins = sum(sum(v is None or v < 0 for v in d[g]) for g in d)
    ties = sum(sum(v == 0 for v in d[g]) for g in d)
    loss = sum(sum(v is not None and v > 0 for v in d[g]) for g in d)
    print(f"{t:45s} W/T/L {wins}/{ties}/{loss}")
    print("    " + "  ".join(f"{g}:{','.join('-' if v is None else f'{v:+d}' for v in d[g])}"
                              for g in graphs if g in d))
if invalid:
    print("KaPoCE unfinished (shown as -, counted as wins):", ", ".join(invalid))
