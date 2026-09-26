"""Root star-packing and P3-packing bounds of the KaPoCE branch-and-bound
(its `lbounds` program, dense implementation) on the SNAP graphs small enough
for it; one CSV row per graph and bound, resumable.

  python experiments/kapoce_root_snap.py OUT.csv GRAPH...
"""
import csv
import os
import platform
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as D  # noqa: E402
from ablation import kapoce_root_bounds  # noqa: E402


def cpu():
    try:
        out = subprocess.run(["lscpu"], capture_output=True, text=True).stdout
        return next(l.split(":", 1)[1].strip() for l in out.splitlines() if "Model name" in l)
    except Exception:
        return platform.processor()


def main(out, graphs):
    done = set()
    if os.path.exists(out):
        with open(out) as fh:
            done = {r["graph"] for r in csv.DictReader(fh)}
    for name in graphs:
        if name in done:
            continue
        g = D.load(name)
        res = kapoce_root_bounds(g, timeout=3600)
        rows = [dict(graph=name, n=g.n, m=g.m, bound=k, lb=v, time=t, cpu=cpu())
                for k, (v, t) in (res or {"timeout": (-1, 3600.0)}).items()]
        new = not os.path.exists(out)
        with open(out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["graph", "n", "m", "bound", "lb", "time", "cpu"])
            if new:
                w.writeheader()
            w.writerows(rows)
        print(rows, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
