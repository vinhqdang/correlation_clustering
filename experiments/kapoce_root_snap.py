"""Root star-packing and P3-packing bounds of the KaPoCE branch-and-bound
(its `lbounds` program, dense implementation) on the SNAP graphs small enough
for it; one CSV row per graph and bound, resumable.

  python experiments/kapoce_root_snap.py OUT.csv GRAPH...

The program is taken from KAPOCE_LB_BIN (see experiments/kapoce/build_lbounds.sh);
the time limit per graph is KAPOCE_LB_TIMEOUT seconds (default 14400).
"""
import csv
import os
import platform
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as D  # noqa: E402
from ablation import kapoce_root_bounds  # noqa: E402

TIMEOUT = int(os.environ.get("KAPOCE_LB_TIMEOUT", "14400"))


def cpu():
    try:
        out = subprocess.run(["lscpu"], capture_output=True, text=True).stdout
        return next(l.split(":", 1)[1].strip() for l in out.splitlines() if "Model name" in l)
    except Exception:
        return platform.processor()


def main(out, graphs):
    from ablation import KAPOCE_LB
    if not os.path.exists(KAPOCE_LB):
        sys.exit(f"{KAPOCE_LB} not found; build it with experiments/kapoce/build_lbounds.sh "
                 "and set KAPOCE_LB_BIN")
    done = set()
    if os.path.exists(out):
        with open(out) as fh:
            done = {r["graph"] for r in csv.DictReader(fh)}
    import time
    for name in graphs:
        if name in done:
            print(f"[{time.strftime('%H:%M:%S')}] {name}: already done", flush=True)
            continue
        print(f"[{time.strftime('%H:%M:%S')}] {name}: running (the largest graphs take up to "
              "an hour or more) ...", flush=True)
        g = D.load(name)
        res = kapoce_root_bounds(g, timeout=TIMEOUT)
        rows = [dict(graph=name, n=g.n, m=g.m, bound=k, lb=v, time=t, cpu=cpu(), limit=TIMEOUT)
                for k, (v, t) in (res or {"timeout": (-1, float(TIMEOUT))}).items()]
        new = not os.path.exists(out)
        with open(out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["graph", "n", "m", "bound", "lb", "time", "cpu",
                                               "limit"], extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerows(rows)
        print(f"[{time.strftime('%H:%M:%S')}]     {name}: " +
              ", ".join(f"{r['bound']} {r['lb']} ({r['time']:.0f} s)" for r in rows), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
