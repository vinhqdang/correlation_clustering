"""Multicut dual bound of RAMA on the distance-two support P, as a lower
bound for correlation clustering (Section 2, multicut paragraph).

For every clustering C, cost(C) >= MC(C) + |N2|, where MC is the multicut
objective on the graph (V, P) with weight +1 on the edges and -1 on the pairs
at distance two.  So the RAMA lower bound for this multicut instance plus |N2|
is a lower bound for correlation clustering (in floating point, not
certified).  One CSV row per graph, resumable.

  python experiments/rama/rama_snap.py RAMA_BINARY OUT.csv [GRAPH...]

Without graphs, all SNAP graphs whose support was computed are run, smallest
first.  RAMA_TIMEOUT (seconds, default 7200) limits each run.
"""
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments")]
import datasets as D  # noqa: E402
from ccbench.support import build_support  # noqa: E402

TIMEOUT = int(os.environ.get("RAMA_TIMEOUT", "7200"))
KEYS = ["graph", "n", "m", "pairs", "n2", "mc_lb", "cc_lb", "time", "status", "solver", "gpu"]


def gpu_name():
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=30).stdout.strip()
        return out.splitlines()[0] if out else ""
    except Exception:
        return ""


def run(binary, name):
    g = D.load(name)
    sup = build_support(g)
    cost = np.where(np.arange(sup.npairs) < sup.npos, 1, -1)
    row = dict(graph=name, n=g.n, m=g.m, pairs=sup.npairs, n2=sup.nneg,
               solver=os.path.basename(binary), gpu=gpu_name())
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "mc.txt")
        with open(f, "w") as fh:
            fh.write("MULTICUT\n")
            np.savetxt(fh, np.c_[sup.pu, sup.pv, cost], fmt="%d")
        del sup
        t0 = time.time()
        try:
            p = subprocess.run([binary, f, "--only_lb"], capture_output=True, text=True,
                               timeout=TIMEOUT)
            out = p.stdout + p.stderr
            m = re.findall(r"final lower bound:\s*(-?[0-9.eE+-]+)", out)
            if p.returncode != 0 or not m:
                row["status"] = f"failed (exit {p.returncode}): {out[-300:]!r}"
            else:
                row["mc_lb"] = float(m[-1])
                row["cc_lb"] = row["mc_lb"] + row["n2"]
                row["status"] = "ok"
        except subprocess.TimeoutExpired:
            row["status"] = "timeout"
        row["time"] = round(time.time() - t0, 1)
    return row


def main(binary, out, graphs):
    if not graphs:
        cache = json.load(open(os.path.join(ROOT, "results", "mpc", "instances.json")))
        graphs = sorted((k for k, v in cache.items() if v.get("pairs")),
                        key=lambda k: cache[k]["pairs"])
    done = set()
    if os.path.exists(out):
        done = {r["graph"] for r in csv.DictReader(open(out))}
    for name in graphs:
        if name in done:
            continue
        row = run(binary, name)
        new = not os.path.exists(out)
        with open(out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(row)
        print({k: row.get(k) for k in ("graph", "pairs", "cc_lb", "time", "status")}, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
