"""Ablation of the lower-bound components and scaling measurements.

  python experiments/ablation.py bounds --graphs ca-GrQc ca-HepTh --budget 300 --out results/ablation_lb.csv
  python experiments/ablation.py scaling --out results/scaling.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as D  # noqa: E402
import baselines as B  # noqa: E402

import ccbench as cc  # noqa: E402
from ccbench.support import build_support  # noqa: E402
from ccbench.dual import PackingBound, packing_bound, star_packing, star_packing_ls  # noqa: E402
from ccbench.blockdual import BlockDualBound, add_star_packing  # noqa: E402
from ccbench.certiflip import block_bound, _packing_rate  # noqa: E402
from ccbench.flip import iterated_flip  # noqa: E402

KAPOCE_LB = os.environ.get("KAPOCE_LB_BIN", "/home/user/baselines/kapoce/build_exp/lbounds")


def kapoce_root_bounds(g, timeout=3600):
    if g.n > 25000 or not os.path.exists(KAPOCE_LB):
        return None
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "g.gr")
        B.write_pace(g, f)
        t0 = time.time()
        with open(f) as fin:
            try:
                out = subprocess.run([KAPOCE_LB], stdin=fin, capture_output=True, text=True,
                                     timeout=timeout).stdout
            except subprocess.TimeoutExpired:
                return None
    res = {}
    for line in out.splitlines():
        p = line.split()
        if len(p) == 3:
            res[p[0]] = (int(p[1]), float(p[2]))
    return res


def bounds(graphs, budget, out):
    rows = []
    for name in graphs:
        g = D.load(name)
        lab = iterated_flip(g, cc.pivot(g, 0), 5, rng=0)
        ub = cc.cost(g, lab)

        def rec(method, value, t):
            r = dict(graph=name, n=g.n, m=g.m, method=method, lb=float(value), time=t, ub=ub)
            rows.append(r)
            print(r, flush=True)
        t = time.time(); sup = build_support(g); t_sup = time.time() - t
        rec("support", 0, t_sup)
        t = time.time(); rec("greedy triangle packing", packing_bound(g, sup), time.time() - t)
        t = time.time(); P = PackingBound(g, sup); v = P.run(0.05, 1e-6)
        rec("triangle packing (MWU)", v, time.time() - t)
        bd0 = BlockDualBound(g, sup)
        t = time.time(); v = star_packing(g, sup, pgraph=(bd0.ptr, bd0.idx, bd0.pid),
                                          order="degree-desc", rng=0)[0]
        rec("greedy star packing", v, time.time() - t)
        _, rate = _packing_rate(g, bd0)
        t = time.time()
        v = star_packing_ls(g, sup, pgraph=(bd0.ptr, bd0.idx, bd0.pid),
                            iters=int(rate * 0.3 * budget), seed=0)[0]
        rec("star packing + local search", v, time.time() - t)
        t = time.time(); bd = block_bound(g, sup, time_limit=budget, packing_fraction=0.0)
        rec("block LP (cold start)", bd.bound(), time.time() - t)
        t = time.time(); bd = block_bound(g, sup, time_limit=budget)
        rec("packing + block LP", bd.bound(), time.time() - t)
        t = time.time(); bd = block_bound(g, sup, time_limit=budget, labels=lab)
        rec("packing + block LP + gap blocks", bd.bound(), time.time() - t)
        kb = kapoce_root_bounds(g)
        if kb:
            for k, (v, tt) in kb.items():
                rec(f"KaPoCE root {k} packing", v, tt)
    write(rows, out)


def scaling(out, sizes=(10000, 30000, 100000, 300000, 1000000)):
    rows = []
    for n in sizes:
        g, _ = D.synthetic(f"sparse-{n}-10-2.0-0.7-0")
        rec = {}
        t = time.time(); lab = cc.pivot(g, 0); rec["pivot"] = time.time() - t
        t = time.time(); cc.multilevel(g, lab, 0); rec["louvain"] = time.time() - t
        t = time.time(); iterated_flip(g, lab, 5, rng=0); rec["iterated flip"] = time.time() - t
        t = time.time(); sup = build_support(g); rec["support"] = time.time() - t
        bd = BlockDualBound(g, sup)
        t = time.time(); star_packing_ls(g, sup, pgraph=(bd.ptr, bd.idx, bd.pid), iters=10 ** 6)
        rec["star packing (1e6 iterations)"] = time.time() - t
        for k, v in rec.items():
            r = dict(n=g.n, m=g.m, pairs=sup.npairs, step=k, time=v)
            rows.append(r)
            print(r, flush=True)
    write(rows, out)


def write(rows, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["bounds", "scaling"])
    ap.add_argument("--graphs", nargs="*", default=["ca-GrQc", "ca-HepTh", "BitcoinAlpha+",
                                                     "BitcoinOTC+", "ca-CondMat"])
    ap.add_argument("--budget", type=float, default=300)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.what == "bounds":
        bounds(a.graphs, a.budget, a.out)
    else:
        scaling(a.out)


if __name__ == "__main__":
    main()
