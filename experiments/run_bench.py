"""Benchmark runner.

Usage examples
  python experiments/run_bench.py --suite pace-exact --algos all --out results/pace_exact.csv
  python experiments/run_bench.py --instances ca-GrQc ca-HepTh --algos flip,kapoce --budget 60

Every row of the output CSV is one (instance, algorithm, seed) run with the
objective value and wall-clock time; lower bounds are rows with algo
starting with "lb:".
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import time
import traceback
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as D  # noqa: E402
import baselines as B  # noqa: E402

import ccbench as cc  # noqa: E402
from ccbench.flip import iterated_flip  # noqa: E402
from ccbench.insertion import insertion_search  # noqa: E402

ALGOS = ["pivot", "pivot50", "vote", "louvain", "leiden", "mfp", "ins", "flip", "certiflip",
         "kapoce"]
BOUNDS = ["lb:greedy", "lb:tri-mwu"]


INSTANCE = [None]  # name of the instance being solved (for certificate files)


def run_algo(name, g, seed, budget):
    t0 = time.time()
    extra = {}
    if name == "pivot":
        lab = cc.pivot(g, seed)
    elif name == "pivot50":
        lab = cc.best_pivot(g, 50, seed)
    elif name == "vote":
        lab = cc.local_search(g, np.arange(g.n), seed)
    elif name == "louvain":
        lab = cc.multilevel(g, None, seed)
    elif name == "leiden":
        lab = B.leiden_cpm(g, seed)
    elif name == "mfp":
        lab = B.match_flip_pivot(g, seed)
    elif name == "ins":
        lab = insertion_search(g, cc.pivot(g, seed), seed)
    elif name == "flip":
        lab = iterated_flip(g, cc.pivot(g, seed), 5, rng=seed)
    elif name == "certiflip":
        from ccbench.certiflip import certiflip
        cert = None
        if os.environ.get("CERT_DIR") and INSTANCE[0]:
            os.makedirs(os.environ["CERT_DIR"], exist_ok=True)
            cert = os.path.join(os.environ["CERT_DIR"],
                                INSTANCE[0].replace("/", "_") + f"_s{seed}.npz")
        res = certiflip(g, time_limit=budget, rng=seed, cert_path=cert)
        lab = res.labels
        extra = {"lb": res.lower_bound, "lb_time": res.lb_time}
    elif name == "kapoce":
        lab, _ = B.kapoce(g, budget)
    else:
        raise KeyError(name)
    return lab, time.time() - t0, extra


def run_bound(name, g, budget):
    from ccbench.support import build_support
    from ccbench.dual import PackingBound, packing_bound
    t0 = time.time()
    deg = g.degrees.astype(np.float64)
    if min(g.m + float((deg * (deg - 1) / 2).sum()), g.n * (g.n - 1) / 2) > 3e7:
        raise MemoryError("distance-2 support too large for this machine")
    sup = build_support(g)
    if name == "lb:greedy":
        v = packing_bound(g, sup)
    elif name == "lb:tri-mwu":
        P = PackingBound(g, sup)
        v = P.run(0.05, 1e-6)
    elif name == "lb:block-star":
        from ccbench.certiflip import block_bound
        bd = block_bound(g, sup, time_limit=budget, labels=None)
        v = bd.bound()
    else:
        raise KeyError(name)
    return float(v), time.time() - t0


def job(args):
    inst, algo, seed, budget = args
    try:
        g = D.load(inst)
        INSTANCE[0] = inst
        if algo.startswith("lb:"):
            v, t = run_bound(algo, g, budget)
            return dict(instance=inst, n=g.n, m=g.m, algo=algo, seed=seed, cost="", lb=v, time=t)
        lab, t, extra = run_algo(algo, g, seed, budget)
        c = cc.cost(g, lab)
        row = dict(instance=inst, n=g.n, m=g.m, algo=algo, seed=seed, cost=c,
                   lb=extra.get("lb", ""), time=t)
        return row
    except Exception as exc:  # keep the benchmark going
        traceback.print_exc()
        return dict(instance=inst, n="", m="", algo=algo, seed=seed, cost="", lb="",
                    time="", error=repr(exc))


def suite(name):
    if name == "pace-exact":
        files = sorted(glob.glob(os.path.join(D.ROOT, "pace", "exact", "*.gr")))
        return ["pace-exact/" + os.path.basename(f) for f in files]
    if name == "pace-heur":
        files = sorted(glob.glob(os.path.join(D.ROOT, "pace", "heur", "*.gr")))
        return ["pace-heur/" + os.path.basename(f) for f in files]
    if name == "snap":
        return list(D.REAL)
    if name == "synthetic":
        out = []
        for flip in [0.05, 0.1, 0.2, 0.3]:
            for s in range(3):
                out.append(f"sbm-400-8-{1 - flip:.2f}-{flip / 4:.4f}-{s}")
        for n in [10000, 100000, 1000000]:
            out.append(f"sparse-{n}-10-2.0-0.7-0")
        return out
    raise KeyError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite")
    ap.add_argument("--instances", nargs="*")
    ap.add_argument("--algos", default="all")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--budget", type=float, default=60.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-n", type=int, default=None)
    a = ap.parse_args()
    insts = a.instances or suite(a.suite)
    algos = ALGOS + BOUNDS if a.algos == "all" else a.algos.split(",")
    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            for r in csv.DictReader(fh):
                done.add((r["instance"], r["algo"], int(r["seed"])))
    jobs = []
    for inst in insts:
        for algo in algos:
            seeds = [0] if algo.startswith("lb:") or algo in ("kapoce", "certiflip") else range(a.seeds)
            for s in seeds:
                if (inst, algo, s) not in done:
                    jobs.append((inst, algo, s, a.budget))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fields = ["instance", "n", "m", "algo", "seed", "cost", "lb", "time", "error"]
    new = not os.path.exists(a.out)
    with open(a.out, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if new:
            w.writeheader()
        with Pool(a.workers, maxtasksperchild=20) as pool:
            for row in pool.imap_unordered(job, jobs):
                w.writerow(row)
                fh.flush()
                print(row, flush=True)


if __name__ == "__main__":
    main()
