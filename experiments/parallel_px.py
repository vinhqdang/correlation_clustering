"""Parallel runs combined by partition crossover, and PX as a post-processor.

(a) Four independent pxmem runs (different seeds) are started at the same time
    on four cores.  For k = 1, 2, 4 we compare the best of k runs with the
    partition-crossover combination of the same k runs (same wall-clock time).
(b) The KaPoCE solution (same wall-clock budget, one core) is combined with our
    solutions by partition crossover.

usage: parallel_px.py OUT.csv T SOLDIR GRAPH...
"""
import csv
import itertools
import os
import sys
import time
from multiprocessing import Process

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
import datasets as D  # noqa: E402
import baselines as B  # noqa: E402
from ccbench.memetic import pxmem, partition_crossover  # noqa: E402
from ccbench.reduce import unit  # noqa: E402
from ccbench.localsearch import _compact  # noqa: E402


def run_ours(name, T, seed, path):
    g = D.load(name)
    lab, c = pxmem(g, T, rng=seed)
    np.save(path, lab)


def run_kapoce(name, T, path):
    g = D.load(name)
    lab, _ = B.kapoce(g, T)
    np.save(path, np.asarray(lab))


def px_fold(wg, sols):
    child = _compact(sols[0])[0]
    for s in sols[1:]:
        child = partition_crossover(wg, child, _compact(s)[0])[0]
    return child


def main(out, T, soldir, graphs):
    os.makedirs(soldir, exist_ok=True)
    new = not os.path.exists(out)
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        if new:
            wr.writerow(["graph", "T", "k", "method", "cost"])
        for name in graphs:
            ours = [os.path.join(soldir, f"{name}_ours_{s}.npy") for s in range(4)]
            kap = os.path.join(soldir, f"{name}_kapoce.npy")
            if not all(os.path.exists(p) for p in ours):
                ps = [Process(target=run_ours, args=(name, T, s, ours[s])) for s in range(4)]
                for p in ps:
                    p.start()
                for p in ps:
                    p.join()
            if not os.path.exists(kap):
                run_kapoce(name, T, kap)
            g = D.load(name)
            wg = unit(g)
            sols = [np.load(p) for p in ours]
            costs = [wg.cost(_compact(s)[0]) for s in sols]
            k_sol = _compact(np.load(kap))[0]
            ck = wg.cost(k_sol)
            rows = []
            for k in (1, 2, 4):
                subsets = list(itertools.combinations(range(4), k))
                best = np.mean([min(costs[i] for i in S) for S in subsets])
                px = np.mean([wg.cost(px_fold(wg, [sols[i] for i in S])) for S in subsets])
                rows += [(k, "best-of-k", best), (k, "PX-of-k", px)]
            rows.append((1, "KaPoCE", ck))
            rows.append((1, "PX(KaPoCE, ours)",
                         np.mean([wg.cost(px_fold(wg, [k_sol, s])) for s in sols])))
            rows.append((4, "PX(KaPoCE, 4 x ours)", wg.cost(px_fold(wg, [k_sol] + sols))))
            for k, m, c in rows:
                wr.writerow([name, T, k, m, f"{c:.2f}"])
                print(name, k, m, round(float(c), 2), flush=True)
            fh.flush()


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4:])
