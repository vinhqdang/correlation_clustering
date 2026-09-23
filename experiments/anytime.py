"""Anytime behaviour of the lower bound: bound value over time."""
from __future__ import annotations

import csv
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as D  # noqa: E402

from ccbench.certiflip import block_bound  # noqa: E402
from ccbench.flip import iterated_flip  # noqa: E402
import ccbench as cc  # noqa: E402


def main(graphs, budget, out):
    rows = []
    for name in graphs:
        g = D.load(name)
        lab = iterated_flip(g, cc.pivot(g, 0), 5, rng=0, time_limit=120)
        for variant, kw in [("block LP (cold start)", dict(packing_fraction=0.0)),
                            ("packing + block LP + gap blocks", dict(labels=lab))]:
            t0 = time.time()
            bd = block_bound(g, time_limit=budget, **kw)
            for t, v in bd.history:
                rows.append(dict(graph=name, variant=variant, time=t, lb=v))
            print(name, variant, bd.bound(), time.time() - t0, flush=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["graph", "variant", "time", "lb"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[3:], float(sys.argv[2]), sys.argv[1])
