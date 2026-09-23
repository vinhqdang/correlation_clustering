"""CertiFlip: certified correlation clustering.

Pipeline
  1. Pivot seed (3-approximation in expectation).
  2. Iterated flipping local search with cluster insertions (never increases
     the cost, so the expected guarantee of step 1 is kept).
  3. Anytime lower bound by dual block-coordinate ascent on the triangle +
     star relaxation over the distance-2 support.
  4. Dual-gap guided LNS with exact sub-MIPs (monotone).
The output is a clustering together with a certified lower bound, so the
ratio cost / LB bounds the approximation factor on the instance.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .graph import Graph
from .objective import cost
from .pivot import pivot
from .flip import iterated_flip
from .support import build_support
from .blockdual import BlockDualBound, metis_blocks
from .lns import gap_lns


@dataclass
class CertiFlipResult:
    labels: np.ndarray
    cost: int
    lower_bound: float
    lb_time: float
    history: list = field(default_factory=list)

    @property
    def certified_ratio(self) -> float:
        lb = max(self.lower_bound, 1e-12)
        return self.cost / np.ceil(lb - 1e-6) if self.cost else 1.0


def block_bound(g: Graph, sup=None, time_limit: float = 600.0, block_size: int = 2500,
                max_sweeps: int = 100, stall: float = 1e-4, seed: int = 0,
                block_time: float = 60.0) -> BlockDualBound:
    """Run block sweeps over randomised METIS partitions until the time limit
    or until a sweep improves the bound by less than ``stall`` (relative)."""
    bd = BlockDualBound(g, sup)
    t0 = time.time()
    last = bd.bound()
    for s in range(max_sweeps):
        if time.time() - t0 > time_limit:
            break
        part = metis_blocks(g, block_size, seed + s)
        bd.sweep(part, time_limit=min(block_time, max(1.0, time_limit - (time.time() - t0))))
        cur = bd.bound()
        if s >= 2 and cur - last < stall * max(1.0, abs(cur)):
            break
        last = cur
    return bd


def certiflip(g: Graph, time_limit: float = 300.0, rng=None, lb_fraction: float = 0.5,
              flip_rounds: int = 5, lns_size: int = 40, block_size: int = 2500,
              verbose: bool = False) -> CertiFlipResult:
    rng = np.random.default_rng(rng)
    t0 = time.time()
    hist = []
    lab = pivot(g, rng)
    hist.append(("pivot", time.time() - t0, cost(g, lab)))
    lab = iterated_flip(g, lab, flip_rounds, rng=rng)
    hist.append(("flip", time.time() - t0, cost(g, lab)))
    if verbose:
        print(hist[-1], flush=True)
    t1 = time.time()
    sup = build_support(g)
    remaining = time_limit - (time.time() - t0)
    bd = block_bound(g, sup, time_limit=max(1.0, lb_fraction * remaining),
                     block_size=block_size, seed=int(rng.integers(1 << 30)))
    lb = bd.bound()
    lb_time = time.time() - t1
    hist.append(("bound", time.time() - t0, lb))
    if verbose:
        print(hist[-1], flush=True)
    remaining = time_limit - (time.time() - t0)
    if remaining > 1 and cost(g, lab) > np.ceil(lb - 1e-6):
        lab = gap_lns(g, lab, bd, size=lns_size, iters=1 << 30, time_limit=remaining, rng=rng)
    c = cost(g, lab)
    hist.append(("lns", time.time() - t0, c))
    return CertiFlipResult(lab, c, lb, lb_time, hist)
