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
from .blockdual import (BlockDualBound, metis_blocks, gap_blocks, add_star_packing,
                        add_triangle_packing)
from .dual import star_packing_ls, star_packing_p3x, PackingBound
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
                max_sweeps: int = 100, stall: float = 1e-3, seed: int = 0,
                block_time: float = 120.0, labels: np.ndarray | None = None,
                subgraphs: bool = True, bd: BlockDualBound | None = None,
                packing_fraction: float = 0.3) -> BlockDualBound:
    """Anytime lower bound.

    1. a star packing improved by ruin-and-recreate local search (a feasible
       dual: y = 1 on each packed star inequality), using about
       ``packing_fraction`` of the time;
    2. block sweeps over randomised METIS partitions, which re-optimise the
       stored rows inside each block, until a sweep improves the bound by less
       than ``stall`` (relative);
    3. if a clustering is given, dual-gap guided blocks for the remaining time.
    """
    bd = bd if bd is not None else BlockDualBound(g, sup)
    t0 = time.time()
    kw = dict(subgraphs=subgraphs)
    if packing_fraction > 0:
        # candidate feasible duals; the best one is installed
        pg = (bd.ptr, bd.idx, bd.pid)
        cands = []
        _, rate = _packing_rate(g, bd)
        iters = int(max(1000, rate * packing_fraction * time_limit))
        v, rp, rpairs, rk = star_packing_ls(g, bd.sup, pgraph=pg, iters=iters, seed=seed)
        cands.append((v, "stars", (rp, rpairs, rk)))
        density = 2.0 * g.m / max(1.0, g.n * (g.n - 1.0))
        if density > 0.02:
            v, rp, rpairs, rk = star_packing_p3x(g, bd.sup, pgraph=pg, iters=0, seed=seed)
            cands.append((v, "stars", (rp, rpairs, rk)))
            P = PackingBound(g, bd.sup)
            v = P.run(0.1, 1e-4)
            cands.append((v, "triangles", (P.ta, P.tb, P.tc, P.y)))
        v, kind, data = max(cands, key=lambda c: c[0])
        if kind == "stars":
            add_star_packing(bd, *data)
        else:
            add_triangle_packing(bd, *data)
        bd.history.append((time.time() - bd.t0, bd.bound()))
    if g.n <= block_size:
        bd.solve_block(np.arange(g.n), time_limit=max(1.0, time_limit), max_rounds=1000, **kw)
        bd.history.append((time.time() - bd.t0, bd.bound()))
        return bd
    last = bd.bound()
    for s in range(max_sweeps):
        if time.time() - t0 > time_limit:
            break
        part = metis_blocks(g, block_size, seed + s)
        bd.sweep(part, time_limit=min(block_time, max(1.0, time_limit - (time.time() - t0))),
                 **kw)
        cur = bd.bound()
        if s >= 1 and cur - last < stall * max(1.0, abs(cur)):
            break
        last = cur
    rest = time_limit - (time.time() - t0)
    if labels is not None and rest > 1:
        gap_blocks(bd, labels, block_size, rest, rng=seed,
                   block_kw=dict(time_limit=block_time, **kw))
    return bd


def lp_seed(g: Graph, bd: BlockDualBound, lab: np.ndarray, rng, rounds: int = 10,
            flip_rounds: int = 2):
    """CMSY rounding of the block-LP primal when one block covered the whole
    graph and separation converged, i.e. x is a feasible point of the metric
    LP on the distance-2 support with cost(x) = LB.  Returns (labels, used)."""
    from .lp import lp_pivot
    lb = getattr(bd, "last_block", None)
    if lb is None or not (lb["whole"] and lb["converged"]):
        return lab, False
    x = np.ones(bd.sup.npairs)
    x[lb["gids"]] = lb["x"]
    best, bc = lab, cost(g, lab)
    for _ in range(rounds):
        cand = lp_pivot(g, bd.sup, x, "cmsy", rng, pgraph=(bd.ptr, bd.idx, bd.pid))
        cand = iterated_flip(g, cand, flip_rounds, rng=rng)
        c = cost(g, cand)
        if c < bc:
            best, bc = cand, c
    return best, True


def _packing_rate(g, bd, pilot: int = 20000):
    """Iterations per second of the packing local search on this graph."""
    t = time.time()
    star_packing_ls(g, bd.sup, pgraph=(bd.ptr, bd.idx, bd.pid), iters=pilot, seed=0)
    dt = max(time.time() - t, 1e-3)
    return pilot, pilot / dt


def certiflip(g: Graph, time_limit: float = 300.0, rng=None, lb_fraction: float = 0.5,
              flip_rounds: int = 5, lns_size: int = 40, block_size: int = 2500,
              use_lp_seed: bool = True, max_pairs: float = 3e7,
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
    deg = g.degrees.astype(np.float64)
    est_pairs = min(g.m + float((deg * (deg - 1) / 2).sum()), g.n * (g.n - 1) / 2)
    if est_pairs > max_pairs:
        # support too large for this machine: primal only (disagreement-guided LNS)
        remaining = time_limit - (time.time() - t0)
        if remaining > 1:
            lab = gap_lns(g, lab, None, size=lns_size, iters=1 << 30, time_limit=remaining,
                          rng=rng)
        c = cost(g, lab)
        hist.append(("lns", time.time() - t0, c))
        return CertiFlipResult(lab, c, float("nan"), 0.0, hist)
    sup = build_support(g)
    remaining = time_limit - (time.time() - t0)
    bd = block_bound(g, sup, time_limit=max(1.0, lb_fraction * remaining),
                     block_size=block_size, seed=int(rng.integers(1 << 30)), labels=lab)
    lb = bd.bound()
    lb_time = time.time() - t1
    hist.append(("bound", time.time() - t0, lb))
    if verbose:
        print(hist[-1], flush=True)
    if use_lp_seed:
        lab, used = lp_seed(g, bd, lab, rng)
        if used:
            hist.append(("lp-seed", time.time() - t0, cost(g, lab)))
    remaining = time_limit - (time.time() - t0)
    if remaining > 1 and cost(g, lab) > np.ceil(lb - 1e-6):
        lab = gap_lns(g, lab, bd, size=lns_size, iters=1 << 30, time_limit=remaining, rng=rng)
    c = cost(g, lab)
    hist.append(("lns", time.time() - t0, c))
    return CertiFlipResult(lab, c, lb, lb_time, hist)
