"""Iterated flipping local search (practical version of the 2 - 2/13 scheme of
Cohen-Addad, Lolck, Pilipczuk, Thorup, Yan and Zhang, STOC 2024).

Each round adds beta = 1/2 to the weight of the positive edges cut by the
previous solution, re-optimises the reweighted objective by local search, flips
again, re-optimises, and combines the last three solutions with the 3-way
pivot.  All candidates are compared under the original objective.  Weights
are kept integral by doubling: positive edges cost 2, 3 or 4 and negative
pairs cost kappa = 2.
"""
from __future__ import annotations

import heapq

import numpy as np

from .graph import Graph
from .insertion import insertion_search
from .localsearch import _compact, multilevel
from .objective import cost


def cut_indicator(g: Graph, labels: np.ndarray) -> np.ndarray:
    """1 for CSR slots of positive edges cut by ``labels``."""
    src = np.repeat(np.arange(g.n), g.degrees)
    return (labels[src] != labels[g.indices]).astype(np.int64)


def pivot3(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """3-way pivot: repeatedly take the most frequent label triple among
    unassigned vertices and cluster every unassigned vertex agreeing with it on
    at least two of the three clusterings."""
    n = len(a)
    a = _compact(a)[0]
    b = _compact(b)[0]
    c = _compact(c)[0]
    key = lambda x, y: x * (n + 1) + y
    groups = {}
    for name, k in (("ab", key(a, b)), ("ac", key(a, c)), ("bc", key(b, c))):
        order = np.argsort(k, kind="stable")
        ks = k[order]
        cuts = np.flatnonzero(np.diff(ks)) + 1
        starts = np.concatenate([[0], cuts])
        ends = np.concatenate([cuts, [n]])
        groups[name] = {int(ks[s]): order[s:e] for s, e in zip(starts, ends)}
    tkey = (a.astype(np.int64) * (n + 1) + b) * (n + 1) + c
    uniq, inv, counts = np.unique(tkey, return_inverse=True, return_counts=True)
    rep = np.empty(len(uniq), dtype=np.int64)
    rep[inv] = np.arange(n)
    members = {}
    order = np.argsort(inv, kind="stable")
    bounds = np.concatenate([[0], np.cumsum(counts)])
    for t in range(len(uniq)):
        members[t] = order[bounds[t]:bounds[t + 1]]
    assigned = np.zeros(n, dtype=bool)
    left = counts.astype(np.int64).copy()
    heap = [(-int(left[t]), t) for t in range(len(uniq))]
    heapq.heapify(heap)
    out = np.full(n, -1, dtype=np.int64)
    lab = 0
    while heap:
        negc, t = heapq.heappop(heap)
        cur = int((~assigned[members[t]]).sum())
        if cur == 0:
            continue
        if cur != -negc:
            heapq.heappush(heap, (-cur, t))
            continue
        v = rep[t]
        cand = np.concatenate([groups["ab"][int(key(a[v], b[v]))],
                               groups["ac"][int(key(a[v], c[v]))],
                               groups["bc"][int(key(b[v], c[v]))]])
        cand = cand[~assigned[cand]]
        out[cand] = lab
        assigned[cand] = True
        lab += 1
    return out


def iterated_flip(g: Graph, init: np.ndarray, rounds: int = 5, beta2: int = 1, rng=None,
                  ls: str = "insertion", return_all: bool = False,
                  time_limit: float | None = None):
    """Iterated flipping local search.

    ``beta2`` is 2*beta (1 gives the paper's beta = 1/2).  The local search
    subroutine is ``insertion`` (cluster insertions + multilevel moves) or
    ``multilevel``."""
    rng = np.random.default_rng(rng)
    base = np.full(g.indices.shape[0], 2, dtype=np.int64)  # doubled unit weight
    kappa = 2

    def LS(w, start):
        if ls == "insertion":
            return insertion_search(g, start, rng, weights=w, kappa=kappa)
        return multilevel(g, start, rng, weights=w, kappa=kappa)

    import time as _time
    t0 = _time.time()
    prev = LS(base, init)
    cands = [prev]
    for _ in range(rounds):
        if time_limit is not None and _time.time() - t0 > time_limit:
            break
        w1 = base + beta2 * cut_indicator(g, prev)
        c1 = LS(w1, prev)
        w2 = w1 + beta2 * cut_indicator(g, c1)
        c2 = LS(w2, c1)
        c3 = pivot3(prev, c1, c2)
        c3 = LS(base, c3)
        cands += [c1, c2, c3]
        prev = c2
    costs = [cost(g, x) for x in cands]
    best = cands[int(np.argmin(costs))]
    # final polish on the original objective
    best = LS(base, best)
    if return_all:
        return best, cands, costs
    return best
