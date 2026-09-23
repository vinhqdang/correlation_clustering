"""Pivot (KwikCluster) of Ailon, Charikar and Newman and simple variants."""
from __future__ import annotations

import numpy as np
import numba as nb

from .graph import Graph
from .objective import cost


@nb.njit(cache=True)
def _pivot(n, indptr, indices, order):
    labels = np.full(n, -1, dtype=np.int64)
    c = 0
    for i in range(n):
        u = order[i]
        if labels[u] >= 0:
            continue
        labels[u] = c
        for p in range(indptr[u], indptr[u + 1]):
            v = indices[p]
            if labels[v] < 0:
                labels[v] = c
        c += 1
    return labels


def pivot(g: Graph, rng: np.random.Generator | int | None = None,
          order: np.ndarray | None = None) -> np.ndarray:
    """One run of Pivot with a uniformly random order (3-approx. in expectation)."""
    if order is None:
        rng = np.random.default_rng(rng)
        order = rng.permutation(g.n)
    return _pivot(g.n, g.indptr, g.indices, np.asarray(order, dtype=np.int64))


def best_pivot(g: Graph, runs: int = 50, rng=None):
    """Best of ``runs`` independent Pivot runs."""
    rng = np.random.default_rng(rng)
    best, best_c = None, None
    for _ in range(runs):
        lab = pivot(g, rng)
        c = cost(g, lab)
        if best_c is None or c < best_c:
            best, best_c = lab, c
    return best
