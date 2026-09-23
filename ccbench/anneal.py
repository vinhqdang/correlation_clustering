"""Simulated annealing on vertex moves for correlation clustering.

A move takes a vertex v from its cluster A to cluster B (the cluster of a
random positive neighbour, or a new singleton).  With k_X the number of
positive neighbours of v in X (v excluded), the change of the objective is

    delta = (|B| - 2 k_B) - (|A| - 1 - 2 k_A),

computed in O(deg v).  Moves are accepted by the Metropolis rule with a
geometric temperature schedule; the best clustering seen is returned.
"""
from __future__ import annotations

import time

import numpy as np
import numba as nb

from .graph import Graph
from .objective import cost
from .localsearch import _compact


@nb.njit(cache=True)
def _anneal(n, indptr, indices, labels, iters, t_start, t_end, p_single, seed, cur_cost):
    np.random.seed(seed)
    csize = np.zeros(n, dtype=np.int64)
    for i in range(n):
        csize[labels[i]] += 1
    free = np.empty(n, dtype=np.int64)
    nfree = 0
    for c in range(n - 1, -1, -1):
        if csize[c] == 0:
            free[nfree] = c
            nfree += 1
    best = labels.copy()
    best_cost = cur_cost
    cost_ = cur_cost
    lt0 = np.log(t_start)
    lt1 = np.log(t_end)
    # vertices with at least one neighbour
    cand = np.empty(n, dtype=np.int64)
    nc = 0
    for v in range(n):
        if indptr[v + 1] > indptr[v]:
            cand[nc] = v
            nc += 1
    since_best = 0
    for it in range(iters):
        if it % 1024 == 0:
            T = np.exp(lt0 + (lt1 - lt0) * it / iters)
        v = cand[np.random.randint(nc)]
        a = labels[v]
        d = indptr[v + 1] - indptr[v]
        if np.random.random() < p_single:
            if csize[a] == 1:
                continue
            b = -1
        else:
            u = indices[indptr[v] + np.random.randint(d)]
            b = labels[u]
            if b == a:
                continue
        ka = 0
        kb = 0
        for p in range(indptr[v], indptr[v + 1]):
            lu = labels[indices[p]]
            if lu == a:
                ka += 1
            elif lu == b:
                kb += 1
        sb = 0 if b == -1 else csize[b]
        delta = (sb - 2 * kb) - (csize[a] - 1 - 2 * ka)
        if delta <= 0 or np.random.random() < np.exp(-delta / T):
            if b == -1:
                nfree -= 1
                b = free[nfree]
            csize[a] -= 1
            if csize[a] == 0:
                free[nfree] = a
                nfree += 1
            csize[b] += 1
            labels[v] = b
            cost_ += delta
            if cost_ < best_cost:
                best_cost = cost_
                since_best += 1
                # copy lazily: only when the improvement persists a little
                if since_best >= 1:
                    best[:] = labels
                    since_best = 0
    return best, best_cost


def anneal(g: Graph, labels: np.ndarray, iters: int | None = None, time_limit: float = 60.0,
           t_start: float = 2.0, t_end: float = 0.05, p_single: float = 0.05, rng=None):
    """Simulated annealing from ``labels``.  Returns the best clustering seen."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    c0 = cost(g, lab)
    if iters is None:
        # calibrate the iteration rate
        t = time.time()
        _anneal(g.n, g.indptr, g.indices, lab.copy(), 200000, t_end, t_end, p_single, 1, c0)
        rate = 200000 / max(time.time() - t, 1e-3)
        iters = int(rate * time_limit)
    best, bc = _anneal(g.n, g.indptr, g.indices, lab, iters, t_start, t_end, p_single,
                       int(rng.integers(1 << 30)), c0)
    assert bc == cost(g, best)
    return best
