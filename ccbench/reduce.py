"""Exact critical-clique contraction.

Two vertices are twins if they have the same closed neighbourhood N[v].  The
classes of twins (critical cliques) are cliques of G+, and some optimal
clustering keeps every critical clique inside one cluster (Guo, TCS 2009), so
contracting them loses nothing.  The contracted instance is a weighted graph in
the representation of :mod:`ccbench.localsearch`: node i stands for size[i]
vertices, w[p] = size[i] * size[j] counts the positive pairs between adjacent
nodes, and self_w[i] = size[i] (size[i] - 1) / 2.  Its cost function

    W_total - 2 sum_c W_c + sum_c S_c (S_c - 1) / 2

equals the cost of the expanded clustering on the original graph.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numba as nb

from .graph import Graph


@dataclass
class WGraph:
    n: int
    indptr: np.ndarray
    indices: np.ndarray
    w: np.ndarray
    size: np.ndarray
    self_w: np.ndarray

    @property
    def total_w(self) -> int:
        return int(self.w.sum() // 2 + self.self_w.sum())

    def cost(self, labels: np.ndarray) -> int:
        return int(_wcost(self.n, self.indptr, self.indices, self.w, self.size, self.self_w,
                          np.asarray(labels, dtype=np.int64)))

    def edges(self) -> np.ndarray:
        src = np.repeat(np.arange(self.n), np.diff(self.indptr))
        m = src < self.indices
        return np.stack([src[m], self.indices[m]], axis=1)


@nb.njit(cache=True)
def _wcost(n, indptr, indices, w, size, self_w, labels):
    k = labels.max() + 1 if n else 0
    S = np.zeros(k, dtype=np.int64)
    tot = 0
    inside = 0
    for i in range(n):
        S[labels[i]] += size[i]
        tot += self_w[i]
        inside += self_w[i]
        for p in range(indptr[i], indptr[i + 1]):
            j = indices[p]
            if j > i:
                tot += w[p]
                if labels[j] == labels[i]:
                    inside += w[p]
    c = tot - 2 * inside
    for x in range(k):
        c += S[x] * (S[x] - 1) // 2
    return c


def unit(g: Graph) -> WGraph:
    """The unweighted instance as a WGraph."""
    return WGraph(g.n, g.indptr, g.indices, np.ones(g.indices.shape[0], dtype=np.int64),
                  np.ones(g.n, dtype=np.int64), np.zeros(g.n, dtype=np.int64))


def critical_cliques(g: Graph) -> tuple[np.ndarray, int]:
    """Twin classes (equal closed neighbourhoods).  Returns (class of each vertex, count)."""
    keys = {}
    grp = np.empty(g.n, dtype=np.int64)
    for v in range(g.n):
        nb_ = g.indices[g.indptr[v]:g.indptr[v + 1]]
        key = np.sort(np.append(nb_, v)).astype(np.int64).tobytes()
        grp[v] = keys.setdefault(key, len(keys))
    return grp, len(keys)


def contract(g: Graph) -> tuple[WGraph, np.ndarray]:
    """Critical-clique contraction.  Returns (weighted instance, vertex -> node map)."""
    from .localsearch import _aggregate
    grp, k = critical_cliques(g)
    one = np.ones(g.indices.shape[0], dtype=np.int64)
    ip, ix, w, cn, sz, sw = _aggregate(g.n, g.indptr, g.indices, one, one,
                                       np.ones(g.n, np.int64), np.zeros(g.n, np.int64), grp, k)
    return WGraph(k, ip, ix.astype(np.int32), w, sz, sw), grp
