"""Correlation clustering objective on complete signed graphs."""
from __future__ import annotations

import numpy as np
import numba as nb

from .graph import Graph


@nb.njit(cache=True)
def _cost(n, indptr, indices, labels):
    cut = 0
    for u in range(n):
        lu = labels[u]
        for p in range(indptr[u], indptr[u + 1]):
            v = indices[p]
            if v > u and labels[v] != lu:
                cut += 1
    m = indices.shape[0] // 2
    inside_pos = m - cut
    sizes = np.zeros(labels.max() + 1, dtype=np.int64)
    for u in range(n):
        sizes[labels[u]] += 1
    pairs = 0
    for s in sizes:
        pairs += s * (s - 1) // 2
    return cut, pairs - inside_pos


def disagreements(g: Graph, labels: np.ndarray):
    """Return (positive edges cut, negative pairs inside clusters)."""
    labels = np.ascontiguousarray(normalize_labels(labels))
    return _cost(g.n, g.indptr, g.indices, labels)


def cost(g: Graph, labels: np.ndarray) -> int:
    a, b = disagreements(g, labels)
    return int(a + b)


def normalize_labels(labels: np.ndarray) -> np.ndarray:
    """Relabel clusters as 0..k-1 in order of first appearance."""
    labels = np.asarray(labels)
    _, first, inv = np.unique(labels, return_index=True, return_inverse=True)
    order = np.argsort(np.argsort(first))
    return order[inv].astype(np.int64)


def is_cluster_graph(g: Graph) -> bool:
    """True iff G+ is a disjoint union of cliques (cost 0 instance)."""
    from .graph import components
    comp, k = components(g)
    return cost(g, comp) == 0
