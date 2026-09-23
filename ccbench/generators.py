"""Synthetic instance generators."""
from __future__ import annotations

import numpy as np

from .graph import Graph, from_edges


def planted_partition(n: int, k: int, p_in: float, p_out: float, rng=None,
                      sizes: np.ndarray | None = None):
    """Stochastic block model with k blocks.  Returns (graph, planted labels)."""
    rng = np.random.default_rng(rng)
    if sizes is None:
        labels = rng.integers(0, k, size=n)
    else:
        labels = np.repeat(np.arange(len(sizes)), sizes)
        n = len(labels)
    edges = []
    # intra-block edges
    for c in np.unique(labels):
        mem = np.flatnonzero(labels == c)
        s = len(mem)
        if s < 2:
            continue
        cnt = rng.binomial(s * (s - 1) // 2, p_in)
        e = _sample_pairs(s, cnt, rng)
        edges.append(mem[e])
    # inter-block edges: sample uniformly among all pairs and reject intra
    tot = n * (n - 1) // 2
    cnt = rng.binomial(tot, p_out)
    if cnt:
        e = _sample_pairs(n, int(cnt * 1.0), rng)
        e = e[labels[e[:, 0]] != labels[e[:, 1]]]
        edges.append(e)
    E = np.concatenate(edges) if edges else np.zeros((0, 2), dtype=np.int64)
    return from_edges(n, E), labels


def _sample_pairs(s, cnt, rng):
    """Sample ``cnt`` distinct unordered pairs from s items (approximately uniform)."""
    tot = s * (s - 1) // 2
    cnt = min(cnt, tot)
    if cnt == 0:
        return np.zeros((0, 2), dtype=np.int64)
    if cnt > tot // 3:
        iu = np.triu_indices(s, 1)
        idx = rng.choice(tot, size=cnt, replace=False)
        return np.stack([iu[0][idx], iu[1][idx]], axis=1).astype(np.int64)
    u = rng.integers(0, s, size=int(cnt * 1.2) + 10)
    v = rng.integers(0, s, size=u.shape[0])
    keep = u != v
    a = np.minimum(u[keep], v[keep])
    b = np.maximum(u[keep], v[keep])
    key = np.unique(a.astype(np.int64) * s + b)[:cnt]
    return np.stack([key // s, key % s], axis=1)


def noisy_cliques(n: int, k: int, flip: float, rng=None, powerlaw: float | None = None):
    """Disjoint union of k cliques with each pair flipped independently w.p. ``flip``.

    With ``powerlaw`` set, cluster sizes follow a truncated power law with that
    exponent (heterogeneous cluster sizes).
    """
    rng = np.random.default_rng(rng)
    if powerlaw is None:
        sizes = np.full(k, n // k)
        sizes[: n - sizes.sum()] += 1
    else:
        x = rng.pareto(powerlaw - 1, size=k) + 1
        sizes = np.maximum(1, np.round(x / x.sum() * n)).astype(int)
        sizes[-1] += n - sizes.sum()
        sizes = sizes[sizes > 0]
    # flipping a pair inside a block removes it; between blocks adds it
    # expected positive density inside blocks 1-flip, between blocks flip
    return planted_partition(0, 0, 1.0 - flip, flip, rng=rng, sizes=sizes)


def sparse_planted(n: int, avg_cluster: int, deg_out: float, p_in: float, rng=None):
    """Large sparse planted instance: clusters of mean size ``avg_cluster``
    (geometric sizes), intra density p_in, ~deg_out random inter edges per vertex."""
    rng = np.random.default_rng(rng)
    sizes = []
    tot = 0
    while tot < n:
        s = int(rng.geometric(1.0 / avg_cluster))
        s = min(s, n - tot)
        sizes.append(s)
        tot += s
    sizes = np.array(sizes)
    labels = np.repeat(np.arange(len(sizes)), sizes)
    edges = []
    starts = np.concatenate([[0], np.cumsum(sizes)])
    for c, s in enumerate(sizes):
        if s < 2:
            continue
        cnt = rng.binomial(s * (s - 1) // 2, p_in)
        edges.append(_sample_pairs(s, cnt, rng) + starts[c])
    cnt = int(n * deg_out / 2)
    u = rng.integers(0, n, size=cnt)
    v = rng.integers(0, n, size=cnt)
    e = np.stack([u, v], axis=1)
    e = e[labels[u] != labels[v]]
    edges.append(e)
    return from_edges(n, np.concatenate(edges)), labels
