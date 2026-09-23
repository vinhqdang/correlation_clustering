"""Local search for correlation clustering.

All routines work on a *weighted aggregated* graph: node i stands for a set of
``size[i]`` original vertices, ``w[p]`` is the number of positive original
edges between the two endpoint sets and ``self_w[i]`` the number of positive
original edges inside node i.  With unit sizes and weights this is the input
instance itself.  For a clustering with cluster sizes S_c and internal positive
weight W_c the cost is

    W_total - 2 * sum_c W_c + sum_c S_c (S_c - 1) / 2 ,

so moving node i (size s) from cluster A to cluster B changes the cost by

    delta = s * (S_B - (S_A - s)) + 2 * (w(i, A - i) - w(i, B)).

Every accepted move strictly decreases the cost, hence local search can only
improve a starting clustering (and preserves any approximation guarantee of
the initial solution).
"""
from __future__ import annotations

import numpy as np
import numba as nb

from .graph import Graph


@nb.njit(cache=True)
def _local_move(n, indptr, indices, w, cnt, kappa, size, labels, order, max_moves):
    """Queue-based best-improvement vertex moves.  Modifies ``labels``.

    ``w`` is the cost of cutting each positive slot, ``cnt`` the number of
    original positive pairs it represents and ``kappa`` the cost of a
    negative pair inside a cluster (unweighted instance: w = cnt = kappa = 1).
    Moving node i (size s) from A to B changes the cost by
        w(i,A-i) - w(i,B) + kappa * (s (S_B - S_A + s) - cnt(i,B) + cnt(i,A-i)).
    Returns the total change in cost (<= 0) and the number of moves."""
    csize = np.zeros(n, dtype=np.int64)
    for i in range(n):
        csize[labels[i]] += size[i]
    free = np.empty(n, dtype=np.int64)
    nfree = 0
    for c in range(n - 1, -1, -1):
        if csize[c] == 0:
            free[nfree] = c
            nfree += 1
    acc = np.zeros(n, dtype=np.int64)
    accn = np.zeros(n, dtype=np.int64)
    touched = np.empty(n, dtype=np.int64)
    queue = np.empty(n, dtype=np.int64)
    inq = np.zeros(n, dtype=np.bool_)
    head = 0
    qlen = 0
    for i in range(n):
        queue[i] = order[i]
        inq[order[i]] = True
    qlen = n
    total = 0
    moves = 0
    while qlen > 0 and moves < max_moves:
        i = queue[head]
        head += 1
        if head == n:
            head = 0
        qlen -= 1
        inq[i] = False
        a = labels[i]
        s = size[i]
        nt = 0
        for p in range(indptr[i], indptr[i + 1]):
            j = indices[p]
            if j == i:
                continue
            c = labels[j]
            if acc[c] == 0 and accn[c] == 0:
                touched[nt] = c
                nt += 1
            acc[c] += w[p]
            accn[c] += cnt[p]
        base = acc[a] + kappa * (accn[a] - s * (csize[a] - s))
        best = 0
        bestc = -1
        # move to a new singleton cluster
        if csize[a] > s:
            if base < best:
                best = base
                bestc = -2
        for t in range(nt):
            c = touched[t]
            if c == a:
                continue
            d = base - acc[c] + kappa * (s * csize[c] - accn[c])
            if d < best:
                best = d
                bestc = c
        for t in range(nt):
            acc[touched[t]] = 0
            accn[touched[t]] = 0
        if bestc == -1:
            continue
        if bestc == -2:
            nfree -= 1
            bestc = free[nfree]
        csize[a] -= s
        if csize[a] == 0:
            free[nfree] = a
            nfree += 1
        csize[bestc] += s
        labels[i] = bestc
        total += best
        moves += 1
        for p in range(indptr[i], indptr[i + 1]):
            j = indices[p]
            if not inq[j] and labels[j] != bestc:
                tail = head + qlen
                if tail >= n:
                    tail -= n
                queue[tail] = j
                qlen += 1
                inq[j] = True
    return total, moves


@nb.njit(cache=True)
def _aggregate(n, indptr, indices, w, cnt, size, self_w, labels, k):
    """Contract clusters (labels in 0..k-1) into super nodes."""
    nsize = np.zeros(k, dtype=np.int64)
    nself = np.zeros(k, dtype=np.int64)
    for i in range(n):
        nsize[labels[i]] += size[i]
        nself[labels[i]] += self_w[i]
    # bucket nodes by cluster
    start = np.zeros(k + 1, dtype=np.int64)
    for i in range(n):
        start[labels[i] + 1] += 1
    for c in range(k):
        start[c + 1] += start[c]
    fill = start[:-1].copy()
    members = np.empty(n, dtype=np.int64)
    for i in range(n):
        members[fill[labels[i]]] = i
        fill[labels[i]] += 1
    acc = np.zeros(k, dtype=np.int64)
    accn = np.zeros(k, dtype=np.int64)
    touched = np.empty(k, dtype=np.int64)
    out_ptr = np.zeros(k + 1, dtype=np.int64)
    cap = indices.shape[0]
    out_idx = np.empty(cap, dtype=np.int32)
    out_w = np.empty(cap, dtype=np.int64)
    out_n = np.empty(cap, dtype=np.int64)
    pos = 0
    for c in range(k):
        nt = 0
        for t in range(start[c], start[c + 1]):
            i = members[t]
            for p in range(indptr[i], indptr[i + 1]):
                d = labels[indices[p]]
                if d == c:
                    if indices[p] > i:
                        nself[c] += cnt[p]
                    continue
                if accn[d] == 0 and acc[d] == 0:
                    touched[nt] = d
                    nt += 1
                acc[d] += w[p]
                accn[d] += cnt[p]
        touched[:nt].sort()
        for t in range(nt):
            d = touched[t]
            out_idx[pos] = d
            out_w[pos] = acc[d]
            out_n[pos] = accn[d]
            acc[d] = 0
            accn[d] = 0
            pos += 1
        out_ptr[c + 1] = pos
    return out_ptr, out_idx[:pos].copy(), out_w[:pos].copy(), out_n[:pos].copy(), nsize, nself


def _compact(labels):
    _, inv = np.unique(labels, return_inverse=True)
    return inv.astype(np.int64), int(inv.max()) + 1 if len(inv) else 0


def local_search(g: Graph, labels: np.ndarray, rng=None, max_moves: int = 1 << 62,
                 weights: np.ndarray | None = None, kappa: int = 1):
    """Single-level vertex-move local search starting from ``labels``.

    ``weights`` (aligned with ``g.indices``) are integer costs of cutting
    positive edges; negative pairs inside clusters cost ``kappa``."""
    rng = np.random.default_rng(rng)
    lab, _ = _compact(np.asarray(labels))
    cnt = np.ones(g.indices.shape[0], dtype=np.int64)
    w = cnt if weights is None else np.asarray(weights, dtype=np.int64)
    size = np.ones(g.n, dtype=np.int64)
    order = rng.permutation(g.n).astype(np.int64)
    _local_move(g.n, g.indptr, g.indices, w, cnt, int(kappa), size, lab, order, max_moves)
    return lab


def weighted_cost(g: Graph, labels: np.ndarray, weights: np.ndarray | None, kappa: int = 1):
    """sum of weights of cut positive edges + kappa * negative pairs inside."""
    from .objective import disagreements
    if weights is None:
        a, b = disagreements(g, labels)
        return int(a + kappa * b)
    lab = np.asarray(labels)
    src = np.repeat(np.arange(g.n), g.degrees)
    cut = lab[src] != lab[g.indices]
    wc = int(np.asarray(weights)[cut].sum()) // 2
    _, b = disagreements(g, labels)
    return wc + kappa * int(b)


def multilevel(g: Graph, labels: np.ndarray | None = None, rng=None,
               max_levels: int = 50, vcycles: int = 3,
               weights: np.ndarray | None = None, kappa: int = 1) -> np.ndarray:
    """Louvain-style multilevel local search for the (weighted) CC objective.

    Level 0 starts from ``labels`` (singletons if None).  After local moves
    converge, clusters are contracted and moves are repeated on the contracted
    graph, where moving a node merges whole clusters.  ``vcycles`` repeats the
    whole procedure starting from the current clustering.  The result is never
    worse than the starting clustering.
    """
    rng = np.random.default_rng(rng)
    if labels is None:
        lab = np.arange(g.n, dtype=np.int64)
    else:
        lab, _ = _compact(np.asarray(labels))
    best_c = weighted_cost(g, lab, weights, kappa)
    for _ in range(max(1, vcycles)):
        new = _one_cycle(g, lab.copy(), rng, max_levels, weights, kappa)
        c = weighted_cost(g, new, weights, kappa)
        if c < best_c:
            lab, best_c = new, c
        else:
            break
    return lab


def _one_cycle(g, lab, rng, max_levels, weights, kappa):
    indptr, indices = g.indptr, g.indices
    cnt0 = np.ones(indices.shape[0], dtype=np.int64)
    w0 = cnt0 if weights is None else np.asarray(weights, dtype=np.int64)
    w, cnt = w0, cnt0
    size = np.ones(g.n, dtype=np.int64)
    self_w = np.zeros(g.n, dtype=np.int64)
    n = g.n
    maps = []
    cur = lab
    for _lvl in range(max_levels):
        order = rng.permutation(n).astype(np.int64)
        _local_move(n, indptr, indices, w, cnt, int(kappa), size, cur, order, 1 << 62)
        cur, k = _compact(cur)
        maps.append(cur)
        if k == n:
            break
        indptr, indices, w, cnt, size, self_w = _aggregate(n, indptr, indices, w, cnt, size,
                                                           self_w, cur, k)
        n = k
        cur = np.arange(n, dtype=np.int64)
    out = maps[-1]
    for mp in reversed(maps[:-1]):
        out = out[mp]
    out, _ = _compact(np.asarray(out, dtype=np.int64))
    s0 = np.ones(g.n, dtype=np.int64)
    _local_move(g.n, g.indptr, g.indices, w0, cnt0, int(kappa), s0, out,
                rng.permutation(g.n).astype(np.int64), 1 << 62)
    return out
