"""Distance-2 support of a complete signed instance.

Lemma (diameter-2).  In every optimal clustering each cluster has diameter at
most 2 in G+.  Hence all pairs at G+-distance >= 3 may be fixed as "separated"
without loss, and every LP/dual object only needs the pair set

    P = E+  union  N2,     N2 = {uv not in E+ : N(u) and N(v) intersect}.

Pairs are indexed 0..m-1 for E+ (u < v) and m..m+|N2|-1 for N2.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numba as nb

from .graph import Graph


@dataclass
class Support:
    g: Graph
    eid: np.ndarray      # CSR-aligned id of each E+ slot (length 2m)
    n2ptr: np.ndarray    # CSR of N2 pairs (both directions)
    n2idx: np.ndarray
    n2id: np.ndarray     # id of each N2 slot
    pu: np.ndarray       # endpoints of every pair id (u < v)
    pv: np.ndarray
    npos: int            # = m

    @property
    def npairs(self) -> int:
        return int(self.pu.shape[0])

    @property
    def nneg(self) -> int:
        return self.npairs - self.npos


@nb.njit(cache=True)
def _edge_ids(n, indptr, indices):
    eid = np.empty(indices.shape[0], dtype=np.int64)
    m = indices.shape[0] // 2
    pu = np.empty(m, dtype=np.int32)
    pv = np.empty(m, dtype=np.int32)
    c = 0
    for u in range(n):
        for p in range(indptr[u], indptr[u + 1]):
            v = indices[p]
            if u < v:
                eid[p] = c
                pu[c] = u
                pv[c] = v
                c += 1
    # fill reverse direction by binary search
    for u in range(n):
        for p in range(indptr[u], indptr[u + 1]):
            v = indices[p]
            if u > v:
                lo = indptr[v]
                hi = indptr[v + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if indices[mid] < u:
                        lo = mid + 1
                    else:
                        hi = mid
                eid[p] = eid[lo]
    return eid, pu, pv


@nb.njit(cache=True)
def _count_n2(n, indptr, indices):
    mark = np.full(n, -1, dtype=np.int64)
    cnt = np.zeros(n, dtype=np.int64)
    for u in range(n):
        for p in range(indptr[u], indptr[u + 1]):
            mark[indices[p]] = u
        mark[u] = u
        for p in range(indptr[u], indptr[u + 1]):
            w = indices[p]
            for q in range(indptr[w], indptr[w + 1]):
                v = indices[q]
                if mark[v] != u:
                    mark[v] = u
                    cnt[u] += 1
    return cnt


@nb.njit(cache=True)
def _fill_n2(n, indptr, indices, n2ptr, m):
    mark = np.full(n, -1, dtype=np.int64)
    n2idx = np.empty(n2ptr[n], dtype=np.int32)
    for u in range(n):
        for p in range(indptr[u], indptr[u + 1]):
            mark[indices[p]] = u
        mark[u] = u
        pos = n2ptr[u]
        for p in range(indptr[u], indptr[u + 1]):
            w = indices[p]
            for q in range(indptr[w], indptr[w + 1]):
                v = indices[q]
                if mark[v] != u:
                    mark[v] = u
                    n2idx[pos] = v
                    pos += 1
        n2idx[n2ptr[u]:n2ptr[u + 1]].sort()
    # ids
    n2id = np.empty(n2ptr[n], dtype=np.int64)
    tot = n2ptr[n] // 2
    pu = np.empty(tot, dtype=np.int32)
    pv = np.empty(tot, dtype=np.int32)
    c = 0
    for u in range(n):
        for p in range(n2ptr[u], n2ptr[u + 1]):
            v = n2idx[p]
            if u < v:
                n2id[p] = m + c
                pu[c] = u
                pv[c] = v
                c += 1
    for u in range(n):
        for p in range(n2ptr[u], n2ptr[u + 1]):
            v = n2idx[p]
            if u > v:
                lo = n2ptr[v]
                hi = n2ptr[v + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if n2idx[mid] < u:
                        lo = mid + 1
                    else:
                        hi = mid
                n2id[p] = n2id[lo]
    return n2idx, n2id, pu, pv


def build_support(g: Graph) -> Support:
    eid, pu1, pv1 = _edge_ids(g.n, g.indptr, g.indices)
    cnt = _count_n2(g.n, g.indptr, g.indices)
    n2ptr = np.zeros(g.n + 1, dtype=np.int64)
    n2ptr[1:] = np.cumsum(cnt)
    n2idx, n2id, pu2, pv2 = _fill_n2(g.n, g.indptr, g.indices, n2ptr, g.m)
    return Support(g, eid, n2ptr, n2idx, n2id,
                   np.concatenate([pu1, pu2]), np.concatenate([pv1, pv2]), g.m)


@nb.njit(cache=True)
def _find(ptr, idx, u, v):
    lo = ptr[u]
    hi = ptr[u + 1]
    while lo < hi:
        mid = (lo + hi) >> 1
        if idx[mid] < v:
            lo = mid + 1
        else:
            hi = mid
    if lo < ptr[u + 1] and idx[lo] == v:
        return lo
    return -1


@nb.njit(cache=True)
def pair_id(indptr, indices, eid, n2ptr, n2idx, n2id, u, v):
    """Id of pair uv in P, or -1 if uv is a far pair (distance >= 3)."""
    p = _find(indptr, indices, u, v)
    if p >= 0:
        return eid[p]
    p = _find(n2ptr, n2idx, u, v)
    if p >= 0:
        return n2id[p]
    return -1


@nb.njit(cache=True)
def _count_wedge_triangles(n, indptr, indices):
    """Number of triangles of G+ and number of open wedges (bad triangles)."""
    tri = 0
    wedges = 0
    mark = np.full(n, -1, dtype=np.int64)
    for w in range(n):
        d = indptr[w + 1] - indptr[w]
        wedges += d * (d - 1) // 2
    for u in range(n):
        for p in range(indptr[u], indptr[u + 1]):
            mark[indices[p]] = u
        for p in range(indptr[u], indptr[u + 1]):
            v = indices[p]
            if v <= u:
                continue
            for q in range(indptr[v], indptr[v + 1]):
                x = indices[q]
                if x > v and mark[x] == u:
                    tri += 1
    return tri, wedges - 3 * tri


def wedge_stats(g: Graph):
    """(number of triangles, number of open wedges = bad triangles)."""
    return _count_wedge_triangles(g.n, g.indptr, g.indices)


@nb.njit(cache=True)
def _enumerate_triangles(n, indptr, indices, eid, n2ptr, n2idx, n2id, ntri):
    """All triangles with >= 2 positive edges.

    Returns arrays (a, b, c) of pair ids.  For a bad triangle u-w-v (uw, wv
    positive, uv negative) the negative pair is stored in c.  Positive
    triangles are listed once.
    """
    a = np.empty(ntri, dtype=np.int64)
    b = np.empty(ntri, dtype=np.int64)
    c = np.empty(ntri, dtype=np.int64)
    t = 0
    for w in range(n):
        s = indptr[w]
        e = indptr[w + 1]
        for p in range(s, e):
            u = indices[p]
            for q in range(p + 1, e):
                v = indices[q]
                r = _find(indptr, indices, u, v)
                if r >= 0:
                    # positive triangle; list once with w the smallest vertex
                    if w < u:
                        a[t] = eid[p]
                        b[t] = eid[q]
                        c[t] = eid[r]
                        t += 1
                else:
                    r = _find(n2ptr, n2idx, u, v)
                    a[t] = eid[p]
                    b[t] = eid[q]
                    c[t] = n2id[r]
                    t += 1
    return a[:t], b[:t], c[:t]


def wedge_triangles(sup: Support):
    tri, wed = wedge_stats(sup.g)
    g = sup.g
    return _enumerate_triangles(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr,
                                sup.n2idx, sup.n2id, tri + wed)
