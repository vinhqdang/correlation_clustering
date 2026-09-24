"""Exact cluster-union neighbourhood.

Given a clustering C, choose a set U of clusters that is connected in the
cluster graph and covers at most ``max_size`` vertices, and replace the
clustering of their union X by an optimal clustering of the induced instance
G+[X].  Pairs between X and the rest stay separated, so their cost is a
constant, and the move never increases the cost: the result is optimal among
all clusterings that agree with C outside X and keep X separated from the
rest.  The sub-instances are solved exactly by the cutting-plane ILP on the
distance-two support (:class:`ccbench.lp.SparseLP`).
"""
from __future__ import annotations

import time

import numpy as np

from .graph import Graph, from_edges
from .objective import cost
from .localsearch import _compact


def induced(g: Graph, V: np.ndarray) -> Graph:
    pos = -np.ones(g.n, dtype=np.int64)
    pos[V] = np.arange(len(V))
    src = np.repeat(V, g.degrees[V])
    dst = np.concatenate([g.indices[g.indptr[v]:g.indptr[v + 1]] for v in V]) if len(V) else \
        np.empty(0, dtype=np.int64)
    pd = pos[dst]
    ps = pos[src]
    m = pd > ps
    return from_edges(len(V), np.stack([ps[m], pd[m]], axis=1).astype(np.int64))


def solve_exact(h: Graph, cur: np.ndarray | None = None, time_limit: float = 10.0):
    """Optimal clustering of a small instance, or None if it was not proven in
    time.  The LP relaxation with triangle and star inequalities is solved
    first; when its bound shows that the current clustering ``cur`` is optimal,
    ``cur`` is returned.  Otherwise the same model (with all cuts found so far)
    is solved as a MIP, starting from ``cur``."""
    import highspy
    from .lp import SparseLP
    if h.m == 0:
        return np.arange(h.n, dtype=np.int64), 0
    t0 = time.time()
    lp = SparseLP(h)
    if cur is not None:
        c0 = cost(h, cur)
        v = lp.solve_with_stars(time_limit=time_limit)
        if v > c0 - 1 + 1e-6:
            return cur, c0
    N = lp.sup.npairs
    lp.h.changeColsIntegrality(N, np.arange(N, dtype=np.int32),
                               np.array([highspy.HighsVarType.kInteger] * N))
    lp.integral = True
    if cur is not None:
        s = lp.sup
        sol = highspy.HighsSolution()
        sol.col_value = list((cur[s.pu] != cur[s.pv]).astype(float))
        sol.value_valid = True
        lp.h.setSolution(sol)
    val = lp.solve(time_limit=max(0.5, time_limit - (time.time() - t0)))
    if lp.dual_bound() < val - 1e-6 or not lp.converged:
        return None
    lab = lp.clustering()
    return lab, cost(h, lab)


class ClusterGraph:
    """Clusters of a labelling with their sizes, members and the number of
    positive edges between every pair of adjacent clusters."""

    def __init__(self, g: Graph, lab: np.ndarray):
        self.lab = lab
        k = int(lab.max()) + 1
        self.size = np.bincount(lab, minlength=k)
        order = np.argsort(lab, kind="stable")
        self.start = np.zeros(k + 1, dtype=np.int64)
        self.start[1:] = np.cumsum(self.size)
        self.members = order
        src = np.repeat(np.arange(g.n), g.degrees)
        a, b = lab[src], lab[g.indices]
        m = a != b
        key = a[m].astype(np.int64) * k + b[m]
        uk, cnt = np.unique(key, return_counts=True)
        ca, cb = uk // k, uk % k
        self.ptr = np.zeros(k + 1, dtype=np.int64)
        np.add.at(self.ptr, ca + 1, 1)
        self.ptr = np.cumsum(self.ptr)
        self.nbr = cb
        self.w = cnt
        self.cut = np.bincount(ca, weights=cnt, minlength=k)  # cut edges per cluster

    def cluster_nodes(self, c):
        return self.members[self.start[c]:self.start[c + 1]]

    def grow(self, seed: int, max_size: int, rng) -> list[int]:
        """Greedy union: repeatedly add the adjacent cluster with the most
        edges into the union (random tie-breaking) while it fits."""
        U = [seed]
        inU = {seed}
        tot = int(self.size[seed])
        score = {}
        for p in range(self.ptr[seed], self.ptr[seed + 1]):
            score[int(self.nbr[p])] = score.get(int(self.nbr[p]), 0) + int(self.w[p])
        while score:
            cand = [(s + rng.random(), c) for c, s in score.items()
                    if tot + self.size[c] <= max_size]
            if not cand:
                break
            _, c = max(cand)
            U.append(c)
            inU.add(c)
            tot += int(self.size[c])
            del score[c]
            for p in range(self.ptr[c], self.ptr[c + 1]):
                d = int(self.nbr[p])
                if d not in inU:
                    score[d] = score.get(d, 0) + int(self.w[p])
        return U


def exact_union_search(g: Graph, labels: np.ndarray, time_limit: float = 60.0,
                       max_size: int = 50, sub_time: float = 5.0, rng=None,
                       max_density: float = 4.0,
                       verbose: bool = False, stats: dict | None = None):
    """Iterated exact cluster-union moves (never increases the cost)."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    t0 = time.time()
    cg = ClusterGraph(g, lab)
    seen = set()
    st = {"tried": 0, "improved": 0, "gain": 0, "failed": 0}
    while time.time() - t0 < time_limit:
        cand = np.flatnonzero(cg.cut > 0)
        if len(cand) == 0:
            break
        seed = int(cand[rng.integers(len(cand))])
        U = cg.grow(seed, max_size, rng)
        if len(U) < 2:
            continue
        key = frozenset(U)
        if key in seen:
            if len(seen) > 50 * len(cg.size):
                break
            continue
        seen.add(key)
        X = np.concatenate([cg.cluster_nodes(c) for c in U])
        h = induced(g, X)
        if h.m > max_density * h.n:
            st["dense"] = st.get("dense", 0) + 1
            continue
        cur = cost(h, _compact(lab[X])[0])
        st["tried"] += 1
        res = solve_exact(h, _compact(lab[X])[0],
                          min(sub_time, max(0.5, time_limit - (time.time() - t0))))
        if res is None:
            st["failed"] += 1
            continue
        sub, val = res
        if val < cur:
            lab[X] = lab.max() + 1 + sub
            lab = _compact(lab)[0].copy()
            cg = ClusterGraph(g, lab)
            seen.clear()
            st["improved"] += 1
            st["gain"] += cur - val
            if verbose:
                print(f"  union of {len(U)} clusters / {len(X)} vertices: -{cur - val} "
                      f"-> {cost(g, lab)} ({time.time() - t0:.1f}s)", flush=True)
    if stats is not None:
        stats.update(st)
    return lab


def hub_union_search(g: Graph, labels: np.ndarray, time_limit: float = 60.0,
                     max_size: int = 100, min_deg: int = 15, sub_time: float = 5.0,
                     rng=None, verbose: bool = False, stats: dict | None = None):
    """Exact re-optimisation of the union of the clusters of a vertex and of all
    its neighbours, for vertices of degree >= ``min_deg`` in random order
    (cluster-union neighbourhood centred at hubs; never increases the cost)."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    t0 = time.time()
    hubs = np.flatnonzero(g.degrees >= min_deg)
    rng.shuffle(hubs)
    st = {"tried": 0, "improved": 0, "gain": 0, "failed": 0, "big": 0}
    order = np.argsort(lab, kind="stable")
    start = np.zeros(lab.max() + 2, dtype=np.int64)
    start[1:] = np.cumsum(np.bincount(lab, minlength=lab.max() + 1))
    for v in hubs:
        if time.time() - t0 > time_limit:
            break
        cl = np.unique(lab[np.append(g.indices[g.indptr[v]:g.indptr[v + 1]], v)])
        size = int((start[cl + 1] - start[cl]).sum())
        if size > max_size:
            st["big"] += 1
            continue
        X = np.concatenate([order[start[c]:start[c + 1]] for c in cl])
        h = induced(g, X)
        cur_lab = _compact(lab[X])[0]
        cur = cost(h, cur_lab)
        st["tried"] += 1
        res = solve_exact(h, cur_lab, min(sub_time, max(0.5, time_limit - (time.time() - t0))))
        if res is None:
            st["failed"] += 1
            continue
        sub, val = res
        if val < cur:
            lab[X] = lab.max() + 1 + sub
            lab = _compact(lab)[0].copy()
            order = np.argsort(lab, kind="stable")
            start = np.zeros(lab.max() + 2, dtype=np.int64)
            start[1:] = np.cumsum(np.bincount(lab, minlength=lab.max() + 1))
            st["improved"] += 1
            st["gain"] += cur - val
            if verbose:
                print(f"  hub {v} (deg {g.degrees[v]}): {len(X)} vertices -{cur - val} "
                      f"({time.time() - t0:.1f}s)", flush=True)
    if stats is not None:
        stats.update(st)
    return lab
