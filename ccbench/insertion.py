"""Cluster-insertion local search.

A cluster insertion C + S removes the vertices of S from their clusters and
makes S a new cluster.  Vertex moves and cluster merges are special cases.
A clustering that no insertion improves is a 2-approximation (sum the local
optimality inequalities over the optimal clusters; Cohen-Addad et al., STOC
2024).  Since checking all 2^n insertions is intractable we generate, for every
root r, one candidate S(r) by "cleaning" the closed neighbourhood N+[r]:
a vertex v of N+[r] belongs to S iff joining S is cheaper for v than staying
in what remains of its current cluster (evaluated with exact counts, repeated
for a few sweeps).  The candidate is applied iff the exact change of the
(weighted) objective is negative.
"""
from __future__ import annotations

import numpy as np
import numba as nb

from .graph import Graph
from .localsearch import _compact, multilevel, weighted_cost


@nb.njit(cache=True)
def _insertion_round(n, indptr, indices, w, kappa, labels, order, sweeps, max_root_deg):
    csize = np.zeros(n, dtype=np.int64)
    for i in range(n):
        csize[labels[i]] += 1
    free = np.empty(n, dtype=np.int64)
    nfree = 0
    for c in range(n - 1, -1, -1):
        if csize[c] == 0:
            free[nfree] = c
            nfree += 1
    inS = np.zeros(n, dtype=np.bool_)
    cntS = np.zeros(n, dtype=np.int64)
    U = np.empty(n, dtype=np.int64)
    touched = np.empty(n, dtype=np.int64)
    total = 0
    applied = 0
    for k in range(order.shape[0]):
        r = order[k]
        deg = indptr[r + 1] - indptr[r]
        if deg == 0 or deg > max_root_deg:
            continue
        # candidate = N+[r]
        nu = 0
        U[nu] = r
        nu += 1
        for p in range(indptr[r], indptr[r + 1]):
            U[nu] = indices[p]
            nu += 1
        sizeS = 0
        for i in range(nu):
            v = U[i]
            inS[v] = True
            cntS[labels[v]] += 1
            sizeS += 1
        for _s in range(sweeps):
            changed = False
            for i in range(1, nu):
                v = U[i]
                lv = labels[v]
                wS = 0
                nS = 0
                wO = 0
                nO = 0
                for p in range(indptr[v], indptr[v + 1]):
                    u = indices[p]
                    if inS[u]:
                        wS += w[p]
                        nS += 1
                    elif labels[u] == lv:
                        wO += w[p]
                        nO += 1
                if inS[v]:
                    szS = sizeS - 1
                    szO = csize[lv] - cntS[lv]
                else:
                    szS = sizeS
                    szO = csize[lv] - cntS[lv] - 1
                d = (wO - wS) + kappa * ((szS - nS) - (szO - nO))
                if d < 0 and not inS[v]:
                    inS[v] = True
                    cntS[lv] += 1
                    sizeS += 1
                    changed = True
                elif d > 0 and inS[v]:
                    inS[v] = False
                    cntS[lv] -= 1
                    sizeS -= 1
                    changed = True
            if not changed:
                break
        # exact change of the objective
        dcut = 0
        dnin = 0
        nt = 0
        for i in range(nu):
            v = U[i]
            if not inS[v]:
                continue
            lv = labels[v]
            for p in range(indptr[v], indptr[v + 1]):
                u = indices[p]
                old_cut = labels[u] != lv
                if inS[u]:
                    if u < v:
                        continue
                    new_cut = False
                else:
                    new_cut = True
                if new_cut != old_cut:
                    if new_cut:
                        dcut += w[p]
                        dnin -= 1
                    else:
                        dcut -= w[p]
                        dnin += 1
        dpairs = sizeS * (sizeS - 1) // 2
        for i in range(nu):
            v = U[i]
            c = labels[v]
            if cntS[c] > 0:
                a = cntS[c]
                s = csize[c]
                dpairs += (s - a) * (s - a - 1) // 2 - s * (s - 1) // 2
                touched[nt] = c
                nt += 1
                cntS[c] = -a  # mark as processed (restored below)
        for t in range(nt):
            cntS[touched[t]] = -cntS[touched[t]]
        delta = dcut + kappa * (dpairs - dnin)
        if delta < 0:
            for t in range(nt):
                c = touched[t]
                csize[c] -= cntS[c]
                if csize[c] == 0:
                    free[nfree] = c
                    nfree += 1
            nfree -= 1
            L = free[nfree]
            for i in range(nu):
                v = U[i]
                if inS[v]:
                    labels[v] = L
            csize[L] = sizeS
            total += delta
            applied += 1
        for t in range(nt):
            cntS[touched[t]] = 0
        for i in range(nu):
            inS[U[i]] = False
            cntS[labels[U[i]]] = 0
    return total, applied


def insertion_search(g: Graph, labels: np.ndarray, rng=None, weights=None, kappa: int = 1,
                     sweeps: int = 3, max_rounds: int = 20, max_root_deg: int = 1 << 30,
                     multilevel_between: bool = True) -> np.ndarray:
    """Alternate cluster-insertion rounds with multilevel vertex/merge moves.

    Never increases the (weighted) objective."""
    rng = np.random.default_rng(rng)
    lab, _ = _compact(np.asarray(labels))
    lab = lab.copy()
    w = (np.ones(g.indices.shape[0], dtype=np.int64) if weights is None
         else np.asarray(weights, dtype=np.int64))
    cur = weighted_cost(g, lab, weights, kappa)
    for _ in range(max_rounds):
        order = rng.permutation(g.n).astype(np.int64)
        delta, applied = _insertion_round(g.n, g.indptr, g.indices, w, int(kappa), lab,
                                          order, sweeps, max_root_deg)
        if multilevel_between:
            lab = multilevel(g, lab, rng, weights=weights, kappa=kappa)
        lab, _ = _compact(lab)
        lab = lab.copy()
        new = weighted_cost(g, lab, weights, kappa)
        if new >= cur:
            break
        cur = new
    return lab
