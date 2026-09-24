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
        # calibrate the iteration rate (after compiling)
        _anneal(g.n, g.indptr, g.indices, lab.copy(), 10, t_end, t_end, p_single, 1, c0)
        t = time.time()
        _anneal(g.n, g.indptr, g.indices, lab.copy(), 200000, t_end, t_end, p_single, 1, c0)
        rate = 200000 / max(time.time() - t, 1e-3)
        iters = int(rate * time_limit)
    best, bc = _anneal(g.n, g.indptr, g.indices, lab, iters, t_start, t_end, p_single,
                       int(rng.integers(1 << 30)), c0)
    assert bc == cost(g, best)
    return best


# --------------------------------------------------------------------------
# weighted annealing (nodes of size s, edge weights = positive-pair counts)
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _anneal_w(n, indptr, indices, w, size, labels, iters, t_start, t_end, p_single, p_best,
              seed, cur_cost, nodes=np.empty(0, dtype=np.int64), p_swap=0.0):
    """Annealing where node i stands for size[i] vertices and w[p] counts the
    positive pairs between the endpoint sets.  Moving i (size s) from A to B
    changes the cost by s (S_B - S_A + s) + 2 (w(i, A - i) - w(i, B)).
    Proposals: move to a new singleton (p_single), best move (p_best), swap
    with a random neighbour in another cluster (p_swap), else move to the
    cluster of a random neighbour.  If ``nodes`` is non-empty only these nodes
    are moved.  The best clustering is recovered with an undo journal, so
    recording a new best costs O(1)."""
    np.random.seed(seed)
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
    touched = np.empty(n, dtype=np.int64)
    best = labels.copy()
    stored = True          # best[] holds the best state (else: labels minus journal)
    jv = np.empty(n, dtype=np.int64)
    jl = np.empty(n, dtype=np.int64)
    jn = 0
    best_cost = cur_cost
    cost_ = cur_cost
    lt0 = np.log(t_start)
    lt1 = np.log(t_end)
    cand = np.empty(n, dtype=np.int64)
    nc = 0
    if nodes.shape[0] > 0:
        for v in nodes:
            if indptr[v + 1] > indptr[v]:
                cand[nc] = v
                nc += 1
    else:
        for v in range(n):
            if indptr[v + 1] > indptr[v]:
                cand[nc] = v
                nc += 1
    if nc == 0:
        return best, best_cost
    T = t_start
    for it in range(iters):
        if it % 1024 == 0:
            T = np.exp(lt0 + (lt1 - lt0) * it / iters)
        v = cand[np.random.randint(nc)]
        a = labels[v]
        s = size[v]
        r = np.random.random()
        if r < p_swap:
            u = indices[indptr[v] + np.random.randint(indptr[v + 1] - indptr[v])]
            b = labels[u]
            if b == a:
                continue
            su = size[u]
            wva = 0
            wvb = 0
            wuv = 0
            for p in range(indptr[v], indptr[v + 1]):
                x = indices[p]
                lx = labels[x]
                if lx == a:
                    wva += w[p]
                elif lx == b:
                    wvb += w[p]
                    if x == u:
                        wuv = w[p]
            wua = 0
            wub = 0
            for p in range(indptr[u], indptr[u + 1]):
                lx = labels[indices[p]]
                if lx == a:
                    wua += w[p]
                elif lx == b:
                    wub += w[p]
            d1 = s * (csize[b] - csize[a] + s) + 2 * (wva - wvb)
            sa2 = csize[a] - s
            sb2 = csize[b] + s
            d2 = su * (sa2 - sb2 + su) + 2 * ((wub + wuv) - (wua - wuv))
            delta = d1 + d2
            if delta <= 0 or np.random.random() < np.exp(-delta / T):
                if jn + 2 > n:
                    if not stored:
                        best[:] = labels
                        for q in range(jn - 1, -1, -1):
                            best[jv[q]] = jl[q]
                        stored = True
                    jn = 0
                jv[jn] = v
                jl[jn] = a
                jv[jn + 1] = u
                jl[jn + 1] = b
                jn += 2
                csize[a] += su - s
                csize[b] += s - su
                labels[v] = b
                labels[u] = a
                cost_ += delta
                if cost_ < best_cost:
                    best_cost = cost_
                    stored = False
                    jn = 0
            continue
        # neighbour cluster weights
        nt = 0
        for p in range(indptr[v], indptr[v + 1]):
            c = labels[indices[p]]
            if acc[c] == 0:
                touched[nt] = c
                nt += 1
            acc[c] += w[p]
        wa = acc[a]
        base = -s * (csize[a] - s) + 2 * wa
        if r < p_swap + p_single:
            if csize[a] == s:
                for t in range(nt):
                    acc[touched[t]] = 0
                continue
            b = -1
            delta = base
        elif r < p_swap + p_single + p_best:
            b = -1
            delta = base if csize[a] > s else 1 << 40
            for t in range(nt):
                c = touched[t]
                if c == a:
                    continue
                d = s * csize[c] + base - 2 * acc[c]
                if d < delta:
                    delta = d
                    b = c
            if b == -1 and delta >= (1 << 40):
                for t in range(nt):
                    acc[touched[t]] = 0
                continue
        else:
            u = indices[indptr[v] + np.random.randint(indptr[v + 1] - indptr[v])]
            b = labels[u]
            if b == a:
                for t in range(nt):
                    acc[touched[t]] = 0
                continue
            delta = s * csize[b] + base - 2 * acc[b]
        for t in range(nt):
            acc[touched[t]] = 0
        if delta <= 0 or np.random.random() < np.exp(-delta / T):
            if b == -1:
                nfree -= 1
                b = free[nfree]
            if jn + 1 > n:
                if not stored:
                    best[:] = labels
                    for q in range(jn - 1, -1, -1):
                        best[jv[q]] = jl[q]
                    stored = True
                jn = 0
            jv[jn] = v
            jl[jn] = a
            jn += 1
            csize[a] -= s
            if csize[a] == 0:
                free[nfree] = a
                nfree += 1
            csize[b] += s
            labels[v] = b
            cost_ += delta
            if cost_ < best_cost:
                best_cost = cost_
                stored = False
                jn = 0
    if not stored:
        best[:] = labels
        for q in range(jn - 1, -1, -1):
            best[jv[q]] = jl[q]
    return best, best_cost


def _rate(fn, *args, pilot=200000):
    t = time.time()
    fn(*args[:6], pilot, *args[7:])
    return pilot / max(time.time() - t, 1e-3)


_RATE = {}


def anneal_w(wg, labels: np.ndarray, time_limit: float = 60.0, t_start: float = 0.6,
             t_end: float = 0.03, p_single: float = 0.03, p_best: float = 0.3, rng=None,
             nodes: np.ndarray | None = None, max_sweeps: float | None = None,
             p_swap: float = 0.0, last: bool = False):
    """Node-level annealing with greedy-biased proposals on a weighted instance
    (:class:`ccbench.reduce.WGraph`).  With ``nodes`` only these nodes move; the
    number of proposals is rate * time_limit, capped at max_sweeps * len(nodes).
    Returns the best clustering seen, or with ``last`` the final state."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    c0 = wg.cost(lab)
    args = (wg.n, wg.indptr, wg.indices, wg.w, wg.size)
    key = (id(wg.indptr), wg.n)
    if key not in _RATE:
        _anneal_w(*args, lab.copy(), 10, t_end, t_end, p_single, p_best, 1, c0)  # compile
        t = time.time()
        _anneal_w(*args, lab.copy(), 100000, t_end, t_end, p_single, p_best, 1, c0)
        _RATE[key] = 100000 / max(time.time() - t, 1e-3)
    iters = int(_RATE[key] * time_limit)
    nd = np.empty(0, dtype=np.int64) if nodes is None else np.asarray(nodes, dtype=np.int64)
    if max_sweeps is not None and nodes is not None:
        iters = min(iters, int(max_sweeps * max(1, len(nd))))
    best, bc = _anneal_w(*args, lab, iters, t_start, t_end, p_single, p_best,
                         int(rng.integers(1 << 30)), c0, nd, p_swap)
    return lab if last else best


def anneal2(g: Graph, labels: np.ndarray, time_limit: float = 60.0, t_start: float = 0.6,
            t_end: float = 0.03, p_single: float = 0.03, p_best: float = 0.3, rng=None):
    """Vertex-level annealing with greedy-biased proposals."""
    from .reduce import unit
    return anneal_w(unit(g), labels, time_limit, t_start, t_end, p_single, p_best, rng)


def ml_anneal(g: Graph, labels: np.ndarray, time_limit: float = 60.0, cycles: int = 4,
              t_start: float = 0.6, t_end: float = 0.03, coarse_share: float = 0.25,
              p_best: float = 0.3, rng=None):
    """Multilevel annealing V-cycles: vertex level, then the contracted graph of
    the current clusters (moves = cluster merges / cluster moves), repeated."""
    from .localsearch import _aggregate
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    t0 = time.time()
    per = time_limit / cycles
    for k in range(cycles):
        left = time_limit - (time.time() - t0)
        if left <= 0.5:
            break
        lab = anneal2(g, lab, time_limit=min(per, left) * (1 - coarse_share), t_start=t_start,
                      t_end=t_end, p_best=p_best, rng=rng)
        grp, kk = _compact(lab)
        one = np.ones(g.indices.shape[0], dtype=np.int64)
        ip, ix, w, cn, sz, sw = _aggregate(g.n, g.indptr, g.indices, one, one,
                                           np.ones(g.n, np.int64), np.zeros(g.n, np.int64),
                                           grp, kk)
        init = np.arange(kk, dtype=np.int64)
        c0 = cost(g, lab)
        t = time.time()
        _anneal_w(kk, ip, ix, w, sz, init.copy(), 20000, t_end, t_end, 0.0, p_best, 1, c0)
        rate = 20000 / max(time.time() - t, 1e-3)
        it = int(rate * min(per, max(0.0, time_limit - (time.time() - t0))) * coarse_share)
        cl, cc_ = _anneal_w(kk, ip, ix, w, sz, init, it, t_start * 2, t_end, 0.0, p_best,
                            int(rng.integers(1 << 30)), c0)
        new = cl[grp]
        if cost(g, new) <= cost(g, lab):
            lab = _compact(new)[0].copy()
    return lab


# --------------------------------------------------------------------------
# localized iterated annealing with rollback
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _mv(v, b, s, labels, csize, free, nfree):
    """Move node v (size s) to cluster b (-1: a new cluster).  Returns
    (b, nfree, new_b, emptied_a) so that the move can be undone exactly."""
    a = labels[v]
    newb = False
    if b == -1:
        nfree -= 1
        b = free[nfree]
        newb = True
    csize[a] -= s
    emptied = False
    if csize[a] == 0:
        free[nfree] = a
        nfree += 1
        emptied = True
    csize[b] += s
    labels[v] = b
    return b, nfree, newb, emptied


@nb.njit(cache=True)
def _local_ils(n, indptr, indices, w, size, labels, steps, region_size, sweeps, t_hot, t_cold,
               p_single, p_best, seed, cur_cost, deg_bias=0.0):
    """Localized iterated annealing.  Every step grows a BFS region of at most
    ``region_size`` nodes around a random node, anneals the nodes of the region
    from t_hot to t_cold (sweeps * |region| proposals, the rest of the graph
    fixed), finishes with greedy best moves, and keeps the result only if the
    cost did not increase; otherwise all moves are undone exactly.  The cost
    never increases.  With probability ``deg_bias`` the centre is drawn with
    probability proportional to its degree (the endpoint of a random edge slot)."""
    np.random.seed(seed)
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
    touched = np.empty(n, dtype=np.int64)
    mark = np.zeros(n, dtype=np.int64)
    region = np.empty(region_size, dtype=np.int64)
    cap = sweeps * region_size + 4 * region_size + 8
    jv = np.empty(cap, dtype=np.int64)
    ja = np.empty(cap, dtype=np.int64)
    jb = np.empty(cap, dtype=np.int64)
    jnew = np.empty(cap, dtype=np.bool_)
    jemp = np.empty(cap, dtype=np.bool_)
    cand = np.empty(n, dtype=np.int64)
    nc = 0
    for v in range(n):
        if indptr[v + 1] > indptr[v]:
            cand[nc] = v
            nc += 1
    if nc == 0:
        return labels, cur_cost, 0
    lr = np.log(t_cold / t_hot)
    accepted = 0
    for step in range(steps):
        stamp = step + 1
        if deg_bias > 0.0 and np.random.random() < deg_bias:
            v0 = indices[np.random.randint(indptr[n])]
        else:
            v0 = cand[np.random.randint(nc)]
        region[0] = v0
        mark[v0] = stamp
        nr = 1
        head = 0
        while head < nr and nr < region_size:
            u = region[head]
            head += 1
            deg = indptr[u + 1] - indptr[u]
            off = np.random.randint(deg) if deg > 0 else 0
            for k in range(deg):
                x = indices[indptr[u] + (off + k) % deg]
                if mark[x] != stamp:
                    mark[x] = stamp
                    region[nr] = x
                    nr += 1
                    if nr >= region_size:
                        break
        jn = 0
        dtot = 0
        iters = sweeps * nr
        for it in range(iters + 2 * nr):
            greedy = it >= iters
            if greedy:
                v = region[(it - iters) % nr]
            else:
                v = region[np.random.randint(nr)]
            if indptr[v + 1] == indptr[v]:
                continue
            a = labels[v]
            s = size[v]
            nt = 0
            for p in range(indptr[v], indptr[v + 1]):
                c = labels[indices[p]]
                if acc[c] == 0:
                    touched[nt] = c
                    nt += 1
                acc[c] += w[p]
            base = -s * (csize[a] - s) + 2 * acc[a]
            r = np.random.random()
            b = -2
            delta = 0
            if greedy or (r >= p_single and r < p_single + p_best):
                b = -2
                delta = 1 << 40
                if csize[a] > s:
                    b = -1
                    delta = base
                for t in range(nt):
                    c = touched[t]
                    if c == a:
                        continue
                    d = s * csize[c] + base - 2 * acc[c]
                    if d < delta:
                        delta = d
                        b = c
            elif r < p_single:
                if csize[a] > s:
                    b = -1
                    delta = base
            else:
                u = indices[indptr[v] + np.random.randint(indptr[v + 1] - indptr[v])]
                if labels[u] != a:
                    b = labels[u]
                    delta = s * csize[b] + base - 2 * acc[b]
            for t in range(nt):
                acc[touched[t]] = 0
            if b == -2 or jn >= cap:
                continue
            if greedy:
                ok = delta < 0
            else:
                T = t_hot * np.exp(lr * it / iters)
                ok = delta <= 0 or np.random.random() < np.exp(-delta / T)
            if ok:
                bb, nfree, nw, em = _mv(v, b, s, labels, csize, free, nfree)
                jv[jn] = v
                ja[jn] = a
                jb[jn] = bb
                jnew[jn] = nw
                jemp[jn] = em
                jn += 1
                dtot += delta
        if dtot > 0:
            for q in range(jn - 1, -1, -1):
                v = jv[q]
                a = ja[q]
                b = jb[q]
                s = size[v]
                if jemp[q]:
                    nfree -= 1
                csize[a] += s
                csize[b] -= s
                labels[v] = a
                if jnew[q]:
                    free[nfree] = b
                    nfree += 1
        else:
            cur_cost += dtot
            if dtot < 0:
                accepted += 1
    return labels, cur_cost, accepted


def local_ils(wg, labels: np.ndarray, time_limit: float = 10.0, region_size: int = 30,
              sweeps: int = 20, t_hot: float = 1.0, t_cold: float = 0.05,
              p_single: float = 0.03, p_best: float = 0.3, rng=None, chunk: int = 2000,
              deg_bias: float = 0.0):
    """Localized iterated annealing with exact rollback (never increases the cost)."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    c = wg.cost(lab)
    t0 = time.time()
    while time.time() - t0 < time_limit:
        lab, c, _ = _local_ils(wg.n, wg.indptr, wg.indices, wg.w, wg.size, lab, chunk,
                               region_size, sweeps, t_hot, t_cold, p_single, p_best,
                               int(rng.integers(1 << 30)), c, deg_bias)
    return lab
