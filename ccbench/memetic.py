"""Memetic flipping search.

A population of clusterings produced by the iterated flipping local search is
evolved with
  * overlay recombination: vertices that two parents put together in both are
    contracted into one node; the multilevel local search then optimises the
    contracted instance starting from the better parent, so it can only move
    whole agreement groups (a strong, structure-preserving crossover);
  * flip mutation: the positive pairs cut by the child get weight 3/2 (the flip
    of Cohen-Addad et al.), the weighted objective is re-optimised and the
    result is polished on the original objective;
  * cluster-insertion local search on every offspring.
Every offspring is a local optimum of the original objective; the best
solution never gets worse, so the expected-cost guarantee of the Pivot seeds
(3 OPT) is kept.
"""
from __future__ import annotations

import time

import numpy as np

from .graph import Graph
from .objective import cost
from .pivot import pivot
from .flip import iterated_flip, cut_indicator
from .insertion import insertion_search
from .localsearch import _compact, _aggregate, _local_move, multilevel


def _ml_weighted(n, indptr, indices, w, cnt, size, self_w, init, rng, max_levels=50):
    """Multilevel local search on an already weighted/contracted graph."""
    maps = []
    cur = _compact(init)[0].copy()
    ip, ix, ww, cc, sz, sw = indptr, indices, w, cnt, size, self_w
    nn = n
    for _ in range(max_levels):
        _local_move(nn, ip, ix, ww, cc, 1, sz, cur, rng.permutation(nn).astype(np.int64),
                    1 << 62)
        cur, k = _compact(cur)
        maps.append(cur)
        if k == nn:
            break
        ip, ix, ww, cc, sz, sw = _aggregate(nn, ip, ix, ww, cc, sz, sw, cur, k)
        nn = k
        cur = np.arange(nn, dtype=np.int64)
    out = maps[-1]
    for mp in reversed(maps[:-1]):
        out = out[mp]
    out = _compact(np.asarray(out, dtype=np.int64))[0].copy()
    _local_move(n, indptr, indices, w, cnt, 1, size, out, rng.permutation(n).astype(np.int64),
                1 << 62)
    return out


def overlay_combine(g: Graph, a: np.ndarray, b: np.ndarray, rng) -> np.ndarray:
    """Overlay recombination of two clusterings (see module docstring)."""
    n = g.n
    a = _compact(a)[0]
    b = _compact(b)[0]
    grp, k = _compact(a.astype(np.int64) * (n + 1) + b)
    one = np.ones(g.indices.shape[0], dtype=np.int64)
    ip, ix, w, cn, sz, sw = _aggregate(n, g.indptr, g.indices, one, one, np.ones(n, np.int64),
                                       np.zeros(n, np.int64), grp, k)
    first = np.full(k, -1, dtype=np.int64)
    first[grp[::-1]] = np.arange(n)[::-1]
    start = a if cost(g, a) <= cost(g, b) else b
    init = start[first]
    lab_c = _ml_weighted(k, ip, ix, w, cn, sz, sw, init, rng)
    return lab_c[grp]


def flip_mutation(g: Graph, lab: np.ndarray, rng) -> np.ndarray:
    base = np.full(g.indices.shape[0], 2, dtype=np.int64)
    w = base + cut_indicator(g, lab)
    m = insertion_search(g, lab, rng, weights=w, kappa=2)
    return insertion_search(g, m, rng)


def ruin_mutation(g: Graph, lab: np.ndarray, rng, size: int = 50) -> np.ndarray:
    """Dissolve a random connected region into singletons and repair."""
    lab = _compact(lab)[0].copy()
    v0 = int(rng.integers(g.n))
    seen = {v0}
    frontier = [v0]
    while frontier and len(seen) < size:
        u = frontier.pop(int(rng.integers(len(frontier))))
        for w in g.indices[g.indptr[u]:g.indptr[u + 1]]:
            w = int(w)
            if w not in seen:
                seen.add(w)
                frontier.append(w)
                if len(seen) >= size:
                    break
    nodes = np.fromiter(seen, dtype=np.int64)
    lab[nodes] = lab.max() + 1 + np.arange(len(nodes))
    return insertion_search(g, lab, rng)


def _key(lab):
    return hash(_compact(lab)[0].tobytes())


def _distance(a, b, sample):
    """Fraction of sampled positive edges on which two clusterings disagree."""
    return float(np.mean((a[sample[0]] == a[sample[1]]) != (b[sample[0]] == b[sample[1]])))


def memetic(g: Graph, time_limit: float = 300.0, rng=None, pop_size: int = 10,
            init_rounds: int = 2, init_fraction: float = 0.25, init: list | None = None,
            verbose: bool = False, history: list | None = None,
            p_cross: float = 0.6, p_ruin: float = 0.2):
    rng = np.random.default_rng(rng)
    t0 = time.time()
    pop, costs, keys = [], [], []
    E = g.edges()
    samp = E[rng.choice(len(E), size=min(len(E), 20000), replace=False)].T if len(E) else None

    def record():
        if history is not None:
            history.append((time.time() - t0, min(costs)))

    def add(lab):
        lab = _compact(lab)[0].copy()
        k = _key(lab)
        if k in keys:
            return False
        pop.append(lab)
        costs.append(cost(g, lab))
        keys.append(k)
        record()
        return True

    for lab in init or []:
        add(lab)
    tries = 0
    while len(pop) < pop_size and tries < 3 * pop_size:
        tries += 1
        if len(pop) >= 2 and time.time() - t0 > init_fraction * time_limit:
            break
        s = int(rng.integers(1 << 30))
        add(iterated_flip(g, pivot(g, s), init_rounds, rng=s,
                          time_limit=max(1.0, init_fraction * time_limit / pop_size)))
    gen = 0
    while time.time() - t0 < time_limit:
        gen += 1
        P = len(pop)
        r = rng.random()
        if P >= 2 and r < p_cross:
            i1 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
            i2 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
            if i1 == i2:
                i2 = (i1 + 1 + int(rng.integers(P - 1))) % P
            child = insertion_search(g, overlay_combine(g, pop[i1], pop[i2], rng), rng)
        else:
            i1 = min(rng.choice(P, min(2, P), replace=False), key=lambda i: costs[i])
            if r < p_cross + p_ruin:
                child = ruin_mutation(g, pop[i1], rng, size=int(rng.integers(20, 200)))
            else:
                child = flip_mutation(g, pop[i1], rng)
        child = _compact(child)[0].copy()
        c = cost(g, child)
        k = _key(child)
        if k not in keys:
            d = [_distance(child, q, samp) for q in pop]
            for j in np.argsort(d)[:3]:
                if c < costs[j] or (c == costs[j] and d[j] > 0):
                    pop[j], costs[j], keys[j] = child, c, k
                    break
            else:
                worst = int(np.argmax(costs))
                if c < costs[worst]:
                    pop[worst], costs[worst], keys[worst] = child, c, k
        record()
        if verbose and gen % 20 == 0:
            print(f"gen {gen} best {min(costs)} worst {max(costs)} t={time.time() - t0:.0f}",
                  flush=True)
    best = int(np.argmin(costs))
    return pop[best], costs[best]


def partition_crossover(wg, a: np.ndarray, b: np.ndarray):
    """Exact decomposition crossover.  Link every A-cluster with every B-cluster
    it meets; in each connected component X of this bipartite graph both parents
    partition the same node set, and pairs across components are separated by
    both, so cost = sum_X cost_X(parent) + const.  Taking the better parent on
    every component gives a child of cost <= min(cost(a), cost(b)).
    Returns (child, number of components where b was taken)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    a = _compact(a)[0]
    b = _compact(b)[0]
    ka, kb = int(a.max()) + 1, int(b.max()) + 1
    M = coo_matrix((np.ones(wg.n), (a, b + ka)), shape=(ka + kb, ka + kb))
    ncomp, comp = connected_components(M, directed=False)
    cx = comp[a]  # component of every node
    src = np.repeat(np.arange(wg.n), np.diff(wg.indptr))
    dst = wg.indices
    m = (src < dst) & (cx[src] == cx[dst])
    ca = np.zeros(ncomp)
    cb = np.zeros(ncomp)
    # cut positive pairs inside a component
    np.add.at(ca, cx[src[m]], wg.w[m] * (a[src[m]] != a[dst[m]]))
    np.add.at(cb, cx[src[m]], wg.w[m] * (b[src[m]] != b[dst[m]]))
    # negative pairs joined: sum over clusters of C(S,2) - internal positive weight
    for lab, k, acc in ((a, ka, ca), (b, kb, cb)):
        S = np.bincount(lab, weights=wg.size, minlength=k)
        W = np.bincount(lab, weights=wg.self_w, minlength=k)
        inside = m & (lab[src] == lab[dst])
        W += np.bincount(lab[src[inside]], weights=wg.w[inside], minlength=k)
        rep = np.zeros(k, dtype=np.int64)
        rep[lab] = cx
        np.add.at(acc, rep, S * (S - 1) / 2 - W)
    take_b = cb < ca
    child = np.where(take_b[cx], b + ka, a)
    return _compact(child)[0].copy(), int(take_b.sum())


def _diff_region(wg, a, b, hops: int = 1):
    """Nodes incident to a positive pair on which a and b disagree, plus their
    neighbours up to ``hops``; None if that is (almost) everything."""
    src = np.repeat(np.arange(wg.n), np.diff(wg.indptr))
    dst = wg.indices
    dis = (a[src] == a[dst]) != (b[src] == b[dst])
    mark = np.zeros(wg.n, dtype=bool)
    mark[src[dis]] = True
    for _ in range(hops):
        mark[dst[mark[src]]] = True
    if mark.sum() > 0.7 * wg.n or not mark.any():
        return None
    return np.flatnonzero(mark)


def overlay_combine_w(wg, a: np.ndarray, b: np.ndarray, rng, start=None) -> np.ndarray:
    """Overlay recombination on a weighted instance (see module docstring);
    the contracted search starts from ``start`` (default: the better parent)."""
    n = wg.n
    a = _compact(a)[0]
    b = _compact(b)[0]
    grp, k = _compact(a.astype(np.int64) * (n + 1) + b)
    ip, ix, w, cn, sz, sw = _aggregate(n, wg.indptr, wg.indices, wg.w, wg.w, wg.size,
                                       wg.self_w, grp, k)
    first = np.full(k, -1, dtype=np.int64)
    first[grp[::-1]] = np.arange(n)[::-1]
    if start is None:
        start = a if wg.cost(a) <= wg.cost(b) else b
    lab_c = _ml_weighted(k, ip, ix, w, cn, sz, sw, np.asarray(start)[first], rng)
    return lab_c[grp]


def memetic_w(wg, time_limit: float = 600.0, rng=None, init: list | None = None,
              pop_size: int = 6, init_share: float = 0.4, child_share: float = 0.04,
              t_child: float = 0.3, verbose: bool = False, history: list | None = None,
              t0: float | None = None, px: bool = True, local: bool = True,
              max_sweeps: float = 3000.0):
    """Memetic search with annealing as improvement operator on a weighted
    instance.  ``init`` are starting clusterings of the nodes; each is annealed
    for init_share * time_limit / pop_size seconds to form the population.
    Offspring: overlay recombination of two parents followed by a short, cooler
    annealing run; crowding replacement keeps the population diverse."""
    from .anneal import anneal_w
    rng = np.random.default_rng(rng)
    t0 = time.time() if t0 is None else t0
    E = wg.edges()
    samp = E[rng.choice(len(E), size=min(len(E), 20000), replace=False)].T
    pop, costs, keys = [], [], []

    def rec():
        if history is not None:
            history.append((time.time() - t0, min(costs)))

    t_init = init_share * time_limit / pop_size
    for i, lab in enumerate(init):
        lab = anneal_w(wg, lab, time_limit=max(1.0, t_init), rng=rng)
        lab = _compact(lab)[0].copy()
        pop.append(lab)
        costs.append(wg.cost(lab))
        keys.append(_key(lab))
        rec()
        if verbose:
            print(f"init {i}: {costs[-1]}  t={time.time() - t0:.0f}", flush=True)
    gen = 0
    while time.time() - t0 < time_limit:
        gen += 1
        P = len(pop)
        i1 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        i2 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        if i1 == i2:
            i2 = (i1 + 1 + int(rng.integers(P - 1))) % P
        start = partition_crossover(wg, pop[i1], pop[i2])[0] if px else None
        child = overlay_combine_w(wg, pop[i1], pop[i2], rng, start)
        left = time_limit - (time.time() - t0)
        region = _diff_region(wg, pop[i1], pop[i2]) if local else None
        child = anneal_w(wg, child, time_limit=max(0.5, min(child_share * time_limit, left)),
                         t_start=t_child, rng=rng, nodes=region, max_sweeps=max_sweeps)
        child = _compact(child)[0].copy()
        c = wg.cost(child)
        k = _key(child)
        if k not in keys:
            d = [_distance(child, q, samp) for q in pop]
            for j in np.argsort(d)[:2]:
                if c < costs[j] or (c == costs[j] and d[j] > 0):
                    pop[j], costs[j], keys[j] = child, c, k
                    break
            else:
                worst = int(np.argmax(costs))
                if c < costs[worst]:
                    pop[worst], costs[worst], keys[worst] = child, c, k
        rec()
        if verbose:
            print(f"gen {gen}: child {c} best {min(costs)} worst {max(costs)} "
                  f"t={time.time() - t0:.0f}", flush=True)
    best = int(np.argmin(costs))
    return pop[best], costs[best]


def memetic_twin(g: Graph, time_limit: float = 600.0, rng=None, pop_size: int = 6,
                 verbose: bool = False, history: list | None = None, **kw):
    """memetic_w on the critical-clique contraction of g.  Seeds: Pivot (which
    never splits a critical clique) followed by one round of flipping local
    search, projected onto the contracted nodes."""
    from .reduce import contract
    rng = np.random.default_rng(rng)
    t0 = time.time()
    wg, grp = contract(g)
    first = np.full(wg.n, -1, dtype=np.int64)
    first[grp[::-1]] = np.arange(g.n)[::-1]
    t_seed = 0.2 * 0.4 * time_limit / pop_size
    init = []
    for _ in range(pop_size):
        s = int(rng.integers(1 << 30))
        init.append(iterated_flip(g, pivot(g, s), 1, rng=s, time_limit=t_seed)[first])
    lab, c = memetic_w(wg, time_limit, rng, init, pop_size=pop_size, verbose=verbose,
                       history=history, t0=t0, **kw)
    out = lab[grp]
    assert cost(g, out) == c
    return out, c


def memetic_sa(g: Graph, time_limit: float = 600.0, rng=None, pop_size: int = 6,
               init_share: float = 0.4, child_share: float = 0.04, t_child: float = 0.3,
               verbose: bool = False, history: list | None = None):
    """Memetic search whose improvement operator is simulated annealing.

    Initial population: independent (Pivot -> flip -> annealing) runs.
    Offspring: overlay recombination of two parents followed by a short, cooler
    annealing run; crowding replacement keeps the population diverse."""
    from .anneal import anneal2
    rng = np.random.default_rng(rng)
    t0 = time.time()
    E = g.edges()
    samp = E[rng.choice(len(E), size=min(len(E), 20000), replace=False)].T
    pop, costs, keys = [], [], []

    def rec():
        if history is not None:
            history.append((time.time() - t0, min(costs)))

    t_init = init_share * time_limit / pop_size
    for i in range(pop_size):
        s = int(rng.integers(1 << 30))
        lab = iterated_flip(g, pivot(g, s), 1, rng=s, time_limit=0.2 * t_init)
        lab = anneal2(g, lab, time_limit=max(1.0, t_init - 0.0), rng=s)
        lab = _compact(lab)[0].copy()
        pop.append(lab)
        costs.append(cost(g, lab))
        keys.append(_key(lab))
        rec()
        if verbose:
            print(f"init {i}: {costs[-1]}  t={time.time() - t0:.0f}", flush=True)
    gen = 0
    while time.time() - t0 < time_limit:
        gen += 1
        P = len(pop)
        i1 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        i2 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        if i1 == i2:
            i2 = (i1 + 1 + int(rng.integers(P - 1))) % P
        child = overlay_combine(g, pop[i1], pop[i2], rng)
        left = time_limit - (time.time() - t0)
        child = anneal2(g, child, time_limit=max(0.5, min(child_share * time_limit, left)),
                        t_start=t_child, rng=rng)
        child = _compact(child)[0].copy()
        c = cost(g, child)
        k = _key(child)
        if k not in keys:
            d = [_distance(child, q, samp) for q in pop]
            for j in np.argsort(d)[:2]:
                if c < costs[j] or (c == costs[j] and d[j] > 0):
                    pop[j], costs[j], keys[j] = child, c, k
                    break
            else:
                worst = int(np.argmax(costs))
                if c < costs[worst]:
                    pop[worst], costs[worst], keys[worst] = child, c, k
        rec()
        if verbose:
            print(f"gen {gen}: child {c} best {min(costs)} worst {max(costs)} "
                  f"t={time.time() - t0:.0f}", flush=True)
    best = int(np.argmin(costs))
    return pop[best], costs[best]
