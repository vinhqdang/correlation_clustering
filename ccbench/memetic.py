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
              max_sweeps: float = 3000.0, p_swap: float = 0.0, p_self: float = 0.0,
              px_accept: bool = False, temps=(0.3, 0.6, 1.2, 2.0), explore: float = 0.2,
              adaptive: bool = False):
    """Memetic search with annealing as improvement operator on a weighted
    instance.  ``init`` are starting clusterings of the nodes; each is annealed
    for init_share * time_limit / pop_size seconds to form the population.
    Offspring: overlay recombination of two parents followed by a short, cooler
    annealing run; crowding replacement keeps the population diverse.
    With ``px_accept`` the offspring annealing runs to its final state and is
    accepted region-wise (partition crossover with its starting point).  With
    probability ``p_self`` a generation is instead a PX-annealing step on one
    member (reheating temperature chosen from ``temps`` by a bandit); with
    ``adaptive`` the choice between the two operators is also made by an
    epsilon-greedy bandit on the gain per second."""
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
        lab = anneal_w(wg, lab, time_limit=max(1.0, t_init), rng=rng, p_swap=p_swap)
        lab = _compact(lab)[0].copy()
        pop.append(lab)
        costs.append(wg.cost(lab))
        keys.append(_key(lab))
        rec()
        if verbose:
            print(f"init {i}: {costs[-1]}  t={time.time() - t0:.0f}", flush=True)
    gen = 0
    score = np.full(len(temps), np.inf)
    op_score = np.full(2, np.inf)  # 0: crossover, 1: self step
    while time.time() - t0 < time_limit:
        gen += 1
        P = len(pop)
        if adaptive:
            self_step = (rng.random() < 0.5) if rng.random() < explore else \
                bool(np.argmax(op_score))
        else:
            self_step = rng.random() < p_self
        ts = time.time()
        if self_step:
            i = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
            k = int(rng.integers(len(temps))) if rng.random() < explore else int(np.argmax(score))
            ts = time.time()
            left = time_limit - (time.time() - t0)
            y = anneal_w(wg, pop[i], time_limit=max(0.5, min(child_share * time_limit, left)),
                         t_start=temps[k], rng=rng, p_swap=p_swap, last=True)
            child = partition_crossover(wg, pop[i], y)[0]
            c = wg.cost(child)
            gain = (costs[i] - c) / max(time.time() - ts, 1e-3)
            score[k] = gain if not np.isfinite(score[k]) else 0.6 * score[k] + 0.4 * gain
            op_score[1] = gain if not np.isfinite(op_score[1]) else 0.7 * op_score[1] + 0.3 * gain
            if c < costs[i]:
                pop[i], costs[i], keys[i] = child, c, _key(child)
            rec()
            if verbose:
                print(f"gen {gen}: self T={temps[k]} {c} best {min(costs)} "
                      f"t={time.time() - t0:.0f}", flush=True)
            continue
        i1 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        i2 = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
        if i1 == i2:
            i2 = (i1 + 1 + int(rng.integers(P - 1))) % P
        start = partition_crossover(wg, pop[i1], pop[i2])[0] if px else None
        child = overlay_combine_w(wg, pop[i1], pop[i2], rng, start)
        left = time_limit - (time.time() - t0)
        region = _diff_region(wg, pop[i1], pop[i2]) if local else None
        y = anneal_w(wg, child, time_limit=max(0.5, min(child_share * time_limit, left)),
                     t_start=t_child, rng=rng, nodes=region, max_sweeps=max_sweeps,
                     p_swap=p_swap, last=px_accept)
        child = partition_crossover(wg, child, y)[0] if px_accept else _compact(y)[0].copy()
        c = wg.cost(child)
        gain = max(0, min(costs[i1], costs[i2]) - c) / max(time.time() - ts, 1e-3)
        op_score[0] = gain if not np.isfinite(op_score[0]) else 0.7 * op_score[0] + 0.3 * gain
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
    search, projected onto the contracted nodes.  No operator increases the
    cost of the best member, so E[cost] <= 3 OPT (the Pivot guarantee)."""
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
        p = pivot(g, s)
        f = iterated_flip(g, p, 1, rng=s, time_limit=t_seed)[first]
        # Pivot never splits a critical clique, so its projection is exact; the
        # projection of the flip result may cost more, keep the better one
        init.append(f if wg.cost(f) <= cost(g, p) else p[first])
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
    score = np.full(len(temps), np.inf)
    op_score = np.full(2, np.inf)  # 0: crossover, 1: self step
    while time.time() - t0 < time_limit:
        gen += 1
        P = len(pop)
        if adaptive:
            self_step = (rng.random() < 0.5) if rng.random() < explore else \
                bool(np.argmax(op_score))
        else:
            self_step = rng.random() < p_self
        ts = time.time()
        if self_step:
            i = min(rng.choice(P, 2, replace=False), key=lambda i: costs[i])
            k = int(rng.integers(len(temps))) if rng.random() < explore else int(np.argmax(score))
            ts = time.time()
            left = time_limit - (time.time() - t0)
            y = anneal_w(wg, pop[i], time_limit=max(0.5, min(child_share * time_limit, left)),
                         t_start=temps[k], rng=rng, p_swap=p_swap, last=True)
            child = partition_crossover(wg, pop[i], y)[0]
            c = wg.cost(child)
            gain = (costs[i] - c) / max(time.time() - ts, 1e-3)
            score[k] = gain if not np.isfinite(score[k]) else 0.6 * score[k] + 0.4 * gain
            op_score[1] = gain if not np.isfinite(op_score[1]) else 0.7 * op_score[1] + 0.3 * gain
            if c < costs[i]:
                pop[i], costs[i], keys[i] = child, c, _key(child)
            rec()
            if verbose:
                print(f"gen {gen}: self T={temps[k]} {c} best {min(costs)} "
                      f"t={time.time() - t0:.0f}", flush=True)
            continue
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


def px_anneal(g: Graph, time_limit: float = 600.0, rng=None, step: float | None = None,
              t_reheat: float | None = None, temps=(0.3, 0.6, 1.2, 2.0), explore: float = 0.2,
              init_share: float = 0.1, p_swap: float = 0.0, verbose: bool = False,
              history: list | None = None):
    """PX-annealing: iterated annealing with region-wise acceptance.

    On the critical-clique contraction, starting from an annealed Pivot
    clustering x, repeat: anneal a copy of x from a reheating temperature down
    to the final state y, then x <- partition_crossover(x, y), which keeps,
    independently in every region where x and y differ, the better of the two.
    Improvements found anywhere are kept even when y is worse overall, and
    the cost of x never increases.  The reheating temperature is ``t_reheat``
    if given, else chosen from ``temps`` by an epsilon-greedy bandit on the
    recent gain per second."""
    from .reduce import contract
    from .anneal import anneal_w
    rng = np.random.default_rng(rng)
    t0 = time.time()
    step = float(np.clip(time_limit / 40, 3.0, 20.0)) if step is None else step
    temps = [float(t_reheat)] if t_reheat is not None else list(temps)
    score = np.full(len(temps), np.inf)  # optimistic start: try every arm once
    wg, grp = contract(g)
    first = np.full(wg.n, -1, dtype=np.int64)
    first[grp[::-1]] = np.arange(g.n)[::-1]
    x = pivot(g, int(rng.integers(1 << 30)))[first]
    x = anneal_w(wg, x, max(1.0, init_share * time_limit), rng=rng, p_swap=p_swap)
    c = wg.cost(x)
    if history is not None:
        history.append((time.time() - t0, c))
    while time.time() - t0 < time_limit:
        left = time_limit - (time.time() - t0)
        k = int(rng.integers(len(temps))) if rng.random() < explore else int(np.argmax(score))
        ts = time.time()
        y = anneal_w(wg, x, max(0.5, min(step, left)), t_start=temps[k], rng=rng,
                     p_swap=p_swap, last=True)
        x, nb = partition_crossover(wg, x, y)
        c_new = wg.cost(x)
        gain = (c - c_new) / max(time.time() - ts, 1e-3)
        score[k] = gain if not np.isfinite(score[k]) else 0.6 * score[k] + 0.4 * gain
        c = c_new
        if history is not None:
            history.append((time.time() - t0, c))
        if verbose:
            print(f"t={time.time() - t0:.0f} T={temps[k]} cost {c} regions {nb}", flush=True)
    out = x[grp]
    assert cost(g, out) == c
    return out, c


def pxmem(g: Graph, time_limit: float = 600.0, rng=None, history: list | None = None,
          verbose: bool = False, **kw):
    """PX-memetic annealing (the default configuration): memetic search on the
    critical-clique contraction with a population of 4, region-wise acceptance
    of offspring by partition crossover, PX-annealing self-improvement steps,
    and bandit selection of the operator and of the reheating temperature."""
    opts = dict(pop_size=4, init_share=0.15, px_accept=True, adaptive=True)
    opts.update(kw)
    return memetic_twin(g, time_limit, rng=rng, history=history, verbose=verbose, **opts)
