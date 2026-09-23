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
