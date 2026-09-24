import numpy as np

import ccbench as cc
from ccbench.generators import planted_partition
from ccbench.reduce import contract, unit
from ccbench.pivot import pivot
from ccbench.anneal import anneal_w
from ccbench.memetic import partition_crossover, memetic_twin

from test_core import brute_force_opt, small_instances


def _project(g, grp, k, lab):
    first = np.full(k, -1, dtype=np.int64)
    first[grp[::-1]] = np.arange(g.n)[::-1]
    return lab[first]


def test_contraction_cost_matches():
    g, _ = planted_partition(0, 0, 1.0, 0.01, rng=3, sizes=np.full(8, 12))
    wg, grp = contract(g)
    assert wg.n < g.n
    rng = np.random.default_rng(0)
    for s in range(5):
        lw = rng.integers(0, 10, size=wg.n)
        assert wg.cost(lw) == cc.cost(g, lw[grp])
    lab = pivot(g, 1)
    assert (_project(g, grp, wg.n, lab)[grp] == lab).all()  # Pivot keeps twins together
    assert unit(g).cost(lab) == cc.cost(g, lab)


def test_contraction_keeps_optimum():
    for g in small_instances(12, 8):
        wg, grp = contract(g)
        opt = brute_force_opt(g)
        best = min(wg.cost(np.array(p)) for p in _partitions(wg.n))
        assert best == opt


def _partitions(n):
    def rec(i, lab, k):
        if i == n:
            yield list(lab)
            return
        for c in range(k + 1):
            lab.append(c)
            yield from rec(i + 1, lab, max(k, c + 1))
            lab.pop()
    yield from rec(0, [], 0)


def test_partition_crossover_not_worse():
    g, _ = planted_partition(0, 0, 0.6, 0.02, rng=5, sizes=np.full(20, 15))
    wg, grp = contract(g)
    for s in range(4):
        a = _project(g, grp, wg.n, pivot(g, 2 * s))
        b = _project(g, grp, wg.n, pivot(g, 2 * s + 1))
        c, _ = partition_crossover(wg, a, b)
        assert wg.cost(c) <= min(wg.cost(a), wg.cost(b))
        c2 = anneal_w(wg, c, time_limit=0.2, rng=s)
        assert wg.cost(c2) <= wg.cost(c)


def test_memetic_twin_runs():
    g, _ = planted_partition(0, 0, 0.8, 0.02, rng=7, sizes=np.full(10, 10))
    lab, c = memetic_twin(g, time_limit=6, rng=0, pop_size=3)
    assert c == cc.cost(g, lab)
    assert c <= cc.cost(g, pivot(g, 0))


def test_anneal_journal_and_swaps_consistent():
    from ccbench.anneal import _anneal_w
    g, _ = planted_partition(0, 0, 0.7, 0.03, rng=11, sizes=np.full(15, 12))
    wg, grp = contract(g)
    lab = _project(g, grp, wg.n, pivot(g, 3))
    c0 = wg.cost(lab)
    for ps in (0.0, 0.3):
        for iters in (50, 5000, 200000):
            best, bc = _anneal_w(wg.n, wg.indptr, wg.indices, wg.w, wg.size, lab.copy(), iters,
                                 1.0, 0.05, 0.03, 0.3, 7, c0, np.empty(0, np.int64), ps)
            assert bc == wg.cost(best) and bc <= c0
