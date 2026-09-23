import itertools

import numpy as np
import pytest

import ccbench as cc
from ccbench.generators import planted_partition, sparse_planted
from ccbench.support import build_support, wedge_stats, pair_id
from ccbench.dual import LagrangianBound, PackingBound, packing_bound
from ccbench.lp import SparseLP, lp_pivot, repair_metric, lp_cost
from ccbench.insertion import insertion_search
from ccbench.flip import iterated_flip, pivot3
from ccbench.localsearch import weighted_cost


def brute_force_opt(g):
    """Exact optimum by enumerating set partitions (n <= 9)."""
    n = g.n
    best = None

    def partitions(i, lab, k):
        nonlocal best
        if i == n:
            c = cc.cost(g, np.array(lab))
            if best is None or c < best:
                best = c
            return
        for j in range(k + 1):
            lab.append(j)
            partitions(i + 1, lab, max(k, j + 1))
            lab.pop()

    partitions(0, [], 0)
    return best


def small_instances(count=12, n=8):
    rng = np.random.default_rng(0)
    for _ in range(count):
        p = rng.uniform(0.2, 0.7)
        iu = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < p]
        yield cc.from_edges(n, np.array(iu, dtype=np.int64).reshape(-1, 2))


def naive_cost(g, lab):
    E = set(map(tuple, g.edges().tolist()))
    c = 0
    for i, j in itertools.combinations(range(g.n), 2):
        pos = (i, j) in E
        if pos and lab[i] != lab[j]:
            c += 1
        if not pos and lab[i] == lab[j]:
            c += 1
    return c


def test_cost_matches_naive():
    rng = np.random.default_rng(1)
    for g in small_instances(5, 12):
        lab = rng.integers(0, 4, g.n)
        assert cc.cost(g, lab) == naive_cost(g, lab)


def test_weighted_cost_unit():
    g, _ = planted_partition(0, 0, 0.6, 0.1, rng=2, sizes=np.full(5, 8))
    lab = cc.pivot(g, 0)
    assert weighted_cost(g, lab, None) == cc.cost(g, lab)
    w = np.full(g.indices.shape[0], 2)
    assert weighted_cost(g, lab, w, 2) == 2 * cc.cost(g, lab)


def test_bounds_and_heuristics_bracket_opt():
    for g in small_instances():
        opt = brute_force_opt(g)
        sup = build_support(g)
        lp = SparseLP(g, sup)
        v = lp.solve()
        assert v <= opt + 1e-6
        ilp = SparseLP(g, sup, integral=True)
        vi = ilp.solve()
        assert abs(vi - opt) < 1e-6
        assert cc.cost(g, ilp.clustering()) == opt
        L = LagrangianBound(g, sup)
        assert L.run(50) <= opt + 1e-6
        P = PackingBound(g, sup)
        assert P.run(0.1) <= v + 1e-6
        assert P.upper >= v - 1e-6 or P.upper == float("inf") or wedge_stats(g)[1] == 0
        assert packing_bound(g, sup) <= opt
        for s in range(3):
            assert cc.cost(g, cc.multilevel(g, cc.pivot(g, s), s)) >= opt
            assert cc.cost(g, insertion_search(g, cc.pivot(g, s), s)) >= opt
            assert cc.cost(g, iterated_flip(g, cc.pivot(g, s), 2, rng=s)) >= opt


def test_local_search_monotone():
    g, _ = sparse_planted(3000, 6, 3.0, 0.6, rng=3)
    for s in range(3):
        p = cc.pivot(g, s)
        c0 = cc.cost(g, p)
        c1 = cc.cost(g, cc.local_search(g, p, s))
        c2 = cc.cost(g, cc.multilevel(g, p, s))
        c3 = cc.cost(g, insertion_search(g, p, s))
        assert c1 <= c0 and c2 <= c0 and c3 <= c0


def test_diameter_two_lemma_on_ilp_solutions():
    """Every optimal cluster has G+-diameter <= 2 (all intra pairs lie in P)."""
    for g in small_instances(8, 9):
        opt = brute_force_opt(g)
        sup = build_support(g)
        ilp = SparseLP(g, sup, integral=True)
        assert abs(ilp.solve() - opt) < 1e-6


def test_lp_solution_is_metric_after_repair():
    g, _ = planted_partition(0, 0, 0.7, 0.05, rng=4, sizes=np.full(6, 10))
    sup = build_support(g)
    lp = SparseLP(g, sup)
    v = lp.solve()
    d, deficit, far = repair_metric(g, sup, lp.x)
    # an LP-feasible point is already a metric: repair changes nothing
    assert far == 0
    assert abs(lp_cost(sup, d) - v) < 1e-6
    for s in range(5):
        lab = lp_pivot(g, sup, lp.x, "cmsy", s)
        assert cc.cost(g, lab) >= v - 1e-6


def test_pair_ids_roundtrip():
    g, _ = planted_partition(0, 0, 0.5, 0.05, rng=5, sizes=np.full(4, 10))
    sup = build_support(g)
    for e in range(sup.npairs):
        u, v = int(sup.pu[e]), int(sup.pv[e])
        assert pair_id(g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id, u, v) == e
        assert pair_id(g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id, v, u) == e


def test_pivot3_partition():
    rng = np.random.default_rng(0)
    a, b, c = (rng.integers(0, 5, 200) for _ in range(3))
    lab = pivot3(a, b, c)
    assert lab.min() == 0 and (lab >= 0).all()
    same = pivot3(a, a, a)
    # identical inputs reproduce the clustering
    assert len(np.unique(same)) == len(np.unique(a))
