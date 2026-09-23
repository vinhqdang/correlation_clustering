import numpy as np

import ccbench as cc
from ccbench.generators import planted_partition, sparse_planted
from ccbench.support import build_support
from ccbench.lp import SparseLP
from ccbench.blockdual import BlockDualBound, metis_blocks
from ccbench.lns import gap_map, reoptimize_block, gap_lns
from ccbench.certiflip import certiflip

from test_core import brute_force_opt, small_instances


def test_block_bound_whole_graph_equals_lp():
    g, _ = planted_partition(0, 0, 0.7, 0.05, rng=1, sizes=np.full(6, 10))
    sup = build_support(g)
    lp = SparseLP(g, sup)
    v = lp.solve()
    bd = BlockDualBound(g, sup)
    bd.solve_block(np.arange(g.n), stars=False)
    assert abs(bd.bound() - v) < 1e-5


def test_block_bound_valid_and_monotone():
    for g in small_instances(10, 9):
        opt = brute_force_opt(g)
        bd = BlockDualBound(g)
        last = bd.bound()
        rng = np.random.default_rng(0)
        for s in range(4):
            bd.sweep(rng.integers(0, 2, g.n))
            assert bd.bound() >= last - 1e-7
            last = bd.bound()
        bd.solve_block(np.arange(g.n))
        assert bd.bound() <= opt + 1e-6


def test_gap_decomposition_is_exact():
    g, _ = sparse_planted(800, 6, 2.0, 0.7, rng=3)
    bd = BlockDualBound(g)
    for s in range(3):
        bd.sweep(metis_blocks(g, 300, s))
    lab = cc.multilevel(g, cc.pivot(g, 0), 0)
    vg, pair_part, row_part = gap_map(bd, lab)
    total = cc.cost(g, lab) - bd.bound()
    assert abs(pair_part + row_part - total) < 1e-5
    assert abs(vg.sum() - total) < 1e-5
    assert (vg >= -1e-9).all()


def test_reoptimize_block_is_optimal_on_whole_graph():
    for g in small_instances(6, 8):
        opt = brute_force_opt(g)
        lab = cc.pivot(g, 0)
        new = reoptimize_block(g, lab, np.arange(g.n), 30)
        assert cc.cost(g, new) == opt


def test_lns_and_certiflip_monotone():
    g, _ = sparse_planted(1500, 6, 2.0, 0.7, rng=4)
    lab = cc.multilevel(g, cc.pivot(g, 0), 0)
    new = gap_lns(g, lab, None, size=25, iters=20, rng=0)
    assert cc.cost(g, new) <= cc.cost(g, lab)
    r = certiflip(g, 20, rng=0)
    assert r.lower_bound <= r.cost
    assert cc.cost(g, r.labels) == r.cost


def test_star_packing_valid_and_installed_exactly():
    from ccbench.dual import star_packing_ls
    from ccbench.blockdual import add_star_packing
    for g in small_instances(10, 9):
        opt = brute_force_opt(g)
        bd = BlockDualBound(g)
        v, rp, rpairs, rk = star_packing_ls(g, bd.sup, pgraph=(bd.ptr, bd.idx, bd.pid),
                                            iters=2000, seed=1)
        assert v <= opt
        # pair-disjointness
        assert len(np.unique(rpairs)) == len(rpairs)
        add_star_packing(bd, rp, rpairs, rk)
        assert abs(bd.bound() - v) < 1e-9
        bd.solve_block(np.arange(g.n))
        assert v - 1e-6 <= bd.bound() <= opt + 1e-6


def test_star_packing_ls_monotone_in_iterations():
    from ccbench.dual import star_packing_ls
    g, _ = sparse_planted(2000, 6, 3.0, 0.6, rng=7)
    v0 = star_packing_ls(g, iters=0)[0]
    v1 = star_packing_ls(g, iters=20000)[0]
    assert v1 >= v0
