import os
import sys

import numpy as np
import pytest

import ccbench as cc
from ccbench.generators import planted_partition
from ccbench.certiflip import certiflip
from ccbench.lns import separate_far

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))
import check_certificate as C  # noqa: E402


def _cert(tmp_path, seed=0):
    g, _ = planted_partition(0, 0, 0.7, 0.03, rng=seed, sizes=np.full(12, 9))
    path = str(tmp_path / "c.npz")
    res = certiflip(g, time_limit=8, rng=seed, cert_path=path)
    return g, res, dict(np.load(path))


def test_certificate_matches_bound(tmp_path):
    g, res, cert = _cert(tmp_path)
    rep = C.check(g, cert)
    assert abs(rep["lb_exact"] - res.lower_bound) < 1e-4
    assert rep["certified"] <= res.cost


def test_certificate_rejects_invalid_row(tmp_path):
    g, _, cert = _cert(tmp_path, seed=1)
    cert["b"] = cert["b"].copy()
    cert["b"][0] += 1.0
    with pytest.raises(ValueError):
        C.check(g, cert)


def test_certificate_rejects_far_pair(tmp_path):
    g, _, cert = _cert(tmp_path, seed=2)
    lab = cc.pivot(g, 0)
    # a non-adjacent pair without common neighbour
    far = None
    for u in range(g.n):
        nu = set(g.indices[g.indptr[u]:g.indptr[u + 1]].tolist())
        for v in range(u + 1, g.n):
            nv = set(g.indices[g.indptr[v]:g.indptr[v + 1]].tolist())
            if v not in nu and not (nu & nv):
                far = (u, v)
                break
        if far:
            break
    assert far is not None and lab is not None
    cert["u"] = cert["u"].copy()
    cert["v"] = cert["v"].copy()
    cert["u"][0], cert["v"][0] = far
    with pytest.raises(ValueError):
        C.check(g, cert)


def test_separate_far_improves_and_separates():
    g, _ = planted_partition(0, 0, 0.6, 0.05, rng=4, sizes=np.full(10, 8))
    lab = np.zeros(g.n, dtype=np.int64)  # one cluster: joins far pairs
    out = separate_far(g, lab)
    assert cc.cost(g, out) < cc.cost(g, lab)
    for c in np.unique(out):
        S = np.flatnonzero(out == c)
        for v in S:
            k = np.isin(g.indices[g.indptr[v]:g.indptr[v + 1]], S).sum()
            assert 2 * k >= len(S) - 1
