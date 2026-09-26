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


def _write_pace(g, path):
    with open(path, "w") as fh:
        fh.write(f"p cep {g.n} {g.m}\n")
        for u, v in g.edges():
            fh.write(f"{u + 1} {v + 1}\n")


def _cert(tmp_path, seed=0):
    """A certificate for a random instance, and the instance as the checker
    reads it from the file (the checker never sees the solver's graph)."""
    g, _ = planted_partition(0, 0, 0.7, 0.03, rng=seed, sizes=np.full(12, 9))
    path = str(tmp_path / "c.npz")
    res = certiflip(g, time_limit=8, rng=seed, cert_path=path)
    gr = str(tmp_path / "g.gr")
    _write_pace(g, gr)
    return C.read_instance(gr, "pace"), res, dict(np.load(path)), g


def test_reader_matches_solver_graph(tmp_path):
    inst, _, _, g = _cert(tmp_path)
    assert (inst.n, inst.m) == (g.n, g.m)
    assert np.array_equal(inst.indptr, g.indptr)
    assert np.array_equal(inst.indices, g.indices)


def test_certificate_matches_bound(tmp_path):
    inst, res, cert, _ = _cert(tmp_path)
    rep = C.check(inst, cert)
    assert abs(rep["lb_exact"] - res.lower_bound) < 1e-4
    assert rep["certified"] <= res.cost


def test_certificate_rejects_other_instance(tmp_path):
    inst, _, cert, _ = _cert(tmp_path)
    cert["edge_sha256"] = np.array("0" * 64)
    with pytest.raises(ValueError):
        C.check(inst, cert)
    cert.pop("edge_sha256")
    cert["n"] = np.array(inst.n + 1)
    with pytest.raises(ValueError):
        C.check(inst, cert)


def test_certificate_rejects_invalid_row(tmp_path):
    g, _, cert, _ = _cert(tmp_path, seed=1)
    cert["b"] = cert["b"].copy()
    cert["b"][0] += 1.0
    with pytest.raises(ValueError):
        C.check(g, cert)


def test_certificate_rejects_far_pair(tmp_path):
    g, _, cert, _ = _cert(tmp_path, seed=2)
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
    assert far is not None
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
