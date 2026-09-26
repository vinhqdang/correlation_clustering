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


def test_clustering_cost_matches_solver(tmp_path):
    g, _ = planted_partition(0, 0, 0.6, 0.1, rng=3, sizes=np.full(6, 7))
    gr = str(tmp_path / "g.gr")
    _write_pace(g, gr)
    inst = C.read_instance(gr, "pace")
    rng = np.random.default_rng(3)
    for _ in range(5):
        lab = rng.integers(0, 6, g.n)
        assert C.clustering_cost(inst, lab) == int(cc.cost(g, lab))


def test_certificate_rejects_negative_multiplier(tmp_path):
    g, _, cert, _ = _cert(tmp_path, seed=1)
    cert["y"] = cert["y"].copy()
    cert["y"][0] = -1e-9
    with pytest.raises(ValueError):
        C.check(g, cert)


def test_star_formula_is_implied_by_enumeration(tmp_path):
    """Every star row the closed form accepts is valid by exhaustive
    enumeration, and with R complete the two agree exactly."""
    g, _ = planted_partition(0, 0, 0.5, 0.08, rng=5, sizes=np.full(6, 8))
    gr = str(tmp_path / "g.gr")
    _write_pace(g, gr)
    inst = C.read_instance(gr, "pace")
    close = lambda a, c: bool(C._close(inst.indptr, inst.indices, min(a, c), max(a, c)))
    rng = np.random.default_rng(0)
    checked = 0
    for _ in range(400):
        v = int(rng.integers(g.n))
        cand = [t for t in range(g.n) if t != v and close(v, t)]
        if len(cand) < 2:
            continue
        T = rng.choice(cand, size=min(len(cand), int(rng.integers(2, 9))), replace=False).tolist()
        pairs = [(a, c) for i, a in enumerate(T) for c in T[i + 1:] if close(a, c)]
        complete = rng.random() < 0.5
        R = pairs if complete else [p for p in pairs if rng.random() < 0.5]
        us = [min(v, t) for t in T] + [min(p) for p in R]
        vs = [max(v, t) for t in T] + [max(p) for p in R]
        val = [1] * len(T) + [-1] * len(R)
        for b in range(len(T) - len(R) - 3, len(T) - len(R) + 1):
            ok_star = C._star_valid(us, vs, val, b, close)
            ok_enum = C._brute_rows(inst.indptr, inst.indices, np.array([0, len(us)]),
                                    np.array(us), np.array(vs), np.array(val),
                                    np.array([b]), np.array([0]), 9)[0] == 1
            assert ok_enum or not ok_star
            if complete:
                assert ok_enum == ok_star
            checked += 1
    assert checked > 100
