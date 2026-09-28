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
    inst = C.read_instance(gr, "pace")
    cert = dict(np.load(path))
    cert.update({"n": np.array(inst.n), "instance_m": np.array(inst.m),
                 "edge_sha256": np.array(inst.edge_sha256),
                 "raw_sha256": np.array(inst.raw_sha256)})
    return inst, res, cert, g


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
    cert["edge_sha256"] = np.array(inst.edge_sha256)
    cert["n"] = np.array(inst.n + 1)
    with pytest.raises(ValueError):
        C.check(inst, cert)
    # a certificate without the instance hash is rejected
    cert["n"] = np.array(inst.n)
    cert.pop("edge_sha256")
    with pytest.raises(ValueError):
        C.check(inst, cert)


def test_empty_row_and_large_star(tmp_path):
    """An empty row is decided without looping; a star row with many leaves
    does not exhaust the recursion limit."""
    inst, _, cert, _ = _cert(tmp_path, seed=1)
    for key in ("ptr", "b", "y"):
        cert[key] = cert[key].copy()
    cert["ptr"] = np.concatenate([cert["ptr"], cert["ptr"][-1:]])
    cert["b"] = np.concatenate([cert["b"], [0.0]])
    cert["y"] = np.concatenate([cert["y"], [1.0]])
    C.check(inst, cert)
    cert["b"][-1] = 1.0
    with pytest.raises(ValueError):
        C.check(inst, cert)
    adj = {t: set() for t in range(1500)}
    adj[0].add(1)
    adj[1].add(0)
    assert C._max_independent(adj) == 1499


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


def test_certificate_rejects_huge_multiplier(tmp_path):
    """A multiplier too large for int64 must be rejected, not wrap around to
    a negative value (a zero row with b = -1 would then add to the bound)."""
    g, _, cert, _ = _cert(tmp_path, seed=1)
    for key in ("y", "b"):
        cert[key] = cert[key].copy()
    cert["y"][0] = 1e10
    with pytest.raises(ValueError):
        C.check(g, cert)
    cert["y"][0] = 0.0
    cert["b"][0] = -2.0 ** 60
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


def test_directory_check_uses_the_instance_table(tmp_path):
    """--dir takes the parser settings from INSTANCES.tsv: a certificate that
    declares another format or sign column is rejected, and an instance file
    missing from the manifest or the table is rejected."""
    import hashlib
    inst, _, cert, _ = _cert(tmp_path)
    raw = tmp_path / "raw"
    (raw / "pace").mkdir(parents=True)
    os.replace(tmp_path / "g.gr", raw / "pace" / "g.gr")
    sha = hashlib.sha256((raw / "pace" / "g.gr").read_bytes()).hexdigest()
    (raw / "MANIFEST.sha256").write_text(f"{sha}  pace/g.gr\n")
    (raw / "INSTANCES.tsv").write_text("# test\npace/*\tpace\t-1\n")
    certs = tmp_path / "certs"
    certs.mkdir()
    cert["instance_file"] = np.array("pace/g.gr")
    np.savez_compressed(certs / "a_ok.npz", **cert)
    np.savez_compressed(certs / "b_sign.npz", **{**cert, "sign_col": np.array(2)})
    np.savez_compressed(certs / "c_fmt.npz", **{**cert, "instance_format": np.array("snap")})
    np.savez_compressed(certs / "d_file.npz", **{**cert, "instance_file": np.array("pace/h.gr")})
    out = str(tmp_path / "out.csv")
    assert not C.check_dir(str(raw), str(certs), out)
    import csv
    status = {r["certificate"]: r["status"] for r in csv.DictReader(open(out))}
    assert status["a_ok.npz"] == "ok"
    assert all(status[k].startswith("rejected") for k in ("b_sign.npz", "c_fmt.npz", "d_file.npz"))


def _tiny(tmp_path, n, edges):
    gr = str(tmp_path / "t.gr")
    with open(gr, "w") as fh:
        fh.write(f"p cep {n} {len(edges)}\n")
        for a, b in edges:
            fh.write(f"{a + 1} {b + 1}\n")
    inst = C.read_instance(gr, "pace")
    ident = {"n": np.array(inst.n), "instance_m": np.array(inst.m),
             "edge_sha256": np.array(inst.edge_sha256), "raw_sha256": np.array(inst.raw_sha256)}
    return inst, ident


def test_certificate_rejects_int64_wrap_in_row(tmp_path):
    """An invalid row whose left-hand side wraps around in int64 (4097 copies
    of -2^51 on one edge, b = 0) must be rejected, not accepted as valid."""
    inst, ident = _tiny(tmp_path, 3, [(0, 1), (1, 2)])
    k = 4097
    cert = dict(ident, ptr=np.array([0, k]), u=np.zeros(k, dtype=np.int64),
                v=np.ones(k, dtype=np.int64), val=np.full(k, -2.0 ** 51), b=np.array([0.0]),
                y=np.array([0.0]))
    with pytest.raises(ValueError, match="same pair twice"):
        C.check(inst, cert)


def test_certificate_rejects_huge_row_sum(tmp_path):
    """Distinct pairs whose coefficients sum to 2^61 or more are rejected."""
    n = 51
    inst, ident = _tiny(tmp_path, n, [(0, t) for t in range(1, n)])
    pairs = [(0, t) for t in range(1, n)] + [(a, c) for a in range(1, n) for c in range(a + 1, n)]
    u = np.array([a for a, _ in pairs], dtype=np.int64)
    v = np.array([c for _, c in pairs], dtype=np.int64)
    cert = dict(ident, ptr=np.array([0, len(pairs)]), u=u, v=v,
                val=np.full(len(pairs), -2.0 ** 51), b=np.array([0.0]), y=np.array([0.0]))
    with pytest.raises(ValueError, match="too large for exact int64 enumeration"):
        C.check(inst, cert)


@pytest.mark.parametrize("change", [
    ("val", lambda c: np.where(np.arange(len(c["val"])) == 0, np.nan, c["val"])),
    ("val", lambda c: np.where(np.arange(len(c["val"])) == 0, np.inf, c["val"])),
    ("b", lambda c: np.where(np.arange(len(c["b"])) == 0, -np.inf, c["b"])),
    ("b", lambda c: np.where(np.arange(len(c["b"])) == 0, np.nan, c["b"])),
    ("y", lambda c: np.where(np.arange(len(c["y"])) == 0, np.nan, c["y"])),
    ("v", lambda c: c["v"][:-1]),
    ("val", lambda c: np.concatenate([c["val"], [1.0]])),
    ("val", lambda c: c["val"][:1]),
    ("b", lambda c: c["b"][:-1]),
    ("u", lambda c: np.stack([c["u"], c["u"]])),
    ("u", lambda c: c["u"].astype(np.uint64) + np.uint64(2 ** 63)),
    ("v", lambda c: c["v"].astype(np.float64)),
    ("val", lambda c: c["val"].astype(np.complex128)),
    ("val", lambda c: np.where(np.arange(len(c["val"])) == 0, 2.0 ** 60, c["val"])),
    ("ptr", lambda c: np.concatenate([c["ptr"][:1], c["ptr"][2:]])),
    ("ptr", lambda c: c["ptr"] + 1),
])
def test_certificate_rejects_malformed_arrays(tmp_path, change):
    """NaN, infinity, arrays of different lengths or shapes, wrapped or huge
    integers and non-real types are rejected, never evaluated."""
    g, _, cert, _ = _cert(tmp_path, seed=1)
    key, f = change
    cert[key] = f(cert)
    with pytest.raises(ValueError):
        C.check(g, cert)


def test_certificate_without_array_is_rejected(tmp_path):
    g, _, cert, _ = _cert(tmp_path, seed=1)
    cert.pop("val")
    with pytest.raises(ValueError):
        C.check(g, cert)
