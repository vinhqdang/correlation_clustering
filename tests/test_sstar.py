"""The sparse star packing (experiments/sstar): valid, disjoint stars whose
value is the reported bound, and a certificate the checker accepts."""
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments")]


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    b = str(tmp_path_factory.mktemp("sstar") / "sstar")
    r = subprocess.run(["g++", "-O2", "-std=c++17", "-o", b,
                        os.path.join(ROOT, "experiments", "sstar", "sstar.cc")],
                       capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("no C++ compiler: " + r.stderr[-200:])
    return b


def random_graph(n, p, seed):
    rng = np.random.default_rng(seed)
    iu, ju = np.triu_indices(n, 1)
    keep = rng.random(len(iu)) < p
    return iu[keep], ju[keep]


def run(binary, tmp_path, n, eu, ev, T=1.0, seed=0, env=None):
    gr = tmp_path / "g.gr"
    with open(gr, "w") as fh:
        fh.write(f"p cep {n} {len(eu)}\n")
        for a, b in zip(eu, ev):
            fh.write(f"{a + 1} {b + 1}\n")
    st = tmp_path / "stars.txt"
    with open(gr) as fin:
        p = subprocess.run([binary, str(st), str(T), str(seed)], stdin=fin,
                           capture_output=True, text=True, timeout=120,
                           env=dict(os.environ, **(env or {})))
    assert p.returncode == 0, p.stderr
    value = int([l for l in p.stdout.splitlines() if l.startswith("star ")][0].split()[1])
    stars = [list(map(int, l.split())) for l in open(st) if l.strip()]
    return value, stars, gr


@pytest.mark.parametrize("n,p,seed", [(30, 0.2, 0), (60, 0.1, 1), (80, 0.3, 2)])
def test_stars_are_valid_and_disjoint(binary, tmp_path, n, p, seed):
    eu, ev = random_graph(n, p, seed)
    E = set(zip(eu.tolist(), ev.tolist()))
    edge = lambda a, b: (min(a, b), max(a, b)) in E
    value, stars, _ = run(binary, tmp_path, n, eu, ev, seed=seed)
    used = set()
    total = 0
    for st in stars:
        c, L = st[0], st[1:]
        assert len(L) >= 2 and len(set(L)) == len(L) and c not in L
        pairs = [(c, l) for l in L] + [(L[i], L[j]) for i in range(len(L))
                                       for j in range(i + 1, len(L))]
        for a, b in pairs[:len(L)]:
            assert edge(a, b)
        for a, b in pairs[len(L):]:
            assert not edge(a, b)
        for a, b in pairs:
            k = (min(a, b), max(a, b))
            assert k not in used
            used.add(k)
        total += len(L) - 1
    assert total == value


def test_leaf_cap(binary, tmp_path):
    eu, ev = random_graph(80, 0.3, 4)
    value, stars, _ = run(binary, tmp_path, 80, eu, ev, env={"SSTAR_KMAX": "3"})
    assert stars and max(len(st) - 1 for st in stars) <= 3
    assert value == sum(len(st) - 2 for st in stars)


def test_certificate_is_accepted(binary, tmp_path):
    import check_certificate as C
    from colab_run_lp import star_rows
    n = 70
    eu, ev = random_graph(n, 0.15, 3)
    value, stars, gr = run(binary, tmp_path, n, eu, ev)
    u, v, vals, ptr, b = star_rows(stars)
    inst = C.read_instance(str(gr), "pace")
    cert = tmp_path / "c.npz"
    np.savez_compressed(cert, n=n, ptr=ptr, u=u, v=v, val=vals, b=b, y=np.ones(len(b)),
                        bound=float(value), instance_m=np.array(inst.m),
                        edge_sha256=np.array(inst.edge_sha256),
                        raw_sha256=np.array(inst.raw_sha256))
    rep = C.check_file(str(gr), str(cert))
    assert rep["certified"] == value
