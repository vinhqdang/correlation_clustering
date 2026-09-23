"""External and simple baselines."""
from __future__ import annotations

import os
import subprocess
import tempfile
import time

import numpy as np

import ccbench as cc
from ccbench.graph import Graph, from_edges, components

KAPOCE = os.environ.get("KAPOCE_BIN", "/home/user/baselines/kapoce/build/ClusterEditing")


def write_pace(g: Graph, path: str) -> None:
    e = g.edges() + 1
    with open(path, "w") as fh:
        fh.write(f"p cep {g.n} {g.m}\n")
        np.savetxt(fh, e, fmt="%d")


def apply_edits(g: Graph, edits: np.ndarray) -> np.ndarray:
    """Clustering obtained by toggling the given pairs (0-indexed) in G+."""
    E = g.edges().astype(np.int64)
    key = set((E[:, 0] * g.n + E[:, 1]).tolist())
    for u, v in edits:
        a, b = min(u, v), max(u, v)
        k = a * g.n + b
        if k in key:
            key.remove(k)
        else:
            key.add(k)
    keys = np.fromiter(key, dtype=np.int64, count=len(key))
    h = from_edges(g.n, np.stack([keys // g.n, keys % g.n], axis=1))
    return components(h)[0]


def kapoce(g: Graph, time_limit: float = 60.0):
    """KaPoCE heuristic (PACE 2021 winner).  Stopped with SIGTERM at the limit."""
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "in.gr")
        write_pace(g, inp)
        t0 = time.time()
        with open(inp) as fin:
            p = subprocess.Popen([KAPOCE], stdin=fin, stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL, text=True)
            try:
                out, _ = p.communicate(timeout=time_limit)
            except subprocess.TimeoutExpired:
                p.terminate()
                out, _ = p.communicate()
        elapsed = time.time() - t0
    rows = []
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            rows.append((int(parts[0]) - 1, int(parts[1]) - 1))
    lab = apply_edits(g, rows)
    return lab, elapsed


def leiden_cpm(g: Graph, rng=None, iterations: int = -1):
    """Leiden with the constant Potts model at resolution 1/2.

    Maximising sum_c (e_c - n_c(n_c-1)/4) ... equals minimising the CC objective
    (Veldt, Gleich and Wirth: LambdaCC with lambda = 1/2)."""
    import igraph as ig
    import leidenalg as la
    seed = int(np.random.default_rng(rng).integers(1 << 30))
    G = ig.Graph(n=g.n, edges=g.edges().tolist())
    part = la.find_partition(G, la.CPMVertexPartition, resolution_parameter=0.5,
                             n_iterations=iterations, seed=seed)
    return np.asarray(part.membership, dtype=np.int64)


def match_flip_pivot(g: Graph, rng=None):
    """MatchFlipPivot (Veldt, ICML 2022): maximal pair-disjoint set of open
    wedges, flip all their pairs, run Pivot on the flipped graph."""
    from ccbench.support import build_support
    rng = np.random.default_rng(rng)
    sup = build_support(g)
    order = rng.permutation(g.n)
    used = np.zeros(sup.npairs, dtype=bool)
    flip = []
    ip, ix = g.indptr, g.indices
    adj = [set(ix[ip[v]:ip[v + 1]].tolist()) for v in range(g.n)] if g.n < 200000 else None
    from ccbench.support import pair_id
    for w in order:
        nb = ix[ip[w]:ip[w + 1]]
        for i in range(len(nb)):
            e1 = sup.eid[ip[w] + i]
            if used[e1]:
                continue
            for j in range(i + 1, len(nb)):
                e2 = sup.eid[ip[w] + j]
                if used[e2] or used[e1]:
                    continue
                u, v = int(nb[i]), int(nb[j])
                if v in adj[u]:
                    continue
                e3 = pair_id(ip, ix, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id, u, v)
                if used[e3]:
                    continue
                used[e1] = used[e2] = used[e3] = True
                flip += [e1, e2, e3]
    flip = np.array(flip, dtype=np.int64)
    fl = np.zeros(sup.npairs, dtype=bool)
    fl[flip] = True
    pos = np.arange(sup.npairs) < sup.npos
    keep = pos ^ fl
    h = from_edges(g.n, np.stack([sup.pu[keep], sup.pv[keep]], axis=1))
    return cc.pivot(h, rng)
