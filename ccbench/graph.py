"""Graph representation for complete signed instances.

An instance of (unweighted, complete-graph) correlation clustering is given by
the positive graph G+ = (V, E+); every pair not in E+ is negative.  We store
G+ in CSR form with sorted neighbour lists and no self-loops.
"""
from __future__ import annotations

import gzip
import os
from dataclasses import dataclass

import numpy as np
import numba as nb


@dataclass
class Graph:
    n: int
    indptr: np.ndarray   # int64, length n+1
    indices: np.ndarray  # int32, length 2m, sorted within each row

    @property
    def m(self) -> int:
        return int(self.indices.shape[0] // 2)

    @property
    def degrees(self) -> np.ndarray:
        return np.diff(self.indptr)

    def neighbors(self, v: int) -> np.ndarray:
        return self.indices[self.indptr[v]:self.indptr[v + 1]]

    def edges(self) -> np.ndarray:
        """Return E+ as an (m, 2) array with u < v."""
        src = np.repeat(np.arange(self.n, dtype=np.int32), self.degrees)
        mask = src < self.indices
        return np.stack([src[mask], self.indices[mask]], axis=1)

    def __repr__(self) -> str:
        return f"Graph(n={self.n}, m={self.m})"


def from_edges(n: int, edges: np.ndarray) -> Graph:
    """Build a Graph from an (k, 2) integer array of undirected edges.

    Self-loops and duplicate edges are removed.
    """
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    if edges.size:
        if edges.min() < 0 or edges.max() >= n:
            raise ValueError("edge endpoint out of range")
    u = np.minimum(edges[:, 0], edges[:, 1])
    v = np.maximum(edges[:, 0], edges[:, 1])
    keep = u != v
    u, v = u[keep], v[keep]
    key = np.unique(u * n + v)
    u, v = key // n, key % n
    src = np.concatenate([u, v])
    dst = np.concatenate([v, u])
    order = np.lexsort((dst, src))
    src, dst = src[order], dst[order]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(indptr, src + 1, 1)
    indptr = np.cumsum(indptr)
    return Graph(n, indptr, dst.astype(np.int32))


def relabel_edges(edges: np.ndarray):
    """Map arbitrary vertex ids to 0..n-1.  Returns (n, edges, original_ids)."""
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    ids, inv = np.unique(edges.ravel(), return_inverse=True)
    return len(ids), inv.reshape(-1, 2), ids


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_edgelist(path: str, comment: str = "#", sign_col: int | None = None,
                  keep_sign: int = 1) -> Graph:
    """Read a whitespace/comma separated edge list (SNAP style).

    If ``sign_col`` is given, only rows whose value in that column has the
    sign of ``keep_sign`` are kept (used for signed networks, where the
    positive edges define G+).
    """
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(comment) or line.startswith("%"):
                continue
            parts = line.replace(",", " ").split()
            if sign_col is not None:
                s = float(parts[sign_col])
                if s * keep_sign <= 0:
                    continue
            rows.append((int(parts[0]), int(parts[1])))
    n, e, _ = relabel_edges(np.array(rows, dtype=np.int64))
    return from_edges(n, e)


def read_pace(path: str) -> Graph:
    """Read a PACE 2021 cluster-editing instance (``p cep n m`` header, 1-indexed)."""
    n = None
    rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("c") or not line.strip():
                continue
            parts = line.split()
            if parts[0] == "p":
                n = int(parts[2])
                continue
            rows.append((int(parts[0]) - 1, int(parts[1]) - 1))
    if n is None:
        raise ValueError(f"{path}: missing 'p cep' header")
    return from_edges(n, np.array(rows, dtype=np.int64).reshape(-1, 2))


def write_edgelist(g: Graph, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savetxt(path, g.edges(), fmt="%d")


@nb.njit(cache=True)
def _components(n, indptr, indices):
    comp = np.full(n, -1, dtype=np.int64)
    stack = np.empty(n, dtype=np.int64)
    c = 0
    for s in range(n):
        if comp[s] >= 0:
            continue
        top = 0
        stack[top] = s
        top += 1
        comp[s] = c
        while top > 0:
            top -= 1
            u = stack[top]
            for p in range(indptr[u], indptr[u + 1]):
                w = indices[p]
                if comp[w] < 0:
                    comp[w] = c
                    stack[top] = w
                    top += 1
        c += 1
    return comp, c


def components(g: Graph):
    """Connected components of G+.  Returns (labels, count)."""
    return _components(g.n, g.indptr, g.indices)


def induced_subgraph(g: Graph, nodes: np.ndarray) -> Graph:
    nodes = np.asarray(nodes, dtype=np.int64)
    pos = np.full(g.n, -1, dtype=np.int64)
    pos[nodes] = np.arange(len(nodes))
    e = g.edges().astype(np.int64)
    keep = (pos[e[:, 0]] >= 0) & (pos[e[:, 1]] >= 0)
    e = pos[e[keep]]
    return from_edges(len(nodes), e)
