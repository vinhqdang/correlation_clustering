"""Thick one-hot chain with guests (locks).  See notes/wave8_hardness.md Sec. 5.
thin chain: w8c_oh.oh_instance(m).  Thick: each thin vertex v -> proto-cluster {(v,i): i<p}; each thin member -> p^2 Latin members.
Pads: a perfect matching of the proto-cluster, taken out of the hard glue clique, owned by pad members.
Guests: permanent guest z (forced into its host proto-cluster) and closable lock guest y with anchor x."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_thick import Builder
from w8c_oh import oh_instance

def one_factorization(p):
    """round-robin 1-factorization of K_p (p even): list of p-1 perfect matchings"""
    assert p % 2 == 0
    pts = list(range(p - 1)) + ['inf']
    res = []
    for r in range(p - 1):
        M = []
        for i in range(p // 2):
            a = (r + i) % (p - 1)
            b = (r - i) % (p - 1)
            if i == 0: M.append((a, p - 1))
            else: M.append((min(a, b), max(a, b)))
        res.append(sorted(set(tuple(sorted(x)) for x in M)))
    # sanity
    allp = set()
    for M in res:
        vs = sorted(v for e in M for v in e); assert vs == list(range(p)), (p, M)
        for e in M: assert e not in allp; allp.add(e)
    assert len(allp) == p * (p - 1) // 2
    return res

def cluster_sets(m):
    """thin OH: for each position a, the clusters (as frozensets of thin names) of the unique ZE solution."""
    n, E, members, names, idx = oh_instance(m)
    sols = ze_solve(n, E, members, all_solutions=True, limit=1000)
    out = {}
    for lab in sols:
        a = None
        for j in range(m):
            x, y, z = idx[('s', 2 * j)], idx[('s', 2 * j + 1)], idx[('s', 2 * j + 2)]
            if lab[x] == lab[y] == lab[z]: a = j
        d = {}
        for v, l in enumerate(lab): d.setdefault(l, set()).add(names[v])
        out[a] = [frozenset(c) for c in d.values()]
    return out

class Thick:
    def __init__(self, m, p):
        self.m, self.p = m, p
        self.n0, self.E0, self.mem0, self.names0, self.idx0 = oh_instance(m)
        self.B = Builder()
        self.fact = one_factorization(p)
        self.used = {}           # thin vertex -> number of matchings used
        self.pads = {}           # thin vertex -> set of frozenset pairs used as pads
        for (a, b, c) in self.mem0:
            A, Bn, C = self.names0[a], self.names0[b], self.names0[c]
            for i in range(p):
                for j in range(p):
                    self.B.member((A, i), (Bn, j), (C, (i + j) % p))
        self.guest_members = []
    def take_matching(self, thin_v):
        k = self.used.get(thin_v, 0)
        assert k < len(self.fact), ('out of matchings at', thin_v)
        self.used[thin_v] = k + 1
        M = self.fact[k]
        for (i, j) in M: self.pads.setdefault(thin_v, set()).add((i, j))
        return M
    def finish(self):
        p = self.p
        for v in self.names0:
            padset = self.pads.get(v, set())
            # glue connectivity
            import networkx as nx
            g = nx.Graph(); g.add_nodes_from(range(p))
            for i in range(p):
                for k in range(i + 1, p):
                    if (i, k) not in padset:
                        g.add_edge(i, k); self.B.hard((v, i), (v, k))
            assert nx.is_connected(g), ('glue disconnected at', v, len(padset))
        return self.B.finish()
