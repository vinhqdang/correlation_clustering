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
        soft0 = set()
        for t in self.mem0: soft0 |= pairs_of(t)
        for (u, v) in self.E0:
            if (u, v) in soft0: continue
            for i in range(p):
                for j in range(p):
                    self.B.hard((self.names0[u], i), (self.names0[v], j))
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


def host_lists(m, host_thin):
    """from the unique ZE solution of each position a: partner sets of the host vertex."""
    cs = cluster_sets(m)
    per = {}
    for a, clusters in cs.items():
        for c in clusters:
            if host_thin in c: per[a] = frozenset(c - {host_thin})
    return per

def add_permanent_guest(T, host, zname, per):
    """z forced into proto-cluster of thin vertex `host`; soft to all clusters the host can be in."""
    p = T.p
    M = T.take_matching(host)
    for (i, k) in M:
        T.B.member((host, i), zname, (host, k))
    nbrs = set()
    for a, c in per.items(): nbrs |= set(c)
    for nv in sorted(nbrs, key=str):
        Mn = T.take_matching(nv)
        for (i, k) in Mn:
            T.B.member(zname, (nv, i), (nv, k))
    return nbrs

def add_lock_guest(T, yname, xname, zname, host, per, forbidden_pos):
    p = T.p
    allowed = set([host])
    for a, c in per.items():
        if a != forbidden_pos: allowed |= set(c)
    for nv in sorted(allowed, key=str):
        Mn = T.take_matching(nv)
        for (i, k) in Mn:
            T.B.member(yname, (nv, i), (nv, k))
    T.B.member(xname, yname, zname)
    return allowed

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test1':
    m = int(sys.argv[2]); j = int(sys.argv[3]); p = int(sys.argv[4]); side = sys.argv[5]   # side 'L' host s_{2j}, 'R' host s_{2j+2}
    T = Thick(m, p)
    host = ('s', 2 * j) if side == 'L' else ('s', 2 * j + 2)
    per = host_lists(m, host)
    print('host', host, 'allowed partners (union over a != j):', sorted(set().union(*[set(c) for a, c in per.items() if a != j]), key=str))
    print('cluster of host in config j (forbidden):', sorted(per[j], key=str))
    nb = add_permanent_guest(T, host, 'z', per)
    al = add_lock_guest(T, 'y', 'x', 'z', host, per, j)
    n, E, members = T.finish()
    print('n', n, '|H|', len(members), check_instance(n, E, members), 'c =', max(sum(1 for t in members if v in t) for v in range(n)), flush=True)
    idx = T.B.idx
    import time; t0 = time.time()
    sols = ze_solve(n, E, members, all_solutions=True, limit=200)
    def pos(lab):
        for q in range(m):
            a, b, c = T.idx0[('s', 2 * q)], T.idx0[('s', 2 * q + 1)], T.idx0[('s', 2 * q + 2)]
            if lab[idx[(('s', 2 * q), 0)]] == lab[idx[(('s', 2 * q + 1), 0)]] == lab[idx[(('s', 2 * q + 2), 0)]]: return q
    from collections import Counter
    print('unlocked #ZE', len(sols), 'positions', sorted(Counter(pos(l) for l in sols).items()), round(time.time() - t0, 1), flush=True)
    sols2 = ze_solve(n, E, members, singletons=[idx['x']], all_solutions=True, limit=200)
    print('locked  #ZE', len(sols2), 'positions', sorted(Counter(pos(l) for l in sols2).items()), flush=True)


def candidates(m):
    cs = cluster_sets(m)
    names = [('s', k) for k in range(2 * m + 1)] + [('e', j) for j in range(m + 1)]
    def C(a, z):
        for c in cs[a]:
            if z in c: return c
    out = {}
    for a in range(m):
        out[a] = {}
        for z in names:
            U = set()
            for b in range(m):
                if b != a: U |= set(C(b, z)) - {z}
            F = set(C(a, z)) - {z} - U
            if F: out[a][z] = (U, F)
    return cs, out

def build_pfg(m, p, host_rule=None, positions=None):
    T = Thick(m, p)
    cs, cand = candidates(m)
    host_of = {}
    positions = list(range(m)) if positions is None else positions
    for a in positions:
        pref = [('s', 2 * a + 2), ('s', 2 * a)] if a % 2 == 0 else [('s', 2 * a), ('s', 2 * a + 2)]
        if host_rule: pref = host_rule(a) + pref
        for z in pref:
            if z in cand[a]:
                host_of[a] = z; break
        else:
            raise RuntimeError(('no host', a, list(cand[a])))
    # permanent guests
    perm = {}
    for h in sorted(set(host_of.values()), key=str):
        per = {a: frozenset(c - {h}) for a, cl_ in cs.items() for c in cl_ if h in c}
        zname = ('zg', h)
        add_permanent_guest(T, h, zname, per)
        perm[h] = zname
    anchors = {}
    for a in positions:
        h = host_of[a]
        per = {b: frozenset(c - {h}) for b, cl_ in cs.items() for c in cl_ if h in c}
        yname, xname = ('y', a), ('x', a)
        add_lock_guest(T, yname, xname, perm[h], h, per, a)
        anchors[a] = xname
    return T, host_of, anchors

def position_of(T, lab, m):
    idx = T.B.idx
    for q in range(m):
        if lab[idx[(('s', 2 * q), 0)]] == lab[idx[(('s', 2 * q + 1), 0)]] == lab[idx[(('s', 2 * q + 2), 0)]]: return q

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'pfg':
    m, p = int(sys.argv[2]), int(sys.argv[3])
    import time
    t0 = time.time()
    T, host_of, anchors = build_pfg(m, p)
    n, E, members = T.finish()
    print('hosts', host_of)
    print('pad matchings used per thin vertex (max):', max(T.used.values()), 'of', len(T.fact))
    print('n', n, '|H|', len(members), check_instance(n, E, members), 'c =', max(sum(1 for t in members if v in t) for v in range(n)), round(time.time() - t0, 1), flush=True)
    idx = T.B.idx
    A = [idx[anchors[a]] for a in range(m)]
    res = ze_solve(n, E, members, singletons=A)
    print('ALL locked feasible?', res is not None, flush=True)
    for i in range(m):
        S = [A[a] for a in range(m) if a != i]
        lab = ze_solve(n, E, members, singletons=S)
        print(' all but', i, 'feasible', lab is not None, 'M position', position_of(T, lab, m) if lab else None, flush=True)
