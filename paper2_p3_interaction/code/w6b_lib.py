"""Wave 6 (prover B) library. Written from scratch; imports nothing from other code/ scripts.
Conventions: vertices 0..n-1, adj = n x n 0/1 numpy array, H = list of triples (x,y,z) with centre y
(induced P3: xy, yz edges, xz non-edge).  soft pair = pair inside a member of H, hard = any other pair.
A partition is a label vector; x(p,q) = 1 iff same label ('final adjacency').  Edited pair = x != G.
"""
import itertools, random
import numpy as np

_PART_CACHE = {}

def all_partitions(n):
    """matrix (Bell(n), n) of restricted growth strings."""
    if n in _PART_CACHE:
        return _PART_CACHE[n]
    out = []
    def rec(i, lab, m):
        if i == n:
            out.append(tuple(lab)); return
        for c in range(m + 1):
            lab.append(c); rec(i + 1, lab, max(m, c + 1)); lab.pop()
    rec(0, [], 0)
    arr = np.array(out, dtype=np.int8)
    _PART_CACHE[n] = arr
    return arr

def pair_list(n):
    return list(itertools.combinations(range(n), 2))

def costs_all(adj):
    """cost of every partition (vector)."""
    n = len(adj)
    P = all_partitions(n)
    pl = pair_list(n)
    ia = np.array([p[0] for p in pl]); ib = np.array([p[1] for p in pl])
    same = (P[:, ia] == P[:, ib])
    g = np.array([adj[a][b] for a, b in pl], dtype=bool)
    return (same != g[None, :]).sum(axis=1)

def mk(n, edges):
    a = np.zeros((n, n), dtype=np.int8)
    for x, y in edges:
        a[x, y] = a[y, x] = 1
    return a

def pairs3(t):
    return {frozenset(p) for p in itertools.combinations(t, 2)}

def induced_p3s(adj):
    n = len(adj); res = []
    for b in range(n):
        nb = [x for x in range(n) if adj[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not adj[a][c]:
                res.append((a, b, c))
    return res

def packed_set(H):
    s = set()
    for t in H:
        s |= pairs3(t)
    return s

def is_packing(adj, H):
    seen = set()
    for (a, b, c) in H:
        if not (adj[a][b] and adj[b][c] and not adj[a][c]):
            return False
        p = pairs3((a, b, c))
        if p & seen:
            return False
        seen |= p
    return True

def pair_maximal(adj, H):
    pk = packed_set(H)
    return all(pairs3(t) & pk for t in induced_p3s(adj))

def random_pm_packing(adj, rng, maxc=None, drop=0.0):
    """random modification-disjoint packing, extended to pair-maximal.  maxc: optional cap on members per vertex
    (the extension step ignores the cap only if needed for maximality)."""
    ts = induced_p3s(adj); rng.shuffle(ts)
    H = []; used = set(); cnt = {}
    def ok(t, cap):
        if pairs3(t) & used: return False
        if cap is not None and any(cnt.get(v, 0) >= cap for v in t): return False
        return True
    for t in ts:
        if rng.random() < drop: continue
        if ok(t, maxc):
            H.append(t); used |= pairs3(t)
            for v in t: cnt[v] = cnt.get(v, 0) + 1
    for t in ts:
        if ok(t, None):
            H.append(t); used |= pairs3(t)
            for v in t: cnt[v] = cnt.get(v, 0) + 1
    return H

class Inst:
    """instance (G,H) with derived structure."""
    def __init__(self, adj, H):
        self.adj = np.array(adj, dtype=np.int8); self.n = len(adj); self.H = [tuple(t) for t in H]
        n = self.n
        self.soft = packed_set(self.H)
        self.Q = [set() for _ in range(n)]
        for p in self.soft:
            a, b = tuple(p); self.Q[a].add(b); self.Q[b].add(a)
        self.Nh = [set(y for y in range(n) if self.adj[v][y] and frozenset((v, y)) not in self.soft) for v in range(n)]
        self.cnt = [sum(1 for t in self.H if v in t) for v in range(n)]
        self.c = max(self.cnt) if n else 0
        self.D = 2 * self.c
        self.member_of = [[i for i, t in enumerate(self.H) if v in t] for v in range(n)]
        # intersection graph I(H)
        m = len(self.H)
        self.I = [set() for _ in range(m)]
        for i in range(m):
            for j in range(i + 1, m):
                if set(self.H[i]) & set(self.H[j]):
                    self.I[i].add(j); self.I[j].add(i)
        self.C = None; self.P = None; self.cost = None
    def is_soft(self, a, b): return frozenset((a, b)) in self.soft
    def hard_pairs(self):
        return [(a, b) for a, b in pair_list(self.n) if not self.is_soft(a, b)]
    def solve(self):
        if self.cost is None:
            self.P = all_partitions(self.n)
            self.cost = costs_all(self.adj)
        return int(self.cost.min())
    def opt_partitions(self, extra=0):
        self.solve()
        m = self.cost.min() + extra
        return self.P[self.cost <= m] if extra else self.P[self.cost == m]
    def edits(self, lab):
        """(hard edited pairs U, all edited pairs S) for a label vector."""
        U = []; S = []
        for a, b in pair_list(self.n):
            e = (lab[a] == lab[b]) != bool(self.adj[a][b])
            if e:
                S.append((a, b))
                if not self.is_soft(a, b): U.append((a, b))
        return U, S
    # I(H) distances between members
    def member_dist(self):
        m = len(self.H); INF = 10 ** 6
        d = [[INF] * m for _ in range(m)]
        for s in range(m):
            d[s][s] = 0; fr = [s]
            while fr:
                nx = []
                for u in fr:
                    for v in self.I[u]:
                        if d[s][v] == INF:
                            d[s][v] = d[s][u] + 1; nx.append(v)
                fr = nx
        return d
    def vertex_dist(self, a, b):
        """I(H)-distance between the closest members of a and b; INF if a or b in no member / different components;
        0 if they share a member."""
        INF = 10 ** 6
        if not self.member_of[a] or not self.member_of[b]: return INF
        if not hasattr(self, '_md'): self._md = self.member_dist()
        return min(self._md[i][j] for i in self.member_of[a] for j in self.member_of[b])

def delta(inst, z, b):
    """Delta(z,b) = {y != z,b : zy, by hard, exactly one is an edge of G}."""
    out = set()
    for y in range(inst.n):
        if y in (z, b): continue
        if inst.is_soft(z, y) or inst.is_soft(b, y): continue
        if bool(inst.adj[z][y]) != bool(inst.adj[b][y]): out.add(y)
    return out

def common_hard(inst, z, b):
    return inst.Nh[z] & inst.Nh[b]

def rand_graph(n, rng, kind=None):
    kind = kind or rng.choice(['gnp', 'noisy', 'sparse', 'hub', 'bridges'])
    A = np.zeros((n, n), dtype=np.int8)
    def put(a, b):
        if a != b: A[a, b] = A[b, a] = 1
    if kind == 'gnp':
        p = rng.uniform(0.25, 0.6)
        for a, b in pair_list(n):
            if rng.random() < p: put(a, b)
    elif kind == 'sparse':
        p = rng.uniform(0.12, 0.3)
        for a, b in pair_list(n):
            if rng.random() < p: put(a, b)
    elif kind == 'noisy':
        k = rng.randint(2, max(2, n // 2)); lab = [rng.randrange(k) for _ in range(n)]
        for a, b in pair_list(n):
            if (lab[a] == lab[b]) != (rng.random() < 0.12): put(a, b)
    elif kind == 'hub':
        h = rng.randrange(n)
        for v in range(n):
            if v != h and rng.random() < 0.7: put(h, v)
        for a, b in pair_list(n):
            if rng.random() < 0.15: put(a, b)
    else:  # bridges: several small cliques joined by few edges
        k = rng.randint(2, 4); lab = [rng.randrange(k) for _ in range(n)]
        for a, b in pair_list(n):
            if lab[a] == lab[b]: put(a, b)
        for _ in range(rng.randint(1, 3)):
            a, b = rng.sample(range(n), 2); put(a, b)
    return A
