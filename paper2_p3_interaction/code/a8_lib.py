"""Audit 8: independent brute-force library (written from scratch; numpy only).
Vertices 0..n-1, graph = symmetric 0/1 numpy matrix, member = (a,b,c) with centre b (ab, bc edges, ac non-edge).
All partitions of n <= 10 vertices (Bell(10)=115975) are enumerated; cost(pi) = #pairs with [edge] != [same block]."""
import itertools, random
import numpy as np
from functools import lru_cache

def prs(n): return list(itertools.combinations(range(n), 2))

@lru_cache(maxsize=None)
def ptab(n):
    res = []
    def rec(lab, m):
        if len(lab) == n: res.append(tuple(lab)); return
        for c in range(m + 1):
            lab.append(c); rec(lab, max(m, c + 1)); lab.pop()
    rec([], 0)
    lab = np.array(res, dtype=np.int8)
    P = prs(n)
    same = np.stack([lab[:, i] == lab[:, j] for i, j in P], 1)
    return lab, same, P, {p: k for k, p in enumerate(P)}

def key(a, b): return (a, b) if a < b else (b, a)
def mpairs(t):
    a, b, c = t; return [key(a, b), key(b, c), key(a, c)]

class Inst:
    def __init__(self, A, H):
        self.A = np.array(A, dtype=np.int8); self.n = len(A); self.H = list(H)
        self.lab, self.same, self.P, self.pidx = ptab(self.n)
        self.e = np.array([self.A[i, j] for i, j in self.P], bool)
        self.diff = self.same != self.e          # edited pairs, per partition
        self.cost = self.diff.sum(1)
        self.soft = {}                            # pair -> member index
        for i, t in enumerate(self.H):
            for p in mpairs(t): self.soft[p] = i
        self.hardmask = np.array([p not in self.soft for p in self.P])
        self.softmask = ~self.hardmask
    def valid(self):
        A = self.A; seen = set()
        for (a, b, c) in self.H:
            if not (A[a, b] and A[b, c] and not A[a, c]): return False
            for p in mpairs((a, b, c)):
                if p in seen: return False
                seen.add(p)
        return True
    def is_pair_maximal(self):
        n = self.n; A = self.A
        for b in range(n):
            nb = [x for x in range(n) if A[b, x]]
            for a, c in itertools.combinations(nb, 2):
                if not A[a, c] and not any(p in self.soft for p in (key(a, b), key(b, c), key(a, c))): return False
        return True
    def Q(self, v):
        return {x for (a, b) in self.soft for x in (a, b) if v in (a, b) and x != v}
    def hard_nbrs(self, v):
        return [x for x in range(self.n) if x != v and self.A[v, x] and key(v, x) not in self.soft]
    def opt(self): return int(self.cost.min())
    def zero_excess(self): return np.nonzero(self.cost == len(self.H))[0]
    def state(self, k):
        """state vector of partition k (only meaningful when every member has exactly one edit)"""
        s = []
        for (a, b, c) in self.H:
            tab, tbc, tac = (self.same[k, self.pidx[key(a, b)]], self.same[k, self.pidx[key(b, c)]], self.same[k, self.pidx[key(a, c)]])
            s.append('L' if (not tab and tbc and not tac) else 'R' if (tab and not tbc and not tac) else 'M' if (tab and tbc and tac) else '?')
        return ''.join(s)
    def together(self, k, u, v): return bool(self.same[k, self.pidx[key(u, v)]])

def induced_p3s(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b, x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a, c]: out.append((a, b, c))
    return out

def greedy_pack(A, rng, maximal=True):
    L = induced_p3s(A); rng.shuffle(L); used = set(); H = []
    for t in L:
        ps = mpairs(t)
        if not any(p in used for p in ps):
            H.append(t); used.update(ps)
    return H

def gen_planted(n, rng, nedit=None, bigcluster=False):
    k = rng.randint(1, max(1, n // 2))
    lab = [rng.randrange(k) for _ in range(n)]
    if bigcluster:
        for v in range(n - 2): lab[v] = 0
    A = np.zeros((n, n), np.int8)
    for i, j in prs(n):
        if lab[i] == lab[j]: A[i, j] = A[j, i] = 1
    for _ in range(nedit if nedit is not None else rng.randint(1, 4)):
        i, j = rng.sample(range(n), 2); A[i, j] ^= 1; A[j, i] = A[i, j]
    return A

def gen_random(n, rng):
    p = rng.choice([0.25, 0.35, 0.5]); A = np.zeros((n, n), np.int8)
    for i, j in prs(n):
        if rng.random() < p: A[i, j] = A[j, i] = 1
    return A

def make_inst(n, rng, kind=None):
    kind = kind or rng.choice(['planted', 'planted', 'random'])
    A = gen_planted(n, rng) if kind == 'planted' else gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: return None
    return Inst(A, H)

def comps_I(H):
    m = len(H); adj = {i: set() for i in range(m)}
    for i in range(m):
        for j in range(i + 1, m):
            if set(H[i]) & set(H[j]): adj[i].add(j); adj[j].add(i)
    return adj
