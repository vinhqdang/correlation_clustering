"""Audit-11 library (written from scratch; numpy + itertools only; no import of w*/a10* files).
Graph = (n, set of sorted-tuple edges).  Member (a,b,c): centre b, ab,bc in E, ac not in E.
Ground truth: ALL set partitions of V (restricted growth strings), n <= 10 (Bell(10)=115975).
"""
import itertools, random
import numpy as np

def K(a, b): return (a, b) if a < b else (b, a)

class Gr:
    def __init__(s, n, edges):
        s.n = n
        s.E = set(K(a, b) for a, b in edges)
        s.adj = [set() for _ in range(n)]
        for a, b in s.E: s.adj[a].add(b); s.adj[b].add(a)
        s.pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
        s.pidx = {p: i for i, p in enumerate(s.pairs)}
    def e(s, a, b): return K(a, b) in s.E
    def p3(s):
        out = []
        for b in range(s.n):
            nb = sorted(s.adj[b])
            for i in range(len(nb)):
                for j in range(i + 1, len(nb)):
                    if not s.e(nb[i], nb[j]): out.append((nb[i], b, nb[j]))
        return out

def pr(t):
    a, b, c = t
    return [K(a, b), K(b, c), K(a, c)]

def is_packing(g, H):
    seen = set()
    for t in H:
        a, b, c = t
        if len({a, b, c}) < 3 or not (g.e(a, b) and g.e(b, c) and not g.e(a, c)): return False
        for p in pr(t):
            if p in seen: return False
            seen.add(p)
    return True

def soft_owner(H):
    o = {}
    for i, t in enumerate(H):
        for p in pr(t): o[p] = i
    return o

def pair_maximal(g, H):
    o = soft_owner(H)
    return all(any(p in o for p in pr(t)) for t in g.p3())

def maximalise(g, H, rng=None):
    H = list(H); o = soft_owner(H)
    P = g.p3()
    if rng: rng.shuffle(P)
    for t in P:
        if not any(p in o for p in pr(t)):
            H.append(t)
            for p in pr(t): o[p] = len(H) - 1
    return H

def random_packing(g, rng):
    """random greedy pair-maximal packing"""
    return maximalise(g, [], rng)

def I_graph(H):
    return {i: {j for j in range(len(H)) if j != i and set(H[i]) & set(H[j])} for i in range(len(H))}

def tw_exact(adj):
    """exact treewidth by DP over subsets (<= ~16 nodes), adj dict node->set"""
    nodes = list(adj); k = len(nodes)
    if k == 0: return 0
    idx = {v: i for i, v in enumerate(nodes)}
    A = [0] * k
    for u in adj:
        for v in adj[u]: A[idx[u]] |= 1 << idx[v]
    full = (1 << k) - 1
    from functools import lru_cache
    def q(S, v):
        seen = 1 << v; stack = [v]; out = 0
        while stack:
            x = stack.pop()
            nb = A[x] & ~seen; seen |= nb
            for y in range(k):
                if nb >> y & 1:
                    if S >> y & 1: stack.append(y)
                    else: out |= 1 << y
        return bin(out).count('1')
    @lru_cache(None)
    def f(S):
        if S == full: return 0
        return min(max(f(S | 1 << v), q(S, v)) for v in range(k) if not S >> v & 1)
    return f(0)

# ---------- partitions ----------
_PART = {}
def partitions(n):
    if n in _PART: return _PART[n]
    out = []
    def rec(pref, m):
        if len(pref) == n: out.append(pref[:]); return
        for l in range(m + 1):
            pref.append(l); rec(pref, max(m, l + 1)); pref.pop()
    rec([], 0)
    _PART[n] = np.array(out, dtype=np.int8)
    return _PART[n]

_DIFF = {}
def diff_matrix(g):
    """bool array (numpartitions x numpairs): True iff pair is edited by that partition"""
    key = (g.n, frozenset(g.E))
    if key in _DIFF: return _DIFF[key]
    P = partitions(g.n)
    ia = np.array([a for a, b in g.pairs]); ib = np.array([b for a, b in g.pairs])
    same = P[:, ia] == P[:, ib]
    adj = np.array([p in g.E for p in g.pairs])
    D = same != adj[None, :]
    if len(_DIFF) > 4: _DIFF.clear()
    _DIFF[key] = D
    return D

def opt_and_optima(g):
    D = diff_matrix(g)
    cost = D.sum(1)
    o = int(cost.min())
    rows = np.nonzero(cost == o)[0]
    return o, [frozenset(g.pairs[i] for i in np.nonzero(D[r])[0]) for r in rows]

def all_graphs(n, connected=True):
    """all graphs on n labelled-up-to-iso vertices via nauty geng"""
    import subprocess
    out = subprocess.run(['nauty-geng', '-q'] + (['-c'] if connected else []) + [str(n)], capture_output=True, text=True).stdout.split()
    import networkx as nx
    for s in out:
        G = nx.from_graph6_bytes(s.encode())
        yield Gr(n, G.edges())

def max_packings(g):
    """all maximum modification-disjoint packings (DFS with bound)"""
    P = g.p3(); pp = [pr(t) for t in P]
    best = [0, []]
    def rec(i, cur, used):
        if len(cur) + (len(P) - i) < best[0]: return
        if i == len(P):
            if len(cur) > best[0]: best[0] = len(cur); best[1] = [tuple(cur)]
            elif len(cur) == best[0]: best[1].append(tuple(cur))
            return
        if not any(p in used for p in pp[i]):
            for p in pp[i]: used.add(p)
            cur.append(P[i]); rec(i + 1, cur, used); cur.pop()
            for p in pp[i]: used.discard(p)
        rec(i + 1, cur, used)
    rec(0, [], set())
    return best[0], best[1]

def opt_ilp(g):
    """exact OPT for larger n via scipy milp (own formulation)"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = g.n; pid = {p: i for i, p in enumerate(g.pairs)}; m = len(pid)
    c = np.array([-1.0 if p in g.E else 1.0 for p in g.pairs])
    rows = []
    for a, b, d in itertools.combinations(range(n), 3):
        x, y, z = pid[(a, b)], pid[(b, d)], pid[(a, d)]
        for (u, v, w) in ((x, y, z), (x, z, y), (y, z, x)):
            r = np.zeros(m); r[u] = 1; r[v] = 1; r[w] = -1; rows.append(r)
    res = milp(c, constraints=LinearConstraint(np.array(rows), -np.inf, 1), integrality=np.ones(m), bounds=Bounds(0, 1))
    return int(round(res.fun)) + len(g.E), np.round(res.x).astype(int)

def nu_mod_ilp(g):
    from scipy.optimize import milp, LinearConstraint, Bounds
    P = g.p3()
    if not P: return 0
    pairs = {}
    for i, t in enumerate(P):
        for p in pr(t): pairs.setdefault(p, []).append(i)
    A = np.zeros((len(pairs), len(P)))
    for r, (p, l) in enumerate(pairs.items()):
        for i in l: A[r, i] = 1
    res = milp(-np.ones(len(P)), constraints=LinearConstraint(A, -np.inf, 1), integrality=np.ones(len(P)), bounds=Bounds(0, 1))
    return int(round(-res.fun))

def hard_of(H, S):
    o = soft_owner(H)
    return [p for p in S if p not in o]
