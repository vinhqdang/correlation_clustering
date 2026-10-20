"""Wave 8 library (repacking).  Own code; imports numpy, networkx, scipy only (+ w7_opt.solve for exact optima, n<=14).
Graph: n, edge set E of frozenset-free tuples (u<v).  Packing member = (a,b,c) centre b, ab,bc in E, ac not in E.
"""
import itertools, sys
import numpy as np
import networkx as nx
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from w7_opt import solve as _solve

def key(a, b): return (a, b) if a < b else (b, a)

class G8:
    def __init__(self, n, edges):
        self.n = n
        self.E = set(key(a, b) for a, b in edges)
        self.adj = [set() for _ in range(n)]
        for a, b in self.E:
            self.adj[a].add(b); self.adj[b].add(a)
        self._p3 = None
    def isE(self, a, b): return key(a, b) in self.E
    def p3s(self):
        """all induced P3 as (a,b,c) centre b, a<c"""
        if self._p3 is None:
            res = []
            for b in range(self.n):
                nb = sorted(self.adj[b])
                for i in range(len(nb)):
                    for j in range(i + 1, len(nb)):
                        a, c = nb[i], nb[j]
                        if not self.isE(a, c): res.append((a, b, c))
            self._p3 = res
        return self._p3
    def opt(self, enum=True, limit=200000):
        return _solve(self.n, sorted(self.E), enum=enum, limit=limit)

def mp(t):
    a, b, c = t
    return [key(a, b), key(b, c), key(a, c)]

def valid_packing(g, H):
    seen = set()
    for t in H:
        a, b, c = t
        if not (g.isE(a, b) and g.isE(b, c) and not g.isE(a, c)): return False
        for p in mp(t):
            if p in seen: return False
            seen.add(p)
    return True

def soft_pairs(H):
    s = {}
    for i, t in enumerate(H):
        for p in mp(t): s[p] = i
    return s

def pair_maximal(g, H):
    s = soft_pairs(H)
    for t in g.p3s():
        if not any(p in s for p in mp(t)): return False
    return True

def maximalize(g, H):
    H = list(H); s = soft_pairs(H)
    for t in g.p3s():
        if not any(p in s for p in mp(t)):
            H.append(t)
            for p in mp(t): s[p] = len(H) - 1
    return H

def Igraph(H):
    I = nx.Graph(); I.add_nodes_from(range(len(H)))
    for i in range(len(H)):
        for j in range(i + 1, len(H)):
            if set(H[i]) & set(H[j]): I.add_edge(i, j)
    return I

def tw_exact(I, cap=20):
    """exact treewidth by DP over vertex subsets of elimination orderings (<= ~18 nodes); else min-fill upper bound"""
    nodes = list(I.nodes)
    k = len(nodes)
    if k == 0: return 0
    if k > cap:
        from networkx.algorithms.approximation import treewidth_min_fill_in
        return treewidth_min_fill_in(I)[0]
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * k
    for u, v in I.edges: adj[idx[u]] |= 1 << idx[v]; adj[idx[v]] |= 1 << idx[u]
    from functools import lru_cache
    full = (1 << k) - 1
    def q(S, v):
        # number of vertices outside S u {v} reachable from v through S
        seen = 1 << v; stack = [v]; cnt = 0
        reach_out = 0
        while stack:
            x = stack.pop()
            nb = adj[x] & ~seen
            seen |= nb
            for y in range(k):
                if nb >> y & 1:
                    if S >> y & 1: stack.append(y)
                    else: reach_out |= 1 << y
        return bin(reach_out).count('1')
    @lru_cache(None)
    def f(S):
        if S == full: return 0
        best = 99
        for v in range(k):
            if not S >> v & 1:
                c = max(f(S | 1 << v), q(S, v))
                if c < best: best = c
        return best
    return f(0)

def enum_packings(g, minsize=0, maxcount=10**6):
    """enumerate all modification-disjoint packings (as sorted tuples of triples) of size >= minsize.  DFS in P3 index order."""
    P = g.p3s()
    pp = [mp(t) for t in P]
    out = []
    def rec(i, cur, used):
        if len(out) >= maxcount: return
        if i == len(P):
            if len(cur) >= minsize: out.append(tuple(cur))
            return
        # prune: remaining can add at most len(P)-i
        if len(cur) + (len(P) - i) < minsize: return
        # take
        if not any(p in used for p in pp[i]):
            for p in pp[i]: used.add(p)
            cur.append(P[i]); rec(i + 1, cur, used); cur.pop()
            for p in pp[i]: used.discard(p)
        rec(i + 1, cur, used)
    rec(0, [], set())
    return out

def nu_mod(g):
    """max modification-disjoint packing by ILP (scipy HiGHS); returns (size, packing)"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    P = g.p3s()
    if not P: return 0, []
    pairs = {}
    for i, t in enumerate(P):
        for p in mp(t): pairs.setdefault(p, []).append(i)
    A = np.zeros((len(pairs), len(P)))
    for r, (p, l) in enumerate(pairs.items()):
        for i in l: A[r, i] = 1
    res = milp(-np.ones(len(P)), constraints=LinearConstraint(A, -np.inf, 1), integrality=np.ones(len(P)), bounds=Bounds(0, 1))
    x = np.round(res.x).astype(int)
    return int(x.sum()), [P[i] for i in range(len(P)) if x[i]]

def opt_ilp(g, extra=None):
    """exact OPT by ILP on pair variables x_e=1 iff together (transitivity); returns (opt, partition list)"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = g.n
    pid = {}
    for a in range(n):
        for b in range(a + 1, n): pid[(a, b)] = len(pid)
    m = len(pid)
    c = np.array([-1.0 if p in g.E else 1.0 for p in pid])   # cost: nonedge together +1, edge together -1 (+|E| const)
    rows = []
    for a, b, d in itertools.combinations(range(n), 3):
        ab, bd, ad = pid[(a, b)], pid[(b, d)], pid[(a, d)]
        for (x, y, z) in ((ab, bd, ad), (ab, ad, bd), (bd, ad, ab)):
            r = np.zeros(m); r[x] = 1; r[y] = 1; r[z] = -1; rows.append(r)   # x+y-z<=1
    A = np.array(rows)
    res = milp(c, constraints=LinearConstraint(A, -np.inf, 1), integrality=np.ones(m), bounds=Bounds(0, 1))
    x = np.round(res.x).astype(int)
    opt = int(round(res.fun)) + len(g.E)
    # partition
    comp = list(range(n))
    def find(u):
        while comp[u] != u: comp[u] = comp[comp[u]]; u = comp[u]
        return u
    for (a, b), i in pid.items():
        if x[i]: comp[find(a)] = find(b)
    cl = {}
    for v in range(n): cl.setdefault(find(v), []).append(v)
    return opt, list(cl.values())

def edits_of_partition(g, part):
    lab = {}
    for i, C in enumerate(part):
        for v in C: lab[v] = i
    S = set()
    for a in range(g.n):
        for b in range(a + 1, g.n):
            tog = lab[a] == lab[b]
            if tog != g.isE(a, b): S.add((a, b))
    return S

def part_from_mask(masks, n):
    return [[v for v in range(n) if C >> v & 1] for C in masks]

def hard_edits(H, S):
    s = soft_pairs(H)
    return [p for p in S if p not in s]

def member_dist(H, a, b, I=None):
    """min I(H)-distance between a member containing a and a member containing b; None-> a or b outside V(H) (returns -1); inf if disconnected"""
    A = [i for i, t in enumerate(H) if a in t]; B = [i for i, t in enumerate(H) if b in t]
    if not A or not B: return -1
    if I is None: I = Igraph(H)
    best = 10**9
    for x in A:
        d = nx.single_source_shortest_path_length(I, x)
        for y in B:
            if y in d: best = min(best, d[y])
    return best
