"""Fast exact cluster editing by subset DP (numba) with optional cluster-validity mask.
Graph: adjacency bitmasks adj[v]. cost = |E| + sum_clusters (nonedges_in - edges_in)."""
import itertools, random
import numpy as np
from numba import njit

@njit(cache=True)
def _inner(n, adj, allowed):
    N = 1 << n
    # cc[S] = nonedges(S) - edges(S)
    cc = np.zeros(N, np.int32)
    for S in range(1, N):
        v = 0
        while not (S >> v) & 1: v += 1
        T = S & (S - 1)
        k = 0; m = T
        while m:
            low = m & -m
            u = 0
            while not (low >> u) & 1: u += 1
            k += 1; m ^= low
        e = 0; m = T & adj[v]
        while m:
            m &= m - 1; e += 1
        cc[S] = cc[T] + (k - e) - e
    INF = 10**8
    f = np.full(N, INF, np.int32); f[0] = 0
    for S in range(1, N):
        low = S & -S; rest = S ^ low; sub = rest; best = INF
        while True:
            C = sub | low
            if allowed[C]:
                r = f[S ^ C]
                if r < INF:
                    c = r + cc[C]
                    if c < best: best = c
            if sub == 0: break
            sub = (sub - 1) & rest
        f[S] = best
    return f[N - 1]

def adjmasks(n, E):
    adj = np.zeros(n, np.int64)
    for e in E:
        a, b = tuple(e); adj[a] |= 1 << b; adj[b] |= 1 << a
    return adj

def opt_cost(n, E, allowed=None):
    adj = adjmasks(n, E)
    if allowed is None: allowed = np.ones(1 << n, np.bool_)
    r = _inner(n, adj, allowed)
    return len(E) + int(r)

def fs(*a): return frozenset(a)
def pairs(P):
    a, b, c = P
    return [fs(a, b), fs(b, c), fs(a, c)]
def induced_p3s(n, E):
    out = []
    for b in range(n):
        for a, c in itertools.combinations(range(n), 2):
            if a != b and c != b and fs(a, b) in E and fs(b, c) in E and fs(a, c) not in E:
                out.append((a, b, c))
    return out
def maximal_packing(n, E, rng, prefer=None):
    cand = induced_p3s(n, E); rng.shuffle(cand)
    if prefer is not None: cand.sort(key=lambda P: sum(x in prefer for x in P))
    used, H = set(), []
    for P in cand:
        ps = pairs(P)
        if not any(q in used for q in ps): H.append(P); used.update(ps)
    return H
def rclasses(n, E, H):
    Vh = {x for P in H for x in P}; R = [v for v in range(n) if v not in Vh]
    seen = set(); cls = []
    for r in R:
        if r in seen: continue
        comp = [r]; seen.add(r); st = [r]
        while st:
            x = st.pop()
            for y in R:
                if y not in seen and fs(x, y) in E: seen.add(y); comp.append(y); st.append(y)
        cls.append(sorted(comp))
    return Vh, R, cls
def mask(S):
    m = 0
    for v in S: m |= 1 << v
    return m
def rand_graph(rng, n, p): return {fs(a, b) for a, b in itertools.combinations(range(n), 2) if rng.random() < p}
def structured_graph(rng, n):
    """core random graph + attached cliques (twin classes) with disjoint attachment sets"""
    k = rng.randint(4, max(4, n - 2)); core = list(range(k)); p = rng.choice([.3, .5, .7])
    E = rand_graph(rng, k, p); nxt = k; pool = core[:]; rng.shuffle(pool)
    while nxt < n:
        s = min(rng.randint(1, 3), n - nxt); X = list(range(nxt, nxt + s)); nxt += s
        a = rng.randint(0, min(len(pool), 4)); A = [pool.pop() for _ in range(a)]
        for x, y in itertools.combinations(X, 2): E.add(fs(x, y))
        for x in X:
            for z in A: E.add(fs(x, z))
    return E
