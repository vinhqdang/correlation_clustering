"""Audit 5 library, written from scratch (no import of a4_lib / w5a_*). Graph = set of frozenset edges / adjacency bitmasks.
packed (soft) pair = two vertices in a common member of H; hard = other pair."""
import itertools, random
def partitions(n):
    def rec(i, lab, m):
        if i == n: yield tuple(lab); return
        for c in range(m + 1):
            lab.append(c); yield from rec(i + 1, lab, max(m, c + 1)); lab.pop()
    yield from rec(0, [], 0)
def cost(adj, lab):
    n = len(adj); s = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (lab[i] == lab[j]) != bool(adj[i][j]): s += 1
    return s
def optima(adj):
    best = 10**9; out = []
    for lab in partitions(len(adj)):
        c = cost(adj, lab)
        if c < best: best = c; out = [lab]
        elif c == best: out.append(lab)
    return best, out
def mk(n, edges):
    a = [[0]*n for _ in range(n)]
    for x, y in edges: a[x][y] = a[y][x] = 1
    return a
def p3s(adj):
    n = len(adj); r = []
    for b in range(n):
        nb = [x for x in range(n) if adj[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not adj[a][c]: r.append((a, b, c))
    return r
def pairs3(t): return {frozenset(p) for p in itertools.combinations(t, 2)}
def is_packing(adj, H):
    seen = set()
    for (a, b, c) in H:
        if not (adj[a][b] and adj[b][c] and not adj[a][c]): return False
        p = pairs3((a, b, c))
        if p & seen: return False
        seen |= p
    return True
def packed(H): 
    s = set()
    for t in H: s |= pairs3(t)
    return s
def pair_maximal(adj, H):
    pk = packed(H)
    return all(pairs3(t) & pk for t in p3s(adj))
def rand_pm(adj, rng, stop=0.0):
    ts = p3s(adj); rng.shuffle(ts); H = []; used = set()
    for t in ts:
        if rng.random() < stop: continue
        p = pairs3(t)
        if not p & used: H.append(t); used |= p
    for t in ts:                       # extend to pair-maximal
        p = pairs3(t)
        if not p & used: H.append(t); used |= p
    return H
def pieces(adj, pk, lab):
    n = len(adj); par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j] and lab[i] == lab[j] and frozenset((i, j)) not in pk: par[f(i)] = f(j)
    return [f(i) for i in range(n)]
def U_set(adj, pk, lab):
    n = len(adj)
    return [(i, j) for i in range(n) for j in range(i+1, n)
            if (lab[i] == lab[j]) != bool(adj[i][j]) and frozenset((i, j)) not in pk]
