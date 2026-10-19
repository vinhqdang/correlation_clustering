"""Independent audit library (written by the proof auditor, not the authors).
Graphs on n vertices: edge set as python int bitmask over C(n,2) pair indices.
All partitions enumerated via restricted-growth strings; per partition a bitmask of
'same cluster' pairs.  Solution S = E xor M.
"""
import itertools, random
import numpy as np

_cache = {}

def pidx(n):
    d = {}
    k = 0
    for u in range(n):
        for v in range(u + 1, n):
            d[(u, v)] = k; d[(v, u)] = k; k += 1
    return d

def all_partitions(n):
    """returns (masks uint64 array, labels list-of-tuples)"""
    if n in _cache: return _cache[n]
    P = pidx(n)
    masks = []; labs = []
    def rec(i, lab, m):
        if i == n:
            masks.append(m); labs.append(tuple(lab)); return
        for c in range(m + 1):
            lab.append(c); rec(i + 1, lab, max(m, c + 1)) ; lab.pop()
    # simple recursion
    def go(lab):
        if len(lab) == n:
            mk = 0
            for u in range(n):
                for v in range(u + 1, n):
                    if lab[u] == lab[v]: mk |= 1 << P[(u, v)]
            masks.append(mk); labs.append(tuple(lab)); return
        mx = max(lab) + 1 if lab else 0
        for c in range(mx + 1):
            go(lab + [c])
    go([])
    res = (np.array(masks, dtype=np.uint64), labs)
    _cache[n] = res
    return res

def popcount(a): return np.bitwise_count(a)

def bit(n, u, v): return 1 << pidx(n)[(u, v)]

def adj_mask(n, edges):
    P = pidx(n); m = 0
    for u, v in edges: m |= 1 << P[(u, v)]
    return m

def has(n, E, u, v): return (E >> pidx(n)[(u, v)]) & 1

def is_cluster(n, E):
    P = pidx(n)
    for b in range(n):
        for a in range(n):
            if a == b: continue
            if not (E >> P[(a, b)]) & 1: continue
            for c in range(a + 1, n):
                if c == b: continue
                if (E >> P[(b, c)]) & 1 and not (E >> P[(a, c)]) & 1: return False
    return True

def induced_p3s(n, E):
    P = pidx(n); out = []
    for b in range(n):
        for a in range(n):
            for c in range(a + 1, n):
                if a == b or c == b: continue
                if (E >> P[(a, b)]) & 1 and (E >> P[(b, c)]) & 1 and not (E >> P[(a, c)]) & 1:
                    out.append((a, b, c))
    return out

def pairmask(n, p):
    a, b, c = p; P = pidx(n)
    return (1 << P[(a, b)]) | (1 << P[(b, c)]) | (1 << P[(a, c)])

def rand_graph(rng, n, kind):
    edges = []
    if kind == 'gnp':
        p = rng.choice([0.25, 0.4, 0.55, 0.7])
        edges = [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < p]
    elif kind == 'hub':      # few hubs, adjacent to most, sparse elsewhere
        hubs = rng.sample(range(n), rng.choice([1, 2]))
        for u in range(n):
            for v in range(u + 1, n):
                if u in hubs or v in hubs:
                    if rng.random() < 0.8: edges.append((u, v))
                elif rng.random() < 0.1: edges.append((u, v))
    elif kind == 'cliques':  # noisy cluster graph
        lab = [rng.randrange(max(2, n // 3)) for _ in range(n)]
        for u in range(n):
            for v in range(u + 1, n):
                same = lab[u] == lab[v]
                if (same and rng.random() < 0.85) or ((not same) and rng.random() < 0.12): edges.append((u, v))
    elif kind == 'path':     # path / cycle-like sparse
        for u in range(n - 1):
            if rng.random() < 0.9: edges.append((u, u + 1))
        for _ in range(rng.randrange(3)):
            u, v = rng.sample(range(n), 2)
            edges.append((min(u, v), max(u, v)))
        edges = list(set(edges))
    return adj_mask(n, edges)

def rand_packing(rng, n, E, mode):
    """random modification-disjoint packing.  mode: 'max' (random greedy maximal),
    'part' (random greedy, stops early), 'overlap' (prefers P3s that share vertices)."""
    cand = induced_p3s(n, E)
    rng.shuffle(cand)
    used = 0; H = []
    lim = len(cand) if mode == 'max' else rng.randrange(0, len(cand) + 1) if cand else 0
    if mode == 'overlap':
        # sort by number of vertices already covered, descending, repeatedly
        pool = cand[:]
        while pool:
            cov = {}
            for P in H:
                for x in P: cov[x] = cov.get(x, 0) + 1
            pool.sort(key=lambda P: -sum(cov.get(x, 0) for x in P) + rng.random())
            P = pool.pop(0)
            m = pairmask(n, P)
            if used & m: continue
            H.append(P); used |= m
            if len(H) >= lim and rng.random() < 0.3: break
        return H
    for P in cand:
        if len(H) >= lim: break
        m = pairmask(n, P)
        if used & m: continue
        H.append(P); used |= m
    return H

def check_packing(n, E, H):
    seen = 0
    for (a, b, c) in H:
        assert has(n, E, a, b) and has(n, E, b, c) and not has(n, E, a, c)
        m = pairmask(n, (a, b, c))
        assert not (seen & m)
        seen |= m
    return seen

def opt_and_costs(n, E):
    masks, labs = all_partitions(n)
    S = masks ^ np.uint64(E)
    return S, popcount(S).astype(np.int64)
