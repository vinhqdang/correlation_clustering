"""Audit 6 library (written from scratch).  G = symmetric 0/1 numpy matrix, P3 = (a,b,c) with b the centre."""
import itertools, random
import numpy as np
from collections import deque

def pairs(n): return list(itertools.combinations(range(n), 2))

def induced_p3s(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a][c]: out.append((a, b, c))
    return out

def p3pairs(t):
    a, b, c = t
    return {tuple(sorted(p)) for p in ((a, b), (b, c), (a, c))}

def greedy_pair_maximal(A, rng, start=()):
    """pair-maximal packing of induced P3s, optionally extending a partial packing `start`."""
    H = list(start); used = set()
    for t in H: used |= p3pairs(t)
    ts = induced_p3s(A); rng.shuffle(ts)
    for t in ts:
        ps = p3pairs(t)
        if not (ps & used): H.append(t); used |= ps
    return H

def soft_matrix(n, H):
    S = np.zeros((n, n), bool)
    for t in H:
        for x, y in itertools.combinations(t, 2): S[x, y] = S[y, x] = True
    return S

def hard_components(A, S):
    n = len(A); hg = A.astype(bool) & ~S; comp = [-1] * n; k = 0
    for s in range(n):
        if comp[s] < 0:
            comp[s] = k; st = [s]
            while st:
                v = st.pop()
                for u in range(n):
                    if hg[v, u] and comp[u] < 0: comp[u] = k; st.append(u)
            k += 1
    return comp, hg

def partitions(n):
    """all set partitions of range(n) as restricted-growth label arrays, shape (Bell(n), n)."""
    cur = np.zeros((1, 1), np.int8)
    for m in range(1, n):
        mx = cur.max(1) + 1; rows = []
        for k in range(m + 1):
            sel = mx >= k
            blk = cur[sel]
            rows.append(np.hstack([blk, np.full((len(blk), 1), k, np.int8)]))
        cur = np.vstack(rows)
    return cur

def part_costs(A, P):
    n = len(A); c = np.zeros(len(P), np.int32)
    for i, j in pairs(n): c += ((P[:, i] == P[:, j]) != bool(A[i][j]))
    return c

def hard_edit_masks(A, S, P):
    """int64 bitmask over hard pairs (index in `hp`) of the hard pairs edited by each partition."""
    n = len(A); hp = [(i, j) for i, j in pairs(n) if not S[i, j]]
    m = np.zeros(len(P), np.int64)
    for k, (i, j) in enumerate(hp):
        ed = ((P[:, i] == P[:, j]) != bool(A[i][j]))
        m |= ed.astype(np.int64) << k
    return hp, m

def bfs_ball(A, allowed, src, R):
    dist = {src: 0}; q = deque([src])
    while q:
        x = q.popleft()
        if dist[x] == R: continue
        for y in allowed:
            if A[x][y] and y not in dist: dist[y] = dist[x] + 1; q.append(y)
    return dist

def rand_graph(n, rng):
    kind = rng.choice(['gnp', 'planted', 'hubby', 'sparse', 'hub2'])
    A = np.zeros((n, n), int)
    def add(i, j): A[i][j] = A[j][i] = 1
    if kind == 'gnp':
        p = rng.choice([.25, .4, .55, .7])
        for i, j in pairs(n):
            if rng.random() < p: add(i, j)
    elif kind == 'sparse':
        for i, j in pairs(n):
            if rng.random() < .2: add(i, j)
    elif kind == 'planted':
        k = rng.randint(2, 3); lab = [rng.randrange(k) for _ in range(n)]; q = rng.choice([.05, .1, .2])
        for i, j in pairs(n):
            e = lab[i] == lab[j]; e = (not e) if rng.random() < q else e
            if e: add(i, j)
    elif kind == 'hubby':
        k = rng.randint(2, 3); lab = [rng.randrange(k) for _ in range(n)]
        for i, j in pairs(n):
            if lab[i] == lab[j] and rng.random() < .85: add(i, j)
            elif rng.random() < .08: add(i, j)
    else:  # big clique-ish hard component plus pendant / low-degree attachments (BIG components containing low-degree vertices)
        m = rng.randint(max(4, n // 2), n - 1)
        for i, j in itertools.combinations(range(m), 2):
            if rng.random() < .9: add(i, j)
        for v in range(m, n):
            for u in rng.sample(range(n), rng.randint(1, 2)):
                if u != v: add(u, v)
    return A
