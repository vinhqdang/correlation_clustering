"""Wave 5 prover B library (own code). Graph = n x n 0/1 numpy matrix; P3 = (a,b,c), b centre.
soft = pair inside a member of H ("packed"); hard = every other pair.  hard edge = hard pair in E."""
import itertools, random
import numpy as np

def prs(n): return list(itertools.combinations(range(n), 2))
def key(a, b): return (a, b) if a < b else (b, a)
def p3_list(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a][c]: out.append((a, b, c))
    return out
def p3_pairs(t):
    a, b, c = t; return {key(a, b), key(b, c), key(a, c)}
def used_pairs(H):
    s = set()
    for t in H: s |= p3_pairs(t)
    return s
def pair_maximal_packing(A, rng):
    ts = p3_list(A); rng.shuffle(ts); H = []; used = set()
    for t in ts:
        ps = p3_pairs(t)
        if not (ps & used): H.append(t); used |= ps
    return H
def soft_matrix(n, H):
    pk = np.zeros((n, n), bool)
    for t in H:
        for x, y in itertools.combinations(t, 2): pk[x, y] = pk[y, x] = True
    return pk
def hard_comps(A, pk):
    n = len(A); hg = (A.astype(bool) & ~pk); comp = [-1] * n; k = 0
    for s in range(n):
        if comp[s] < 0:
            st = [s]; comp[s] = k
            while st:
                v = st.pop()
                for u in range(n):
                    if hg[v, u] and comp[u] < 0: comp[u] = k; st.append(u)
            k += 1
    return comp, hg

_PART = {}
def all_partitions(n):
    """all set partitions as label arrays (restricted growth strings), shape (Bell(n), n), int8"""
    if n in _PART: return _PART[n]
    cur = np.zeros((1, 1), np.int8)
    for m in range(1, n):
        mx = cur.max(1) + 1
        rows = []
        for k in range(m + 1):
            sel = mx >= k
            if k == 0: sel = np.ones(len(cur), bool)
            else: sel = mx >= k  # label k allowed iff k <= max+1 i.e. mx>=k
            blk = cur[sel]
            rows.append(np.hstack([blk, np.full((len(blk), 1), k, np.int8)]))
        cur = np.vstack(rows)
    _PART[n] = cur; return cur
def costs(A, P):
    n = len(A); c = np.zeros(len(P), np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            same = P[:, i] == P[:, j]
            c += (same != bool(A[i][j]))
    return c
def rand_graph(n, rng):
    kind = rng.choice(['gnp', 'planted', 'hubby', 'sparse'])
    A = np.zeros((n, n), int)
    if kind == 'gnp':
        p = rng.choice([.25, .4, .55, .7])
        for i, j in prs(n):
            if rng.random() < p: A[i][j] = A[j][i] = 1
    elif kind == 'sparse':
        for i, j in prs(n):
            if rng.random() < .2: A[i][j] = A[j][i] = 1
    elif kind == 'planted':
        k = rng.randint(2, 3); lab = [rng.randrange(k) for _ in range(n)]; q = rng.choice([.05, .1, .2])
        for i, j in prs(n):
            e = (lab[i] == lab[j]); e = (not e) if rng.random() < q else e
            if e: A[i][j] = A[j][i] = 1
    else:  # a few large near-cliques with attachments: large hard degrees
        k = rng.randint(2, 3); lab = [rng.randrange(k) for _ in range(n)]
        for i, j in prs(n):
            if lab[i] == lab[j] and rng.random() < .85: A[i][j] = A[j][i] = 1
            elif rng.random() < .08: A[i][j] = A[j][i] = 1
    return A
