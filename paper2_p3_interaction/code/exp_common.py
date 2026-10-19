"""Shared helpers for the P3-packing / interaction-graph experiments.

Packings are sets of induced P3s (a,b,c) (edges ab, bc, non-edge ac) whose three
vertex pairs are pairwise disjoint over the packing (modification-disjoint).
All randomness is seeded.  Tree-width numbers are heuristic UPPER bounds
(elimination orderings); lower bounds are stated separately.
"""
import gzip, heapq, os
import numpy as np
from numba import njit

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# ---------------------------------------------------------------- readers
def read_pace(path):
    edges, n = [], 0
    for line in open(path):
        t = line.split()
        if not t or t[0] == 'c':
            continue
        if t[0] == 'p':
            n = int(t[2]); continue
        u, v = int(t[0]) - 1, int(t[1]) - 1
        if u != v:
            edges.append((u, v))
    return n, edges

def read_snap(path, signcol=-1):
    ids, edges = {}, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            if line[0] in '#%':
                continue
            t = line.replace(',', ' ').split()
            if len(t) < 2:
                continue
            if signcol >= 0 and float(t[signcol]) <= 0:
                continue
            u, v = t[0], t[1]
            if u == v:
                continue
            for x in (u, v):
                if x not in ids:
                    ids[x] = len(ids)
            edges.append((ids[u], ids[v]))
    return len(ids), edges

def adjacency(n, edges):
    A = np.zeros((n, n), dtype=np.uint8)
    for u, v in edges:
        A[u, v] = 1; A[v, u] = 1
    return A

def csr(A):
    n = A.shape[0]
    deg = A.sum(1).astype(np.int64)
    indptr = np.zeros(n + 1, dtype=np.int64); indptr[1:] = np.cumsum(deg)
    indices = np.nonzero(A)[1].astype(np.int32)   # row-major order -> grouped by row
    return indptr, indices

# ---------------------------------------------------------------- packers
@njit(cache=True)
def greedy_pack(A, indptr, indices, seed, cap):
    """One pass over centres in random order; for a centre b scan neighbour pairs
    in random order and take every induced P3 whose 3 pairs are unused (and whose
    vertices all carry < cap packed P3s).  cap >= n gives an unrestricted packing,
    which is maximal; cap < n need not be."""
    np.random.seed(seed)
    n = A.shape[0]
    used = np.zeros((n, n), dtype=np.uint8)
    load = np.zeros(n, dtype=np.int32)
    out = np.empty((n * 4 + 16, 3), dtype=np.int32)
    cnt = 0
    cen = np.random.permutation(n)
    for b in cen:
        d = indptr[b + 1] - indptr[b]
        if d < 2:
            continue
        nb = indices[indptr[b]:indptr[b + 1]].copy()
        nb = nb[np.random.permutation(d)]
        for i in range(d):
            a = nb[i]
            if used[a, b] or load[a] >= cap or load[b] >= cap:
                continue
            for j in range(i + 1, d):
                c = nb[j]
                if A[a, c] or used[b, c] or used[a, c] or load[c] >= cap:
                    continue
                used[a, b] = used[b, a] = 1
                used[b, c] = used[c, b] = 1
                used[a, c] = used[c, a] = 1
                load[a] += 1; load[b] += 1; load[c] += 1
                if cnt == out.shape[0]:
                    out2 = np.empty((out.shape[0] * 2, 3), dtype=np.int32)
                    out2[:cnt] = out
                    out = out2
                out[cnt, 0] = a; out[cnt, 1] = b; out[cnt, 2] = c
                cnt += 1
                break
            if load[b] >= cap:
                break
    return out[:cnt].copy()

@njit(cache=True)
def local_search(A, H0, iters, seed, tries):
    """(1,k)-swap plateau search.  Repeat: pick a random packed P3 Q, free its 3
    pairs, collect every induced P3 that is entirely free (every such P3 must use
    a freed pair, by maximality of the current packing), then re-pack greedily in
    random order (best of `tries` random orders).  The result has >= 1 P3, so
    |H| never decreases; k >= 2 is a strict improvement, k = 1 a plateau move
    (a different P3 may replace Q).  Packing stays maximal."""
    np.random.seed(seed)
    n = A.shape[0]
    cap = H0.shape[0] * 3 + iters * 2 + 16
    P = np.empty((cap, 3), dtype=np.int32)
    m = H0.shape[0]
    P[:m] = H0
    owner = -np.ones((n, n), dtype=np.int32)
    for i in range(m):
        a, b, c = P[i, 0], P[i, 1], P[i, 2]
        owner[a, b] = owner[b, a] = i
        owner[b, c] = owner[c, b] = i
        owner[a, c] = owner[c, a] = i
    cand = np.empty((6 * n + 6, 3), dtype=np.int32)
    mark = np.zeros((n, n), dtype=np.uint8)
    best = np.empty((8, 3), dtype=np.int32)
    cur = np.empty((8, 3), dtype=np.int32)
    gains = 0
    for it in range(iters):
        if m == 0:
            break
        q = np.random.randint(m)
        qa, qb, qc = P[q, 0], P[q, 1], P[q, 2]
        owner[qa, qb] = owner[qb, qa] = -1
        owner[qb, qc] = owner[qc, qb] = -1
        owner[qa, qc] = owner[qc, qa] = -1
        nc = 0
        # pairs of Q: (qa,qb) edge, (qb,qc) edge, (qa,qc) non-edge
        for t in range(3):
            if t == 0:
                u, v = qa, qb
            elif t == 1:
                u, v = qb, qc
            else:
                u, v = qa, qc
            if A[u, v]:
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if A[u, w] and not A[v, w] and owner[u, w] < 0 and owner[v, w] < 0:
                        cand[nc, 0] = w; cand[nc, 1] = u; cand[nc, 2] = v; nc += 1
                    if A[v, w] and not A[u, w] and owner[v, w] < 0 and owner[u, w] < 0:
                        cand[nc, 0] = u; cand[nc, 1] = v; cand[nc, 2] = w; nc += 1
            else:
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if A[u, w] and A[v, w] and owner[u, w] < 0 and owner[v, w] < 0:
                        cand[nc, 0] = u; cand[nc, 1] = w; cand[nc, 2] = v; nc += 1
        bk = 0
        for tr in range(tries):
            perm = np.random.permutation(nc)
            k = 0
            for idx in perm:
                a, b, c = cand[idx, 0], cand[idx, 1], cand[idx, 2]
                if mark[a, b] or mark[b, c] or mark[a, c]:
                    continue
                mark[a, b] = mark[b, a] = 1; mark[b, c] = mark[c, b] = 1; mark[a, c] = mark[c, a] = 1
                if k < 8:
                    cur[k, 0] = a; cur[k, 1] = b; cur[k, 2] = c
                k += 1
            for idx in range(k if k < 8 else 8):
                a, b, c = cur[idx, 0], cur[idx, 1], cur[idx, 2]
                mark[a, b] = mark[b, a] = 0; mark[b, c] = mark[c, b] = 0; mark[a, c] = mark[c, a] = 0
            if k > 8:  # cannot happen for a P3 repack (<=3 pairs freed -> <=3 P3s)
                k = 8
            if k > bk:
                bk = k
                best[:k] = cur[:k]
                if bk >= 2 and tr >= 1:
                    break
        # (a freed Q has 3 pairs, so at most 3 disjoint new P3s: k<=3 <8)
        # write result: first new P3 overwrites slot q, others appended
        for j in range(bk):
            if j == 0:
                s = q
            else:
                s = m; m += 1
            P[s, 0], P[s, 1], P[s, 2] = best[j, 0], best[j, 1], best[j, 2]
            a, b, c = best[j, 0], best[j, 1], best[j, 2]
            owner[a, b] = owner[b, a] = s
            owner[b, c] = owner[c, b] = s
            owner[a, c] = owner[c, a] = s
        if bk == 0:  # should not happen (Q itself is a candidate); restore
            owner[qa, qb] = owner[qb, qa] = q
            owner[qb, qc] = owner[qc, qb] = q
            owner[qa, qc] = owner[qc, qa] = q
        gains += max(bk - 1, 0)
    return P[:m].copy(), gains

def verify_packing(A, H):
    """Raise if H is not a modification-disjoint induced P3 packing."""
    seen = set()
    for a, b, c in H:
        assert A[a, b] and A[b, c] and not A[a, c] and a != c
        for p in ((a, b), (b, c), (a, c)):
            k = (min(p), max(p))
            assert k not in seen, 'pair reused'
            seen.add(k)
    return True

def is_maximal(A, H):
    """No induced P3 of G has all three pairs outside the packing."""
    n = A.shape[0]
    used = np.zeros((n, n), dtype=bool)
    for a, b, c in H:
        for x, y in ((a, b), (b, c), (a, c)):
            used[x, y] = used[y, x] = True
    for b in range(n):
        nb = np.nonzero(A[b])[0]
        for i, a in enumerate(nb):
            if used[a, b]:
                continue
            for c in nb[i + 1:]:
                if not A[a, c] and not used[b, c] and not used[a, c]:
                    return False
    return True

# ---------------------------------------------------------------- I(H)
def interaction_stats(H, n, tw_cap=100, edge_cap=4_000_000):
    """Statistics of the intersection graph I(H) (nodes = P3s, adjacent iff they
    share a vertex).  Returns dict.  tw_ub is an elimination-ordering UPPER bound
    (None if it exceeded tw_cap or I(H) was too large to build); tw_lb is
    max(c-1, degeneracy) (a LOWER bound, degeneracy only when I(H) was built)."""
    p = len(H)
    lists = [[] for _ in range(n)]
    for i, t in enumerate(H):
        for v in t:
            lists[v].append(i)
    load = np.array([len(l) for l in lists])
    c = int(load.max()) if p else 0
    # components via union-find over vertices
    par = list(range(p))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for l in lists:
        for j in l[1:]:
            ra, rb = find(l[0]), find(j)
            if ra != rb:
                par[ra] = rb
    comps = len({find(i) for i in range(p)}) if p else 0
    est = int((load.astype(np.int64) ** 2).sum() // 2)
    res = dict(p=p, c=c, comps=comps, tw_lb=max(c - 1, 0), tw_ub=None, maxdeg=None,
               avgdeg=None, tw_ub_minfill=None, tw_ub_mindeg=None)
    # max degree without building edges
    nbrs_size = np.zeros(p, dtype=np.int64)
    if est <= edge_cap:
        adj = [set() for _ in range(p)]
        for l in lists:
            for x in l:
                adj[x].update(l)
        for i in range(p):
            adj[i].discard(i)
        degs = np.array([len(a) for a in adj])
        res['maxdeg'] = int(degs.max()) if p else 0
        res['avgdeg'] = float(degs.mean()) if p else 0.0
        res['tw_lb'] = max(res['tw_lb'], degeneracy(adj))
        if res['tw_lb'] <= tw_cap:
            w1, _ = elimination(adj, tw_cap, 'mindeg')
            w2 = None
            if p <= 3000:
                w2, _ = elimination(adj, tw_cap, 'minfill')
            res['tw_ub_mindeg'] = w1; res['tw_ub_minfill'] = w2
            ws = [w for w in (w1, w2) if w is not None]
            res['tw_ub'] = min(ws) if ws else None
    else:
        # degree of P3 i = |union of the vertex lists| - 1 (no edge set)
        md = 0
        for i, t in enumerate(H):
            s = set(lists[t[0]]) | set(lists[t[1]]) | set(lists[t[2]])
            md = max(md, len(s) - 1)
        res['maxdeg'] = md
    return res

def degeneracy(adj):
    p = len(adj)
    deg = [len(a) for a in adj]
    alive = [True] * p
    h = [(d, i) for i, d in enumerate(deg)]
    heapq.heapify(h)
    best = 0
    while h:
        d, i = heapq.heappop(h)
        if not alive[i] or d != deg[i]:
            continue
        best = max(best, d); alive[i] = False
        for j in adj[i]:
            if alive[j]:
                deg[j] -= 1
                heapq.heappush(h, (deg[j], j))
    return best

def elimination(adj, cap, rule):
    """Greedy elimination ordering on a copy of adj (list of sets).  Returns
    (width, order); width None if it would exceed cap.  rule in {'mindeg','minfill'}.
    width = max degree at elimination = upper bound on treewidth."""
    p = len(adj)
    g = [set(a) for a in adj]
    alive = set(range(p))
    order, width = [], 0
    if rule == 'mindeg':
        h = [(len(g[i]), i) for i in range(p)]
        heapq.heapify(h)
        while h:
            d, i = heapq.heappop(h)
            if i not in alive or d != len(g[i]):
                continue
            if d > cap:
                return None, order
            width = max(width, d)
            nb = list(g[i])
            for x in nb:
                g[x].discard(i)
            for a in range(len(nb)):
                ga = g[nb[a]]
                for b in range(a + 1, len(nb)):
                    ga.add(nb[b]); g[nb[b]].add(nb[a])
            alive.discard(i); order.append(i); g[i] = set()
            for x in nb:
                heapq.heappush(h, (len(g[x]), x))
        return width, order
    # min-fill with lazy recomputation
    def fill(i):
        nb = list(g[i]); f = 0
        for a in range(len(nb)):
            ga = g[nb[a]]
            for b in range(a + 1, len(nb)):
                if nb[b] not in ga:
                    f += 1
        return f
    score = {i: fill(i) for i in alive}
    while alive:
        i = min(alive, key=lambda x: (score[x], len(g[x]), x))
        d = len(g[i])
        if d > cap:
            return None, order
        width = max(width, d)
        nb = list(g[i])
        for x in nb:
            g[x].discard(i)
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                g[nb[a]].add(nb[b]); g[nb[b]].add(nb[a])
        alive.discard(i); order.append(i); g[i] = set()
        affected = set(nb)
        for x in nb:
            affected |= g[x]
        for x in affected:
            if x in alive:
                score[x] = fill(x)
        del score[i]
    return width, order
