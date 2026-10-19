"""Audit 7 library (own).  Graph = numpy 0/1 matrix, P3 = (a,b,c) with b the centre."""
import itertools, random
import numpy as np

def prs(n): return list(itertools.combinations(range(n), 2))

def ind_p3s(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b][x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a][c]: out.append((a, b, c))
    return out

def p3p(t):
    a, b, c = t
    return {tuple(sorted(p)) for p in ((a, b), (b, c), (a, c))}

def greedy_pm(A, rng, start=()):
    H = list(start); used = set()
    for t in H: used |= p3p(t)
    ts = ind_p3s(A); rng.shuffle(ts)
    for t in ts:
        ps = p3p(t)
        if not (ps & used): H.append(t); used |= ps
    return H

def soft_of(n, H):
    S = np.zeros((n, n), bool)
    for t in H:
        for x, y in itertools.combinations(t, 2): S[x, y] = S[y, x] = True
    return S

def is_pair_maximal(A, H):
    S = soft_of(len(A), H)
    for (a, b, c) in ind_p3s(A):
        if not (S[a, b] or S[b, c] or S[a, c]): return False
    return True

def is_packing(A, H):
    used = set()
    for t in H:
        a, b, c = t
        if not (A[a][b] and A[b][c] and not A[a][c]): return False
        ps = p3p(t)
        if ps & used: return False
        used |= ps
    return True

_PC = {}
def partitions(n):
    if n in _PC: return _PC[n]
    cur = np.zeros((1, 1), np.int8)
    for m in range(1, n):
        mx = cur.max(1) + 1; rows = []
        for k in range(m + 1):
            blk = cur[mx >= k]
            rows.append(np.hstack([blk, np.full((len(blk), 1), k, np.int8)]))
        cur = np.vstack(rows)
    _PC[n] = cur
    return cur

def same_matrix(P):
    n = P.shape[1]; pl = prs(n)
    X = np.zeros((len(P), len(pl)), bool)
    for k, (i, j) in enumerate(pl): X[:, k] = P[:, i] == P[:, j]
    return X, pl

def costs(A, P):
    X, pl = same_matrix(P)
    e = np.array([bool(A[i][j]) for i, j in pl])
    return (X != e).sum(1), X, pl

def opt_subset(A):
    """exact cluster editing optimum by subset DP (n <= ~12)."""
    n = len(A); N = 1 << n
    adj = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    cost = [0] * N   # cost of making subset one cluster: nonedges inside + (edges leaving counted separately)
    inner_e = [0] * N
    for m in range(1, N):
        low = (m & -m).bit_length() - 1; r = m & (m - 1)
        inner_e[m] = inner_e[r] + bin(adj[low] & r).count("1")
    tot_e = sum(bin(a).count("1") for a in adj) // 2
    # cost(partition)= |E| + sum_C (C(|C|,2) - 2 e(C))
    f = [0] * N
    for m in range(N):
        s = bin(m).count("1"); f[m] = s * (s - 1) // 2 - 2 * inner_e[m]
    INF = 10 ** 9; dp = [INF] * N; dp[0] = 0
    for m in range(1, N):
        low = m & -m; rest = m ^ low; sub = rest; best = INF
        while True:
            c = sub | low
            v = dp[m ^ c] + f[c]
            if v < best: best = v
            if sub == 0: break
            sub = (sub - 1) & rest
        dp[m] = best
    return tot_e + dp[N - 1]

def gen_graph(n, rng):
    mode = rng.choice(["rand", "clus", "clus", "sparse", "hub"])
    A = np.zeros((n, n), int)
    if mode == "rand":
        p = rng.uniform(.2, .6)
        for i, j in prs(n):
            if rng.random() < p: A[i][j] = A[j][i] = 1
    elif mode == "sparse":
        p = rng.uniform(.1, .3)
        for i, j in prs(n):
            if rng.random() < p: A[i][j] = A[j][i] = 1
    elif mode == "clus":
        k = rng.randint(2, 4); lab = [rng.randrange(k) for _ in range(n)]; q = rng.uniform(.05, .25)
        for i, j in prs(n):
            e = lab[i] == lab[j]
            if rng.random() < q: e = not e
            if e: A[i][j] = A[j][i] = 1
    else:
        h = rng.randrange(n)
        for j in range(n):
            if j != h and rng.random() < .7: A[h][j] = A[j][h] = 1
        for i, j in prs(n):
            if i != h and j != h and rng.random() < .2: A[i][j] = A[j][i] = 1
    return A

def gen_instance(n, rng, partial=True):
    while True:
        A = gen_graph(n, rng)
        if not ind_p3s(A): continue
        ts = ind_p3s(A); start = []
        if partial and rng.random() < .5:
            used = set(); rng.shuffle(ts)
            for t in ts[:rng.randint(0, 2)]:
                if not (p3p(t) & used): start.append(t); used |= p3p(t)
        H = greedy_pm(A, rng, start)
        return A, H

def delta_gamma(A, S):
    """Delta[z][b], Gamma_h[z][b] as sets (both hard)."""
    n = len(A); Dl = {}; Gm = {}
    for z, b in prs(n):
        d = set(); g = set()
        for y in range(n):
            if y in (z, b) or S[z, y] or S[b, y]: continue
            ez, eb = A[z][y], A[b][y]
            if ez != eb: d.add(y)
            elif ez and eb: g.add(y)
        Dl[z, b] = d; Gm[z, b] = g
    return Dl, Gm
