"""Audit 4 library (written from scratch; imports nothing from the authors' or earlier auditors' scripts).
Graph = n x n 0/1 numpy matrix.  P3 = (a,b,c): b centre, ab,bc edges, ac non-edge, a<c.
Packed pair = two vertices lying in a common member of H ('soft').  Every other pair = 'hard' (edge or non-edge).
"""
import itertools, random
import numpy as np

# ---------------------------------------------------------------- basics
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
def packed_matrix(n, H):
    pk = np.zeros((n, n), bool)
    for t in H:
        for x, y in itertools.combinations(t, 2): pk[x, y] = pk[y, x] = True
    return pk
def is_packing(A, H):
    seen = set()
    for t in H:
        a, b, c = t
        if not (A[a][b] and A[b][c] and not A[a][c]): return False
        ps = p3_pairs(t)
        if ps & seen: return False
        seen |= ps
    return True
def used_pairs(H):
    s = set()
    for t in H: s |= p3_pairs(t)
    return s
def greedy_extend(A, H, rng, prefer_overlap=False):
    """add pair-disjoint induced P3s until none is left (pair-maximal)"""
    H = list(H); used = used_pairs(H); ts = p3_list(A); rng.shuffle(ts)
    if prefer_overlap:
        cnt = {}
        for t in H:
            for v in t: cnt[v] = cnt.get(v, 0) + 1
        ts.sort(key=lambda t: -sum(cnt.get(v, 0) for v in t) + rng.random())
    for t in ts:
        ps = p3_pairs(t)
        if not (ps & used):
            H.append(t); used |= ps
    return H
def rand_packing(A, rng, stop=0.0, prefer_overlap=False):
    ts = p3_list(A); rng.shuffle(ts); H = []; used = set(); cnt = {}
    if prefer_overlap: ts.sort(key=lambda t: -sum(cnt.get(v, 0) for v in t) + rng.random())
    for t in ts:
        if rng.random() < stop: continue
        ps = p3_pairs(t)
        if not (ps & used):
            H.append(t); used |= ps
            for v in t: cnt[v] = cnt.get(v, 0) + 1
    return H
def is_pair_maximal(A, H):
    used = used_pairs(H)
    return all(p3_pairs(t) & used for t in p3_list(A))
def cmax(n, H):
    cnt = [0] * n
    for t in H:
        for v in t: cnt[v] += 1
    return max(cnt) if cnt else 0
def soft_degree(pk): return int(pk.sum(1).max()) if len(pk) else 0

# ---------------------------------------------------------------- treewidth (exact subset DP, <= ~15 vertices)
def tw_exact(adjm):
    n = len(adjm)
    if n == 0: return -1
    arr = [0] * (1 << n)
    def Q(S, v):
        reach = 1 << v; frontier = reach; out = 0
        while frontier:
            nxt = 0; f = frontier
            while f:
                b = f & -f; f ^= b; nxt |= adjm[b.bit_length() - 1]
            nxt &= ~reach; out |= nxt & ~S; nxt &= S; reach |= nxt; frontier = nxt
        return bin(out).count('1')
    for S in range(1, 1 << n):
        best = 10 ** 9; f = S
        while f:
            b = f & -f; f ^= b; v = b.bit_length() - 1; T = S ^ b
            val = max(arr[T], Q(T, v))
            if val < best: best = val
        arr[S] = best
    return arr[(1 << n) - 1]
def ih_adj(H):
    m = len(H); adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if set(H[i]) & set(H[j]): adj[i] |= 1 << j; adj[j] |= 1 << i
    return adj
def soft_adj(n, H):
    pk = packed_matrix(n, H); return [sum(1 << j for j in range(n) if pk[i, j]) for i in range(n)]

# ---------------------------------------------------------------- partitions (n <= 9)
_PT = {}
def part_table(n):
    if n in _PT: return _PT[n]
    res = []
    def rec(i, lab, m):
        if i == n: res.append(tuple(lab)); return
        for c in range(m + 1):
            lab.append(c); rec(i + 1, lab, max(m, c + 1)); lab.pop()
    rec(0, [], 0)
    lab = np.array(res, dtype=np.int8)
    P = prs(n)
    same = np.stack([lab[:, i] == lab[:, j] for i, j in P], 1)
    _PT[n] = (lab, same, P); return _PT[n]
def all_costs(A):
    n = len(A); lab, same, P = part_table(n)
    e = np.array([A[i][j] for i, j in P], bool)
    return lab, same, (same ^ e).sum(1)

# ---------------------------------------------------------------- piece structure of a solution
def pieces_of(A, pk, same_fn):
    """same_fn(i,j) -> bool, i.e. i,j in same cluster.  returns (piece id list, U, F, I)"""
    n = len(A); par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    U = []; F = []; I = []
    for i, j in prs(n):
        s = same_fn(i, j); e = bool(A[i][j])
        if s != e and not pk[i, j]:
            U.append((i, j)); (F if e else I).append((i, j))
        if e and s and not pk[i, j]: par[f(i)] = f(j)       # kept hard edge
    return [f(i) for i in range(n)], U, F, I

def sm_check(A, pk, lab, Ds):
    """check SM (i),(ii),(iii) (+ my sharper bound) on partition lab (labels).  returns (violations, stats)"""
    n = len(A)
    pid, U, F, I = pieces_of(A, pk, lambda i, j: lab[i] == lab[j])
    Iset = set(I); u = len(U); viol = []; stat = {'multi': 0, 'maxC': 0, 'u': u}
    for c in set(lab):
        C = [v for v in range(n) if lab[v] == c]
        ps = {}
        for v in C: ps.setdefault(pid[v], []).append(v)
        if len(ps) < 2: continue
        stat['multi'] += 1; stat['maxC'] = max(stat['maxC'], len(C))
        plist = list(ps.values())
        for A_, B_ in itertools.combinations(plist, 2):
            ins = 0
            for a in A_:
                for b in B_:
                    if pk[a, b]: continue
                    if key(a, b) in Iset: ins += 1
                    else: viol.append(('i', a, b))
            mn, mx = min(len(A_), len(B_)), max(len(A_), len(B_))
            if mn * (mx - Ds) > ins or ins > u: viol.append(('ii-count', len(A_), len(B_), ins))
            if mx > Ds + u: viol.append(('ii-size', mx, Ds, u))
        l_ = u
        if len(C) > max(2 * Ds + l_, 2 * l_): viol.append(('iii', len(C), Ds, l_))
        if len(C) > max(Ds + l_ + 1, 2 * Ds, 2 * l_): viol.append(('iii-sharp', len(C), Ds, l_))
    return viol, stat

# ---------------------------------------------------------------- HUB
def max_indep(nodes, bad):
    """max independent set; bad(x,y) true if x,y may not both be chosen"""
    nodes = list(nodes); best = 0
    def rec(i, cur):
        nonlocal best
        if len(cur) + len(nodes) - i <= best: return
        if i == len(nodes): best = max(best, len(cur)); return
        x = nodes[i]
        if all(not bad(x, y) for y in cur): rec(i + 1, cur + [x])
        rec(i + 1, cur)
    rec(0, []); return best
def hub_stats(A, pk, H):
    """for each u: (|X| max hard-independent set, alpha(G[N_h(u)]), |N_h|)"""
    n = len(A); out = []
    for u in range(n):
        Nh = [x for x in range(n) if A[u][x] and not pk[u, x]]
        X = max_indep(Nh, lambda x, y: A[x][y] or pk[x, y])
        al = max_indep(Nh, lambda x, y: A[x][y])
        out.append((X, al, len(Nh)))
    return out

# ---------------------------------------------------------------- MAX
def rand_M(A, H, rng):
    """random maximal set of induced P3s, pairs pairwise disjoint and disjoint from pairs(H)"""
    used = used_pairs(H); ts = p3_list(A); rng.shuffle(ts); M = []
    for t in ts:
        ps = p3_pairs(t)
        if not (ps & used): M.append(t); used |= ps
    return M

# ---------------------------------------------------------------- generators
def gen_graph(rng, n, fam=None):
    fam = fam or rng.choice(['gnp', 'noisy', 'hub', 'sparse', 'fano'])
    if fam == 'fano' and n < 7: fam = 'gnp'
    A = np.zeros((n, n), np.int8)
    def add(i, j, v=1): A[i, j] = A[j, i] = v
    if fam == 'gnp':
        p = rng.uniform(.2, .7)
        for i, j in prs(n):
            if rng.random() < p: add(i, j)
    elif fam == 'sparse':
        p = rng.uniform(.1, .3)
        for i, j in prs(n):
            if rng.random() < p: add(i, j)
    elif fam == 'noisy':
        k = rng.randint(2, 3); g = [rng.randrange(k) for _ in range(n)]
        for i, j in prs(n):
            if rng.random() < (.85 if g[i] == g[j] else .08): add(i, j)
    elif fam == 'hub':
        hubs = rng.sample(range(n), rng.randint(1, 2))
        for h in hubs:
            for x in range(n):
                if x != h and rng.random() < .75: add(h, x)
        for i, j in prs(n):
            if i not in hubs and j not in hubs and rng.random() < .18: add(i, j)
    elif fam == 'fano':
        # soft-clique gadget: lines of Fano plane (7) or AG(2,3) (9), each line gets exactly one non-edge
        if n >= 9:
            pts = [(x, y) for x in range(3) for y in range(3)]; lines = set()
            for p, q in itertools.combinations(range(9), 2):
                a, b = pts[p], pts[q]; dx, dy = (b[0] - a[0]) % 3, (b[1] - a[1]) % 3
                lines.add(frozenset(i for i, r in enumerate(pts) if ((r[0] - a[0]) * dy - (r[1] - a[1]) * dx) % 3 == 0))
            m = 9
        else:
            lines = [frozenset(((i + d) % 7) for d in (0, 1, 3)) for i in range(7)]; m = 7
        for L in lines:
            L = sorted(L); skip = rng.choice([(0, 1), (0, 2), (1, 2)])
            for x, y in itertools.combinations(range(3), 2):
                if (x, y) != skip: add(L[x], L[y])
        for v in range(m, n):
            for x in range(n):
                if x != v and rng.random() < .4: add(v, x)
    return A
def noise(A, rng, k):
    A = A.copy(); n = len(A)
    for _ in range(k):
        i, j = rng.sample(range(n), 2); A[i, j] = A[j, i] = 1 - A[i, j]
    return A

# ---------------------------------------------------------------- ILP (HiGHS via scipy)
_TRI = {}
def _tri(n):
    if n in _TRI: return _TRI[n]
    from scipy.sparse import lil_matrix
    P = prs(n); idx = {p: k for k, p in enumerate(P)}; rows = []
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = idx[(i, j)], idx[(j, k)], idx[(i, k)]
        rows += [(a, b, c), (a, c, b), (b, c, a)]
    M = lil_matrix((len(rows), len(P)))
    for r, (a, b, c) in enumerate(rows): M[r, a] = 1; M[r, b] = 1; M[r, c] = -1
    _TRI[n] = (M.tocsr(), idx); return _TRI[n]
def ilp(A, obj=None, cost_le=None, fixed=None):
    """z_p = 1 iff pair p in the same cluster.  minimise obj.z (default: cost)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = len(A); P = prs(n); T, idx = _tri(n); m = len(P)
    e = np.array([A[i][j] for i, j in P], float); cvec = 1 - 2 * e; const = e.sum()
    cons = [LinearConstraint(T, -np.inf, 1)]
    if cost_le is not None: cons.append(LinearConstraint(cvec.reshape(1, -1), -np.inf, cost_le - const))
    lb = np.zeros(m); ub = np.ones(m)
    for p, v in (fixed or {}).items(): lb[idx[key(*p)]] = ub[idx[key(*p)]] = v
    r = milp(cvec if obj is None else obj, constraints=cons, integrality=np.ones(m), bounds=Bounds(lb, ub))
    if r.status != 0: return None
    z = np.round(r.x).astype(int); cost = int(round(const + cvec @ z))
    return z, cost, idx
def labels_from_z(n, z, idx):
    lab = [-1] * n; c = 0
    for v in range(n):
        if lab[v] < 0:
            lab[v] = c
            for u in range(v + 1, n):
                if z[idx[(v, u)]]: lab[u] = c
            c += 1
    return lab
def opt_ilp(A):
    r = ilp(A); return r[1]
