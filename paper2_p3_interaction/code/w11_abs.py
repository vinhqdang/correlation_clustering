"""w11_abs.py -- abstract states-only chain model (soft-triangle-free, c<=2): path of members, consecutive link relations over D5 (all enumerated link
configurations), unit/2-list locks available at every member; minimum number of locks with h(J) > l via CEGAR.  Tests the weighted Helly / bounded increment
question in the cleanest model."""
import sys, itertools, random, time
import numpy as np
from ortools.sat.python import cp_model
from w11_maj2 import centre_variants
from w11_maj5 import rel5

COST = [0, 0, 0, 1, 2]
LISTS = [(0,), (1,), (2,), (0, 3), (1, 3), (3, 4), (4,)] if False else [(0,), (1,), (2,), (0, 3), (1, 3)]

_link_cache = {}
def link_matrices():
    """all 5x5 boolean relations of a link (A,B sharing q) over D5, one per (centre choice A, centre B, cross bits): returns list"""
    if 'm' in _link_cache: return _link_cache['m']
    out = []; seen = set()
    cross = [(1,3),(1,4),(2,3),(2,4)]
    for ma in centre_variants([0,1,2]):
        for mb in centre_variants([0,3,4]):
            for bits in itertools.product((0,1), repeat=4):
                hard = {tuple(sorted(p)): b for p, b in zip(cross, bits)}
                R = rel5([ma, mb], 5, hard, ())
                M = np.zeros((5,5), dtype=bool)
                for a, b in R.tolist(): M[a, b] = True
                key = M.tobytes()
                if key in seen: continue
                seen.add(key); out.append(M)
    _link_cache['m'] = out
    return out

def hval(mats, lockmasks):
    """min cost over assignments of members (D5) satisfying link relations mats[i] between i,i+1 and unary masks (list of 5-bool arrays); inf if none"""
    N = len(mats) + 1
    INF = 10**6
    f = np.array([COST[s] if lockmasks[0][s] else INF for s in range(5)], dtype=np.int64)
    for i in range(N - 1):
        g = np.full(5, INF, dtype=np.int64)
        for t in range(5):
            if not lockmasks[i + 1][t]: continue
            best = INF
            for s in range(5):
                if mats[i][s, t] and f[s] < INF: best = min(best, f[s])
            if best < INF: g[t] = best + COST[t]
        f = g
    return int(f.min())

def hassign(mats, lockmasks, bound):
    """one assignment with cost<=bound (or None)"""
    N = len(mats) + 1
    INF = 10**6
    F = [None] * N; P = [None] * N
    f = np.array([COST[s] if lockmasks[0][s] else INF for s in range(5)], dtype=np.int64); F[0] = f
    for i in range(N - 1):
        g = np.full(5, INF, dtype=np.int64); p = [-1] * 5
        for t in range(5):
            if not lockmasks[i + 1][t]: continue
            for s in range(5):
                if mats[i][s, t] and f[s] < INF and f[s] + COST[t] < g[t]:
                    g[t] = f[s] + COST[t]; p[t] = s
        f = g; F[i + 1] = f; P[i + 1] = p
    t = int(np.argmin(f))
    if f[t] > bound: return None
    x = [t]
    for i in range(N - 1, 0, -1):
        t = P[i][t]; x.append(t)
    return x[::-1]

def min_lock_set(mats, l, maxsize=12, rnd=None):
    """minimum number of locks (member, list) with h > l ; None if even all locks give h <= l"""
    N = len(mats) + 1
    full = [np.ones(5, dtype=bool) for _ in range(N)]
    allmask = []
    pool = [(i, lst) for i in range(N) for lst in LISTS]
    def masks_of(sel):
        ms = [np.ones(5, dtype=bool) for _ in range(N)]
        for (i, lst) in sel:
            m = np.zeros(5, dtype=bool); m[list(lst)] = True
            ms[i] = ms[i] & m
        return ms
    if hval(mats, masks_of(pool)) <= l: return None
    cons = []; sel = set()
    for it in range(2000):
        ms = masks_of(sel)
        x = hassign(mats, ms, l)
        if x is None: return len(sel), sorted(sel)
        # assignment x (cost<=l) satisfies sel; require a lock in pool violated by x
        cons.append([k for k, (i, lst) in enumerate(pool) if x[i] not in lst])
        m = cp_model.CpModel(); v = [m.NewBoolVar('v%d' % k) for k in range(len(pool))]
        for c in cons: m.AddBoolOr([v[k] for k in c])
        m.Minimize(sum(v))
        s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 60
        st = s.Solve(m)
        if st != cp_model.OPTIMAL: raise RuntimeError('master')
        sel = set(pool[k] for k in range(len(pool)) if s.Value(v[k]))
        if len(sel) > maxsize: return len(sel), None
    raise RuntimeError('no convergence')

if __name__ == '__main__':
    N = int(sys.argv[1]); l = int(sys.argv[2]); samples = int(sys.argv[3]); seed = int(sys.argv[4])
    rnd = random.Random(seed)
    mats_pool = link_matrices(); print('distinct link relations over D5:', len(mats_pool))
    from collections import Counter
    res = Counter(); best = []
    t0 = time.time()
    for s in range(samples):
        while True:
            mats = [rnd.choice(mats_pool) for _ in range(N - 1)]
            if hval(mats, [np.ones(5, dtype=bool) for _ in range(N)]) == 0: break
        r = min_lock_set(mats, l)
        k = None if r is None else r[0]
        res[k] += 1
        if k is not None and k >= 4: best.append((k, r[1]))
    print('N', N, 'l', l, sorted(res.items(), key=str), round(time.time() - t0), 's')
    best.sort(key=lambda x: -x[0])
    for b in best[:3]: print(b)
