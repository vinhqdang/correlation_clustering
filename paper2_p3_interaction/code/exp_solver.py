"""Exact solvers for CEaMP on small instances.

  brute_opt(A)        optimum cluster-editing cost by enumeration of all set partitions
  ceamp_brute(A,H,l)  CEaMP decision: OPT <= |H| + l
  dp_l0(A,H)          l = 0 decision via the CSP view (Observation 1/2 of analysis_01):
                      variables x_i in {0,1,2} (which pair of P_i is edited), one
                      constraint per vertex triple; solved by dynamic programming over
                      an elimination ordering (= a tree decomposition whose bags are
                      {v} u N+(v)) of the constraint graph, a subgraph of I(H).
Run as a script: validation against brute force on random small instances.
    python3 exp_solver.py [n_instances]   (seed fixed)
"""
import itertools, random, sys, time
import numpy as np
from exp_common import elimination, greedy_pack, csr, verify_packing

# ------------------------------------------------------------ brute force
_PART = {}
def partition_matrix(n):
    """rows = set partitions (restricted growth strings); cols = vertex pairs;
    entry 1 iff the pair lies in one cluster."""
    if n in _PART:
        return _PART[n]
    parts = []
    def rec(rg, mx):
        if len(rg) == n:
            parts.append(tuple(rg)); return
        for k in range(mx + 2):
            rg.append(k); rec(rg, max(mx, k)); rg.pop()
    rec([0], 0)
    R = np.array(parts)
    iu, ju = np.triu_indices(n, 1)
    M = (R[:, iu] == R[:, ju]).astype(np.int8)
    _PART[n] = (M, iu, ju)
    return _PART[n]

def brute_opt(A):
    n = A.shape[0]
    M, iu, ju = partition_matrix(n)
    adj = A[iu, ju].astype(np.int8)
    return int((M != adj[None, :]).sum(1).min())

def ceamp_brute(A, H, l):
    return brute_opt(A) <= len(H) + l

# ------------------------------------------------------------ CSP / DP
def build_csp(A, H):
    """Return (feasible_trivially, constraints) with constraints = list of (scope, table).
    Pair (u,v) of H[i] has index j: 0=ab, 1=bc, 2=ac.  After choosing x, pair status
    s' = s XOR [x_i == j].  A triple is bad iff exactly two of its pairs are edges.
    Triples with no packed pair are not enumerated here; they are bad iff G has an
    induced P3 outside the packing's pair set, which the caller checks via 'maximal'
    (here: checked explicitly when n is small, see unpacked_p3)."""
    n = A.shape[0]
    owner = {}
    for i, (a, b, c) in enumerate(H):
        owner[(min(a, b), max(a, b))] = (i, 0)
        owner[(min(b, c), max(b, c))] = (i, 1)
        owner[(min(a, c), max(a, c))] = (i, 2)
    # M[u,w] = 1 if the pair can be an edge (edge in G, or packed)
    M = A.astype(bool).copy()
    for (u, v) in owner:
        M[u, v] = M[v, u] = True
    triples = set()
    for (u, v) in owner:
        for w in np.nonzero(M[u] | M[v])[0]:
            if w != u and w != v:
                triples.add(tuple(sorted((u, v, int(w)))))
    cons = {}
    cache = {}
    for t in triples:
        u, v, w = t
        prs = [(u, v), (u, w), (v, w)]
        info = []
        scope = []
        for pr in prs:
            o = owner.get(pr)
            info.append((int(A[pr]), o))
            if o is not None and o[0] not in scope:
                scope.append(o[0])
        scope = tuple(sorted(scope))
        key = (tuple((s, None if o is None else (scope.index(o[0]), o[1])) for s, o in info), len(scope))
        tab = cache.get(key)
        if tab is None:
            tab = np.ones((3,) * len(scope), dtype=bool)
            for x in itertools.product(range(3), repeat=len(scope)):
                ne = 0
                for s, o in info:
                    if o is not None and x[scope.index(o[0])] == o[1]:
                        s ^= 1
                    ne += s
                if ne == 2:
                    tab[x] = False
            cache[key] = tab
        if tab.all():
            continue
        if scope in cons:
            cons[scope] = cons[scope] & tab
        else:
            cons[scope] = tab.copy() if tab.base is not None or True else tab
    return cons

def unpacked_p3(A, H):
    """True iff G has an induced P3 with none of its three pairs packed."""
    n = A.shape[0]
    used = np.zeros((n, n), dtype=bool)
    for a, b, c in H:
        for x, y in ((a, b), (b, c), (a, c)):
            used[x, y] = used[y, x] = True
    for b in range(n):
        nb = np.nonzero(A[b])[0]
        for i, a in enumerate(nb):
            for c in nb[i + 1:]:
                if not A[a, c] and not used[a, b] and not used[b, c] and not used[a, c]:
                    return True
    return False

def constraint_graph(p, cons):
    adj = [set() for _ in range(p)]
    for scope in cons:
        for a in scope:
            for b in scope:
                if a != b:
                    adj[a].add(b)
    return adj

def dp_l0(A, H, max_width=14, check_unpacked=True):
    """Decide l = 0.  Returns (answer, info) where answer is True/False or None if
    the elimination width of the constraint graph exceeds max_width.  info has the
    width used and the number of non-trivial constraints."""
    p = len(H)
    if check_unpacked and unpacked_p3(A, H):
        return False, dict(width=0, ncons=0, reason='unpacked P3')
    cons = build_csp(A, H)
    adj = constraint_graph(p, cons)
    wd, order = elimination(adj, max_width, 'mindeg')
    wf, order_f = elimination(adj, max_width, 'minfill') if p <= 3000 else (None, [])
    if wd is None and wf is None:
        return None, dict(width=None, ncons=len(cons))
    if wf is not None and (wd is None or wf < wd):
        wd, order = wf, order_f
    factors = [(s, t) for s, t in cons.items()]
    pos = {v: k for k, v in enumerate(order)}
    # bucket elimination, boolean tables
    buckets = [[] for _ in range(p)]
    for s, t in factors:
        if s:
            buckets[min(pos[v] for v in s)].append((s, t))
    for k, v in enumerate(order):
        fs = buckets[k]
        if not fs:
            continue
        U = sorted(set().union(*[set(s) for s, _ in fs]))
        tab = np.ones((3,) * len(U), dtype=bool)
        for s, t in fs:
            shape = [3 if u in s else 1 for u in U]
            tab = tab & t.reshape(shape)
        ax = U.index(v)
        res = tab.any(axis=ax)
        rest = tuple(u for u in U if u != v)
        if not res.any():
            return False, dict(width=wd, ncons=len(cons))
        if rest:
            buckets[min(pos[u] for u in rest)].append((rest, res))
    return True, dict(width=wd, ncons=len(cons))

# ------------------------------------------------------------ validation
def random_instance(rng, n):
    kind = rng.choice(['gnp', 'planted', 'planted'])
    A = np.zeros((n, n), dtype=np.uint8)
    if kind == 'gnp':
        p = rng.choice([0.2, 0.35, 0.5, 0.7])
        for u in range(n):
            for v in range(u + 1, n):
                if rng.random() < p:
                    A[u, v] = A[v, u] = 1
    else:
        k = rng.randint(1, max(1, n // 2))
        lab = [rng.randrange(k) for _ in range(n)]
        for u in range(n):
            for v in range(u + 1, n):
                if lab[u] == lab[v]:
                    A[u, v] = A[v, u] = 1
        for _ in range(rng.randint(1, 4)):
            u, v = rng.sample(range(n), 2)
            A[u, v] = A[v, u] = 1 - A[u, v]
    return A

def random_packing(A, rng, maximal):
    n = A.shape[0]
    cand = []
    for b in range(n):
        nb = [x for x in range(n) if A[b, x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a, c]:
                cand.append((a, b, c))
    rng.shuffle(cand)
    used, H = set(), []
    for a, b, c in cand:
        ps = [frozenset(x) for x in ((a, b), (b, c), (a, c))]
        if any(q in used for q in ps):
            continue
        if not maximal and rng.random() < 0.3:
            continue
        used.update(ps); H.append((a, b, c))
    return H

def validate(N, seed=12345):
    rng = random.Random(seed)
    stats = dict(total=0, yes=0, no=0, mism_l0=0, mism_l=0, maximal=0, nonmax=0, tested_l=0)
    t0 = time.time()
    for it in range(N):
        n = rng.randint(4, 9 if it % 5 == 0 else 8)
        A = random_instance(rng, n)
        maximal = rng.random() < 0.7
        H = random_packing(A, rng, maximal)
        if not H and maximal:
            continue
        verify_packing(A, H)
        opt = brute_opt(A)
        truth = opt <= len(H)
        ans, info = dp_l0(A, H, max_width=20)
        stats['total'] += 1; stats['yes' if truth else 'no'] += 1
        stats['maximal' if maximal else 'nonmax'] += 1
        if ans != truth:
            stats['mism_l0'] += 1
            print('MISMATCH', n, A.tolist(), H, opt, ans)
        # brute CEaMP for l in 0..2 must be monotone and agree with opt
        for l in range(3):
            stats['tested_l'] += 1
            if ceamp_brute(A, H, l) != (opt - len(H) <= l):
                stats['mism_l'] += 1
    stats['seconds'] = round(time.time() - t0, 1)
    return stats

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    print(validate(N))
