"""Audit: independent implementation of the (U enumeration + clique-scope CSP bucket elimination) algorithm of
Theorem A / B1 (domain 7, also domain 8), and of the Theorem B search tree."""
import itertools, random, sys
import numpy as np
from audit_lib import *

INF = 10 ** 9

def build(n, E, H):
    Pm = pidx(n)
    bits = []
    for (a, b, c) in H:
        bits.append([1 << Pm[(a, b)], 1 << Pm[(b, c)], 1 << Pm[(a, c)]])
    own = {}
    for i, bs in enumerate(bits):
        for b in bs: own[b] = i
    return Pm, bits, own

def subset_mask(bits_i, s):
    m = 0
    for j in range(3):
        if (s >> j) & 1: m |= bits_i[j]
    return m

def is_p3_triple(n, Ep, tri):
    u, v, w = tri
    P = pidx(n)
    e = [(Ep >> P[(u, v)]) & 1, (Ep >> P[(v, w)]) & 1, (Ep >> P[(u, w)]) & 1]
    return sum(e) == 2

def min_fill_order(nodes, adj):
    adj = {v: set(adj[v]) for v in nodes}
    order = []; width = -1
    rem = set(nodes)
    while rem:
        best = None
        for v in rem:
            nb = adj[v]
            fill = sum(1 for a, b in itertools.combinations(nb, 2) if b not in adj[a])
            key = (fill, len(nb))
            if best is None or key < best[0]: best = (key, v)
        v = best[1]
        nb = adj[v]
        width = max(width, len(nb))
        for a, b in itertools.combinations(nb, 2):
            adj[a].add(b); adj[b].add(a)
        for a in nb: adj[a].discard(v)
        rem.discard(v); order.append(v); del adj[v]
    return order, width

def solve_csp(n, E, H, U, domain, Iadj, order, width):
    """min sum |x_i| s.t. every triple non-P3 in (E xor U) xor x; returns (min, maxscope)."""
    Pm, bits, own = build(n, E, H)
    p = len(H)
    Ep = E ^ U
    dom = list(domain)
    factors = []
    # unary costs
    for i in range(p):
        factors.append(((i,), {(s,): bin(s).count('1') for s in dom}))
    # triple constraints
    for tri in itertools.combinations(range(n), 3):
        u, v, w = tri
        pr = [1 << Pm[(u, v)], 1 << Pm[(v, w)], 1 << Pm[(u, w)]]
        sc = sorted({own[b] for b in pr if b in own})
        # scope must be a clique of I(H)  (Obs 2) -- asserted here
        for a, b in itertools.combinations(sc, 2): assert b in Iadj[a], 'scope not a clique of I(H)'
        tab = {}
        for asg in itertools.product(dom, repeat=len(sc)):
            flip = 0
            for i, s in zip(sc, asg): flip |= subset_mask(bits[i], s)
            Ef = Ep ^ flip
            e = [(Ef >> Pm[(u, v)]) & 1, (Ef >> Pm[(v, w)]) & 1, (Ef >> Pm[(u, w)]) & 1]
            tab[asg] = 0 if sum(e) != 2 else INF
        factors.append((tuple(sc), tab))
    maxscope = 0
    for v in order:
        inv = [f for f in factors if v in f[0]]
        factors = [f for f in factors if v not in f[0]]
        sc = sorted({x for f in inv for x in f[0]} - {v})
        maxscope = max(maxscope, len(sc) + 1)
        tab = {}
        for asg in itertools.product(dom, repeat=len(sc)):
            a = dict(zip(sc, asg)); best = INF
            for val in dom:
                a[v] = val; c = 0
                for f in inv:
                    c += f[1][tuple(a[x] for x in f[0])]
                    if c >= INF: break
                best = min(best, c)
            tab[asg] = best
        factors.append((tuple(sc), tab))
    total = sum(f[1][()] for f in factors)
    return min(total, INF), maxscope

def dp_value(n, E, H, l, domain=range(1, 8)):
    """min over |U|<=l of |U| + CSP min (INF if none).  Also returns max scope seen vs width bound."""
    Pm, bits, own = build(n, E, H)
    p = len(H)
    packed = 0
    for bs in bits:
        for b in bs: packed |= b
    unp = [k for k in range(n * (n - 1) // 2) if not (packed >> k) & 1]
    I = {i: set() for i in range(p)}
    for i in range(p):
        for j in range(i + 1, p):
            if set(H[i]) & set(H[j]): I[i].add(j); I[j].add(i)
    order, width = min_fill_order(list(range(p)), I)
    best = INF; worst_scope = 0; bound_ok = True
    for sz in range(0, l + 1):
        for Ucomb in itertools.combinations(unp, sz):
            U = 0
            for k in Ucomb: U |= 1 << k
            val, ms = solve_csp(n, E, H, U, domain, I, order, width)
            worst_scope = max(worst_scope, ms)
            if ms > width + 1: bound_ok = False
            if val < INF: best = min(best, val + sz)
    return best, width, worst_scope, bound_ok

def thmB_search(n, Ex, packed, budget, stats):
    """Theorem B search tree on current graph Ex; packed pairs frozen. True iff a completion with <= budget flips exists."""
    stats[0] += 1
    p3 = induced_p3s(n, Ex)
    if not p3: return True
    if budget == 0: return False
    a, b, c = p3[0]
    P = pidx(n)
    for (x, y) in ((a, b), (b, c), (a, c)):
        k = P[(x, y)]
        if (packed >> k) & 1: continue
        if thmB_search(n, Ex ^ (1 << k), packed, budget - 1, stats): return True
    return False

def run(seed, trials):
    rng = random.Random(seed)
    cnt = dict(inst=0, A_checked=0, A8_checked=0, B_checked=0, bound_viol=0, maxw=0, A_mismatch=0, B_mismatch=0,
               B_nodes_exceed=0, A_yes=0, A_no=0)
    for t in range(trials):
        n = rng.randint(5, 8)
        E = rand_graph(rng, n, rng.choice(['gnp', 'hub', 'cliques', 'path']))
        H = rand_packing(rng, n, E, rng.choice(['max', 'part', 'overlap']))
        if not H or len(H) > 7: continue
        pk = check_packing(n, E, H)
        p = len(H); cnt['inst'] += 1
        S, cost = opt_and_costs(n, E)
        opt = int(cost.min())
        for l in (0, 1, 2):
            v7, width, ws, ok = dp_value(n, E, H, l, range(1, 8))
            cnt['maxw'] = max(cnt['maxw'], width)
            if not ok: cnt['bound_viol'] += 1
            if opt <= p + l:
                if v7 != opt: cnt['A_mismatch'] += 1; print('A mismatch', n, E, H, opt, v7)
                cnt['A_yes'] += 1
            else:
                if v7 <= p + l: cnt['A_mismatch'] += 1; print('A false accept', n, E, H, opt, v7)
                cnt['A_no'] += 1
            cnt['A_checked'] += 1
        if t % 4 == 0:
            v8, *_ = dp_value(n, E, H, 1, range(0, 8))
            v7b, *_ = dp_value(n, E, H, 1, range(1, 8))
            assert v8 == v7b, 'domain 8 vs 7 differ'
            cnt['A8_checked'] += 1
        # Theorem B: random complete x
        packed = pk
        for _ in range(6):
            xs = [rng.randint(1, 7) for _ in H]
            _, bits, _ = build(n, E, H)
            xm = 0
            for i, s in enumerate(xs): xm |= subset_mask(bits[i], s)
            Sarr = S
            m = (Sarr & np.uint64(packed)) == np.uint64(xm)
            ex = INF
            if m.any(): ex = int(popcount(Sarr[m] & ~np.uint64(packed)).min())
            Ex = E ^ xm
            got = INF
            for b in range(0, 5):
                st = [0]
                if thmB_search(n, Ex, packed, b, st):
                    got = b;
                    if st[0] > (3 ** (b + 1) - 1) // 2 + 0: cnt['B_nodes_exceed'] += 1
                    break
            if ex <= 4 or got < INF:
                if ex != got: cnt['B_mismatch'] += 1; print('B mismatch', ex, got)
            cnt['B_checked'] += 1
    return cnt

if __name__ == '__main__':
    print(run(int(sys.argv[1]), int(sys.argv[2])))
