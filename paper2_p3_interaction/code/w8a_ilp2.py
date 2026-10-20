"""Joint ILP 2: choose a packing H* (|H*| >= size) and an OPTIMAL partition S minimising the number of hard edits (a,b), both ends covered by H*,
whose soft-graph distance d_soft(a,b) > r.  (d_soft<=2 <=> T1-local / adjacent members; d_I(H*) <= d_soft-1 <= 2 d_I + 1.)
r=0 gives the OUT-count of w8a_ilp.py."""
import itertools, sys, time, collections, random
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from w8a_lib import *

def joint2(g, opt, size, r, tl=60):
    n = g.n; P = g.p3s()
    pid = {}
    for a in range(n):
        for b in range(a + 1, n): pid[(a, b)] = len(pid)
    m = len(pid); nt = len(P)
    cnt = 0
    def newvars(k):
        nonlocal cnt
        s = cnt; cnt += k; return s
    oy = newvars(m); ox = newvars(nt); ov = newvars(n); ob = newvars(m)
    # reach[k][pair-ordered (a,b)] for k=2..r  (k=1 is the soft indicator, expressed via x)
    cover = {}
    for i, t in enumerate(P):
        for p in mp(t): cover.setdefault(p, []).append(i)
    reach = {}
    for k in range(2, r + 1):
        for a in range(n):
            for b in range(n):
                if a != b: reach[(k, a, b)] = newvars(1)
    # w variables for products: w[k,a,u,b] <= reach[k-1,a,u], <= soft(u,b)
    wv = {}
    for k in range(2, r + 1):
        for a in range(n):
            for b in range(n):
                if a == b: continue
                for u in range(n):
                    if u == a or u == b: continue
                    wv[(k, a, u, b)] = newvars(1)
    N = cnt
    rows = []; lo = []; hi = []
    def add(coefs, l, h):
        rows.append(coefs); lo.append(l); hi.append(h)
    for a, b, d in itertools.combinations(range(n), 3):
        ab, bd, ad = pid[(a, b)], pid[(b, d)], pid[(a, d)]
        add([(oy + ab, 1), (oy + bd, 1), (oy + ad, -1)], -np.inf, 1)
        add([(oy + ab, 1), (oy + ad, 1), (oy + bd, -1)], -np.inf, 1)
        add([(oy + bd, 1), (oy + ad, 1), (oy + ab, -1)], -np.inf, 1)
    add([(oy + i, -1 if p in g.E else 1) for p, i in pid.items()], opt - len(g.E), opt - len(g.E))
    for p, l in cover.items(): add([(ox + i, 1) for i in l], -np.inf, 1)
    add([(ox + i, 1) for i in range(nt)], size, np.inf)
    for i, t in enumerate(P):
        for a in t: add([(ov + a, 1), (ox + i, -1)], 0, np.inf)
    def soft(a, b): return [(ox + j, 1) for j in cover.get(key(a, b), [])]
    # reach definitions: reach[k,a,b] <= sum_u w[k,a,u,b] + reach[k-1,a,b] (monotone) ; w <= reach[k-1,a,u]; w <= soft(u,b)
    for k in range(2, r + 1):
        for a in range(n):
            for b in range(n):
                if a == b: continue
                prev = soft(a, b) if k == 2 else [(reach[(k - 1, a, b)], 1)]
                terms = [(wv[(k, a, u, b)], -1) for u in range(n) if u != a and u != b]
                add([(reach[(k, a, b)], 1)] + terms + [(q, -c) for q, c in prev], -np.inf, 0)
                for u in range(n):
                    if u == a or u == b: continue
                    w = wv[(k, a, u, b)]
                    pr = soft(a, u) if k == 2 else [(reach[(k - 1, a, u)], 1)]
                    add([(w, 1)] + [(q, -c) for q, c in pr], -np.inf, 0)
                    add([(w, 1)] + [(q, -c) for q, c in soft(u, b)], -np.inf, 0)
    for p, i in pid.items():
        a, b = p
        zc = [(oy + i, -1)] if p in g.E else [(oy + i, 1)]
        zconst = 1 if p in g.E else 0
        # bad >= z + v_a + v_b - soft - loc - 2 ; loc = reach[r,a,b] (r>=2) else 0
        coefs = [(ob + i, 1)] + [(q, -c) for q, c in zc] + [(ov + a, -1), (ov + b, -1)] + soft(a, b)
        if r >= 2: coefs.append((reach[(r, a, b)], 1))
        add(coefs, zconst - 2, np.inf)
    A = lil_matrix((len(rows), N))
    for ri, cs in enumerate(rows):
        for q, c in cs: A[ri, q] += c
    c = np.zeros(N); c[ob:ob + m] = 1
    integ = np.ones(N); integ[ob:ob + m] = 0
    ub = np.ones(N)
    res = milp(c, constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=integ, bounds=Bounds(0, 1), options={'time_limit': tl})
    if res.x is None: return None
    return int(round(res.fun))

if __name__ == '__main__':
    # self-test on FAR with H forced? just compute for FAR graphs: r=0..3
    for m in (3, 4):
        n = 2 * m + 1
        g = G8(n, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)])
        opt, _ = opt_ilp(g); nu, _ = nu_mod(g)
        print('FAR', m, [joint2(g, opt, nu, r) for r in (0, 2, 3)])
    # K2 v 3K1 (bad case n=5): r=0,2
    g = G8(5, [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,4)])
    opt,_=opt_ilp(g); nu,_=nu_mod(g); print('K2+3K1', opt, nu, [joint2(g,opt,nu,r) for r in (0,2,3)])
