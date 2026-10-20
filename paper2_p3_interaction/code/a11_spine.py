"""Audit 11: spine family Sp(F) from paper/sections/app_spine.tex, built from the paper text (not from w6/w8 code);
checks Prop 3.2: H' = H - X_0 + Y + Z is a packing of size |H|+1, OPT = |H|+1, tw(I(H')) <= 2, plus local-hub claim."""
import itertools, sys
from a11_lib import *

def orient(F, nv):
    """orientation of cubic graph F with indegree <= 2: brute force"""
    for bits in itertools.product([0, 1], repeat=len(F)):
        ind = [0] * nv; arcs = []
        for (u, v), b in zip(F, bits):
            t, h = (u, v) if b == 0 else (v, u)
            ind[h] += 1; arcs.append((t, h))
        if max(ind) <= 2: return arcs
    raise Exception

def build(F, nv):
    arcs = orient(F, nv)
    star = nv
    t0 = arcs[0][0]
    e0 = (t0, star)
    rest = arcs[:]
    # ordering search: e0 first, consecutive heads differ, head(e1) != t0
    def rec(order, left):
        if not left: return order
        for i, a in enumerate(left):
            if a[1] == order[-1][1]: continue
            if len(order) == 1 and a[1] == t0: continue
            r = rec(order + [a], left[:i] + left[i + 1:])
            if r: return r
        return None
    order = rec([e0], rest)
    N = len(order) - 1
    idx = lambda L, i: 4 * i + L            # a=0,b=1,c=2,u=3
    z = 4 * (N + 1); y = z + 1; n = y + 1
    E = set()
    Kp = {}
    for i, (t, h) in enumerate(order):
        Kp.setdefault(t, set()).add(idx(0, i))
        for L in (1, 2, 3): Kp.setdefault(h, set()).add(idx(L, i))
    for v, S in Kp.items():
        for a, b in itertools.combinations(sorted(S), 2): E.add(K(a, b))
    for i in range(N + 1): E.add(K(idx(0, i), idx(1, i)))
    for i in range(N): E.add(K(idx(1, i), idx(3, i + 1)))
    E.add(K(idx(1, N), z))
    E.add(K(y, idx(3, 0))); E.add(K(y, idx(2, 0)))
    g = Gr(n, E)
    P = [(idx(0, i), idx(1, i), idx(2, i)) for i in range(N + 1)]
    X = [(idx(3, i), idx(1, i), idx(3, i + 1) if i < N else z) for i in range(N + 1)]
    H = []
    for i in range(N + 1): H += [P[i], X[i]]
    return g, H, dict(N=N, idx=idx, z=z, y=y, P=P, X=X, order=order)

for name, F, nv in (('triangle', [(0, 1), (1, 2), (2, 0)], 3), ('K4', [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], 4),
                    ('K33', [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)], 6)):
    g, H, d = build(F, nv)
    N, idx, y = d['N'], d['idx'], d['y']
    print(f'== Sp({name}): n={g.n} N={N} |H|={len(H)} valid={is_packing(g, H)} pairmax={pair_maximal(g, H)} tw(I(H))={tw_exact(I_graph(H)) if len(H)<=16 else "n/a"}')
    # explicit partition of paper: cliques K+_v, z alone, y into K+_*  -> cost
    # rebuild cliques as components of 'clique edges': use partition from edit count directly
    # cost via ILP check
    opt, x = opt_ilp(g) if g.n <= 31 else (None, None)
    print('   OPT(ILP) =', opt, ' |H|+1 =', len(H) + 1)
    # H' = H - X_0 + Y + Z
    u0, b0, u1 = d['X'][0]
    Y = (y, u0, b0)
    Z = (b0, u1, idx(2, 1))
    Hp = [t for t in H if t != d['X'][0]] + [Y, Z]
    print('   Y,Z induced P3 and H\' modification-disjoint packing:', is_packing(g, Hp), ' |H\'| =', len(Hp), ' tw(I(H\')) =', tw_exact(I_graph(Hp)) if len(Hp) <= 17 else 'n/a')
    # local hub claim: Ball(y) cap V(H) within X_0 u P_0 which share b_0
    ball = set(); fr = {y}
    for _ in range(2):
        fr = {w for v in fr for w in g.adj[v]} - ball - {y}; ball |= fr
    print('   Ball(y,2) =', sorted(ball), ' = {u0,c0,b0}:', ball == {u0, idx(2, 0), b0}, ' X0&P0 share b0:', set(d['X'][0]) & set(d['P'][0]) == {b0})
    # hard degrees <= 7
    own = soft_owner(H)
    hd = [sum(1 for w in g.adj[v] if K(v, w) not in own) for v in range(g.n)]
    print('   max hard degree', max(hd))
