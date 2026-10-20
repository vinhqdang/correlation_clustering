"""Theorem PL check (semantics, not the DP itself).
Given (G,H) and a tree decomposition (T,beta) of I(H) (networkx min-fill), define
  X_t = union of V(P), P in beta_t ;  a,b co-bagged iff some X_t contains both.
  U_loc = hard pairs {a,b} with a,b in V(H) co-bagged, or {r,v}, r outside V(H), v in V(H), such that every partner candidate of r
          (the vertices of V(H) within distance 2 of r in G) is pairwise co-bagged  ('local hub').
CLASS-OPT := min over U subset U_loc, |U|<=ell, of the cheapest solution whose hard edits are exactly U  (ILP).
Promise holds iff CLASS-OPT = OPT.  Also: 'visibility lemma' : every vertex triple can be verified from bag-local data (labels of in-bag vertices + owner states)."""
import itertools, sys, random, collections, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from networkx.algorithms.approximation import treewidth_min_fill_in
from w8a_lib import *

def decomposition(H):
    I = Igraph(H)
    if I.number_of_edges() == 0:
        return [frozenset([i]) for i in range(len(H))] or [frozenset()]
    w, T = treewidth_min_fill_in(I)
    return [frozenset(b) for b in T.nodes] + [frozenset([i]) for i in range(len(H))]

def Xsets(H, bags):
    return [frozenset(v for i in b for v in H[i]) for b in bags]

def cobagged(X, a, b):
    return any(a in x and b in x for x in X)

def uloc(g, H, bags):
    X = Xsets(H, bags); n = g.n
    soft = soft_pairs(H); VH = set(v for t in H for v in t)
    cand = []
    for a in range(n):
        for b in range(a + 1, n):
            p = (a, b)
            if p in soft: continue
            if a in VH and b in VH and cobagged(X, a, b): cand.append(p)
    # local hubs
    for r in range(n):
        if r in VH: continue
        ball = set()
        for x in g.adj[r]:
            ball.add(x); ball |= g.adj[x]
        ball = (ball - {r}) & VH
        if all(cobagged(X, u, v) for u, v in itertools.combinations(sorted(ball), 2)) and all(any(u in x for x in X) for u in ball):
            for v in ball:
                p = key(r, v)
                if p not in soft: cand.append(p)
    return sorted(set(cand))

def best_with_U(g, H, U):
    """min cost of a solution whose hard edits are exactly U (hard = not in soft pairs); None if infeasible"""
    n = g.n; soft = soft_pairs(H)
    pid = {}
    for a in range(n):
        for b in range(a + 1, n): pid[(a, b)] = len(pid)
    m = len(pid)
    rows = []
    for a, b, d in itertools.combinations(range(n), 3):
        ab, bd, ad = pid[(a, b)], pid[(b, d)], pid[(a, d)]
        for (x, y, z) in ((ab, bd, ad), (ab, ad, bd), (bd, ad, ab)):
            r = np.zeros(m); r[x] = 1; r[y] = 1; r[z] = -1; rows.append(r)
    lo = np.zeros(m); hi = np.ones(m)
    Uset = set(U)
    for p, i in pid.items():
        if p in soft: continue
        tog = (p in g.E) != (p in Uset)      # together after editing
        lo[i] = hi[i] = 1 if tog else 0
    c = np.array([-1.0 if p in g.E else 1.0 for p in pid])
    res = milp(c, constraints=LinearConstraint(np.array(rows), -np.inf, 1), integrality=np.ones(m), bounds=Bounds(lo, hi))
    if res.x is None: return None
    return int(round(res.fun)) + len(g.E)

def class_opt(g, H, ell, bags=None):
    bags = decomposition(H) if bags is None else bags
    cand = uloc(g, H, bags)
    best = None; bestU = None
    for k in range(0, ell + 1):
        for U in itertools.combinations(cand, k):
            v = best_with_U(g, H, U)
            if v is not None and (best is None or v < best): best, bestU = v, U
    return best, bestU, cand

def visibility_check(g, H, bags, rng, trials=200):
    """random S in the class (random soft states + random U subset of U_loc); local verification vs global cluster-graph test"""
    X = Xsets(H, bags); cand = uloc(g, H, bags); soft = soft_pairs(H); n = g.n
    VH = set(v for t in H for v in t)
    ballH = {}
    for r in range(n):
        if r in VH: continue
        bb = set()
        for x in g.adj[r]: bb.add(x); bb |= g.adj[x]
        ballH[r] = (bb - {r}) & VH
    def bag_of_members(ms):
        for bi, b in enumerate(bags):
            if set(ms) <= set(b): return bi
        return None
    bad = 0
    for _ in range(trials):
        U = set(p for p in cand if rng.random() < 0.15)
        S = set(U)
        for i, t in enumerate(H):
            S |= set(p for p in mp(t) if rng.random() < 0.4) or {rng.choice(mp(t))}
        # global validity
        adj = [[(key(a, b) in g.E) != (key(a, b) in S) if a != b else False for b in range(n)] for a in range(n)]
        glob = all(not (adj[a][b] and adj[b][c] and not adj[a][c]) and not (adj[a][b] and adj[a][c] and not adj[b][c]) and not (adj[a][c] and adj[b][c] and not adj[a][b])
                   for a, b, c in itertools.combinations(range(n), 3))
        # local verification: statuses from bag-visible data only
        labels = {v: set() for v in range(n)}                       # E_v for covered v
        for (a, b) in U:
            if a in VH and b in VH: labels[a].add(b); labels[b].add(a)
            elif a in VH: labels[a].add(b)
            else: labels[b].add(a)
        ok = True
        for tr in itertools.combinations(range(n), 3):
            prs = [key(tr[0], tr[1]), key(tr[1], tr[2]), key(tr[0], tr[2])]
            owners = set(soft[p] for p in prs if p in soft)
            if owners:
                bi = bag_of_members(owners)
                assert bi is not None, 'Helly failure'
                vis = X[bi]
            else:
                us = [p for p in prs if p in U]
                if not us: continue                       # unedited all-hard triple: fine by pair-maximality (checked separately)
                a, b = us[0]
                hub = a if a not in VH else (b if b not in VH else None)
                need = set([x for x in (a, b) if x in VH])
                if hub is not None: need |= ballH[hub]
                vis = None
                for xs in X:
                    if need <= xs: vis = xs; break
                if vis is None: ok = False; continue
            st = []
            for p in prs:
                if p in soft: s_ = (p in S)
                else:
                    x, y = p
                    if x not in VH and y not in ballH.get(x, set()): e = False        # static: no edit possible at (r,c), c outside the hub's candidate ball
                    elif y not in VH and x not in ballH.get(y, set()): e = False
                    else:
                        seen = (x in vis and x in VH) or (y in vis and y in VH)
                        assert seen, 'invisible hard pair'
                        e = (y in labels[x] if x in VH else False) or (x in labels[y] if y in VH else False)
                    s_ = e
                st.append(((p in g.E) != s_))
            if sum(st) == 2: ok = False
        if ok != glob: bad += 1
    return bad

if __name__ == '__main__':
    rng = random.Random(5)
    for m in (3, 4, 5, 6):
        n = 2 * m + 1
        g = G8(n, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)])
        opt, _ = opt_ilp(g)
        H0 = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
        Hs = H0[:m - 1] + [(1, 2 * m - 1, 2 * m - 2)]
        for name, H in (('H=path tiling', H0), ('H*=cycle repack', Hs)):
            co, U, cand = class_opt(g, H, 1)
            print('FAR m=%d %-16s OPT=%d CLASS-OPT=%s  best U=%s  |U_loc|=%d  promise holds: %s' % (m, name, opt, co, U, len(cand), co == opt))
    # visibility lemma on random instances
    tot = bad = 0
    for it in range(40):
        n = rng.randint(7, 10); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < 0.3]
        g = G8(n, E)
        if not g.p3s(): continue
        H = []; used = set(); P = g.p3s()[:]; rng.shuffle(P)
        for t in P:
            if not any(p in used for p in mp(t)): H.append(t); used.update(mp(t))
        assert pair_maximal(g, H)
        bags = decomposition(H)
        bad += visibility_check(g, H, bags, rng, 60); tot += 60
    print('visibility lemma: %d mismatches in %d random (S in class) checks' % (bad, tot))

def visibility_on_optima(trials=60, seed=9):
    """for optimal solutions S whose hard edits lie in U_loc: the local verification must ACCEPT (S is globally valid); count how many were tested"""
    rng = random.Random(seed); tested = bad = 0
    for it in range(trials):
        n = rng.randint(7, 10); E = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < 0.3]
        g = G8(n, E)
        if not g.p3s(): continue
        H = []; used = set(); P = g.p3s()[:]; rng.shuffle(P)
        for t in P:
            if not any(p in used for p in mp(t)): H.append(t); used.update(mp(t))
        bags = decomposition(H); cand = set(uloc(g, H, bags))
        opt, sols = g.opt(limit=200)
        for s in sols:
            S = edits_of_partition(g, part_from_mask(s, n))
            U = set(hard_edits(H, S))
            if not U <= cand: continue
            tested += 1
            # reuse local check with this exact S: temporarily emulate by calling internal logic
            X = Xsets(H, bags); soft = soft_pairs(H); VH = set(v for t in H for v in t)
            ballH = {}
            for r in range(n):
                if r in VH: continue
                bb = set()
                for x in g.adj[r]: bb.add(x); bb |= g.adj[x]
                ballH[r] = (bb - {r}) & VH
            labels = {v: set() for v in range(n)}
            for (a, b) in U:
                if a in VH and b in VH: labels[a].add(b); labels[b].add(a)
                elif a in VH: labels[a].add(b)
                else: labels[b].add(a)
            ok = True
            for tr in itertools.combinations(range(n), 3):
                prs = [key(tr[0], tr[1]), key(tr[1], tr[2]), key(tr[0], tr[2])]
                owners = set(soft[p] for p in prs if p in soft)
                if owners:
                    bi = next((i for i, b in enumerate(bags) if owners <= set(b)), None)
                    assert bi is not None; vis = X[bi]
                else:
                    us = [p for p in prs if p in U]
                    if not us: continue
                    a, b = us[0]; hub = a if a not in VH else (b if b not in VH else None)
                    need = set(x for x in (a, b) if x in VH)
                    if hub is not None: need |= ballH[hub]
                    vis = next((xs for xs in X if need <= xs), None)
                    if vis is None: ok = False; continue
                st = []
                for p in prs:
                    if p in soft: e = p in S
                    else:
                        x, y = p
                        if x not in VH and y not in ballH.get(x, set()): e = False
                        elif y not in VH and x not in ballH.get(y, set()): e = False
                        else:
                            assert (x in vis and x in VH) or (y in vis and y in VH)
                            e = (y in labels[x] if x in VH else False) or (x in labels[y] if y in VH else False)
                    st.append((p in g.E) != e)
                if sum(st) == 2: ok = False
            if not ok: bad += 1
    return tested, bad
if __name__ == '__main__' and True:
    print('optimal solutions in the class: (tested, rejected by local verification) =', visibility_on_optima())
