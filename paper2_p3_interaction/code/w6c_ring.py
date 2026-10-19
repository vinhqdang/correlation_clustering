"""Rotation-symmetric ring search for a lock-core of size N.
Ring of N periods; period j has vertices s_j, t_j, p_j, x_j (x_j = clean anchor).  Members: A_j on {x_j, s_j, t_j} (anchor member),
F_j on {t_{j-1}, p_j, s_j} (free member).  All hard pairs not involving anchors are random but rotation invariant.
Target: OPT = |H| (tight), every set of N-1 anchors can be singletons in some optimum, the full set cannot."""
import sys, itertools, random, numpy as np
from w6c_lib import *
N = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed)
def vid(kind, j):
    return {'s': 0, 't': 1, 'p': 2, 'x': 3}[kind] * N + (j % N)
n = 4 * N
def build(aorder, forder, orbit_bits, orbit_list):
    # member triples
    members = []
    for j in range(N):
        s, t, p, x = vid('s', j), vid('t', j), vid('p', j), vid('x', j)
        tprev = vid('t', j - 1)
        A = [None] * 3; verts = [x, s, t]
        # aorder: permutation of verts, x must be at end
        A = tuple(verts[i] for i in aorder)
        F = tuple([tprev, p, vid('s', j)][i] for i in forder)
        members += [A, F]
    return members
def orbits():
    # orbits of pairs among non-anchor vertices under rotation
    seen = set(); orbs = []
    non = [v for v in range(n) if v < 3 * N]
    for u, v in itertools.combinations(non, 2):
        if key(u, v) in seen: continue
        orb = set()
        for r in range(N):
            def rot(a): return (a // N) * N + ((a % N + r) % N)
            orb.add(key(rot(u), rot(v)))
        seen |= orb; orbs.append(sorted(orb))
    return orbs
ORB = orbits()
t0 = __import__('time').time(); found = 0; tested = 0; tight = 0
import time
while time.time() - t0 < T:
    aorder = rng.choice([(0, 1, 2), (0, 2, 1), (1, 2, 0)[::-1]])
    # x at an end: orders of [x,s,t] indices with index0 (x) at position 0 or 2
    aorder = rng.choice([(0, 1, 2), (0, 2, 1), (2, 1, 0), (1, 2, 0)])
    forder = rng.choice(list(itertools.permutations(range(3))))
    members = build(aorder, forder, None, None)
    ok, _ = check_packing(np.zeros((n, n), np.int8), [])  # dummy
    edges = set(); soft = set()
    for (a, b, c) in members:
        edges |= {key(a, b), key(b, c)}; soft |= pairs_of((a, b, c))
    # random orbit bits for hard non-anchor pairs
    E = set(edges)
    for orb in ORB:
        if any(p in soft for p in orb): continue
        if rng.random() < rng.choice([0.2, 0.5]):
            E |= set(orb)
    A = adj_from_edges(n, E)
    ok, why = check_packing(A, members)
    if not ok: continue
    tested += 1
    opt, _ = ilp_opt(A)
    if opt != len(members): continue
    tight += 1
    xs = [vid('x', j) for j in range(N)]
    def face(S):
        fixed = {}
        for a in S:
            for t in range(n):
                if t != a: fixed[key(a, t)] = 0
        r = ilp_opt(A, fixed=fixed)
        return r is not None and r[0] == opt
    sub = [face([xs[i] for i in range(N) if i != j]) for j in range(N)]
    full = face(xs)
    if all(sub) and not full:
        found += 1; print('RING CORE', N, aorder, forder, sorted(E), flush=True)
        break
print('N', N, 'tested', tested, 'tight', tight, 'found', found, flush=True)
