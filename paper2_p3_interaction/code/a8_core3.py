"""Audit 8: independent check of wave6 Sec 4.4: tight connected graphs have no minimal non-face of size >= 3 (faces of vertex sets: S is a face iff some optimal
partition makes all x in S singletons; ZE solutions = optimal partitions when OPT = max P3-packing number)."""
import sys, random, itertools
import numpy as np, networkx as nx
from networkx.generators.atlas import graph_atlas_g
from a8_lib import *
def packing_number(A):
    n = len(A); L = induced_p3s(A); best = 0
    P = [set(mpairs(t)) for t in L]
    def rec(i, used, cnt):
        nonlocal best
        if cnt + (len(L) - i) <= best: return
        if i == len(L): best = max(best, cnt); return
        if not (P[i] & used): rec(i + 1, used | P[i], cnt + 1)
        rec(i + 1, used, cnt)
    rec(0, frozenset(), 0); return best
def analyse(A):
    n = len(A); I = Inst(A, [])
    opt = I.opt()
    if opt == 0: return None
    if packing_number(A) != opt: return 'notight'
    ks = np.nonzero(I.cost == opt)[0]
    sing = np.array([[all(I.lab[k, v] != I.lab[k, x] for v in range(n) if v != x) for x in range(n)] for k in ks])
    F = {}
    mx = 0
    for r in range(1, n + 1):
        any_ = False
        for S in itertools.combinations(range(n), r):
            f = bool(sing[:, list(S)].all(1).any()); F[S] = f
            if not f and all((tuple(y for y in S if y != x) in F and F[tuple(y for y in S if y != x)]) or len(S) == 1 for x in S):
                mx = max(mx, r)
    return mx
res = {}
for G in graph_atlas_g():
    n = G.number_of_nodes()
    if n < 4 or n > 7 or not nx.is_connected(G): continue
    A = nx.to_numpy_array(G, dtype=np.int8)
    r = analyse(A); k = ('atlas', n, r); res[k] = res.get(k, 0) + 1
print('atlas', sorted(res.items(), key=str))
rng = random.Random(3); res = {}
for it in range(int(sys.argv[1]) if len(sys.argv) > 1 else 1500):
    n = rng.choice([8, 9]); A = gen_planted(n, rng, nedit=rng.randint(2, 6)) if rng.random() < 0.7 else gen_random(n, rng)
    if not nx.is_connected(nx.from_numpy_array(A)): continue
    r = analyse(A); k = (n, r); res[k] = res.get(k, 0) + 1
print('random', sorted(res.items(), key=str))
