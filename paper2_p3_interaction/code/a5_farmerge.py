"""c=1, H = three vertex-disjoint P3s (I(H) edgeless, n=9). Is there an instance where EVERY optimum has a MERGE-type inserted pair
(different pieces) between vertices of DIFFERENT members (infinite I(H)-distance)?  numpy over all set partitions."""
import random, itertools, numpy as np
from a5_lib import *
n = 9; H = [(0,1,2),(3,4,5),(6,7,8)]; mem = lambda v: v//3
P = list(itertools.combinations(range(n), 2))
labs = np.array(list(partitions(n)), dtype=np.int8)
same = np.stack([labs[:, i] == labs[:, j] for i, j in P], 1)
pk = packed(H); soft = np.array([frozenset(p) in pk for p in P])
far = np.array([mem(i) != mem(j) for i, j in P])
found = 0; tot = 0
for sd in range(3000):
    rng = random.Random(sd); p = rng.choice([.15, .25, .4, .6])
    E = [(0,1),(1,2),(3,4),(4,5),(6,7),(7,8)] + [(i, j) for i, j in P if mem(i) != mem(j) and rng.random() < p]
    adj = mk(n, E)
    if not (is_packing(adj, H) and pair_maximal(adj, H)): continue
    tot += 1
    e = np.array([adj[i][j] for i, j in P], bool)
    cost = (same ^ e[None, :]).sum(1); idx = np.nonzero(cost == cost.min())[0]
    allf = True
    for k in idx:
        lab = tuple(labs[k]); pid = pieces(adj, pk, lab)
        ok = any(lab[i] == lab[j] and not adj[i][j] and mem(i) != mem(j) and pid[i] != pid[j] for i, j in P)
        if not ok: allf = False; break
    if allf:
        found += 1
        if found <= 3: print('FORCED far merge: edges', E, 'optima', [tuple(labs[k]) for k in idx])
print('instances', tot, 'forced far merge-type', found)
