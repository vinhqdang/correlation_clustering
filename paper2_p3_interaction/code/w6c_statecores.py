"""State-level cores: minimal pairwise-consistent partial assignments with values in {L,R} (the only states a clean-end lock can force) that do not extend to a zero-excess solution."""
import sys, itertools, random, time, numpy as np, networkx as nx
from w6c_lib import *
from w6c_nflib import packings_of_size
from w6c_majority import states
n = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
tot = 0; withcore = 0; ex = None
while time.time() - t0 < T:
    g = nx.gnm_random_graph(n, rng.randint(n, int(2.3 * n)), seed=rng.randrange(10 ** 9))
    if not nx.is_connected(g): continue
    A = np.array(nx.to_numpy_array(g), dtype=np.int8)
    opt, lab, same, cost, P = brute_opt(A)
    Hs = packings_of_size(A, opt)
    if not Hs: continue
    optidx = np.nonzero(cost == opt)[0]
    for H in rng.sample(Hs, min(2, len(Hs))):
        Sig = {states(lab[i], H) for i in optidx}
        p = len(H)
        if p < 3 or len(Sig) < 2: continue
        tot += 1
        S = np.array(sorted(Sig))
        found = None
        # partial assignments on 3 coordinates with values in {0,1} (L,R) pairwise consistent but not extendable
        for coords in itertools.combinations(range(p), 3):
            sub = S[:, coords]
            trip = {tuple(r) for r in sub}
            pairs = {(i, j): {(r[i], r[j]) for r in sub} for i, j in ((0, 1), (0, 2), (1, 2))}
            for v in itertools.product((0, 1), repeat=3):
                if v in trip: continue
                if all((v[i], v[j]) in pairs[(i, j)] for i, j in ((0, 1), (0, 2), (1, 2))):
                    found = (coords, v); break
            if found: break
        if found:
            withcore += 1
            if ex is None: ex = (sorted(g.edges()), H, found)
print('n', n, 'instances', tot, 'with a pairwise-consistent non-extendable 3-assignment over {L,R}:', withcore)
if ex: print(ex)
