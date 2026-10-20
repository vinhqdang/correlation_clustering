"""Among tight instances that have a pairwise-consistent non-extendable L/R 3-assignment (state-level core), is any minimal non-face of size >= 3 realisable by singleton locks on all-soft vertices?"""
import sys, itertools, random, time, numpy as np, networkx as nx
from w6c_lib import *
from w6c_nflib import packings_of_size, minimal_nonfaces
from w6c_majority import states
n = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
tot = 0; withcore = 0; withlock3 = 0; anchor_counts = []
while time.time() - t0 < T:
    g = nx.gnm_random_graph(n, rng.randint(n, int(2.3 * n)), seed=rng.randrange(10 ** 9))
    if not nx.is_connected(g): continue
    A = np.array(nx.to_numpy_array(g), dtype=np.int8)
    opt, lab, same, cost, P = brute_opt(A)
    Hs = packings_of_size(A, opt)
    if not Hs: continue
    optidx = np.nonzero(cost == opt)[0]
    for H in rng.sample(Hs, min(3, len(Hs))):
        Sig = {states(lab[i], H) for i in optidx}
        p = len(H)
        if p < 3 or len(Sig) < 2: continue
        tot += 1
        S = np.array(sorted(Sig)); found = False
        for coords in itertools.combinations(range(p), 3):
            sub = S[:, coords]; trip = {tuple(r) for r in sub}
            pairs = {(i, j): {(r[i], r[j]) for r in sub} for i, j in ((0, 1), (0, 2), (1, 2))}
            for v in itertools.product((0, 1), repeat=3):
                if v not in trip and all((v[i], v[j]) in pairs[(i, j)] for i, j in ((0, 1), (0, 2), (1, 2))): found = True; break
            if found: break
        if not found: continue
        withcore += 1
        used = set()
        for t in H: used |= pairs_of(t)
        inm = {}
        for i, m in enumerate(H):
            for pos, v in enumerate(m): inm.setdefault(v, []).append(pos)
        anch = [v for v in range(n) if v in inm and all(pos in (0, 2) for pos in inm[v]) and all(key(v, u) in used for u in range(n) if A[v][u])]
        anchor_counts.append(len(anch))
        if len(anch) < 3: continue
        sing = set()
        for i in optidx:
            l = lab[i]; cl = {}
            for v in range(n): cl.setdefault(l[v], []).append(v)
            sing.add(frozenset(v for b in cl.values() if len(b) == 1 for v in b))
        nf = minimal_nonfaces(list(sing), anch)
        if any(len(s) >= 3 for s in nf):
            withlock3 += 1; print('LOCK 3-CORE', sorted(g.edges()), H, [sorted(s) for s in nf if len(s) >= 3], flush=True)
print('n', n, 'instance-packings', tot, 'with state-level L/R 3-core', withcore, 'of which with >=3 anchors', sum(1 for c in anchor_counts if c >= 3), 'lock cores of size >= 3:', withlock3)
