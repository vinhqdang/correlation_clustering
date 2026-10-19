"""Audit 7: is sum_pi c(pi) = 2D|U| attained, and for which D?  random search over all partitions with cost<=|H|+3."""
import sys, random
import numpy as np
import a7_lib as L
rng = random.Random(int(sys.argv[1])); n = int(sys.argv[2]); N = int(sys.argv[3])
best = {}   # D -> max ratio
for it in range(N):
    A = L.gen_graph(n, rng)
    ts = L.ind_p3s(A)
    if not ts: continue
    # random packing; random drop to get non-pair-maximal too? keep pair-maximal
    mode = rng.random()
    H = L.greedy_pm(A, rng)
    if mode < .3 and len(H) > 1: H = H[:-1]
    if not H: continue
    if not L.is_pair_maximal(A, H): continue
    S = L.soft_of(n, H); Dl, Gm = L.delta_gamma(A, S)
    c = max(sum(1 for t in H if v in t) for v in range(n)); D = 2 * c
    P = L.partitions(n); cs, X, pl = L.costs(A, P); idx = {p: k for k, p in enumerate(pl)}
    E = np.array([bool(A[i][j]) for i, j in pl])
    for r in np.where(cs <= len(H) + 3)[0]:
        ed = X[r] != E
        u = sum(1 for k in np.where(ed)[0] if not S[pl[k]])
        if not u: continue
        tot = sum((len(Dl[p]) if X[r][idx[p]] else len(Gm[p])) for p in pl if S[p])
        rt = tot / (2 * D * u)
        if rt > best.get(D, (0,))[0]: best[D] = (rt, A.tolist(), H, X[r].tolist() and [int(x) for x in P[r]])
for D in sorted(best): print("D", D, "max ratio", best[D][0])
if len(best) and 4 in best: print(best[4][1:])
