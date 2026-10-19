"""Wave 5 B: candidate-partner locality (Corollary CAND).  For every OPTIMAL partition with cost <= |H|+l and every edited hard pair
(u,v): u,v have d_h <= B = 2D+3*ex, v is within distance <= 2D+2*nU of u in G[L_B], L_B = {d_h <= B}  (deleted edges: distance 1).
usage: python3 -I w5b_cand.py N SEED NMIN NMAX"""
import sys, random
import numpy as np
from collections import deque
sys.path.insert(0, '.')
from w5b_lib import *
def run(ninst, seed, nmin, nmax):
    rng = random.Random(seed); chk = 0; ins_chk = 0; maxd = 0; fails = []; maxslack = -99
    for it in range(ninst):
        n = rng.randint(nmin, nmax); A = rand_graph(n, rng)
        if A.sum() == 0: continue
        H = pair_maximal_packing(A, rng)
        if not H: continue
        pk = soft_matrix(n, H); D = int(pk.sum(1).max()); comp, hg = hard_comps(A, pk); dh = hg.sum(1)
        P = all_partitions(n); c = costs(A, P); opt = int(c.min()); ex = opt - len(H)
        for ii in np.nonzero(c == opt)[0]:
            lab = P[ii]; same = lab[:, None] == lab[None, :]
            U = [(x, y) for x, y in prs(n) if (same[x, y] != bool(A[x][y])) and not pk[x, y]]
            B = 2 * D + 3 * ex; L = [v for v in range(n) if dh[v] <= B]; Ls = set(L)
            for (u, v) in U:
                chk += 1
                if u not in Ls or v not in Ls: fails.append(('L', it)); continue
                dist = {u: 0}; q = deque([u])
                while q:
                    x = q.popleft()
                    for y in L:
                        if A[x][y] and y not in dist: dist[y] = dist[x] + 1; q.append(y)
                d = dist.get(v, 10**9)
                bound = 1 if A[u][v] else 2 * D + 2 * len(U)
                if not A[u][v]: ins_chk += 1; maxd = max(maxd, d)
                maxslack = max(maxslack, d - bound)
                if d > bound: fails.append(('DIST', it, u, v, d, bound))
    return dict(pairs=chk, inserted=ins_chk, max_dist_inserted=maxd, max_slack=maxslack, fails=fails[:5], nfails=len(fails))
if __name__ == '__main__':
    print(run(*map(int, sys.argv[1:5])))
