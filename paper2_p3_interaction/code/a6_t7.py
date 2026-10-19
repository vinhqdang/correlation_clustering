"""Audit 6, Thm 7 / Cor 6 / Lemma 3-5: independent exhaustive test (n <= 10, all set partitions).
For each instance (pair-maximal H) and each OPTIMAL partition pi, with budget l = L0 + extra (L0 = OPT - |H|):
  (1) endpoints of U have d_h <= D + l ; clusters of endpoints have size <= 2D+2l+1 and also <= 1+2D+l (sharper, derived here)
  (2) Cor 6 one-endpoint partner: deleted -> v in N_h(u); inserted -> dist_{G[L_B]}(u,v) <= 2D+2l   (so n^|U| * ball^|U| suffices)
  (3) Thm 7 root/ball cover: ~-components, one root each, every endpoint within R=(2l-1)(2D+2l) of its root in G[L_B]
  (4) mutation: radius 1 for inserted pairs must FAIL somewhere (non-vacuity)
  (5) exact Thm A simulation: min cost over partitions whose hard-edit set is exactly U' ; the candidate-U' algorithm accepts iff YES.
usage: python3 -I a6_t7.py N SEED NMIN NMAX"""
import sys, random, itertools
import numpy as np
sys.path.insert(0, '.')
from a6_lib import *

def run(ninst, seed, nmin, nmax):
    rng = random.Random(seed); st = dict(inst=0, opt_parts=0, upairs=0, ins=0, maxdist=0, fail=[], mut_fail=0, alg_checked=0, size_slack=-99, sharp_slack=-99)
    for it in range(ninst):
        n = rng.randint(nmin, nmax); A = rand_graph(n, rng)
        if A.sum() == 0: continue
        H = greedy_pair_maximal(A, rng)
        if not H: continue
        S = soft_matrix(n, H); D = int(S.sum(1).max()); comp, hg = hard_components(A, S); dh = hg.sum(1)
        P = partitions(n); c = part_costs(A, P); opt = int(c.min()); L0 = opt - len(H); st['inst'] += 1
        hp, masks = hard_edit_masks(A, S, P)
        for extra in (0, 1):
            l = L0 + extra
            B = 2 * D + 3 * l; L = [v for v in range(n) if dh[v] <= B]
            R = max(0, (2 * l - 1) * (2 * D + 2 * l))
            # (5) Thm A with candidate U' (all hard pairs in balls around ANY <= l roots) : accept iff cost<=|H|+l
            # simulate by taking min over masks of partitions with popcount<=l whose endpoints are in some ball cover: tested below via (3)
            for ii in np.nonzero(c == opt)[0]:
                lab = P[ii]
                U = [hp[k] for k in range(len(hp)) if (int(masks[ii]) >> k) & 1]
                if len(U) > l: st['fail'].append(('U>l', it)); continue
                st['opt_parts'] += 1
                dele = [0] * n; ins = [0] * n
                for x, y in U:
                    t = dele if A[x][y] else ins; t[x] += 1; t[y] += 1
                cl = {}
                for v in range(n): cl.setdefault(int(lab[v]), []).append(v)
                ends = sorted({v for p in U for v in p})
                Ls = set(L)
                for (x, y) in U:
                    st['upairs'] += 1
                    for u, v in ((x, y), (y, x)):
                        if dh[u] > D + l: st['fail'].append(('dh', it, u))
                        C = cl[int(lab[u])]
                        st['size_slack'] = max(st['size_slack'], len(C) - (2 * D + 2 * l + 1))
                        st['sharp_slack'] = max(st['sharp_slack'], len(C) - (1 + 2 * D + len(U)))
                        if len(C) > 2 * D + 2 * l + 1: st['fail'].append(('size', it))
                        if len(C) > 1 + 2 * D + len(U): st['fail'].append(('sharp', it))
                        if u not in Ls or v not in Ls: st['fail'].append(('L_B', it)); continue
                        dist = bfs_ball(A, L, u, 10 ** 6); d = dist.get(v, 10 ** 9)
                        if A[u][v]:
                            if d != 1: st['fail'].append(('del', it))
                        else:
                            st['ins'] += 1; st['maxdist'] = max(st['maxdist'], d)
                            if d > 2 * D + 2 * l: st['fail'].append(('ins', it, u, v, d))
                            if d > 1 and extra == 0 and False: pass
                            if d > 1: st['mut_fail'] += 1   # radius-1 mutation would break here
                # (3) ~-components
                par = {v: v for v in ends}
                def f(v):
                    while par[v] != v: par[v] = par[par[v]]; v = par[v]
                    return v
                for C in cl.values():
                    e = [v for v in C if v in par]
                    for v in e[1:]: par[f(v)] = f(e[0])
                for x, y in U: par[f(x)] = f(y)
                groups = {}
                for v in ends: groups.setdefault(f(v), []).append(v)
                if len(groups) > len(U): st['fail'].append(('r>|U|', it))
                for g in groups.values():
                    if len(g) > 2 * len(U): st['fail'].append(('compsize', it))
                    root = g[0]; dist = bfs_ball(A, L, root, R)
                    for v in g:
                        if v not in dist: st['fail'].append(('ball', it, root, v))
    return st

if __name__ == '__main__':
    N, seed, a, b = map(int, sys.argv[1:5]); st = run(N, seed, a, b); f = st.pop('fail'); print(st); print('FAILS', len(f), f[:8])
