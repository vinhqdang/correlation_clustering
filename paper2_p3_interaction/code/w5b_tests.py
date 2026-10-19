"""Wave 5 B: brute-force tests of the cluster-partition lemmas (all set partitions, n <= 10).
usage: python3 -I w5b_tests.py N_INSTANCES SEED NMIN NMAX"""
import sys, random, itertools
import numpy as np
sys.path.insert(0, '.')
from w5b_lib import *

LMAX = 3
def run(ninst, seed, nmin, nmax):
    rng = random.Random(seed); st = dict(inst=0, sols=0, opt_sols=0, hardpairs_checked=0, nontriv_big=0,
        maxgap_d=-99, maxgap_C=-99, fine_d=-99, fine_C=-99, maxC=0, beta_multi_max=0, beta_checked=0,
        fails=[], idchk=0, he=0, disj=0, multi_comp_alpha=0)
    for it in range(ninst):
        n = rng.randint(nmin, nmax); A = rand_graph(n, rng)
        if A.sum() == 0: continue
        H = pair_maximal_packing(A, rng)
        if not H: continue
        pk = soft_matrix(n, H); D = int(pk.sum(1).max()); comp, hg = hard_comps(A, pk)
        dh = hg.sum(1); nE = int(A.sum() // 2); h = len(H)
        # Lemma H1/H2 (pair-maximality consequences)
        for a in range(n):
            for b in range(n):
                if hg[a, b]:
                    for z in range(n):
                        if z != a and z != b and hg[a, z]:
                            st['he'] += 1
                            if not (A[z][b] or pk[z, b]): st['fails'].append(('HE', it))
        for u, v in prs(n):
            if not A[u][v] and not pk[u, v]:
                st['disj'] += 1
                if (hg[u] & hg[v]).any(): st['fails'].append(('H2', it))
        P = all_partitions(n); c = costs(A, P); opt = c.min()
        st['inst'] += 1
        idx = np.nonzero(c <= h + LMAX)[0]
        for ii in idx:
            lab = P[ii]; cost = int(c[ii]); st['sols'] += 1
            same = lab[:, None] == lab[None, :]
            edited = (same != A.astype(bool)); np.fill_diagonal(edited, False)
            clusters = {}
            for v in range(n): clusters.setdefault(int(lab[v]), []).append(v)
            # identity: cost = |E| + sum_C [C(|C|,2) - 2 e(C)]
            tot = nE
            for C in clusters.values():
                e = sum(A[x][y] for x, y in itertools.combinations(C, 2))
                tot += len(C) * (len(C) - 1) // 2 - 2 * e
            st['idchk'] += 1
            if tot != cost: st['fails'].append(('ID', it))
            soft_ed = int((edited & pk).sum() // 2); hard_ed = int((edited & ~pk).sum() // 2)
            if soft_ed + hard_ed != cost or soft_ed < h: st['fails'].append(('SPLIT', it))
            U = [(x, y) for x, y in prs(n) if edited[x, y] and not pk[x, y]]
            nU = len(U)
            if nU > cost - h: st['fails'].append(('U<=l', it))
            dele = [0] * n; ins = [0] * n
            for x, y in U:
                t = dele if A[x][y] else ins
                t[x] += 1; t[y] += 1
            ends = set(v for p in U for v in p)
            for (x, y) in U:
                for (u, v) in ((x, y), (y, x)):
                    st['hardpairs_checked'] += 1
                    if A[u][v]: bound = D + dele[u] + dele[v] - 1
                    else: bound = D + ins[v] - 1 + dele[u]
                    st['fine_d'] = max(st['fine_d'], dh[u] - bound)
                    st['maxgap_d'] = max(st['maxgap_d'], dh[u] - (D + nU))
                    if dh[u] > bound or dh[u] > D + nU: st['fails'].append(('EX1', it, u, v))
                    C = clusters[int(lab[u])]
                    b2 = 1 + D + dh[u] + ins[u]
                    st['fine_C'] = max(st['fine_C'], len(C) - b2)
                    st['maxgap_C'] = max(st['maxgap_C'], len(C) - (2 * D + 2 * nU + 1))
                    st['maxC'] = max(st['maxC'], len(C))
                    if len(C) > b2 or len(C) > 2 * D + 2 * nU + 1: st['fails'].append(('EX2', it, u))
                    for y in C:
                        if dh[y] > len(C) - 1 + dele[y]: st['fails'].append(('EX3', it, y))
            # normal form of clusters with no endpoint
            for C in clusters.values():
                if ends & set(C): continue
                st['beta_checked'] += 1
                S = set(C)
                for v in C:
                    for u in range(n):
                        if hg[v, u] and u not in S: st['fails'].append(('NFclosed', it))
                cs = set(comp[v] for v in C)
                for x, y in itertools.combinations(C, 2):
                    if comp[x] != comp[y] and not pk[x, y]: st['fails'].append(('NFsoftcomplete', it))
                if len(cs) >= 2:
                    st['beta_multi_max'] = max(st['beta_multi_max'], len(C))
                    if len(C) > 2 * D: st['fails'].append(('NFsize', it))
            # F1 on optima
            if cost == opt:
                st['opt_sols'] += 1
                for C in clusters.values():
                    S = set(C); seen = {C[0]}; stk = [C[0]]
                    while stk:
                        v = stk.pop()
                        for u in C:
                            if A[v][u] and u not in seen: seen.add(u); stk.append(u)
                    if seen != S: st['fails'].append(('F1', it))
            # BIG: vertex with d_h > 2D+3*ex lies in a cluster equal to its hard component
            ex = cost - h
            for v in range(n):
                if dh[v] > 2 * D + 3 * ex:
                    st['nontriv_big'] += 1
                    C = set(clusters[int(lab[v])]); K = set(u for u in range(n) if comp[u] == comp[v])
                    if C != K: st['fails'].append(('BIG', it, v))
    return st
if __name__ == '__main__':
    ninst, seed, nmin, nmax = map(int, sys.argv[1:5])
    st = run(ninst, seed, nmin, nmax)
    f = st.pop('fails'); print(st); print('FAILS:', len(f), f[:10])
