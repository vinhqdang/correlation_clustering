"""Localizability test.  Is there an OPTIMAL solution all of whose inserted hard pairs (a,b) have a common soft neighbour
(W-distance exactly 2, any z in V)?  Exhaustive n<=9 (min over all optimal partitions), ILP n=10..13 (min via side objective).
Also the same with 'hard-edge-distance' alternative: counts pairs that are W^2-nonlocal."""
import sys, random, itertools, numpy as np
from w5a_lib import *
from a4_lib import _tri, ilp

def w2local(pk):
    n = len(pk); L = (pk.astype(int) @ pk.astype(int)) > 0
    np.fill_diagonal(L, False); return L & ~pk

def exh_min(A, pk, H_):
    n = len(A); lab_t, same, P = part_table(n); lab, same_, cost = all_costs(A); opt = cost.min()
    L = w2local(pk); E = np.array([A[i][j] for i, j in P], bool)
    VH = {v for t in H_ for v in t}
    bad = np.array([ (not pk[i, j]) and (not L[i, j]) and i in VH and j in VH for i, j in P])      # hard pair, not W^2-local
    ins = same & ~E[None, :]                                                   # inserted (hard or soft) pairs
    nb = (ins & bad[None, :]).sum(1)
    idx = np.nonzero(cost == opt)[0]
    return opt, int(nb[idx].min()), int(nb[idx].max()), len(idx)

def ilp_min(A, pk, H_):
    n = len(A); P = prs(n); r0 = ilp(A); opt = r0[1]
    L = w2local(pk); E = np.array([A[i][j] for i, j in P], float)
    VH = {v for t in H_ for v in t}
    badv = np.array([1.0 if ((not pk[i, j]) and (not L[i, j]) and not A[i][j] and i in VH and j in VH) else 0.0 for i, j in P])
    # objective: cost*1000 + #bad inserted pairs   (z=1 means same cluster; inserted bad = badv*z)
    cvec = 1 - 2 * E
    r = ilp(A, obj=1000 * cvec + badv, cost_le=opt)
    z, cost, idx = r
    return opt, int(round(badv @ z)), cost

def structured(rng, n, kind):
    A = np.zeros((n, n), np.int8)
    def add(i, j): A[i, j] = A[j, i] = 1
    if kind == 'path':
        for i in range(n - 1): add(i, i + 1)
    elif kind == 'cycle':
        for i in range(n): add(i, (i + 1) % n)
    elif kind == 'pathplus':
        for i in range(n - 1): add(i, i + 1)
        for _ in range(rng.randint(1, 3)):
            i, j = rng.sample(range(n), 2); add(i, j)
    elif kind == 'caterpillar':
        k = n // 2
        for i in range(k - 1): add(i, i + 1)
        for i in range(k, n): add(i, rng.randrange(k))
    elif kind == 'twopaths':
        k = n // 2
        for i in range(k - 1): add(i, i + 1)
        for i in range(k, n - 1): add(i, i + 1)
        for _ in range(2):
            add(rng.randrange(k), rng.randrange(k, n))
    return A

if __name__ == '__main__':
    mode = sys.argv[1]; n = int(sys.argv[2]); N = int(sys.argv[3])
    tot = nonloc = 0; worst = None
    for sd in range(N):
        rng = random.Random(f"w5aloc-{mode}-{n}-{sd}")
        if mode in ('path', 'cycle', 'pathplus', 'caterpillar', 'twopaths'):
            A = structured(rng, n, mode)
        else:
            A = gen_graph(rng, n, mode)
            if rng.random() < .5: A = noise(A, rng, rng.randint(1, 3))
        H = greedy_extend(A, rand_packing(A, rng, stop=rng.choice([0, .3, .6]), prefer_overlap=rng.random() < .5), rng, prefer_overlap=rng.random() < .5)
        if not H or not is_pair_maximal(A, H): continue
        pk = packed_matrix(n, H)
        if n <= 9: opt, mn, mx, k = exh_min(A, pk, H)
        else: opt, mn, _ = ilp_min(A, pk, H); mx = -1
        tot += 1
        if mn > 0:
            nonloc += 1
            if worst is None: worst = (sd, H, A.copy(), opt, mn)
    print(mode, n, 'instances', tot, 'no-W2-local-optimum', nonloc)
    if worst:
        sd, H, A, opt, mn = worst
        print(' first: seed', sd, 'opt', opt, '|H|', len(H), 'min #nonlocal inserted', mn, 'H', H)
        print(' edges', [(i, j) for i, j in prs(len(A)) if A[i][j]])
