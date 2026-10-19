"""Mutation test: each deliberately weakened statement must be falsified by the same machinery (shows tests are not vacuous)."""
import random, collections, itertools
import numpy as np, a4_lib as L
rng = random.Random(77); F = collections.Counter(); N = 0
for it in range(1500):
    n = rng.randint(6, 9); A = L.gen_graph(rng, n)
    if rng.random() < .5: A = L.noise(A, rng, 2)
    H = L.greedy_extend(A, L.rand_packing(A, rng, stop=.3), rng); pk = L.packed_matrix(n, H); h = len(H)
    D = L.soft_degree(pk); c = L.cmax(n, H); lab, same, cost = L.all_costs(A); OPT = int(cost.min()); N += 1
    hs = L.hub_stats(A, pk, H)
    for l in range(4):
        if OPT <= h + l:
            for X, al, nh in hs:
                if X > 2 * l: F['HUB with 2l instead of 2l+1'] += 1
                if al > (c) * (2 * l + 1) : F['HUB alpha with c instead of c+1'] += 1
    if 0 < h <= 12:
        w = L.tw_exact(L.ih_adj(H)); tw = L.tw_exact(L.soft_adj(n, H))
        if tw > 3 * (w + 1) - 2: F['softtw 3(w+1)-2'] += 1
    M = L.rand_M(A, H, rng)
    if len(H) <= 11 and len(H + M) <= 13 and M:
        if L.tw_exact(L.ih_adj(H + M)) > L.tw_exact(L.ih_adj(H)) + len(M) - 1: F['MAX tw w+|M|-1'] += 1
    for r in np.nonzero(cost <= h + 2)[0][:60]:
        v, st = L.sm_check(A, pk, lab[r], max(D - 1, 0))
        for x in v: F['SM with Ds-1 (' + x[0] + ')'] += 1
print(N, 'instances; falsifications of weakened claims:', dict(F))
