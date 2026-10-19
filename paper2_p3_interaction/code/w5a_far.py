"""Do far INSERTED HH pointers (different components of I(H), here c=1 vertex-disjoint P3s) occur in every optimum?
n=9: three disjoint P3 (0,1,2),(3,4,5),(6,7,8) + random hard edges between them; exhaustive over partitions."""
import random, itertools, numpy as np
from w5a_lib import *
n = 9; H = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
comp = {v: i // 3 for i, v in enumerate(range(9))}
lab_t, same, P = part_table(n)
forced = 0; tot = 0; ex = None; cnt_t2 = 0
for sd in range(4000):
    rng = random.Random(sd); A = np.zeros((n, n), np.int8)
    for t in H:
        a, b, c = t; A[a, b] = A[b, a] = A[b, c] = A[c, b] = 1
    p = rng.choice([.15, .25, .4, .6])
    for i, j in prs(n):
        if comp[i] != comp[j] and rng.random() < p: A[i, j] = A[j, i] = 1
    if not is_packing(A, H) or not is_pair_maximal(A, H): continue
    tot += 1
    lab, _, cost = all_costs(A); opt = cost.min()
    E = np.array([A[i][j] for i, j in P], bool)
    farp = np.array([comp[i] != comp[j] for i, j in P])
    ins = (same & ~E[None, :] & farp[None, :]).sum(1)
    idx = np.nonzero(cost == opt)[0]
    if ins[idx].min() > 0:
        forced += 1
        if ex is None: ex = (sd, opt, [(i, j) for i, j in prs(n) if A[i][j]], [lab_t[k].tolist() for k in idx[:2]])
print('instances', tot, 'every optimum has an inserted pair between different P3-components:', forced)
print(ex)
