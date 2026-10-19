"""Lemma 3 on ALL partitions (n=7, any packing, not nec. maximal/optimal): |{a: az,ab hard, [az in E] != [ab in E]}| <= |U| for z,b in one cluster.
Also checks the anchor form M(z,b)=N_h(z) \\ (N[b] u Q(b)) for soft (z,b)."""
import random, itertools
from a5_lib import *
viol = 0; chk = 0; anch = 0
for sd in range(60):
    rng = random.Random(sd); n = 7; p = rng.uniform(.2, .7)
    adj = mk(n, [(i, j) for i in range(n) for j in range(i+1, n) if rng.random() < p])
    H = rand_pm(adj, rng, rng.choice([0, .5])); pk = packed(H)
    for lab in partitions(n):
        U = len(U_set(adj, pk, lab))
        for z, b in itertools.permutations(range(n), 2):
            if lab[z] != lab[b]: continue
            S = [a for a in range(n) if a not in (z, b) and frozenset((a, z)) not in pk and frozenset((a, b)) not in pk and adj[a][z] != adj[a][b]]
            chk += 1
            if len(S) > U: viol += 1
            if frozenset((z, b)) in pk:
                M = [a for a in range(n) if a != z and a != b and adj[a][z] and frozenset((a, z)) not in pk and not adj[a][b] and frozenset((a, b)) not in pk]
                anch += 1
                if len(M) > U: viol += 1
print('checks', chk, 'anchor checks', anch, 'violations', viol)
