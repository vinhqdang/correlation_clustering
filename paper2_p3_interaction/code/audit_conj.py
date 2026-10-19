"""Audit of the C4 example generalisation and of 'Conjecture C1' (no optimum forced to use an unpacked edit with both endpoints in R)."""
import random, sys
import numpy as np
from audit_lib import *

rng = random.Random(int(sys.argv[1])); trials = int(sys.argv[2])
touchR = rr = tot = 0; ex = None
for t in range(trials):
    n = rng.randint(5, 9)
    E = rand_graph(rng, n, rng.choice(['gnp', 'hub', 'cliques', 'path']))
    H = rand_packing(rng, n, E, 'max')
    if not H: continue
    pk = check_packing(n, E, H)
    if not all(pairmask(n, q) & pk for q in induced_p3s(n, E)): continue
    Vh = set().union(*[set(q) for q in H])
    R = [x for x in range(n) if x not in Vh]
    if len(R) < 2: continue
    P = pidx(n)
    S, c = opt_and_costs(n, E); o = int(c.min()); sel = S[c == o]
    unp = ~np.uint64(pk)
    maskR = 0; maskRR = 0
    for r in R:
        for v in range(n):
            if v != r: maskR |= 1 << P[(r, v)]
    for r in R:
        for s in R:
            if r < s: maskRR |= 1 << P[(r, s)]
    tot += 1
    a = ((sel & unp & np.uint64(maskR)) != 0)
    b = ((sel & unp & np.uint64(maskRR)) != 0)
    if a.all(): touchR += 1
    # forced RR: every optimum uses an RR unpacked pair
    if b.all():
        rr += 1; ex = (n, E, H)
    # also: is there an optimum with no unpacked edit with both ends in R? (C1)
    if not (~b).any():
        print('C1 counterexample', n, E, H)
print('maximal-H instances with |R|>=2:', tot, ' every optimum touches R by an unpacked edit:', touchR, ' forced R-R:', rr)
