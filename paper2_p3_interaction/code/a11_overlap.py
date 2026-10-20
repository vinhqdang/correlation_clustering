"""Audit 11: are 'A-local (some witness with two soft pairs)' and 'absorbable (some witness with exactly one soft pair)' exclusive? count."""
import random, sys
from a11_lib import *
rng = random.Random(7); c = dict(edits=0, both=0, only_two=0, only_one=0)
for _ in range(1500):
    n = rng.randint(6, 9); p = rng.choice([0.25, 0.35, 0.5])
    g = Gr(n, [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < p])
    if not g.p3(): continue
    H = random_packing(g, rng); own = soft_owner(H)
    o, opts = opt_and_optima(g)
    for S in opts:
        for e in hard_of(H, S):
            ks = [sum(q in own for q in pr(t)) for t in g.p3() if e in pr(t)]
            c['edits'] += 1
            t2, t1 = 2 in ks, 1 in ks
            c['both'] += t1 and t2; c['only_two'] += t2 and not t1; c['only_one'] += t1 and not t2
print(c)
