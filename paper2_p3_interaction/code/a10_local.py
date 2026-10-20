"""Lemma 2 locality: a hard pair {p,r} flipped; which triple-constraints change, by scope. Claim: only triples containing p,r;
unary changes sit on members through p or r; a binary change (A,B) needs p in A, r in B (or v.v.) and A,B adjacent (sharing the third vertex)."""
import sys, random, itertools
import numpy as np
from a10_lib import *
from a10_closure import apply
def table(A, H, t, ow):
    out = {}
    for st in itertools.product(range(3), repeat=len(ow)):
        B = apply(A, H, dict(zip(ow, st))); x, y, z = t
        out[st] = int(B[x, y]) + int(B[y, z]) + int(B[x, z]) != 2
    return out
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
viol = 0; chg = {1: 0, 2: 0, 0: 0}; checked = 0; far_unary_only = 0; far_pairs = 0
for it in range(N):
    n = rng.randint(6, 9); g = planted(n, rng.randint(2, 6), rng, big=rng.random() < .5)
    if g is None: continue
    A, H = g; o = owner_map(H)
    if soft_triangles(n, H): continue
    mem = lambda v: {i for i, t in enumerate(H) if v in t}
    for p, r in itertools.combinations(range(n), 2):
        if frozenset((p, r)) in o: continue
        B = A.copy(); B[p, r] = B[r, p] = not B[p, r]
        far = not any(set(H[i]) & set(H[j]) for i in mem(p) for j in mem(r))
        far_pairs += far; changed_scopes = []
        for t in itertools.combinations(range(n), 3):
            ow = sorted({o[q] for q in map(frozenset, itertools.combinations(t, 2)) if q in o})
            tA, tB = table(A, H, t, ow), table(B, H, t, ow)
            if tA != tB:
                checked += 1; chg[len(ow)] += 1; changed_scopes.append(ow)
                if not (p in t and r in t): viol += 1
                if not set(ow) <= (mem(p) | mem(r)): viol += 1
                if len(ow) == 2:
                    a, b = ow
                    okk = (p in H[a] and r in H[b]) or (p in H[b] and r in H[a])
                    if not (okk and set(H[a]) & set(H[b])): viol += 1
        if far and any(len(s) > 1 for s in changed_scopes): viol += 1
print("changed triple-constraints", checked, chg, "far pairs", far_pairs, "violations", viol)
