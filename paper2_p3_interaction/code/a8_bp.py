"""Audit 8: Lemma 7b, abstract part: (1) Baker-Pixley step on random relations closed under random uniform majority rules;
(2) 'S is a face iff the forced partial vector extends to Sigma' on real instances (anchors force L/R states)."""
import sys, random, itertools
import numpy as np
from a8_lib import *
rng = random.Random(7)
# ---- (1) random relations
SYMS = 'LRM'
def rule_random():
    return {p: rng.choice(SYMS) for p in itertools.permutations(SYMS)}
def maj(rule, a, b, c):
    out = []
    for x, y, z in zip(a, b, c):
        out.append(x if (x == y or x == z) else y if y == z else rule[(x, y, z)])
    return ''.join(out)
def close(S, rule):
    S = set(S)
    while True:
        new = {maj(rule, a, b, c) for a in S for b in S for c in S} - S
        if not new: return S
        S |= new
        if len(S) > 30: return None
cnt = 0; nontriv = 0
for it in range(3000):
    m = rng.choice([4, 5, 6])
    seeds = [''.join(rng.choice(SYMS) for _ in range(m)) for _ in range(rng.randint(2, 5))]
    R = close(seeds, rule_random())
    if R is None or len(R) > 30: continue
    cnt += 1
    for r in (2, 3):
        for T in itertools.combinations(range(m), r):
            projT = {tuple(v[i] for i in T) for v in R}
            for t in itertools.product(SYMS, repeat=r):
                inproj = t in projT
                binary = all(((t[a], t[b]) in {(v[T[a]], v[T[b]]) for v in R}) for a, b in itertools.combinations(range(r), 2))
                if r == 3 and binary != inproj: print('BAKER-PIXLEY FAIL', R, T, t); sys.exit()
                if r == 3 and not inproj: nontriv += 1
                if r == 3 and binary and not inproj: print('FAIL'); sys.exit()
print('Baker-Pixley step: relations', cnt, 'triple coordinate sets tested, nonmember tuples', nontriv, 'no failure')

# ---- (2) forced states
n_ok = 0; chk = 0; tried = 0; mism = 0
while n_ok < 500 and tried < 100000:
    tried += 1
    n = rng.choice([7, 8, 9])
    A = gen_planted(n, rng, nedit=rng.randint(2, 6))
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid() or I.opt() != len(H): continue
    anchors = [x for x in range(n) if not I.hard_nbrs(x) and any(x in t for t in H)]
    if len(anchors) < 2: continue
    n_ok += 1
    ZE = I.zero_excess(); vecs = [I.state(k) for k in ZE]
    for r in (1, 2, 3):
        for S in itertools.combinations(anchors, r):
            face = any(all(all(I.lab[k, v] != I.lab[k, x] for v in range(n) if v != x) for x in S) for k in ZE)
            forced = {}; bad = False
            for x in S:
                for i, t in enumerate(H):
                    if x in t:
                        if t[1] == x: bad = True
                        else:
                            s = 'L' if t[0] == x else 'R'
                            if forced.setdefault(i, s) != s: bad = True
            ext = (not bad) and any(all(v[i] == s for i, s in forced.items()) for v in vecs)
            chk += 1
            if ext != face: mism += 1; print('MISMATCH', A.tolist(), H, S)
print('forced-state equivalence: instances', n_ok, 'subsets', chk, 'mismatches', mism)
