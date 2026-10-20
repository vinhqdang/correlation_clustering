"""w11_wit.py -- are the pairs of a minimum certificate witnesses (induced P3 with exactly one soft pair, other pairs hard) of a far edit?"""
import sys, random, itertools
from collections import Counter
from w11_lib import *
from w11_gen import *
rnd = random.Random(int(sys.argv[1])); samples = int(sys.argv[2]); m = int(sys.argv[3]); npend = int(sys.argv[4])
st = Counter(); fr = Counter()
for s in range(samples):
    nv, H = random_members_c2(m, rnd)
    n = nv + npend
    k = rnd.choice([2, 3, 4, 6]); cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
    I = make_inst(nv, H, npend, cl)
    if not (I.valid() and I.pair_maximal() and I.cmax() <= 2): continue
    sz, T = min_certificate(I, 0, max_size=30, workers=1)
    if not sz: continue
    st['inst'] += 1
    for e in T:
        # witnesses: induced P3s of G containing e with exactly one soft pair
        found = False; found_tri = False
        for z in range(n):
            if z in e: continue
            ps = [key(e[0], z), key(e[1], z), e]
            soft = [p for p in ps if p in I.owner]
            edges = sum(1 for p in ps if p in I.E)
            if len(soft) == 1:
                if edges == 2: found = True
                if edges == 3: found_tri = True
        fr[('P3-witness' if found else ('triangle-only' if found_tri else 'none'))] += 1
print(st, fr)
