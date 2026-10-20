"""w11_t3.py -- MC0 without pair-maximality: random hard graph (independent pair bits), c<=2, min certificate l=0 must be <= 6."""
import sys, random, itertools, time
from collections import Counter
from w11_lib import *
from w11_gen import *
rnd = random.Random(int(sys.argv[1])); samples = int(sys.argv[2]); m = int(sys.argv[3]); npend = int(sys.argv[4])
sizes = Counter(); t0 = time.time()
for s in range(samples):
    nv, H = random_members_c2(m, rnd)
    n = nv + npend
    E = set()
    soft = set()
    for (a, b, c) in H:
        E.add(key(a, b)); E.add(key(b, c)); soft |= {key(a, b), key(b, c), key(a, c)}
    p = rnd.choice([0.1, 0.2, 0.3])
    for u, v in itertools.combinations(range(n), 2):
        if key(u, v) not in soft and rnd.random() < p: E.add(key(u, v))
    I = Inst(n, E, H)
    if not (I.valid() and I.cmax() <= 2): continue
    sz, T = min_certificate(I, 0, max_size=30, workers=1)
    sizes[sz] += 1
    if sz is not None and sz > 6: print('VIOLATION', sz, nv, H, sorted(E))
print(sorted(sizes.items(), key=str), round(time.time() - t0))
