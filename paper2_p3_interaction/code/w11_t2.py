"""w11_t2.py -- K4 component test: octahedron of 4 members (+ optional extra members/pendants), random hard clusters; Sigma_0 closure under M-majority; min certificate."""
import sys, random, itertools, time
from collections import Counter
from w11_lib import *
from w11_gen import make_inst
from w11_t1 import sigma0, closed_M
rnd = random.Random(int(sys.argv[1])); samples = int(sys.argv[2]); npend = int(sys.argv[3])
mem = [[0,1,2],[0,3,4],[1,3,5],[2,4,5]]
cnt = Counter(); sizes = Counter(); t0 = time.time()
for s in range(samples):
    H = []
    for vs in mem:
        t = vs[:]; rnd.shuffle(t); H.append(tuple(t))
    nv = 6; n = nv + npend
    k = rnd.choice([2, 3, 4]); cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
    I = make_inst(nv, H, npend, cl)
    if not (I.valid() and I.pair_maximal()): cnt['skip'] += 1; continue
    S = sigma0(I)
    if S: cnt['closed' if closed_M(S) else 'NOT closed'] += 1
    else: cnt['empty'] += 1
    sz, T = min_certificate(I, 0, max_size=30, workers=1)
    sizes[sz] += 1
print(cnt, sorted(sizes.items(), key=str), round(time.time() - t0))
