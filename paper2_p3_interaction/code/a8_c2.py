"""Audit 8: Conjecture C2 (c<=2 => ZE set closed under M-majority: majority if two agree, M if all distinct) -- independent search."""
import sys, random, itertools
import numpy as np
from a8_lib import *
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
stats = {}; bad = []
def bump(k, v=1): stats[k] = stats.get(k, 0) + v
def mmaj(a, b, c):
    out = []
    for x, y, z in zip(a, b, c):
        out.append(x if (x == y or x == z) else y if y == z else 'M')
    return ''.join(out)
done = 0; tried = 0
while done < N and tried < 200 * N:
    tried += 1
    n = rng.choice([7, 8, 9, 10])
    A = gen_planted(n, rng, nedit=rng.randint(2, 7), bigcluster=rng.random() < 0.3) if rng.random() < 0.8 else gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid() or I.opt() != len(H): continue
    c = max(sum(1 for t in H if v in t) for v in range(n))
    if c > 2: continue
    ZE = I.zero_excess(); vecs = sorted({I.state(k) for k in ZE})
    if len(vecs) < 3 or len(vecs) > 60: continue
    done += 1; bump('inst'); vs = set(vecs)
    ok = all(mmaj(a, b, cc) in vs for a in vecs for b in vecs for cc in vecs)
    bump('closed' if ok else 'NOT_closed')
    if not ok: bad.append((I.A.tolist(), H, vecs))
print(stats)
for b in bad[:3]: print(b)
