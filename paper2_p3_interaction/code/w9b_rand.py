"""random instances (c<=cmax) with feasible Sigma; record max list core."""
import sys, random, itertools, time
from collections import Counter
from w9b_lib import *
from w9b_sigma import *

def gen(rng, m, cmax, nv, p_hard):
    """m members on nv vertices: pick triples (a,b,c) with disjoint pair sets, c<=cmax; then random hard edges."""
    for _ in range(200):
        H = []; usedp = set(); cnt = [0]*nv; ok = True
        tries = 0
        while len(H) < m and tries < 500:
            tries += 1
            t = tuple(rng.sample(range(nv), 3))
            ps = [key(t[0], t[1]), key(t[1], t[2]), key(t[0], t[2])]
            if any(p in usedp for p in ps): continue
            if any(cnt[v] >= cmax for v in t): continue
            H.append(t); usedp |= set(ps)
            for v in t: cnt[v] += 1
        if len(H) < m: continue
        E = set()
        for (a, b, c) in H: E.add(key(a, b)); E.add(key(b, c))
        for u, v in itertools.combinations(range(nv), 2):
            if key(u, v) not in usedp and rng.random() < p_hard: E.add((u, v))
        I = Inst(nv, E, H)
        if I.valid(): return I
    return None

if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); m = int(sys.argv[3]); nv = int(sys.argv[4]); cmax = int(sys.argv[5])
    rng = random.Random(seed); c1 = Counter(); c2 = Counter(); best = (0, None)
    t0 = time.time()
    for it in range(N):
        p = rng.choice([0.15, 0.25, 0.35, 0.5])
        I = gen(rng, m, cmax, nv, p)
        if I is None: continue
        Sg = sigma(I, drop_far=False)
        if len(Sg) < 2: continue
        k1 = max_core(Sg, singleton_only=True); k2 = max_core(Sg)
        c1[k1] += 1; c2[k2] += 1
        if k1 > best[0]: best = (k1, (I.n, sorted(I.E), I.H))
    print('m', m, 'nv', nv, 'c<=', cmax, 'k1 dist', dict(sorted(c1.items())), 'k2 dist', dict(sorted(c2.items())), 'time', round(time.time()-t0))
    if best[0] >= 3: print('best', best)
