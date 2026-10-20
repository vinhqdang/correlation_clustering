"""hill climbing over (member tree structure fixed, hard clusters mutable, pendants) maximising the minimum certificate size."""
import sys, time, random
from w11_lib import *
from w11_gen import *

def score(I, l, cap):
    sz, T = min_certificate(I, l, max_size=cap)
    return -1 if sz is None else sz

def run(m, npend, l, seed, iters, cap=24, restarts=1, verbose=True):
    rnd = random.Random(seed)
    best_overall = (-1, None)
    for r in range(restarts):
        # init: random member tree, random clusters; retry until some certificate
        while True:
            nv, H = random_members(m, rnd)
            n = nv + npend
            k = rnd.choice([3, 4, 6])
            cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
            I = make_inst(nv, H, npend, cl)
            s = score(I, l, cap)
            if s >= 1: break
        cur = s
        for it in range(iters):
            cl2 = cl[:]
            for _ in range(rnd.choice([1, 1, 2, 3])):
                v = rnd.randrange(n)
                cl2[v] = rnd.randrange(-1, 8) if rnd.random() < 0.7 else rnd.choice(cl2)
            H2 = H
            if rnd.random() < 0.15:
                # mutate centre choices
                H2 = [tuple(t) if rnd.random() < 0.8 else (t[1], t[0], t[2]) for t in H]
                if rnd.random()<0.5:
                    H2 = [ (t[0],t[2],t[1]) if rnd.random()<0.2 else t for t in H2]
            I2 = make_inst(nv, H2, npend, cl2)
            if not (I2.valid() and I2.pair_maximal() and I2.cmax() <= 2): continue
            s2 = score(I2, l, cap)
            if s2 >= cur:
                if s2 > cur and verbose: print('  seed', seed, 'it', it, 'score', s2, flush=True)
                cur = s2; cl = cl2; H = H2
        if cur > best_overall[0]: best_overall = (cur, (nv, H, npend, cl))
    return best_overall

if __name__ == '__main__':
    m, npend, l, seed, iters = map(int, sys.argv[1:6])
    t0 = time.time()
    b = run(m, npend, l, seed, iters, restarts=int(sys.argv[6]) if len(sys.argv) > 6 else 1)
    print('BEST', b[0], b[1], round(time.time() - t0), 's')
