import sys, time, random
from w11_lib import *
from w11_gen import *
from w11_hc import score
def run(m, npend, l, seed, iters, cmax, cap=40):
    rnd = random.Random(seed)
    while True:
        nv, H = random_members_cc(m, rnd, cmax)
        n = nv + npend
        k = rnd.choice([3, 4, 6])
        cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
        I = make_inst(nv, H, npend, cl)
        if not I.valid() or I.cmax() > cmax: continue
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
            H2 = [ (t[1],t[0],t[2]) if rnd.random()<0.3 else ((t[0],t[2],t[1]) if rnd.random()<0.3 else t) for t in H]
        I2 = make_inst(nv, H2, npend, cl2)
        if not (I2.valid() and I2.pair_maximal()): continue
        s2 = score(I2, l, cap)
        if s2 >= cur:
            if s2 > cur: print('  seed', seed, 'it', it, 'score', s2, flush=True)
            cur = s2; cl = cl2; H = H2
    return cur, (nv, H, npend, cl)
if __name__ == '__main__':
    m, npend, l, seed, iters, cmax = map(int, sys.argv[1:7])
    t0 = time.time()
    b = run(m, npend, l, seed, iters, cmax)
    print('BEST', b[0], b[1], round(time.time() - t0), 's')
