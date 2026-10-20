"""w11_inc2.py -- INC test along greedy chains: Phi_0 = empty; Phi_{i+1} = Phi_i u C_i with C_i a minimum increment; record sizes |C_i| until F'(Phi)=F'(N).
Hill-climb instances (c<=2) to maximise max_i |C_i|."""
import sys, time, random
from w11_lib import *
from w11_gen import *

def chain_sizes(I, cap=30, rnd=None):
    Phi = set(); sizes = []
    for step in range(12):
        sz, C, t = min_increment(I, Phi, max_size=cap, workers=2)
        if sz is None: break
        if C is None: return sizes + [sz]
        sizes.append(sz); Phi |= set(C)
    return sizes

if __name__ == '__main__':
    m, npend, seed, iters = map(int, sys.argv[1:5])
    rnd = random.Random(seed)
    while True:
        nv, H = random_members_c2(m, rnd)
        n = nv + npend
        k = rnd.choice([3, 4, 6]); cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
        I = make_inst(nv, H, npend, cl)
        if not (I.valid() and I.pair_maximal() and I.cmax() <= 2): continue
        s = chain_sizes(I)
        if len(s) >= 2: break
    cur = max(s); curS = s
    for it in range(iters):
        cl2 = cl[:]
        for _ in range(rnd.choice([1, 1, 2, 3])):
            v = rnd.randrange(n); cl2[v] = rnd.randrange(-1, 8) if rnd.random() < 0.7 else rnd.choice(cl2)
        I2 = make_inst(nv, H, npend, cl2)
        if not (I2.valid() and I2.pair_maximal()): continue
        s2 = chain_sizes(I2)
        sc = max(s2) if s2 else 0
        if sc >= cur and len(s2) >= 1:
            if sc > cur: print('  it', it, 'max increment', sc, 'chain', s2, flush=True)
            cur = sc; cl = cl2; curS = s2
    print('BEST max increment', cur, 'chain', curS, (nv, H, npend, cl))
