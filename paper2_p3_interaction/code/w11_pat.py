"""w11_pat.py -- K3 chain with arbitrary pendant attachment patterns: min certificate (l=0) as a function of patterns."""
import sys, time, itertools, random
from w11_lib import *

NAMES = 'pqrcab'

def build(N, mid, start, end, hard_cq=True, extra_mid=()):
    """mid: tuple of subsets (as strings of NAMES) -- pendants attached to every block j with that neighbourhood
       start/end: list of subsets for pendants on block 0 / block N-1 (end uses names relative to block N-1)."""
    vid = {}
    def V(j, name):
        if name == 'a' and j > 0:
            return V(j - 1, 'b')
        if (j, name) not in vid:
            vid[(j, name)] = len(vid)
        return vid[(j, name)]
    H = []; E = set()
    for j in range(N):
        p, q, r, c, b, a = (V(j, x) for x in 'pqrcba')
        H += [(q, p, a), (q, r, b), (p, c, r)]
        E |= {key(q, p), key(p, a), key(q, r), key(r, b), key(p, c), key(c, r)}
        if hard_cq: E.add(key(c, q))
    nv = len(vid)
    def pend(j, S):
        nonlocal nv
        w = nv; nv += 1
        for x in S:
            E.add(key(w, V(j, x)))
    for j in range(N):
        for S in mid: pend(j, S)
    for S in start: pend(0, S)
    for S in end: pend(N - 1, S)
    return Inst(nv, E, H)

if __name__ == '__main__':
    N = int(sys.argv[1]); samples = int(sys.argv[2]); seed = int(sys.argv[3])
    rnd = random.Random(seed)
    subsets = [''.join(c) for k in range(0, 4) for c in itertools.combinations(NAMES, k)]
    res = {}
    t0 = time.time()
    tried = 0; valid = 0; cert = 0
    best = []
    for s in range(samples):
        mid = (rnd.choice(subsets),) if rnd.random() < 0.9 else (rnd.choice(subsets), rnd.choice(subsets))
        st = (rnd.choice(subsets),)
        en = (rnd.choice(subsets),)
        I = build(N, mid, st, en)
        tried += 1
        if not (I.valid() and I.pair_maximal()):
            continue
        valid += 1
        sz, T = min_certificate(I, 0, max_size=14)
        if sz is None: continue
        cert += 1
        res[(mid, st, en)] = sz
        best.append((sz, mid, st, en))
    best.sort(reverse=True)
    from collections import Counter
    print('N', N, 'tried', tried, 'valid/pm', valid, 'with certificate', cert, round(time.time() - t0), 's')
    print(Counter(b[0] for b in best))
    for b in best[:15]: print(b)
