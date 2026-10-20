"""w11_t1.py -- test of Theorem MC0 (c<=2, l=0): (i) Sigma_0 (remote pairs dropped) closed under M-majority unless K4 component;
(ii) min certificate <= 6 (<=12 with K4 component)."""
import sys, time, random, itertools
import numpy as np
from collections import Counter
from w11_lib import *
from w11_gen import *

def sigma0(I):
    """ZE state vectors (0=L,1=R,2=M) satisfying all triples that contain no remote pair; local hard pairs fixed at G-status."""
    soft, lochard, remote = I.classify()
    rem = set(remote)
    X = {0: (0, 1, 0), 1: (1, 0, 0), 2: (1, 1, 1)}
    m = len(I.H)
    cons = []   # (list of (member,k) or const, )
    for t in itertools.combinations(range(I.n), 3):
        ps = [key(t[0], t[1]), key(t[1], t[2]), key(t[0], t[2])]
        if any(p in rem for p in ps): continue
        terms = []; const = 0; owners = set()
        for p in ps:
            if p in I.owner:
                i, k = I.owner[p]; terms.append((i, k)); owners.add(i)
            else:
                const += 1 if p in I.E else 0
        cons.append((terms, const, owners))
    by_last = [[] for _ in range(m)]
    static_bad = False
    for terms, const, owners in cons:
        if not owners:
            if const == 2: static_bad = True
            continue
        by_last[max(owners)].append((terms, const))
    if static_bad: return []
    res = []; st = [0] * m
    def ok(k):
        for terms, const in by_last[k]:
            s = const + sum(X[st[i]][kk] for i, kk in terms)
            if s == 2: return False
        return True
    def rec(k):
        if k == m: res.append(tuple(st)); return
        for s in (0, 1, 2):
            st[k] = s
            if ok(k): rec(k + 1)
    rec(0)
    return res

def Mmaj(x, y, z):
    if x == y or x == z: return x
    if y == z: return y
    return 2

def closed_M(S):
    Sset = set(S)
    for u, v, w in itertools.product(S, repeat=3):
        if tuple(Mmaj(a, b, c) for a, b, c in zip(u, v, w)) not in Sset: return False
    return True

def has_k4(I):
    g = I.IH()
    import networkx as nx
    for comp in nx.connected_components(g):
        if len(comp) == 4 and g.subgraph(comp).number_of_edges() == 6: return True
    return False

if __name__ == '__main__':
    m = int(sys.argv[1]); npend = int(sys.argv[2]); samples = int(sys.argv[3]); seed = int(sys.argv[4])
    rnd = random.Random(seed)
    stat = Counter(); sizes = Counter(); t0 = time.time()
    for s in range(samples):
        nv, H = random_members_c2(m, rnd)
        n = nv + npend
        k = rnd.choice([2, 3, 4, 6, 9])
        cl = [rnd.randrange(-1, k) if rnd.random() < 0.8 else -1 for _ in range(n)]
        I = make_inst(nv, H, npend, cl)
        if not (I.valid() and I.pair_maximal() and I.cmax() <= 2): stat['skip'] += 1; continue
        if n > 22: stat['big'] += 1; continue
        k4 = has_k4(I)
        S = sigma0(I)
        if 0 < len(S) <= 400:
            c = closed_M(S)
            stat[('closed' if c else 'NOTclosed', 'K4' if k4 else 'noK4')] += 1
        sz, T = min_certificate(I, 0, max_size=30)
        sizes[(k4, sz)] += 1
        if sz is not None and sz > (12 if k4 else 6): print('VIOLATION', sz, nv, H, cl, flush=True)
    print(stat); print(sorted(sizes.items(), key=str), round(time.time() - t0), 's')
