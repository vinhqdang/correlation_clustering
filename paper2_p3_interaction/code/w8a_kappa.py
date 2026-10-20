"""kappa(G,H): the minimum number of NEW members |H* \\ H| over packings H* with |H*| >= |H| such that some optimum has no 'far' hard edit
(far := both ends covered by H*, I(H*)-distance between the members of the ends >= 2, i.e. not T1-adjacent).  Input H: random pair-maximal packing with s(H)>0."""
import random, sys, time, collections
from w8a_lib import *
from w8a_ilp import rand_graph
def far_count(H, optima, I=None):
    I = Igraph(H) if I is None else I
    V = set(v for t in H for v in t)
    best = None
    for S in optima:
        c = 0
        for (a, b) in hard_edits(H, S):
            if a in V and b in V:
                d = member_dist(H, a, b, I)
                if d >= 2: c += 1
        if best is None or c < best: best = c
        if best == 0: break
    return best
def random_pm_packing(g, rng):
    P = g.p3s()[:]; rng.shuffle(P)
    H = []; used = set()
    for t in P:
        if not any(p in used for p in mp(t)):
            H.append(t); used.update(mp(t))
    return H
def kappa(g, H, optima, maxpk=200000):
    packs = enum_packings(g, minsize=len(H), maxcount=maxpk)
    Hs = set(H); best = None; bestH = None
    for H2 in packs:
        H2 = list(H2)
        new = sum(1 for t in H2 if t not in Hs)
        if best is not None and new >= best: continue
        if not pair_maximal(g, H2): continue
        if far_count(H2, optima) == 0: best, bestH = new, H2
    return best, bestH, len(packs)
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); nlo, nhi = int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed); st = collections.Counter(); ex = []
    t0 = time.time()
    for it in range(N):
        n = rng.randint(nlo, nhi); mode = rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic', 'planted'])
        E = rand_graph(n, mode, rng); g = G8(n, E)
        gg = nx.Graph(list(g.E))
        if len(g.E) < 3 or gg.number_of_nodes() < n or not nx.is_connected(gg): continue
        H = random_pm_packing(g, rng)
        opt, sols = g.opt(limit=3000)
        optima = [edits_of_partition(g, part_from_mask(s, n)) for s in sols]
        l = opt - len(H)
        sH = far_count(H, optima)
        if sH == 0: st[('s(H)=0', min(l, 3))] += 1; continue
        k, H2, npk = kappa(g, H, optima)
        st[('s(H)>0', 'l=%d' % l, 'kappa=%s' % k)] += 1
        if k is None or k >= 2: ex.append((n, sorted(g.E), H, l, k, H2))
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(st.items(), key=str): print('  ', k, v)
    for e in ex[:6]: print('  EX', e)
