"""Random search: for random connected graphs, does there exist a MAXIMUM packing H* (hence pair-maximal) and an optimum S such that all hard edits are
I(H*)-local (distance <= r)?  Record min over (H*,S) of the max hard-edit distance, rho(G).  Report instances with large rho."""
import random, sys, collections, time
from w8a_lib import *
def rho(g, maxpk=3000, slack=0):
    opt, sols = g.opt(limit=5000)
    parts = [edits_of_partition(g, part_from_mask(s, g.n)) for s in sols]
    nu, _ = nu_mod(g)
    packs = enum_packings(g, minsize=nu - slack, maxcount=maxpk)
    best = None; bestH = None
    for H in packs:
        H = list(H)
        if not pair_maximal(g, H): continue
        I = Igraph(H)
        rb = None
        for S in parts:
            r = 0
            for (a, b) in hard_edits(H, S):
                d = member_dist(H, a, b, I)
                r = max(r, 0 if d == -1 else d)
                if rb is not None and r >= rb: break
            if rb is None or r < rb: rb = r
            if rb == 0: break
        if best is None or rb < best: best, bestH = rb, H
        if best == 0: break
    return opt, nu, len(packs), best, bestH
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); random.seed(seed)
    stat = collections.Counter(); worst = []
    t0 = time.time()
    for it in range(N):
        n = random.randint(7, 11)
        p = random.choice([0.2, 0.25, 0.3, 0.4])
        E = [(a, b) for a in range(n) for b in range(a + 1, n) if random.random() < p]
        g = G8(n, E)
        if not nx.is_connected(nx.Graph(E + [(i, i) for i in range(n)])) : continue
        opt, nu, npk, r, H = rho(g)
        stat[(opt - nu, r)] += 1
        if r is not None and r >= 2: worst.append((r, n, sorted(g.E), H))
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(stat.items(), key=lambda x: (x[0][0], -1 if x[0][1] is None else x[0][1])): print('  (OPT-nu, rho) =', k, v)
    for w in worst[:5]: print('  WORST', w)
