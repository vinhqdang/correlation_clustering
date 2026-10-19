"""Smallest instance where EVERY optimum has a MERGE-type inserted pair (a,b in different pieces, both in V(H)) with no common soft neighbour.
(Intra-piece inserted pairs are determined by hard connectivity and are not 'merge pointers'.)  Random n=6,7 and 8 (stated example)."""
import sys, random, itertools
from a5_lib import *
def bad_merge(adj, pk, lab, H):
    n = len(adj); VH = {v for t in H for v in t}; pid = pieces(adj, pk, lab)
    Q = {v: {x for x in range(n) if frozenset((v, x)) in pk} for v in range(n)}
    return [(a, b) for a in range(n) for b in range(a+1, n) if lab[a] == lab[b] and not adj[a][b] and frozenset((a, b)) not in pk
            and a in VH and b in VH and pid[a] != pid[b] and not (Q[a] & Q[b])]
def test(adj, H):
    pk = packed(H); opt, opts = optima(adj)
    return all(bad_merge(adj, pk, l, H) for l in opts), opts
if __name__ == '__main__':
    E = [(0,2),(0,6),(1,4),(1,6),(1,7),(2,3),(2,5),(3,7),(4,5),(4,6),(4,7),(5,7)]
    H = [(0,6,4),(3,7,4),(2,5,4),(0,2,3),(6,1,7)]
    adj = mk(8, E); ok, opts = test(adj, H); pk = packed(H)
    print('stated n=8 example forced merge-type nonlocal pair:', ok, opts, 'pieces', pieces(adj, pk, opts[0]))
    n = int(sys.argv[1]); N = int(sys.argv[2]); found = 0; tot = 0
    for sd in range(N):
        rng = random.Random(f"a5m-{n}-{sd}"); p = rng.uniform(.2, .6)
        adj = mk(n, [(i, j) for i in range(n) for j in range(i+1, n) if rng.random() < p])
        H = rand_pm(adj, rng, rng.choice([0, .3]))
        if not H: continue
        tot += 1; ok, opts = test(adj, H)
        if ok:
            found += 1
            if found <= 2: print('FORCED', [e for e in itertools.combinations(range(n), 2) if adj[e[0]][e[1]]], H, opts)
    print(n, tot, 'forced merge-type', found)
