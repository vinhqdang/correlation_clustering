"""Check the stated n=8 counterexample, and search n=6,7 for a smaller 'forced T2' instance
(every optimum has an inserted pair (a,b) with Q(a) and Q(b) disjoint, i.e. no common soft neighbour)."""
import sys, random, itertools
from a5_lib import *
H_G = []
def nonlocal_pairs(adj, pk, lab):
    n = len(adj); VH = {v for t in H_G for v in t}; Q = {v: {x for x in range(n) if frozenset((v, x)) in pk} for v in range(n)}
    return [(a, b) for a in range(n) for b in range(a+1, n) if lab[a] == lab[b] and not adj[a][b]
            and frozenset((a, b)) not in pk and not (Q[a] & Q[b]) and a in VH and b in VH]
def stated():
    E = [(0,2),(0,6),(1,4),(1,6),(1,7),(2,3),(2,5),(3,7),(4,5),(4,6),(4,7),(5,7)]
    global H_G
    adj = mk(8, E); H = [(0,6,4),(3,7,4),(2,5,4),(0,2,3),(6,1,7)]
    print('is packing', is_packing(adj, H), 'pair-maximal', pair_maximal(adj, H))
    H_G = H
    opt, opts = optima(adj); pk = packed(H)
    print('OPT', opt, '|H|', len(H), 'num optima', len(opts), opts)
    for lab in opts:
        print(lab, 'non-local inserted pairs', nonlocal_pairs(adj, pk, lab), 'U', U_set(adj, pk, lab))
    # is the stated 'common neighbours 4,7 both hard+soft' right?
    print('14 packed?', frozenset((1,4)) in pk, '45 packed?', frozenset((4,5)) in pk, '17', frozenset((1,7)) in pk, '57', frozenset((5,7)) in pk)
def search(n, N, seed):
    global H_G
    found = 0; tot = 0
    for sd in range(N):
        rng = random.Random(f"a5s-{n}-{seed}-{sd}")
        p = rng.uniform(.2, .6)
        adj = mk(n, [(i, j) for i in range(n) for j in range(i+1, n) if rng.random() < p])
        H = rand_pm(adj, rng, rng.choice([0, .3]))
        if not H: continue
        tot += 1
        H_G = H
        pk = packed(H); opt, opts = optima(adj)
        if all(nonlocal_pairs(adj, pk, lab) for lab in opts):
            found += 1
            if found <= 2: print('FORCED n=%d' % n, [tuple(sorted(e)) for e in itertools.combinations(range(n), 2) if adj[e[0]][e[1]]], H, opts)
    print(n, 'instances', tot, 'forced', found)
if __name__ == '__main__':
    stated()
    if len(sys.argv) > 1:
        search(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 0)
