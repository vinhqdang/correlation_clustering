"""Smallest n with a pair-maximal H and unique/all optima forcing an inserted pair (both ends in V(H)) with no common soft neighbour.
Exhaustive over all graphs on n<=5 vertices and ALL pair-maximal packings; explicit n=6 example."""
import itertools
from a5_lib import *
def all_pm_packings(adj):
    ts = p3s(adj); res = []
    def rec(i, H, used):
        if i == len(ts):
            if pair_maximal(adj, H): res.append(list(H))
            return
        p = pairs3(ts[i])
        if not p & used: rec(i+1, H+[ts[i]], used | p)
        rec(i+1, H, used)
    rec(0, [], set()); return res
def forced(adj, H):
    pk = packed(H); n = len(adj); VH = {v for t in H for v in t}
    Q = {v: {x for x in range(n) if frozenset((v, x)) in pk} for v in range(n)}
    opt, opts = optima(adj)
    def bad(lab):
        return [(a, b) for a in range(n) for b in range(a+1, n) if lab[a] == lab[b] and not adj[a][b] and frozenset((a,b)) not in pk
                and a in VH and b in VH and not (Q[a] & Q[b])]
    return all(bad(l) for l in opts), opts, [bad(l) for l in opts]
if __name__ == '__main__':
    for n in (4, 5):
        P = list(itertools.combinations(range(n), 2)); cnt = 0
        for m in range(1 << len(P)):
            adj = mk(n, [P[i] for i in range(len(P)) if m >> i & 1])
            for H in all_pm_packings(adj):
                if H and forced(adj, H)[0]: cnt += 1; print('n=%d FORCED' % n, m, H)
        print('n', n, 'forced configs (labelled graphs x all pair-maximal H):', cnt)
    E = [(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,5),(2,5),(3,5),(4,5)]; H = [(3,5,4),(0,2,1)]
    adj = mk(6, E); print('n=6 example: packing', is_packing(adj, H), 'pair-maximal', pair_maximal(adj, H), forced(adj, H))
