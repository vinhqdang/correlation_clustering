"""Alternating chain of K4-units: odd units are the flip gadget G* (bridge member with anchor x), even units E are K4 + two slots without bridge.
Consecutive units share the slot vertex (a of unit i = b of unit i+1)."""
import sys, itertools, numpy as np
from w6c_lib import *
from w6c_ringK4 import face
GEc = [(0,1),(0,4),(0,5),(1,4),(1,5),(4,5),(3,0),(3,1),(6,0),(6,5)]
Mc = [(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
def build(units, ring=True):
    """units: string over 'O','E'. slot vertices w_i between unit i-1 and i (cyclic if ring)."""
    U = len(units); idx = {}
    def V(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    members = []; edges = set(); xs = []
    for i, t in enumerate(units):
        loc = lambda v: V(('K', i, v))
        b = V(('w', i)); a = V(('w', (i + 1) % U)) if ring else V(('w', i + 1))
        mp = {0: loc(0), 1: loc(1), 4: loc(4), 5: loc(5), 3: a, 6: b}
        for (u, v) in GEc: edges.add(key(mp[u], mp[v]))
        for m in Mc: members.append(tuple(mp[t_] for t_ in m))
        if t == 'O':
            x = loc(7); edges |= {key(mp[3], mp[6]), key(mp[3], x)}
            members.append((mp[6], mp[3], x)); xs.append(x)
    A = adj_from_edges(len(idx), edges)
    return A, members, xs
if __name__ == '__main__':
    for units in sys.argv[1:]:
        A, H, xs = build(units)
        ok, why = check_packing(A, H)
        if not ok: print(units, 'invalid', why); continue
        o = ilp_opt(A)[0]
        print(units, 'n', len(A), '|H|', len(H), 'OPT excess', o - len(H), 'pairmax', is_pair_maximal(A, H), flush=True)
        if o != len(H) or not xs: continue
        N = len(xs)
        fam = [S for r in range(N + 1) for S in itertools.combinations(range(N), r) if face(A, [xs[i] for i in S], o)]
        maxf = [S for S in fam if not any(set(S) < set(T) for T in fam)]
        print('   maximal faces', maxf, flush=True)
