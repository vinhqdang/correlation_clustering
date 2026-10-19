"""Ring of N 'K4-joiner' gadgets (derived from the 5-member flip gadget found by w6c_flip4.py):
clique K_j={0_j,1_j,4_j,5_j} (hard edges 0-4, 1-5 and soft edges 01,05,14,45), joiner vertices w_j (w_j adjacent to 0_j,1_j and to 0_{j+1},5_{j+1}),
anchor x_j with member (w_{j-1}, w_j, x_j), and members (4,1,w_j), (5,0,w_j), (4,5,w_{j-1}), (w_{j-1},0,1) for K_j.
Claims to test: tight (OPT=|H|); odd N: the full set of anchors {x_j} cannot all be singletons, every N-1 subset can (a lock core of size N)."""
import sys, itertools, numpy as np
from w6c_lib import *
def build(N, twist=False):
    idx = {}
    def v(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    members = []; edges = set()
    for j in range(N):
        K0, K1, K4, K5 = (v(('K', j, t)) for t in (0, 1, 4, 5)); w = v(('w', j)); wp = v(('w', (j - 1) % N)); x = v(('x', j))
        for a, b in ((K0, K1), (K0, K4), (K0, K5), (K1, K4), (K1, K5), (K4, K5)): edges.add(key(a, b))
        for a, b in ((w, K0), (w, K1), (wp, K0), (wp, K5), (wp, w), (w, x)): edges.add(key(a, b))
        members += [(wp, w, x), (K4, K1, w), (K5, K0, w), (K4, K5, wp), (wp, K0, K1)]
    n = len(idx)
    A = adj_from_edges(n, edges)
    return A, members, idx
def face(A, S, opt):
    fixed = {}
    n = len(A)
    for a in S:
        for t in range(n):
            if t != a: fixed[key(a, t)] = 0
    r = ilp_opt(A, fixed=fixed)
    return r is not None and r[0] == opt
if __name__ == '__main__':
    for N in map(int, sys.argv[1:]):
        A, H, idx = build(N)
        ok, why = check_packing(A, H)
        opt = ilp_opt(A)[0]
        xs = [idx[('x', j)] for j in range(N)]
        print(f'N={N} n={len(A)} |H|={len(H)} valid packing={ok} {why} pair-maximal={is_pair_maximal(A,H)} OPT={opt}')
        full = face(A, xs, opt); subs = [face(A, [xs[i] for i in range(N) if i != j], opt) for j in range(N)]
        print('   all N anchors singleton possible:', full, '; each (N-1)-subset possible:', subs)
        ih = ih_graph(H)
        print('   tw(I(H)) <= ', tw_heur(ih), ' max members per vertex c =', max(sum(1 for m in H if v in m) for v in range(len(A))))
