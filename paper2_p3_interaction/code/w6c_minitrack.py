"""Flip gadget G* (found by w6c_flip4.py) as a one-gadget 'track'.
G*: K4 {0,1,4,5} with hard edges 0-4,1-5; a=3, b=6, anchor x=7; members (6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1).
Boundary forcing: a helper t adjacent (hard) to 0,1,4,5,6 but not to 3 forces b into the K4 cluster; with x locked (singleton)
the bridge a~b would drag a into the cluster -> impossible at zero excess; releasing x (x joins a) is the only free way out.
Tests: OPT without stub, with private stub (tunnel to nothing), and the optimum partitions."""
import numpy as np, itertools
from w6c_lib import *
def inst(with_t, with_stub):
    E = [(0,1),(0,4),(0,5),(1,4),(1,5),(4,5),(3,0),(3,1),(6,0),(6,5),(3,6),(3,7)]
    n = 8
    if with_t: E += [(8,v) for v in (0,1,4,5,6)]; n = 9
    if with_stub: E += [(7, n)]; n += 1
    # vertex 2 is unused in G*: relabel by keeping it as an isolated vertex
    return adj_from_edges(n, E)
H = [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
for wt, ws in ((False, False), (True, False), (True, True)):
    A = inst(wt, ws)
    ok, why = check_packing(A, H)
    o = ilp_opt(A)[0]
    print(f't={wt} stub={ws}: n={len(A)} valid packing={ok} pair-maximal={is_pair_maximal(A,H)} OPT={o} excess={o-len(H)}')
    if len(A) <= 9:
        opt, labs, sm, cost, P = brute_opt(A)
        sols = set()
        for l, c in zip(labs, cost):
            if c == opt:
                cl = {}
                for v, t in enumerate(l): cl.setdefault(t, []).append(v)
                sols.add(tuple(sorted(tuple(b) for b in cl.values() if len(b) > 1)))
        for s in sorted(sols): print('    ', s)
