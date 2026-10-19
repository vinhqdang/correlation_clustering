"""Verified semantics of a clean tunnel (hard K2 between two far members).
Member P=(x,y,z), P'=(x',y',z').  Local forcing: z has a hard neighbour t (t-y non-edge, hard), so P cannot keep yz at zero excess:
its only zero-excess state is R (delete yz, keep xy), which makes x non-isolated.  Same for P'.
(a) no tunnel:            OPT = |H| = 2
(b) private stubs at x,x': OPT = |H| + 2   (each deviation pays its own cut)
(c) tunnel x-x':          OPT = |H| + 1   (one cut releases both)   <- cost sharing
(d) inactive anchors: members whose default state is L cost nothing even with the tunnel."""
import numpy as np, itertools
from w6c_lib import *
def inst(tunnel, stubs, forced=True):
    # vertices: x0 y1 z2 t3 | x4 y5 z6 t7 | stubs 8,9
    E = [(0,1),(1,2),(4,5),(5,6)]
    if forced: E += [(2,3),(6,7)]
    n = 8
    if tunnel: E += [(0,4)]
    if stubs:
        n = 10; E += [(0,8),(4,9)]
    A = adj_from_edges(n, E)
    H = [(0,1,2),(4,5,6)]
    return A, H
for name, (tun, st, forced) in {'no tunnel (forced R)':(False,False,True),'private stubs':(False,True,True),'shared tunnel':(True,False,True),
                                'shared tunnel, unforced (default L ok)':(True,False,False)}.items():
    A, H = inst(tun, st, forced)
    ok, why = check_packing(A, H); assert ok, why
    opt = ilp_opt(A)[0]; opt2 = brute_opt(A)[0] if len(A) <= 9 else None
    print(f'{name:42s} |H|={len(H)} OPT={opt} (brute {opt2}) excess={opt-len(H)}  pair-maximal={is_pair_maximal(A,H)}')
