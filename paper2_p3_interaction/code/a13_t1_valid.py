"""a13_t1_valid.py: G_N is a valid input; exact OPT by subset DP."""
import sys, time
from a13_lib import *
for N in range(1, 4):
    inst = build_GN(N)
    ok, c = validate(inst)
    H = len(inst['members'])
    adj = IH(inst)
    ne = sum(len(s) for s in adj.values())//2
    print('N=%d n=%d |H|=%d valid(induced, mod-disjoint, pair-maximal)=%s c=%d  |E(I(H))|=%d (expect 4N-1=%d)' % (N, inst['n'], H, ok, c, ne, 4*N-1), flush=True)
    if inst['n'] <= 21:
        t = time.time()
        o = opt_exact(inst)
        print('   OPT(G_N)=%d = |H|+%d   (%.1fs)' % (o, o-H, time.time()-t), flush=True)
