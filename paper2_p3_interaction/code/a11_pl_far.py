"""Audit 11: PL DP on FAR_m with path tiling H0 and exchanged H*; and Sp(triangle) with original H (DP only, no brute force)."""
import sys, time
from a11_lib import *
from a11_pl import *
from a11_far import far

def run(g, H, l, order, brute=True):
    Ia = I_graph(H)
    bags, edges = elimination_decomposition(Ia, order)
    inst = Inst(g, H, l, bags, edges)
    t = time.time()
    val, cf = dp(inst)
    r = dict(dp=val, cov_fail=cf, tw=max(len(b) for b in bags) - 1, secs=round(time.time() - t, 1))
    if brute:
        cb = brute_class_opt(inst)
        r['class'] = cb[0]
    return r

if __name__ == '__main__':
    for m in (3, 4, 5):
        g, H0, Hs = far(m)
        o = opt_ilp(g)[0]
        print('FAR m=%d OPT=%d' % (m, o))
        print('  H0 (path tiling), l=1:', run(g, H0, 1, list(range(m)), brute=g.n <= 11))
        print('  H* (exchanged),   l=1:', run(g, Hs, 1, list(range(m)), brute=g.n <= 11))
