"""Canonical tiling H = {(2i,2i+1,2i+2)} of the path 0..2m plus a chord (i,j) (odd,odd).  Minimum number of NEW members |H* \\ H| over packings H* (|H*|>=|H|)
such that some optimum has no far hard edit (both ends covered, I(H*)-distance >= 2)."""
import sys, time
from w8a_lib import *
from w8a_kappa import far_count
for m in (5, 6, 7):
    n = 2 * m + 1
    H0 = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    out = []
    for i in range(1, n, 2):
        for j in range(i + 4, n, 2):
            g = G8(n, [(k, k + 1) for k in range(n - 1)] + [(i, j)])
            if not valid_packing(g, H0) or not pair_maximal(g, H0): continue
            opt, sols = g.opt(limit=5000)
            optima = [edits_of_partition(g, part_from_mask(s, n)) for s in sols]
            f0 = far_count(H0, optima)
            packs = enum_packings(g, minsize=m, maxcount=200000)
            best = None
            Hs = set(H0)
            for P in packs:
                P = list(P)
                new = sum(1 for t in P if t not in Hs)
                if best is not None and new >= best: continue
                if pair_maximal(g, P) and far_count(P, optima) == 0: best = new
            out.append(((i, j), 'OPT-|H|=%d' % (opt - m), 'far(H0)=%s' % f0, 'exchanges needed=%s' % best))
    print('m=%d' % m)
    for o in out: print('   ', o)
