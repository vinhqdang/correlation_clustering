"""Augmenting chain: path 0..2m+2 (2m+2 edges), H0 = shifted tiling (1,2,3),(3,4,5),...,(2m-1,2m,2m+1) of size m is pair-maximal but nu = m+1.
Smallest r such that an r-exchange improves H0."""
from w8a_lib import *
from w8a_aug import improve
for m in range(2, 7):
    n = 2 * m + 3
    g = G8(n, [(i, i + 1) for i in range(n - 1)])
    H0 = [(2 * i + 1, 2 * i + 2, 2 * i + 3) for i in range(m)]
    nu, _ = nu_mod(g); opt, _ = opt_ilp(g)
    r = None
    for rr in range(1, m + 1):
        if improve(g, H0, rr): r = rr; break
    print('m=%d n=%d |H0|=%d pairmax=%s nu=%d OPT=%d  smallest improving exchange size r = %s' % (m, n, len(H0), pair_maximal(g, H0), nu, opt, r))
