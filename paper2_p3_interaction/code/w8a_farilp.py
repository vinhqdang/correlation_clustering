import time
from w8a_lib import *
from w8a_ilp2 import joint2
for m in (5, 6, 7, 8):
    n = 2 * m + 1
    g = G8(n, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)])
    opt, _ = opt_ilp(g); nu, _ = nu_mod(g); t = time.time()
    print('FAR', m, 'n', n, 'OPT', opt, 'nu', nu, 'min #hard edits (both covered, soft-dist>2) over max packings x optima =', joint2(g, opt, nu, 2, tl=120), '%.0fs' % (time.time() - t), flush=True)
