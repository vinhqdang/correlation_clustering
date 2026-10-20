import sys, itertools, time, collections
from w8a_lib import *
from w8a_rand import rho
def theta(p, q, r):
    # vertices 0 and 1 are the poles; three paths with p,q,r edges
    E = []; n = 2
    for L in (p, q, r):
        prev = 0
        for i in range(L - 1):
            E.append((prev, n)); prev = n; n += 1
        E.append((prev, 1))
    return n, E
res = []
for (p, q, r) in [(1,3,4),(1,4,4),(1,4,5),(1,5,5),(1,3,6),(1,5,6),(2,3,4),(2,4,5),(3,3,3),(3,4,4),(1,6,6),(1,2,7),(1,7,7),(1,5,7)]:
    n, E = theta(p, q, r)
    if n > 15: continue
    g = G8(n, E); t = time.time()
    opt, nu, npk, r_, H = rho(g, maxpk=20000)
    print('theta', (p, q, r), 'n', n, 'OPT', opt, 'nu', nu, '#max packings', npk, 'rho(min over maxpk,optima of max hard-edit dist)', r_, '%.1fs' % (time.time() - t), flush=True)
