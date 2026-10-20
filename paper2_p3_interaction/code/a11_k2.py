"""Audit 11: Prop 5.1: K_2 v 3K_1; and PL DP on Sp(triangle) with the ORIGINAL packing (hub edit)."""
import sys, time
from a11_lib import *
g = Gr(5, [(3, 4)] + [(x, h) for x in (0, 1, 2) for h in (3, 4)])
nu, MP = max_packings(g)
o, opts = opt_and_optima(g)
print('K2 v 3K1: induced P3s =', g.p3(), ' nu =', nu, ' #max packings =', len(MP), ' OPT =', o, ' #optima =', len(opts))
allcov = all({v for t in M for v in t} == set(range(5)) for M in MP)
print('every maximum packing covers all 5 vertices:', allcov)
# for every max packing and every optimum: set of hard edits; is there an optimum with no hard edit / with only uncovered-endpoint hard edits?
res = []
for M in MP:
    own = soft_owner(M)
    hs = [sorted(hard_of(M, S)) for S in opts]
    res.append(min(len(h) for h in hs))
print('min #hard edits over optima, per max packing:', res)
print('optima:', [sorted(S) for S in opts])
if len(sys.argv) > 1:
    from a11_pl_far import run
    from a11_spine import build
