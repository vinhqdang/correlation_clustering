"""w11_tight.py -- an explicit c=2, pair-maximal instance whose minimum l=0 certificate has exactly 4 remote pairs (so the bound 4 of Theorem MC0 is attained)."""
from w11_lib import *
from w11_gen import make_inst
H=[(1, 0, 3), (5, 0, 4), (2, 6, 1), (2, 7, 8)]
cl=[-1, 1, -1, 0, 1, 0, -1, -1, -1, 1, -1, -1, 0]
I=make_inst(9,H,4,cl)
print('n',I.n,'valid',I.valid(),'pair-maximal',I.pair_maximal(),'c',I.cmax(),'soft triangles',len(I.soft_triangles()),'I(H) edges',sorted(I.IH().edges()))
print('edges of G', sorted(I.E))
sz,T=min_certificate(I,0)
print('minimum certificate size',sz,T)
O=FPrime(I,len(I.H))
print('F(empty) feasible at |H|:',O.feasible([]) is not None,' F(N) feasible:',O.feasible(I.remote()) is not None)
import itertools
rem=I.remote()
print('every 3-subset of remote pairs feasible?', all(O.feasible(S) is not None for S in itertools.combinations(rem,3)), 'number of remote pairs',len(rem))
