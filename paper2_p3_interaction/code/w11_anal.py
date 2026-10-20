import sys, itertools
from w11_lib import *
from w11_gen import *
def show(nv,H,npend,cl,l):
    I=make_inst(nv,H,npend,cl)
    print('n',I.n,'members',H,'valid',I.valid(),'pm',I.pair_maximal(),'c',I.cmax())
    g=I.IH(); print('I(H) edges',sorted(g.edges()))
    sz,T=min_certificate(I,l)
    print('min cert l=%d:'%l,sz,T)
    # which pairs: classify as vertex stars
    from collections import defaultdict
    byv=defaultdict(list)
    for (a,b) in T:
        byv[a].append(b); byv[b].append(a)
    print('vertex incidence in T:',{v:sorted(w) for v,w in byv.items()})
    # unnecessary check and sub-certificates per l-1
    O=FPrime(I,len(I.H)+l)
    # which smaller budgets: min cert for l-1..0
    for ll in range(0,l):
        s2,T2=min_certificate(I,ll); print('  min cert at l=%d:'%ll,s2,T2)
    # soft triangles, clusters
    print('clusters',cl)
    return I,T
if __name__=='__main__':
    nv,H,npend,cl=(11, [(2, 4, 0), (0, 1, 3), (1, 6, 5), (2, 8, 7), (3, 9, 10)], 5, [5, 1, 7, 1, 7, 2, 2, 4, 3, -1, 2, -1, 7, 1, 2, 4])
    show(nv,H,npend,cl,1)
