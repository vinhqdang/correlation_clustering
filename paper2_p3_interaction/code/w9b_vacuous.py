"""search pendant adjacency sets U (w adjacent to U) such that the only non-vacuous guard of the pendant is the intended one."""
import itertools, sys
import numpy as np
from w9b_lib import *
from w9b_k3chain import chain
from w9b_octtree import enum_sigma
N=int(sys.argv[1])
I=chain(N,(0,0,0,0))
S=enum_sigma(I,limit=10**6); Sg=np.array(S,dtype=np.int8); K=len(S)
X={0:(0,1,1),1:(1,0,1),2:(0,0,1)}
print('Sigma rows',K,'members',len(I.H))
# block j=1 (middle): members 3,4,5 ; C_1 = member 5 ; its ac pair = vertices of C
j=1; C=I.H[3*j+2]; a,b,c=C; pr=(a,c)
print('C',C,'intended pair',pr)
verts=sorted({v for m in range(3*j,3*j+3) for v in I.H[m]})
print('candidate adjacency vertices',verts)
def guard_masks(U):
    """for pendant adjacent to exactly U: dict soft pair -> (required value, mask of rows satisfying)"""
    res={}
    Us=set(U)
    for x,y in itertools.combinations(range(I.n),2):
        o=I.owner.get(key(x,y))
        if o is None: continue
        s=int(x in Us)+int(y in Us)
        if s==0: continue
        i,k=o
        want=0 if s==1 else 1
        allowed=[st for st in (0,1,2) if X[k][st]==want]
        sat=np.isin(Sg[:,i],allowed)
        res[(x,y)]=(want,sat.all(),sat.sum())
    return res
sols=[]
for r in range(1,len(verts)+1):
    for U in itertools.combinations(verts,r):
        g=guard_masks(U)
        if key(*pr) not in g: continue
        want,vac,cnt=g[key(*pr)]
        if want!=0 or vac: continue
        others=[(p,v) for p,v in g.items() if p!=key(*pr)]
        if all(v[1] for p,v in others):
            sols.append((U,len(others)))
print('adjacency sets U with intended guard non-vacuous and all other guards vacuous:',len(sols))
for U,k in sols[:10]: print(U,k)
