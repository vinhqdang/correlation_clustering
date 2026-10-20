import numpy as np, itertools, io, contextlib
from w7_lib import *
with contextlib.redirect_stdout(io.StringIO()):
    from w7_cores import cores
for N in (5,6,7):
    n=2*N+1
    A=np.zeros((n,n),bool)
    for i in range(2*N): A[i,i+1]=A[i+1,i]=True
    H=[(2*i,2*i+1,2*i+2) for i in range(N)]
    I=Inst(A,H); sg=I.sigma()
    an,co=cores(I,sg)
    # anchors that are lockable alone
    # recompute singleton sets
    sing=[]
    for s in sg:
        B=I.A.copy()
        for (a,b,c),st in zip(I.H,s):
            u,v={'L':(a,b),'R':(b,c),'M':(a,c)}[st]; B[u,v]=not B[u,v]; B[v,u]=not B[v,u]
        sing.append(frozenset(v for v in range(n) if not B[v].any()))
    single=[v for v in an if any(v in x for x in sing)]
    print('N=%d n=%d |Sigma|=%d lockable single anchors %s; cores (min non-faces among anchors): %s'%(N,n,len(sg),single,co))
