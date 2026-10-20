"""LPS23 variable gadget (Sect. 4.1.1): cliques K0,K1,K2 of 5 vertices (F5), P3 v0p-v1q-v2r for every q-p=r-q. Are there soft triangles?"""
import itertools
from a10_lib import *
import numpy as np
n=15
v=lambda j,p:5*j+p
A=np.zeros((n,n),bool)
for j in range(3):
    for p,q in itertools.combinations(range(5),2): A[v(j,p),v(j,q)]=A[v(j,q),v(j,p)]=True
for j in (0,1):
    for p in range(5):
        for q in range(5): A[v(j,p),v(j+1,q)]=A[v(j+1,q),v(j,p)]=True
H=[(v(0,p),v(1,q),v(2,r)) for p in range(5) for q in range(5) for r in range(5) if (q-p)%5==(r-q)%5]
print("members",len(H),"valid(pair-disjoint,induced P3)",valid(A,H))
o=owner_map(H)
cross=[frozenset((v(i,p),v(j,q))) for i in range(3) for j in range(i+1,3) for p in range(5) for q in range(5)]
print("all cross pairs soft:",all(c in o for c in cross), "cmax",cmax(n,H))
st=soft_triangles(n,H); print("soft triangles:",len(st), st[:3])
# a concrete one
p,q,r=st[0]; print([o[frozenset(x)] for x in ((p,q),(q,r),(p,r))],[H[o[frozenset(x)]] for x in ((p,q),(q,r),(p,r))])
