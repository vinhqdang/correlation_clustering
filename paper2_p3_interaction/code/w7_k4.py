import itertools, random, sys
import numpy as np
from w7_lib import *
from w7_rel7 import closed_all_majorities
import io, contextlib
# K4: members 0..3, vertex v_{ij} shared by members i<j  (6 vertices)
vid={}
for (i,j) in itertools.combinations(range(4),2): vid[(i,j)]=len(vid)
mem=[[vid[tuple(sorted((i,j)))] for j in range(4) if j!=i] for i in range(4)]
soft=set()
for m in mem:
    for a,b in itertools.combinations(m,2): soft.add(key(a,b))
hard=[p for p in itertools.combinations(range(6),2) if p not in soft]
print('K4: hard pairs',hard)
seen=set(); tot=0; nd=0; ncl=0
for orders in itertools.product(*[list(itertools.permutations(m)) for m in mem]):
    for eps in itertools.product((0,1),repeat=len(hard)):
        A=np.zeros((6,6),bool)
        for (a,b,c) in orders:
            A[a,b]=A[b,a]=A[b,c]=A[c,b]=True
        for p,e in zip(hard,eps): A[p[0],p[1]]=A[p[1],p[0]]=bool(e)
        I=Inst(A,[tuple(o) for o in orders])
        if not I.valid(): continue
        sg=I.sigma()
        if not sg: continue
        k=frozenset(sg)
        if k in seen: continue
        seen.add(k)
print('K4 tight configurations: distinct Sigma sets',len(seen))
nd=sum(1 for k in seen if not two_decomposable(sorted(k)))
ncl=sum(1 for k in seen if not closed_all_majorities(sorted(k)))
print('non-2-decomposable',nd,'not closed under all majorities',ncl)
