"""Exact formulation of the 'U-only' model: members singly edited, excess = |U| (unpacked edits).  OPT_U - |H| = min{|U| : Sigma(G xor U) != empty}.
Compare with the true OPT (C subset DP) on random soft-triangle-free perturbed instances (n <= 12)."""
import random, itertools, sys
import numpy as np
from w7_lib import *
from w7_dpy import opt
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); LMAX=2
eq=0; less=0; tot=0; hist={}
while tot<N:
    m=rng.randint(3,5); n=rng.randint(int(1.4*m)+2,min(12,3*m))
    I=gen_planted(n,m,rng,cmax=2,maxcl=4,tries=1500)
    if I is None or I.soft_triangles(): continue
    A=I.A.copy()
    hp=[(a,b) for a in range(n) for b in range(a+1,n) if key(a,b) not in I.owner]
    for _ in range(rng.randint(1,2)):
        a,b=rng.choice(hp); A[a,b]=A[b,a]=not A[a,b]
    J=Inst(A,I.H)
    edges=[(a,b) for a in range(n) for b in range(a+1,n) if A[a,b]]
    true,_=opt(n,edges,0)
    ex=true-m
    # U-only value (search up to LMAX)
    uo=None
    for r in range(0,LMAX+1):
        for U in itertools.combinations(hp,r):
            B=A.copy()
            for a,b in U: B[a,b]=B[b,a]=not B[a,b]
            if Inst(B,J.H).sigma(): uo=r;break
        if uo is not None: break
    tot+=1
    key_=(ex,uo); hist[key_]=hist.get(key_,0)+1
    if uo is not None:
        assert uo>=ex
for k in sorted(hist,key=lambda t:(t[0],10 if t[1] is None else t[1])): print('OPT-|H|=%d  U-only value=%s : %d instances'%(k[0],k[1],hist[k]))
