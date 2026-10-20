"""Tunnel-release model, exhaustive over all cut sets of size <= k^2 (k=2): ZE(G - C) feasible  iff  C = one tunnel per cell forming a valid tiling."""
import sys, itertools, random, time
from w10_red import *
k,n=int(sys.argv[1]),int(sys.argv[2]); seed=int(sys.argv[3]); dens=float(sys.argv[4]); maxc=int(sys.argv[5]) if len(sys.argv)>5 else k*k
rng=random.Random(seed)
S=[[set() for _ in range(k)] for _ in range(k)]
for i in range(k):
    for j in range(k):
        for a in range(1,n+1):
            for b in range(1,n+1):
                if rng.random()<dens: S[i][j].add((a,b))
        if not S[i][j]: S[i][j].add((1,1))
W,meta=build_gt(k,n,S,2,thin_mode=True)
I=ThinInst(W)
T=meta['tunnels']
cells=[(i,j) for i in range(k) for j in range(k)]
sols=gt_solvable(k,n,S)
print('S',S,'tunnels',len(T),'GT sols',len(sols),flush=True)
t0=time.time(); tot=0; bad=0; nfeas=0
for r in range(0,maxc+1):
    for C in itertools.combinations(range(len(T)),r):
        # expected
        cs=[T[c][0] for c in C]
        exp=False
        if r==k*k and len(set(cs))==k*k:
            t={T[c][0]:T[c][1] for c in C}
            exp = t in sols
        lab=I.solve(drop_hard=[(T[c][2],T[c][3]) for c in C])
        tot+=1; nfeas+= lab is not None
        if (lab is not None)!=exp:
            bad+=1; print('MISMATCH',[T[c][:2] for c in C],lab is not None,exp,flush=True)
print('cut sets tested',tot,'feasible',nfeas,'mismatches',bad,round(time.time()-t0,1))
