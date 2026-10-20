import random, itertools, numpy as np
from w8b_lib import *
from w8b_branch import gen_chain
from w8b_c1 import gen_c1
rng=random.Random(5); tot=0; ze=0; bad=0
for it in range(600):
    I=gen_chain(rng.choice([2,3]),rng.randint(0,4),rng) if it%2 else gen_c1(2,rng.randint(0,6),rng)
    if I is None: continue
    n=I.n; P=parts(n); pr=I.pairs
    a=np.array([p[0] for p in pr]); b=np.array([p[1] for p in pr])
    same=(P[:,a]==P[:,b]); Ev=np.array([p in I.E for p in pr])
    opt=(same!=Ev[None,:]).sum(1).min()
    zeCost = opt==len(I.H)
    ok=np.ones(len(P),bool)
    for i,p in enumerate(pr):
        if p in I.soft: continue
        ok&= (same[:,i]==Ev[i])
    idx=I.pidx
    for (x,y,z) in I.H:
        ok&= same[:,idx[tuple(sorted((x,y)))]] | same[:,idx[tuple(sorted((y,z)))]]
    att=bool(ok.any()); tot+=1; ze+=zeCost; bad+= (att!=zeCost)
print("instances",tot,"ZE feasible",ze,"mismatch",bad)
