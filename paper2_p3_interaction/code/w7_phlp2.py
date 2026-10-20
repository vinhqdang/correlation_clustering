import random, sys
from w7_phlp import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); k=int(sys.argv[3]); n=int(sys.argv[4]); smax=int(sys.argv[5])
stats={True:0,False:0}; mism=0
for _ in range(N):
    S={}
    for i in range(k):
        for j in range(k):
            allt=[(a,b) for a in range(1,n+1) for b in range(1,n+1)]
            S[(i,j)]=set(rng.sample(allt,rng.randint(1,smax)))
    gt=gt_solve(k,n,S); nid,paths,links=build(k,n,S)
    ph=lpph_solve(nid,paths,links,k*k); stats[gt]+=1
    if gt!=ph: mism+=1; print('MISMATCH',S)
print('k=%d n=%d: YES %d NO %d mismatches %d'%(k,n,stats[True],stats[False],mism))
