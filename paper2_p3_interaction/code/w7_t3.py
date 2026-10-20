"""dense c<=2 instances (n close to 3m/2) with many soft triangles"""
import random, sys, itertools, pickle
from w7_lib import *
from w7_cores import cores
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); CM=int(sys.argv[3])
tot=0; nclM=0; nd=[]; big=[]; ntry=0
while tot<N:
    ntry+=1
    m=rng.randint(4,8); n=rng.randint(int(1.5*m)+1, 2*m+2)
    I=gen_planted(n,m,rng,cmax=CM,maxcl=rng.choice([3,4,5,6]),tries=2000)
    if I is None or not I.soft_triangles(): continue
    sg=I.sigma()
    if len(sg)<3 or len(sg)>120: continue
    tot+=1
    if not closed_mmaj_fast(sg): nclM+=1
    if not two_decomposable(sg): nd.append((I.A.astype(int).tolist(),I.H,sg))
    an,co=cores(I,sg)
    if co and max(len(c) for c in co)>=3: big.append((I.A.astype(int).tolist(),I.H,co))
print('instances',tot,'tries',ntry,'non-Mclosed',nclM,'non-2dec',len(nd),'3-cores',len(big))
pickle.dump((nd,big),open('/tmp/w7_t3_%s_%s.pkl'%(sys.argv[1],CM),'wb'))
