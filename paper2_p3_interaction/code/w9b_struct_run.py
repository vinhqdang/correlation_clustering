import sys, random, time
from w9b_struct import *
seed=int(sys.argv[1]); kind=sys.argv[2]; a=int(sys.argv[3]); b=int(sys.argv[4]); steps=int(sys.argv[5]); cmax=int(sys.argv[6])
rng=random.Random(seed)
if kind=='strip': n,members=strip_struct(a,b)
t0=time.time()
r=anneal(n,members,steps,rng,cmax)
if r is None: print('no feasible init'); sys.exit()
s,k,o,h,I=r
print(kind,a,b,'members',len(members),'n',n,'best k1',k,'score',round(s,2),'time',round(time.time()-t0),flush=True)
if k>=3: print(I.H,sorted(I.E))
