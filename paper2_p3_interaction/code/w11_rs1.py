import sys, time, random
from collections import Counter
from w11_lib import *
from w11_gen import *
m=int(sys.argv[1]); npend=int(sys.argv[2]); samples=int(sys.argv[3]); l=int(sys.argv[4]); seed=int(sys.argv[5]); K=int(sys.argv[6]) if len(sys.argv)>6 else 12
rnd=random.Random(seed)
cnt=Counter(); best=[]
t0=time.time()
for s in range(samples):
    nv,H=random_members(m,rnd)
    n=nv+npend
    k=rnd.choice([2,3,4,6,8])
    cl=[rnd.randrange(-1,k) if rnd.random()<0.8 else -1 for _ in range(n)]
    # singletons: -1 means own cluster
    I=make_inst(nv,H,npend,cl)
    assert I.valid() and I.pair_maximal() and I.cmax()<=2
    sz,T=min_certificate(I,l,max_size=K)
    cnt[sz]+=1
    if sz is not None and sz>=4: best.append((sz,nv+npend,H,cl))
print(m,npend,'l',l,cnt,round(time.time()-t0))
best.sort(key=lambda x:-x[0])
for b in best[:5]: print(b)
