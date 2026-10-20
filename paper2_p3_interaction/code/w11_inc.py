"""w11_inc.py -- increment lemma test: for random instance and random frozen set Phi of remote pairs, size of the smallest C with F'(Phi u C) > F'(Phi)."""
import sys, time, random
from collections import Counter
from w11_lib import *
from w11_gen import *
m=int(sys.argv[1]); npend=int(sys.argv[2]); samples=int(sys.argv[3]); seed=int(sys.argv[4]); cmax=int(sys.argv[5])
rnd=random.Random(seed); cnt=Counter(); worst=[]
t0=time.time()
for s in range(samples):
    if cmax==2: nv,H=random_members_c2(m,rnd)
    else: nv,H=random_members_cc(m,rnd,cmax)
    n=nv+npend
    k=rnd.choice([2,3,4,6]); cl=[rnd.randrange(-1,k) if rnd.random()<0.8 else -1 for _ in range(n)]
    I=make_inst(nv,H,npend,cl)
    if not (I.valid() and I.pair_maximal() and I.cmax()<=cmax) or n>20: continue
    rem=I.remote()
    # random Phi: random subset (size from 0..len/2) possibly structured as star of random vertices
    mode=rnd.choice(['rand','vertex','empty'])
    if mode=='empty': Phi=set()
    elif mode=='rand': Phi=set(rnd.sample(rem,min(len(rem),rnd.randrange(0,max(1,len(rem)//2)))))
    else:
        vs=rnd.sample(range(n),min(n,rnd.randrange(1,5))); Phi=set(p for p in rem if p[0] in vs or p[1] in vs)
    sz,C,t=min_increment(I,Phi,max_size=30,workers=2)
    cnt[sz]+=1
    if sz is not None and sz>=6: worst.append((sz,t-len(H),mode,nv,H,npend,cl,sorted(Phi)))
print(m,npend,cmax,cnt,round(time.time()-t0))
for w in sorted(worst,key=lambda x:-x[0])[:3]: print(w[:4], w[4:7])
