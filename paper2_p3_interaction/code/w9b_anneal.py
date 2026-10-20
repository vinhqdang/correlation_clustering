"""simulated annealing for instances with large list-core (singleton lists k1)"""
import sys, random, itertools, time, math, pickle
import numpy as np
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import enum_sigma

def score(I):
    if not I.valid() or I.cmax()>CMAX: return -1,None
    S=enum_sigma(I,limit=3000)
    if len(S)<2 or len(S)>3000: return -1,None
    Sg=np.array(S,dtype=np.int8)
    k1=max_core(Sg,singleton_only=True,maxk=7)
    return k1*10+min(len(S),30)*0.01, k1

def rand_members(rng,n,m):
    H=[];usedp=set();cnt=[0]*n;tries=0
    while len(H)<m and tries<1000:
        tries+=1
        t=tuple(rng.sample(range(n),3))
        ps=[key(t[0],t[1]),key(t[1],t[2]),key(t[0],t[2])]
        if any(p in usedp for p in ps): continue
        if any(cnt[v]>=CMAX for v in t): continue
        H.append(t); usedp|=set(ps)
        for v in t: cnt[v]+=1
    return H if len(H)==m else None

def make(n,H,hardedges):
    E=set()
    for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
    E|=set(hardedges)
    return Inst(n,E,H)

def run(seed,n,m,steps):
    rng=random.Random(seed)
    # initialise: random members and planted-cluster hard statuses until feasible
    while True:
        H=rand_members(rng,n,m)
        if H is None: continue
        soft=set()
        for (a,b,c) in H: soft|={key(a,b),key(b,c),key(a,c)}
        lab=[rng.randrange(max(2,n//3)) for _ in range(n)]
        he=set(key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft and lab[u]==lab[v] and rng.random()>0.1)
        I=make(n,H,he); s,k=score(I)
        if s>=0: break
    cur=(s,H,he); best=cur; bestk=k
    T0=0.5
    for it in range(steps):
        s0,H0,he0=cur
        H1=list(H0); he1=set(he0)
        mv=rng.random()
        soft0=set()
        for (a,b,c) in H0: soft0|={key(a,b),key(b,c),key(a,c)}
        if mv<0.7:
            # toggle hard pairs
            for _ in range(rng.choice([1,1,2])):
                u,v=rng.sample(range(n),2); p=key(u,v)
                if p in soft0: continue
                if p in he1: he1.remove(p)
                else: he1.add(p)
        else:
            # re-draw one member
            i=rng.randrange(m)
            for _ in range(20):
                t=tuple(rng.sample(range(n),3))
                H2=H1[:i]+[t]+H1[i+1:]
                ps=set(); ok=True
                for (a,b,c) in H2:
                    for p in (key(a,b),key(b,c),key(a,c)):
                        if p in ps: ok=False
                        ps.add(p)
                if ok: H1=H2; break
            soft1=set()
            for (a,b,c) in H1: soft1|={key(a,b),key(b,c),key(a,c)}
            he1={p for p in he1 if p not in soft1}
        I=make(n,H1,he1)
        s,k=score(I)
        if s<0: continue
        T=T0*(1-it/steps)+0.01
        if s>=s0 or rng.random()<math.exp((s-s0)/T):
            cur=(s,H1,he1)
            if s>best[0]: best=cur; bestk=k
    return best,bestk

if __name__=='__main__':
    seed=int(sys.argv[1]); n=int(sys.argv[2]); m=int(sys.argv[3]); steps=int(sys.argv[4]); CMAX=int(sys.argv[5])
    t0=time.time()
    best,k=run(seed,n,m,steps)
    print('seed',seed,'n',n,'m',m,'c<=',CMAX,'best k1',k,'score',round(best[0],2),'time',round(time.time()-t0))
    if k>=3: print(best[1],sorted(best[2]))
