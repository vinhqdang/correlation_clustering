import itertools, sys, random, time, pickle
from collections import Counter
import numpy as np
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import build_tree, enum_sigma
U=[u for u in pickle.load(open('/tmp/w9b_units.pkl','rb')) if len(u[4])>=2]
def run(seed,trials,shape,kc):
    rng=random.Random(seed)
    edges={'chain2':[(1,0)],'chain3':[(1,0),(2,1)],'chain4':[(1,0),(2,1),(3,2)],'star2':[(1,0),(2,0)],'star3':[(1,0),(2,0),(3,0)],'bin7':[(1,0),(2,0),(3,1),(4,1),(5,2),(6,2)]}[shape]
    nu=1+len(edges)
    dist=Counter(); best=None; feas=0
    for t in range(trials):
        units=[rng.choice(U) for _ in range(nu)]
        used=[set() for _ in range(nu)]; ed=[]
        for (c,p) in edges:
            cv=rng.choice([v for v in range(6) if v not in used[c]]); pv=rng.choice([v for v in range(6) if v not in used[p]])
            used[c].add(cv); used[p].add(pv); ed.append((c,cv,p,pv))
        cross=[]
        for _ in range(rng.choice(kc)):
            (c,cv,p,pv)=rng.choice(ed); cross.append(((c,rng.randrange(6)),(p,rng.randrange(6))))
        I,maps=build_tree([(u[0],u[1]) for u in units],ed,cross)
        if not I.valid() or I.cmax()>4: continue
        S=enum_sigma(I)
        if len(S)<2: continue
        feas+=1
        Sg=np.array(S,dtype=np.int8)
        k1=max_core(Sg,singleton_only=True); k2=max_core(Sg)
        dist[(k1,k2)]+=1
        if best is None or (k1,k2)>best[0]: best=((k1,k2),(I.n,sorted(I.E),I.H),len(S))
    return shape,trials,feas,sorted(dist.items()),best
if __name__=='__main__':
    seed=int(sys.argv[1]); trials=int(sys.argv[2]); shape=sys.argv[3]
    t0=time.time()
    r=run(seed,trials,shape,[0,0,1,2])
    print(r[0],'trials',r[1],'feasible',r[2],'dist',r[3],'time',round(time.time()-t0)); print('best',r[4][0] if r[4] else None)
    if r[4] and r[4][0][0]>=3: print(r[4])
