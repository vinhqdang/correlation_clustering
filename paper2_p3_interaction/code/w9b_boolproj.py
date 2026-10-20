"""Boolean sub-domain projections of Sigma onto member triples: bijunctive or not (majority closure)."""
import itertools, sys, random, time
import numpy as np
from collections import Counter
from w9b_lib import *
from w9b_octtree import enum_sigma
SUBD=[(0,1),(0,2),(1,2)]
def nonbij_triples(Sg,m):
    out=[]
    for i,j,k in itertools.combinations(range(m),3):
        for di in SUBD:
            for dj in SUBD:
                for dk in SUBD:
                    sel=np.isin(Sg[:,i],di)&np.isin(Sg[:,j],dj)&np.isin(Sg[:,k],dk)
                    if sel.sum()<3: continue
                    proj={(di.index(a),dj.index(b),dk.index(c)) for a,b,c in Sg[sel][:,[i,j,k]].tolist()}
                    # majority closure test
                    L=list(proj); bad=False
                    for t1 in L:
                        for t2 in L:
                            for t3 in L:
                                r=tuple((x+y+z)>=2 for x,y,z in zip(t1,t2,t3))
                                if tuple(int(v) for v in r) not in proj: bad=True;break
                            if bad:break
                        if bad:break
                    if bad: out.append(((i,j,k),(di,dj,dk),len(proj)))
    return out
if __name__=='__main__':
    from w9b_locality import rand_tree_instance
    seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); cmax=int(sys.argv[4])
    rng=random.Random(seed); cnt=Counter(); ex=[]
    for _ in range(N):
        I=rand_tree_instance(rng,m,cmax,rng.choice([0.0,0.3,0.5]))
        if I is None or not I.valid(): continue
        if not I.soft_triangles(): continue
        S=enum_sigma(I,limit=3000,drop_far=True)
        if len(S)<3 or len(S)>=3000: continue
        Sg=np.array(S,dtype=np.int8)
        nb=nonbij_triples(Sg,len(I.H))
        cnt[bool(nb)]+=1
        if nb: ex.append((I,nb))
    print('instances with soft triangles',sum(cnt.values()),'non-bijunctive Boolean projection exists in',cnt[True])
    for I,nb in ex[:3]: print(I.H,sorted(I.E),nb[:3])
