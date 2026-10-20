"""search all K3-block chains (block config x out/in member x link eps) for unbounded list-core growth (kappa2(N=3) >= 5)."""
import itertools, sys, time, pickle
import numpy as np
from collections import defaultdict, Counter
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_sigma import max_core
from w9b_scope import role_orders, scope_relation
# local ids: p=0,q=1,a=2,r=3,b=4,c=5 ; K3 members A on {p,q,a}, B on {q,r,b}, C on {r,p,c}
hardpairs=[(2,4),(2,5),(4,5),(2,3),(4,0),(5,1)]
priv={0:2,1:4,2:5}     # private vertex of member 0 (A), 1 (B), 2 (C)
def groups():
    G=defaultdict(list)
    for A in role_orders((0,1,2)):
        for B in role_orders((1,3,4)):
            for C in role_orders((3,0,5)):
                sig=(A.index(2),B.index(4),C.index(5))
                for st in itertools.product((0,1),repeat=6):
                    hard={frozenset(p):s for p,s in zip(hardpairs,st)}
                    R=scope_relation([A,B,C],hard,6)
                    if len(R)>=2: G[(R,sig)].append((A,B,C,st))
    return G
def chain_inst(N,cfg,outm,inm,eps):
    A,B,C,st=cfg
    vid={}
    def V(j,l):
        if j>0 and l==priv[inm]: return V(j-1,priv[outm])
        if (j,l) not in vid: vid[(j,l)]=len(vid)
        return vid[(j,l)]
    H=[];E=set();soft=set()
    for j in range(N):
        for t in (A,B,C): H.append(tuple(V(j,x) for x in t))
    for (a,b,c) in H: E.add(key(a,b));E.add(key(b,c)); soft|={key(a,b),key(b,c),key(a,c)}
    for j in range(N):
        for (x,y),s in zip(hardpairs,st):
            u,v=V(j,x),V(j,y)
            if key(u,v) in soft: continue
            if s: E.add(key(u,v))
    # cross eps between out member's other vertices (in block j) and in member's other vertices (block j+1)
    outm_vs=[v for v in (A,B,C)[outm] if v!=priv[outm]]
    inm_vs=[v for v in (A,B,C)[inm] if v!=priv[inm]]
    pairs=[(x,y) for x in outm_vs for y in inm_vs]
    for j in range(N-1):
        for bit,(x,y) in zip(eps,pairs):
            if bit:
                u,v=V(j,x),V(j+1,y)
                if key(u,v) not in soft: E.add(key(u,v))
    return Inst(len(vid),E,H)
if __name__=='__main__':
    t0=time.time()
    G=groups(); print('groups',len(G),round(time.time()-t0),flush=True)
    keys=list(G); good=[]
    cnt=0
    for (R,sig),cfgs in G.items():
        cfg=cfgs[0]
        for outm in range(3):
            for inm in range(3):
                for eps in itertools.product((0,1),repeat=4):
                    I=chain_inst(3,cfg,outm,inm,eps)
                    if not I.valid(): continue
                    S=enum_sigma(I,limit=5000)
                    if len(S)<2 or len(S)>=5000: continue
                    Sg=np.array(S,dtype=np.int8)
                    k2=max_core(Sg,maxk=8)
                    cnt+=1
                    if k2>=5: good.append((sorted(R),sig,cfg,outm,inm,eps,len(S),k2))
    print('chains evaluated',cnt,'with kappa2(N=3)>=5:',len(good),round(time.time()-t0),'s')
    pickle.dump(good,open('/tmp/w9b_k3good.pkl','wb'))
    for g in good[:20]: print(g)
