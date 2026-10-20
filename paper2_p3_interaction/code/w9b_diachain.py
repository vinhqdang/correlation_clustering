"""chains of diamond blocks (4 members, ports A and D) glued D_j.d == A_{j+1}.a ; search for kappa2 growth."""
import itertools, sys, time, pickle
import numpy as np
from collections import defaultdict
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_sigma import max_core
SETS=[(0,1,3),(1,2,4),(0,2,5),(4,5,6)]
n=7
soft=set()
for t in SETS:
    for p in itertools.combinations(sorted(t),2): soft.add(p)
hard=[p for p in itertools.combinations(range(n),2) if p not in soft]
XTt=np.array([[0,1,1],[1,0,1],[0,0,1]])
states=np.array(list(itertools.product((0,1,2),repeat=4)))
statuses=np.array(list(itertools.product((0,1),repeat=len(hard))))
hid={p:i for i,p in enumerate(hard)}
def groups():
    G={}
    for orient in itertools.product(range(6),repeat=4):
        H=[tuple(itertools.permutations(SETS[k]))[orient[k]] for k in range(4)]
        owner={}
        for i,(a,b,c) in enumerate(H):
            for r,p in enumerate((key(a,b),key(b,c),key(a,c))): owner[p]=(i,r)
        ok=np.ones((512,81),bool)
        for x,y,z in itertools.combinations(range(n),3):
            ps=[key(x,y),key(y,z),key(x,z)]
            ss=np.zeros(81,int); hs=np.zeros(512,int); anys=False
            for p in ps:
                o=owner.get(p)
                if o is not None: ss+=XTt[o[1]][states[:,o[0]]]; anys=True
                else: hs+=statuses[:,hid[p]]
            if not anys: ok&=(hs!=2)[:,None]; continue
            ok&=((ss[None,:]+hs[:,None])!=2)
        sig=(H[0].index(3),H[3].index(6))
        for si in range(512):
            sol=np.nonzero(ok[si])[0]
            if len(sol)<2: continue
            fs=frozenset(map(tuple,states[sol].tolist()))
            G.setdefault((fs,sig),(H,tuple(statuses[si])))
    return G
def chain_inst(N,cfg,eps):
    H0,st=cfg
    vid={}
    def V(j,l):
        if j>0 and l==3: return V(j-1,6)
        if (j,l) not in vid: vid[(j,l)]=len(vid)
        return vid[(j,l)]
    H=[];E=set();sf=set()
    for j in range(N):
        for t in H0: H.append(tuple(V(j,x) for x in t))
    for (a,b,c) in H: E.add(key(a,b));E.add(key(b,c)); sf|={key(a,b),key(b,c),key(a,c)}
    for j in range(N):
        for p,s in zip(hard,st):
            u,v=V(j,p[0]),V(j,p[1])
            if key(u,v) in sf: continue
            if s: E.add(key(u,v))
    pairs=[(4,0),(4,1),(5,0),(5,1)]   # D_j others (b,c) vs A_{j+1} others (p,q)
    for j in range(N-1):
        for bit,(x,y) in zip(eps,pairs):
            if bit:
                u,v=V(j,x),V(j+1,y)
                if key(u,v) not in sf: E.add(key(u,v))
    return Inst(len(vid),E,H)
if __name__=='__main__':
    t0=time.time(); G=groups(); print('groups',len(G),round(time.time()-t0),flush=True)
    good=[]; cnt=0
    for (fs,sig),cfg in G.items():
        for eps in itertools.product((0,1),repeat=4):
            I=chain_inst(3,cfg,eps)
            if not I.valid(): continue
            S=enum_sigma(I,limit=3000)
            if len(S)<2 or len(S)>=3000: continue
            Sg=np.array(S,dtype=np.int8)
            k2=max_core(Sg,maxk=8); cnt+=1
            if k2>=5: good.append((sorted(fs),sig,cfg,eps,len(S),k2)); print('GOOD',sig,cfg,eps,len(S),k2,flush=True)
    print('evaluated',cnt,'good',len(good),round(time.time()-t0),flush=True)
    pickle.dump(good,open('/tmp/w9b_diagood.pkl','wb'))
