"""chain of K3 blocks (block relation {LRL, RLR, RRM} on (A,B,C)), B_j glued to A_{j+1} through their private vertices."""
import itertools, sys
import numpy as np
from collections import Counter
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import enum_sigma
# local ids: p=0,q=1,a=2,r=3,b=4,c=5
A0=(1,0,2); B0=(1,3,4); C0=(0,5,3)
HARD0={(2,4):0,(2,5):0,(4,5):0,(2,3):0,(4,0):0,(5,1):1}   # (a,b),(a,c),(b,c),(a,r),(b,p),(c,q)
def chain(N,eps,orient=None,extra_cross=()):
    vid={}
    def V(j,l):
        if l==2 and j>0: return V(j-1,4)    # a_j == b_{j-1}
        if (j,l) not in vid: vid[(j,l)]=len(vid)
        return vid[(j,l)]
    H=[];E=set();soft=set()
    for j in range(N):
        for t in (A0,B0,C0):
            m=tuple(V(j,x) for x in t); H.append(m)
    for (a,b,c) in H: E.add(key(a,b));E.add(key(b,c)); soft|={key(a,b),key(b,c),key(a,c)}
    for j in range(N):
        for (x,y),s in HARD0.items():
            u,v=V(j,x),V(j,y)
            if key(u,v) in soft: continue
            if s: E.add(key(u,v))
    # cross eps between B_j others (q=1,r=3) and A_{j+1} others (p=0,q=1)
    pairs=[(1,0),(1,1),(3,0),(3,1)]
    for j in range(N-1):
        for bit,(x,y) in zip(eps,pairs):
            u,v=V(j,x),V(j+1,y)
            if bit and key(u,v) not in soft: E.add(key(u,v))
    for (j1,x,j2,y) in extra_cross:
        E.add(key(V(j1,x),V(j2,y)))
    return Inst(len(vid),E,H)
if __name__=='__main__':
    N=int(sys.argv[1])
    for eps in itertools.product((0,1),repeat=4):
        I=chain(N,eps)
        if not I.valid(): print(eps,'invalid'); continue
        S=enum_sigma(I)
        if len(S)==0: print(eps,'N',N,'infeasible'); continue
        Sg=np.array(S,dtype=np.int8)
        print(eps,'N',N,'|Sigma|',len(S),'k1',max_core(Sg,singleton_only=True),'k2',max_core(Sg,maxk=9),flush=True)
