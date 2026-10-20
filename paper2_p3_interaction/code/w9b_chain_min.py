"""exact minimum certificate size for the K3-chain-with-pendants instance (csp variant, l=0), by enumeration of guard sets."""
import itertools, sys, time
import numpy as np
from w9b_chain_cert import *
from w9b_octtree import enum_sigma
def guards_of_pendant(I,w):
    """effective guards of pendant w (all pairs of w are far): triple (w,x,y) with soft pair xy, e(w,x)+e(w,y) in {1,2}: far pairs {wx,wy}.
    returns list of (cost-pairs frozenset, member, allowed set)"""
    X={0:(0,1,1),1:(1,0,1),2:(0,0,1)}
    out=[]
    n=I.n
    for x,y in itertools.combinations(range(n),2):
        if x==w or y==w: continue
        o=I.owner.get(key(x,y))
        if o is None: continue
        s=int(I.adj[w][x])+int(I.adj[w][y])
        if s==0: continue
        i,k=o
        # x_e + s != 2  -> if s==1: x_e=0 ; if s==2: x_e=1
        want=0 if s==1 else 1
        allowed=frozenset(st for st in (0,1,2) if X[k][st]==want)
        out.append((frozenset([key(w,x),key(w,y)]),i,allowed))
    return out
def run(N):
    I,T=build(N)
    pend=sorted({p[0] if p[0]>=len(I.mem) - 0 else p[1] for p in T}) if False else None
    nv=I.n
    pendants=[v for v in range(nv) if not I.mem[v]]
    S=enum_sigma(I,limit=10**6,drop_far=True)
    Sg=np.array(S,dtype=np.int8)
    K=len(S); full=(1<<K)-1
    rows=Sg
    def mask(i,allowed):
        m=0
        idx=np.nonzero(np.isin(rows[:,i],list(allowed)))[0]
        for j in idx: m|=1<<int(j)
        return m
    items=[]  # (pairs, mask)
    for w in pendants:
        for pairs,i,allowed in guards_of_pendant(I,w):
            items.append((pairs,mask(i,allowed)))
    # dedupe by pair-set
    d={}
    for pairs,m in items:
        d[pairs]=d.get(pairs,full)&m
    items=list(d.items())
    # minimum number of far pairs in union of chosen items s.t. AND of masks == 0
    # iterative deepening on total number of far pairs
    best=None
    for budget in range(2,2*N+7):
        # choose items greedily by DFS with union size <= budget
        found=[None]
        def rec(start,cur,pairs):
            if found[0]: return
            if cur==0: found[0]=set(pairs); return
            for j in range(start,len(items)):
                pj,mj=items[j]
                np_=pairs|pj
                if len(np_)>budget: continue
                nxt=cur&mj
                if nxt==cur: continue
                rec(j+1,nxt,np_)
                if found[0]: return
        rec(0,full,frozenset())
        if found[0]: best=(budget,found[0]); break
    return I,T,best,len(items)
if __name__=='__main__':
    for N in range(1,int(sys.argv[1])+1):
        t0=time.time()
        I,T,best,ni=run(N)
        print('N',N,'|T| intended',len(T),'| exact minimum certificate size over guard sets:',best[0] if best else None,'(items',ni,')',round(time.time()-t0,1),flush=True)
