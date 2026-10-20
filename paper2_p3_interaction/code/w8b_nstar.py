import itertools, sys, time
from w8b_lib import Inst
from w8b_ilp import ILP
def nstar(edges_oriented, tstar):
    # edges_oriented: list of (tail,head) in desired order, first one must be (t,'*')
    N=len(edges_oriented)-1
    nid={}; cnt=[0]
    def new(name): nid[name]=cnt[0]; cnt[0]+=1; return nid[name]
    for i in range(N+1):
        for s in 'abcu': new((s,i))
    z=new('z'); y=new('y')
    E=set()
    def add(p,q): E.add((min(p,q),max(p,q)))
    Kv={}
    for i,(t,h) in enumerate(edges_oriented):
        Kv.setdefault(t,[]).append(nid[('a',i)])
        Kv.setdefault(h,[]).extend([nid[('b',i)],nid[('c',i)],nid[('u',i)]])
    for v,K in Kv.items():
        for p,q in itertools.combinations(K,2): add(p,q)
    for i in range(N+1): add(nid[('a',i)],nid[('b',i)])
    for i in range(N): add(nid[('b',i)],nid[('u',i+1)])
    add(nid[('b',N)],z)
    add(y,nid[('u',0)]); add(y,nid[('c',0)])
    H=[]
    for i in range(N+1):
        H.append((nid[('a',i)],nid[('b',i)],nid[('c',i)]))
        H.append((nid[('u',i)],nid[('b',i)],nid[('u',i+1)] if i<N else z))
    return Inst(cnt[0],E,H)
if __name__=="__main__":
    K4=[(0,'*'),(0,1),(0,2),(1,3),(1,2),(3,0),(2,3)]
    I=nstar(K4,0); print("n",I.n,"|H|",len(I.H),"valid",I.valid(),"pairmax",I.pair_maximal())
    t=time.time(); o=ILP(I); print("width",o.width,"nonlocal",len(o.N),"c")
    r,st=o.run(1); print("l=1",r,st,time.time()-t)
    r,st=o.run(0); print("l=0",r,st)
