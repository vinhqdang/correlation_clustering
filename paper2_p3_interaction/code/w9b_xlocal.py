"""counterexample under the tree-decomposition locality of wave 8 (X_u cap X_v nonempty): minimal certificate = T plus the link cross pairs."""
import itertools, sys, time
from w9b_lib import *
from w9b_leaf import leaf_cost, min_certificate
from w9b_chain_cert import build
def xlocal_inst(N):
    I,T=build(N)
    m=len(I.H)
    # members: block j -> A=3j, B=3j+1, C=3j+2 ; bags: block bags and link bags
    bags=[]
    for j in range(N):
        bags.append(frozenset({3*j,3*j+1,3*j+2}))
        if j<N-1: bags.append(frozenset({3*j+1,3*(j+1)}))
    mem=I.mem
    X={}
    for v in range(I.n):
        if not mem[v]: X[v]=frozenset(); continue
        X[v]=frozenset(k for k,b in enumerate(bags) if set(mem[v])<=b)
    I.local=lambda u,v,X=X: bool(X[u]&X[v])
    return I,T
if __name__=='__main__':
    for N in range(2,int(sys.argv[1])+1):
        I,T=xlocal_inst(N)
        far=set(I.far_pairs())
        base=len(I.H)
        # certificate by deletion starting from all far pairs but in a fixed order: put pendant pairs LAST so they are deleted last? we test the structured candidate instead
        cross=[p for p in far if all(I.mem[x] for x in p)]
        Tfull=list(T)+cross
        c,_=leaf_cost(I,Tfull); c0,_=leaf_cost(I,[])
        print('N',N,'far pairs',len(far),'among V(H):',len(cross),'| F(T+cross)-|H|=',c-base,'F(empty)-|H|=',c0-base,flush=True)
        # minimality of T+(needed cross pairs): delete greedily cross pairs then T
        cur=list(Tfull)
        for p in list(cross):
            t2=[q for q in cur if q!=p]
            cc,_=leaf_cost(I,t2)
            if cc is None or cc>base: cur=t2
        needT=[p for p in cur if p in set(T)]
        ok=all((lambda cc:(cc is not None and cc<=base))(leaf_cost(I,[q for q in cur if q!=p])[0]) for p in cur)
        print('   greedy-minimal certificate size',len(cur),'(pendant pairs',len(needT),', link/cross pairs',len(cur)-len(needT),') minimal:',ok,flush=True)
