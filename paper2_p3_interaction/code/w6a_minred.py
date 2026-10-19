"""distribution of s(T) = min size of a reducible (tight, R-ok, pair-closed) connected subfamily inside a tight component T, by |T|."""
import sys, random, collections
from w6a_lib import *
from w6a_rule import *
def minred(adj,H,T):
    for s in range(1,len(T)+1):
        for F in conn_subfamilies([H[i] for i in T], s):
            if len(F)!=s: continue
            FF=tuple(T[i] for i in F); W=members_vertices(H,FF)
            if pair_closed(H,FF,W) and any(R_ok(adj,S,W) for S in tight_choices(adj,H,FF,W)): return s
    return None
def main():
    seed,k,a,b=map(int,sys.argv[1:5]); rng=random.Random(seed); D=collections.Counter(); ex={}
    for it in range(k):
        n=rng.randint(a,b); adj=random_instance(n,rng); H=maximal_pair_packing(adj,rng)
        if not H: continue
        for T in I_components(H):
            W=members_vertices(H,T)
            if opt(induced(adj,W))>len(T): continue
            if not any(R_ok(adj,S,W) for S in tight_choices(adj,H,T,W)): D[(len(T),'irr')]+=1; continue
            s=minred(adj,H,T); D[(len(T),s)]+=1
            if s==len(T) and len(T)>=3: ex[len(T)]=(n,sorted(edges_of(adj)),H,T)
    print(sorted(D.items(),key=lambda x:str(x)))
    for m,e in ex.items(): print(m,e)

if __name__=="__main__": main()
