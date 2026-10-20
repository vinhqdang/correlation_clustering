"""pair-maximality of the reduction instance: no induced P3 of G with all three pairs hard (exact check via hard-edge adjacency + hard non-edge structure)."""
import sys, random, itertools
from w10_red import *
def pair_maximal(n, E, members):
    soft=set()
    for t in members: soft|=pairs_of(t)
    hadj=[set() for _ in range(n)]
    for (u,v) in E:
        if (u,v) not in soft: hadj[u].add(v); hadj[v].add(u)
    bad=[]
    for v in range(n):
        nb=sorted(hadj[v])
        for i in range(len(nb)):
            for j in range(i+1,len(nb)):
                u,w=nb[i],nb[j]
                pq=key(u,w)
                if pq in soft: continue          # soft pair
                if pq in E: continue             # hard edge -> triangle
                bad.append((u,v,w))
    return bad
if __name__=='__main__':
    k,n,p=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
    rng=random.Random(1)
    S=[[set((a,b) for a in range(1,n+1) for b in range(1,n+1) if rng.random()<0.6) or {(1,1)} for _ in range(k)] for _ in range(k)]
    W,meta=build_gt(k,n,S,p,thin_mode=False)
    I=Inst(W)
    bad=pair_maximal(I.n,I.E,I.members)
    print('n',I.n,'all-hard induced P3 count',len(bad))
    for u,v,w in bad[:5]: print([W.B.names[x] for x in (u,v,w)])
