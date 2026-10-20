"""Binary relations Rel(N,N') realised by an octahedral K4 unit with two link members (c<=3)."""
import itertools, sys, random, time, pickle
from collections import Counter
import numpy as np
from w9b_lib import *
from w9b_sigma import sigma
from w9b_octchain import unit_members, perms
U=[u for u in pickle.load(open('/tmp/w9b_units.pkl','rb'))]
def build(u, v1, v2, N1roles, N2roles, eps1, eps2, eps12):
    H0,E0=u[0],u[1]
    # vertices 0..5 octahedron; N1 = shared v1 + new 6,7 ; N2 = shared v2 + new 8,9
    n=10
    ords1=list(itertools.permutations((v1,6,7)))[N1roles]
    ords2=list(itertools.permutations((v2,8,9)))[N2roles]
    H=list(H0)+[ords1,ords2]
    E=set(E0)
    for (a,b,c) in (ords1,ords2): E.add(key(a,b)); E.add(key(b,c))
    # eps1: edges between {6,7} and the 5 other oct vertices ; eps2 likewise for {8,9}; eps12 between {6,7} and {8,9}
    others1=[x for x in range(6) if x!=v1]; others2=[x for x in range(6) if x!=v2]
    for k,(a,b) in enumerate(itertools.product((6,7),others1)):
        if eps1>>k&1: E.add(key(a,b))
    for k,(a,b) in enumerate(itertools.product((8,9),others2)):
        if eps2>>k&1: E.add(key(a,b))
    for k,(a,b) in enumerate(itertools.product((6,7),(8,9))):
        if eps12>>k&1: E.add(key(a,b))
    return Inst(n,E,H)
if __name__=='__main__':
    seed=int(sys.argv[1]); trials=int(sys.argv[2])
    rng=random.Random(seed); rels=Counter(); info={}; t0=time.time()
    for t in range(trials):
        u=rng.choice(U)
        v1,v2=rng.sample(range(6),2)
        eps1=0
        for _ in range(rng.choice([0,0,1,2])): eps1|=1<<rng.randrange(10)
        eps2=0
        for _ in range(rng.choice([0,0,1,2])): eps2|=1<<rng.randrange(10)
        eps12=0
        for _ in range(rng.choice([0,0,1])): eps12|=1<<rng.randrange(4)
        I=build(u,v1,v2,rng.randrange(6),rng.randrange(6),eps1,eps2,eps12)
        if not I.valid(): continue
        S=sigma(I,drop_far=False)
        if len(S)==0: continue
        rel=frozenset((int(r[4]),int(r[5])) for r in S)
        rels[rel]+=1
        if rel not in info: info[rel]=(u[2],u[3],v1,v2,eps1,eps2,eps12)
    print('distinct binary relations Rel(N,N\'):',len(rels),'time',round(time.time()-t0))
    pickle.dump((rels,info),open('/tmp/w9b_octlink_%d.pkl'%seed,'wb'))
