"""all diamond (K4 minus an edge in I(H), c=2) gadgets: 4 members A,B,C,D on 7 vertices; every orientation and every status of the 9 hard pairs."""
import itertools, pickle, time
import numpy as np
from w9b_lib import key
# vertices: p0 q1 r2 a3 b4 c5 d6 ; A=(p,q,a) B=(q,r,b) C=(p,r,c) D=(b,c,d)
SETS=[(0,1,3),(1,2,4),(0,2,5),(4,5,6)]
n=7
soft=set()
for t in SETS:
    for p in itertools.combinations(sorted(t),2): soft.add(p)
hard=[p for p in itertools.combinations(range(n),2) if p not in soft]
print('hard pairs',len(hard),hard)
XTt=np.array([[0,1,1],[1,0,1],[0,0,1]])
states=np.array(list(itertools.product((0,1,2),repeat=4)))   # 81 x 4
statuses=np.array(list(itertools.product((0,1),repeat=len(hard))))  # 512 x 9
hid={p:i for i,p in enumerate(hard)}
rels={}
t0=time.time()
orients=list(itertools.product(range(6),repeat=4))
for oi,orient in enumerate(orients):
    H=[tuple(itertools.permutations(SETS[k]))[orient[k]] for k in range(4)]
    owner={}
    for i,(a,b,c) in enumerate(H):
        for r,p in enumerate((key(a,b),key(b,c),key(a,c))): owner[p]=(i,r)
    ok=np.ones((512,81),bool)
    for x,y,z in itertools.combinations(range(n),3):
        ps=[key(x,y),key(y,z),key(x,z)]
        soft_sum=np.zeros(81,int); hs=np.zeros(512,int); anysoft=False
        for p in ps:
            o=owner.get(p)
            if o is not None: soft_sum+=XTt[o[1]][states[:,o[0]]]; anysoft=True
            else:
                # hard pair status contributes; but member edge status for soft pair is handled by XTt (togetherness), hard: edge iff status
                hs+=statuses[:,hid[p]]
        if not anysoft:
            ok&=(hs!=2)[:,None]; continue
        ok&=((soft_sum[None,:]+hs[:,None])!=2)
    # validity of the members as P3 in G is by construction (orientation); hard statuses free
    for si in range(512):
        sol=np.nonzero(ok[si])[0]
        if len(sol)==0: continue
        fs=frozenset(map(tuple,states[sol].tolist()))
        if fs not in rels: rels[fs]=(orient,tuple(statuses[si]))
print('distinct nonempty diamond Sigma:',len(rels),'time',round(time.time()-t0))
pickle.dump(rels,open('/tmp/w9b_diamond_rels.pkl','wb'))
from collections import Counter
print(sorted(Counter(len(r) for r in rels).items()))
