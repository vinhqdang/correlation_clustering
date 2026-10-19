"""Central member C=(c-perm) and three anchored members A_i=(x_i, c_i, w_i) or (x_i, w_i, c_i); x_i clean.
Search the 9 free hard pairs among {c_i, w_i} for a 3-core: any two locks compatible, all three not."""
import itertools, numpy as np
from w6c_lib import *
n=9
lab, same, P = part_table(n); same_f=same.astype(np.float32)
pidx={p:i for i,p in enumerate(P)}
c=[0,1,2]; x=[3,4,5]; w=[6,7,8]
found=0; tested=0
for cperm in itertools.permutations(c):
  for ashape in itertools.product((0,1),repeat=3):
    members=[tuple(cperm)]
    for i in range(3):
        members.append((x[i],c[i],w[i]) if ashape[i]==0 else (x[i],w[i],c[i]))
    soft=set(); fe=set()
    for m in members:
        soft|=pairs_of(m); fe|={key(m[0],m[1]),key(m[1],m[2])}
    hard_free=[p for p in P if p not in soft and not (p[0] in x or p[1] in x)]
    assert len(hard_free)==9,len(hard_free)
    Hm=((np.arange(512)[:,None]>>np.arange(9)[None,:])&1).astype(np.float32)
    E=np.zeros((512,len(P)),np.float32)
    for p in fe: E[:,pidx[p]]=1
    for i,p in enumerate(hard_free): E[:,pidx[p]]=Hm[:,i]
    cost=same_f@(1-2*E.T)+E.sum(1)[None,:]
    k=4
    # states of anchored members: lock state = L if x at pos0 -> edge (x,mid) deleted => x apart from mid
    def lockmask(i):
        m=members[i+1]; return ~same[:,pidx[key(m[0],m[1])]]  # x separated from its neighbour
    # a lock means x singleton: x apart from all  (zero-excess => same as member delete edge)
    sing=[np.ones(len(same),bool) for _ in range(3)]
    for i in range(3):
        for t in range(n):
            if t!=x[i]: sing[i]&=~same[:,pidx[key(x[i],t)]]
    for j in range(512):
        ok=cost[:,j]==k
        if not ok.any(): continue
        tested+=1
        f=lambda S: (ok & np.logical_and.reduce([sing[i] for i in S])).any()
        if all(f(S) for S in itertools.combinations(range(3),2)) and not f((0,1,2)):
            found+=1
            if found<=5: print('3-CORE',cperm,ashape,j)
print('tested nonempty patterns',tested,'3-cores found',found)
