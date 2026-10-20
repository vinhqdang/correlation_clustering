import random, sys
from collections import Counter
from w8b_ilp import ILP
from w8b_branch import gen_chain
rng=random.Random(int(sys.argv[1])); m=int(sys.argv[2]); N=int(sys.argv[3]); c=Counter(); best=(0,None)
for _ in range(N):
    I=gen_chain(m,rng.randint(3,3*m),rng)
    if I is None: continue
    o=ILP(I); B=len(I.H)
    if o.F({},o.N)<=B: continue
    T=list(o.N)
    for q in list(T):
        T2=[x for x in T if x!=q]
        if o.F({},T2)>B: T=T2
    c[len(T)]+=1
    if len(T)>best[0]: best=(len(T),(sorted(I.E),I.H))
print("m",m,"l=0 core sizes",sorted(c.items()),"max",best[0])
