import random, sys
from collections import Counter, defaultdict
from w8b_ilp import ILP
from w8b_csp import ILPcsp
from w8b_branch import gen_chain
seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); mode=sys.argv[4]; rng=random.Random(seed)
res=defaultdict(Counter)
for _ in range(N):
    I=gen_chain(m,rng.randint(2,2+m),rng)
    if I is None: continue
    o=(ILP if mode=='part' else ILPcsp)(I)
    for l in (0,1):
        B=len(I.H)+l
        if o.F({},o.N)<=B: continue
        T=list(o.N)
        for q in list(T):
            T2=[x for x in T if x!=q]
            if o.F({},T2)>B: T=T2
        res[l][len(T)]+=1
print("m",m,mode,{l:sorted(c.items()) for l,c in res.items()})
