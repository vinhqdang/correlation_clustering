import random, sys, itertools
from collections import Counter, defaultdict
from w8b_lib import Inst
from w8b_ilp import ILP
def gen_c1(m,extra,rng):
    # m vertex-disjoint P3s on 3m vertices + random extra edges (hard); valid if members stay induced P3
    n=3*m; E=set(); H=[]
    for i in range(m):
        a,b,c=3*i,3*i+1,3*i+2; E.add((a,b)); E.add((b,c)); H.append((a,b,c))
    cand=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in E and not (i//3==j//3)]
    for p in rng.sample(cand,min(extra,len(cand))): E.add(p)
    I=Inst(n,E,H); return I if I.valid() else None
if __name__=="__main__":
    seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); rng=random.Random(seed)
    res=defaultdict(Counter)
    for _ in range(N):
        I=gen_c1(m,rng.randint(3,4*m),rng)
        if I is None: continue
        o=ILP(I)
        for l in (0,1,2):
            B=len(I.H)+l
            if o.F({},o.N)<=B: continue
            T=list(o.N)
            for q in list(T):
                T2=[x for x in T if x!=q]
                if o.F({},T2)>B: T=T2
            res[l][len(T)]+=1
    print("c=1 m",m,{l:sorted(c.items()) for l,c in res.items()})
