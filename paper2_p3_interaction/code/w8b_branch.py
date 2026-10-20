"""Branching algorithm on far (non-local) hard pairs, brute-force oracle F(T). See notes/wave8_branch.md.
solve(inst,l): returns (answer, stats). stats: max core size, nodes."""
import numpy as np, itertools, sys, random
from w8b_lib import *

class Solver:
    def __init__(self,inst):
        self.o=Oracle(inst); self.I=inst; self.m=len(inst.H)
        o=self.o
        self.npair=len(inst.pairs)
        self.nodes=0; self.maxcore=0; self.cores=[]
    def prep(self,flip):
        o=self.o
        E=o.Ev^flip
        D=(o.same!=E[None,:])         # P x npairs, disagreement
        cost=D[:,o.counted].sum(1)
        return D,cost
    def F(self,D,cost,T):
        if T: ok=~D[:,T].any(1)
        else: ok=np.ones(len(cost),bool)
        c=cost[ok]
        return int(c.min()) if len(c) else 10**9
    def solve(self,flip,l):
        self.nodes+=1
        o=self.o; B=self.m+l
        D,cost=self.prep(flip)
        N=o.N
        if self.F(D,cost,N)<=B: return True
        if l==0: return False
        # minimal core by deletion
        T=list(N)
        for q in list(T):
            T2=[x for x in T if x!=q]
            if self.F(D,cost,T2)>B: T=T2
        self.maxcore=max(self.maxcore,len(T)); self.cores.append(len(T))
        for e in T:
            f2=flip.copy(); f2[e]=~f2[e]
            if self.solve(f2,l-1): return True
        return False
    def run(self,l):
        self.nodes=0; self.maxcore=0; self.cores=[]
        r=self.solve(np.zeros(self.npair,bool),l)
        return r
def exact_opt(inst):
    o=Oracle(inst); E=o.Ev
    D=(o.same!=E[None,:]); return int(D.sum(1).min())

def check(inst,ls=(0,1,2,3)):
    s=Solver(inst); opt=exact_opt(inst); bad=0; mc=0
    for l in ls:
        r=s.run(l); mc=max(mc,s.maxcore)
        if r!=(opt<=len(inst.H)+l): bad+=1
    return bad,opt-len(inst.H),mc,s.o

def gen_chain(m,extra,rng):
    n=2*m+1
    E=set((i,i+1) for i in range(2*m))
    H=[(2*i,2*i+1,2*i+2) for i in range(m)]
    allp=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in E and j-i>1]
    for p in rng.sample(allp,extra): E.add(p)
    I=Inst(n,E,H)
    return I if I.valid() else None

if __name__=="__main__":
    seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); rng=random.Random(seed)
    tot=bad=0; from collections import Counter; mcs=Counter(); exc=Counter()
    for _ in range(N):
        I=gen_chain(m,rng.randint(1,5),rng)
        if I is None: continue
        b,ex,mc,o=check(I)
        tot+=1; bad+=b; mcs[mc]+=1; exc[ex]+=1
    print("instances",tot,"disagreements",bad,"maxcore dist",dict(mcs),"excess dist",dict(exc))
