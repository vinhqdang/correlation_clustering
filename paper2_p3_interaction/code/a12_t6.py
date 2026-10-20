"""Claim (6): core step of Theorem F: for U=empty, r<Ex(empty): minimal site set C with Ex_C > r-1; every nonempty far W with |W|+Ex(W)<=r has an endpoint in V(C)."""
import sys, itertools, random
from a12_lib import *
from a12_t3 import Inst, isP3, gen
S5=['L','R','M','O','X']; EXTRA={'L':0,'R':0,'M':0,'O':1,'X':2}
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
def ex_min(ins,E,keep,limit):
    """min extra cost <=limit over (x,F) s.t. triples in keep are non-P3 (keep=None: whole graph cluster). returns cost or None"""
    cr=sorted(ins.cross&set(ins.hard))
    best=None
    for kF in range(0,limit+1):
        for F in itertools.combinations(cr,kF):
            EF=set(E)^set(F)
            for xs in itertools.product(S5,repeat=len(ins.H)):
                c=kF+sum(EXTRA[s] for s in xs)
                if c>limit or (best is not None and c>=best): continue
                E2=set(EF)
                for t,s in zip(ins.H,xs): E2=apply_state(E2,t,ST[s])
                ok=is_cluster(ins.n,E2) if keep is None else all(not isP3(E2,T) for T in keep)
                if ok: best=c
    return best
fails=0; tested=0; nontriv=0; csizes={}
for _ in range(N):
    n=rng.randint(6,8); E,H=gen(rng,n,3,rng.choice([.25,.35]))
    ins=Inst(n,E,H)
    if any(all(p in ins.hard for p in prs(t)) for t in p3s(n,ins.E0)): continue
    r=int(sys.argv[3])
    ex0=ex_min(ins,ins.E0,None,r)
    if ex0 is not None: continue          # need Ex(empty) > r
    # sites
    sites=[('m',i) for i in range(len(H))]+[('e',i,j) for i in range(len(H)) for j in range(i+1,len(H)) if set(H[i])&set(H[j])]
    crossedge={}
    for i in range(len(H)):
        for j in range(i+1,len(H)):
            for q in set(H[i])&set(H[j]):
                for p in set(H[i])-{q}:
                    for rr in set(H[j])-{q}: crossedge.setdefault(K(p,rr),set()).add(('e',i,j))
    # triple -> associated sites
    def assoc(T):
        ps=[K(a,b) for a,b in itertools.combinations(T,2)]
        A=set()
        for p in ps:
            if p in ins.soft:
                for i,t in enumerate(H):
                    if p in prs(t): A.add(('m',i))
            if p in crossedge: A|=crossedge[p]
        # members touching T's far pair endpoints are NOT added: group triples are attached to the member owning the soft pair
        return A
    chk=[T for T in ins.triples if any(K(a,b) in ins.soft or K(a,b) in crossedge for a,b in itertools.combinations(T,2))]
    hasfar=lambda T:any(K(a,b) in ins.far for a,b in itertools.combinations(T,2))
    def keep_for(C):
        return [T for T in chk if (not hasfar(T)) or (assoc(T)&set(C))]
    t=r-1
    full=set(sites)
    if ex_min(ins,ins.E0,keep_for(full),t) is not None: continue   # Ex_full > t required (should equal Ex>r in this setting)
    C=list(sites)
    for s in list(C):
        C2=[x for x in C if x!=s]
        if ex_min(ins,ins.E0,keep_for(C2),t) is None: C=C2
    VC=set()
    for s in C:
        if s[0]=='m': VC|=set(H[s[1]])
        else: VC|=set(H[s[1]])|set(H[s[2]])
    csizes[len(C)]=csizes.get(len(C),0)+1
    # all W nonempty far with |W|+Ex(W)<=r
    tested+=1
    for kW in range(1,r+1):
        for W in itertools.combinations(ins.far,kW):
            ex=ex_min(ins,ins.E0^set(W),None,r-kW)
            if ex is None: continue
            nontriv+=1
            if not any(set(p)&VC for p in W): fails+=1; print('FAIL',n,sorted(E),H,W,C)
print('tested',tested,'good W found',nontriv,'fails',fails,'|C| histogram',csizes)
