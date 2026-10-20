"""Polynomial algorithm for l=0 on soft-triangle-free instances: build the binary CSP of Theorem A from the vertex triples and decide it by
path consistency (valid because every relation is majority closed). Compare with brute-force enumeration of Sigma."""
import itertools, random, sys
import numpy as np
from w7_lib import *
X={'ab':{'L':0,'R':1,'M':1},'bc':{'L':1,'R':0,'M':1},'ac':{'L':0,'R':0,'M':1}}
def pairname(t,u,v):
    a,b,c=t
    s={u,v}
    if s=={a,b}: return 'ab'
    if s=={b,c}: return 'bc'
    return 'ac'
def build_csp(I):
    m=len(I.H); n=I.n
    dom=[set(ST) for _ in range(m)]
    rel={}   # (i,j) i<j -> set of allowed (si,sj); absent = all
    def xval(owner,p,q,s):
        return X[pairname(I.H[owner],p,q)][s]
    for p,q,r in itertools.combinations(range(n),3):
        prs=[(p,q),(q,r),(p,r)]
        own=[I.owner.get(key(*pr)) for pr in prs]
        owners=sorted({o for o in own if o is not None})
        if len(owners)==3: raise ValueError('soft triangle')
        if len(owners)==0:
            ed=[int(I.A[a,b]) for a,b in prs]
            if sum(ed)==2: return None   # static induced P3 -> infeasible
            continue
        for states in itertools.product(ST,repeat=len(owners)):
            st=dict(zip(owners,states))
            vals=[]
            for (a,b),o in zip(prs,own):
                vals.append(xval(o,a,b,st[o]) if o is not None else int(I.A[a,b]))
            ok = sum(vals)!=2
            if not ok:
                if len(owners)==1: dom[owners[0]].discard(states[0])
                else:
                    i,j=owners
                    rel.setdefault((i,j),set(itertools.product(ST,ST)))
                    rel[(i,j)].discard(states)
    return dom,rel
def path_consistency(dom,rel,m):
    # relation matrices
    R={}
    for i in range(m):
        for j in range(m):
            if i==j: continue
            a,b=(i,j) if i<j else (j,i)
            allowed=rel.get((a,b),set(itertools.product(ST,ST)))
            R[(i,j)]={(s,t) for (s,t) in allowed if s in dom[a] and t in dom[b]} if i<j else {(t,s) for (s,t) in allowed if s in dom[a] and t in dom[b]}
    # make domains consistent with R
    changed=True
    while changed:
        changed=False
        for i in range(m):
            for j in range(m):
                if i==j: continue
                # revise R[i,j] via every k: R_ij &= R_ik o R_kj
                for k in range(m):
                    if k in (i,j): continue
                    new={(s,t) for (s,t) in R[(i,j)] if any((s,u) in R[(i,k)] and (u,t) in R[(k,j)] for u in ST)}
                    if new!=R[(i,j)]:
                        R[(i,j)]=new; R[(j,i)]={(t,s) for (s,t) in new}; changed=True
                        if not new: return False
        # unary: domain of i must have support in each R[i,j]
        for i in range(m):
            for j in range(m):
                if i==j: continue
                sup={s for (s,t) in R[(i,j)]}
                if not sup: return False
    return True
if __name__=='__main__':
    rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
    agree=0; yes=0; tot=0
    while tot<N:
        m=rng.randint(2,6); n=rng.randint(int(1.4*m)+2,3*m+1)
        # mix planted (feasible) and perturbed (maybe infeasible) instances
        I=gen_planted(n,m,rng,cmax=rng.choice([2,3,4]),maxcl=rng.choice([3,4,5]),tries=1500)
        if I is None or I.soft_triangles(): continue
        A=I.A.copy()
        if rng.random()<0.6:
            # flip one random hard pair (keeps packing valid)
            hp=[(a,b) for a in range(n) for b in range(a+1,n) if key(a,b) not in I.owner]
            for _ in range(rng.randint(1,2)):
                a,b=rng.choice(hp); A[a,b]=A[b,a]=not A[a,b]
            I=Inst(A,I.H)
        tot+=1
        sg=I.sigma(); truth=bool(sg)
        c=build_csp(I)
        if c is None: pred=False
        else:
            dom,rel=c
            pred=all(dom[i] for i in range(len(I.H))) and path_consistency(dom,rel,len(I.H))
        agree+= (pred==truth); yes+=truth
    print('instances',tot,'feasible',yes,'path-consistency agrees with brute force on',agree)
