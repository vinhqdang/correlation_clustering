"""Triangle hub: members A,B,C pairwise sharing u=AB, v=AC, w=BC; thirds a3,b3,c3 each carry a branch (chain of members ending in an anchor).
Sigma_tri from brute force on the 6-vertex graph (all 216 orders x 64 hard-pair bits). Branch image families from chain calculus (w7_hub3), generously
closed under intersection with arbitrary unary lists. Look for a 3-core (all locked infeasible, each pair locked feasible, all unlocked feasible)."""
import itertools, io, contextlib, sys
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub3 import node, cfgs, img, D, S
    from w7_lib import Inst
f=lambda s:''.join(sorted(s)) or '0'
def branch_family(rmu,depth=4):
    out={}
    for crho in 'abc':
        for (TLc,TUc),desc in node(crho,depth):
            for eps in itertools.product((0,1),repeat=4):
                R,LP,LC,E=cfgs[(rmu,crho,eps)]
                if desc[0]=='leaf':
                    e=desc[2]
                    if any(E[(p,e)]==1 for p in [x for x in 'abc' if x!=rmu]): continue
                CL=LP & img(R,TLc&LC); CU=LP & img(R,TUc&LC)
                out[(CL,CU)]=1
    return out
fam={r:branch_family(r) for r in 'abc'}
lists=[frozenset(x) for r in range(1,4) for x in itertools.combinations(S,r)]
def closure(F):
    out=set()
    for (L,U) in F:
        for lam in lists+[D]:
            out.add((L&lam,U&lam))
    return out
famc={r:closure(fam[r]) for r in 'abc'}
# also the "no branch" option (unlocked, no lock): pair (D,D) with lists
for r in 'abc':
    if r=='a': famc[r].add((frozenset('L'),D))
    if r=='c': famc[r].add((frozenset('R'),D))
    famc[r].add((D,D))
famc_none=None
print({r:len(famc[r]) for r in 'abc'})
verts=['u','v','w','a3','b3','c3']
idx={x:i for i,x in enumerate(verts)}
mem={'A':('u','v','a3'),'B':('u','w','b3'),'C':('v','w','c3')}
hardpairs=[('a3','b3'),('a3','c3'),('b3','c3'),('a3','w'),('b3','v'),('c3','u')]
stlist=[('L',(0,1)),('R',(1,2)),('M',(0,2))]
sigmas={}  # (orderA,orderB,orderC) -> {eps tuple -> set of (sA,sB,sC)}
for oA in itertools.permutations(mem['A']):
    for oB in itertools.permutations(mem['B']):
        for oC in itertools.permutations(mem['C']):
            H=[tuple(idx[x] for x in o) for o in (oA,oB,oC)]
            for eps in itertools.product((0,1),repeat=6):
                A=np.zeros((6,6),bool)
                for (a,b,c) in H:
                    A[a,b]=A[b,a]=A[b,c]=A[c,b]=True
                for (x,y),e in zip(hardpairs,eps):
                    A[idx[x],idx[y]]=A[idx[y],idx[x]]=bool(e)
                I=Inst(A,H)
                if not I.valid(): continue
                sg=I.sigma()
                if sg:
                    sigmas[(oA,oB,oC,eps)]=set(tuple(s) for s in sg)
print('valid tight triangle configs:',len(sigmas))
# role of third vertex in each member's order: position 0->'a',1->'b',2->'c'
posrole=lambda o,x:'abc'[o.index(x)]
found=[]
seen=set()
for (oA,oB,oC,eps),sig in sigmas.items():
    rA=posrole(oA,'a3'); rB=posrole(oB,'b3'); rC=posrole(oC,'c3')
    key=(rA,rB,rC,frozenset(sig))
    if key in seen: continue
    seen.add(key)
    for (aL,aU) in famc[rA]:
        for (bL,bU) in famc[rB]:
            for (cL,cU) in famc[rC]:
                def feas(sa,sb,sc): return any((x in sa and y in sb and z in sc) for (x,y,z) in sig)
                # states are letters; sig entries are tuples of chars
                if not feas(aU,bU,cU): continue
                if feas(aL,bL,cL): continue
                if feas(aL,bL,cU) and feas(aL,bU,cL) and feas(aU,bL,cL):
                    found.append((rA,rB,rC,sorted(sig),(f(aL),f(aU)),(f(bL),f(bU)),(f(cL),f(cU))))
print('distinct (roles,Sigma) systems:',len(seen),' triangle-hub 3-cores found:',len(found))
for x in found[:5]: print(x)
