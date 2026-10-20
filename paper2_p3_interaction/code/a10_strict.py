"""Lemma S re-check on real 5-vertex graphs (member A parent, member B child sharing q=0; four cross pairs hard).
R = (sA,sB) passing the four binary triples (p,q,r); LamA = sA passing triples (1,2,r); LamB = sB passing (p,3,4).
For X<=Y<={L,R,M}: X'=LamA & img_R(X & LamB); claim X'!=Y' => |X'|<=1."""
import itertools, numpy as np
from a10_lib import *
def edited(A,H,s):
    B=A.copy()
    for (a,b,c),st in zip(H,s):
        x,y=[(a,b),(b,c),(a,c)][st]; B[x,y]=B[y,x]=not B[x,y]
    return B
def p3(B,t):
    x,y,z=t; return int(B[x,y])+int(B[y,z])+int(B[x,z])==2
subsets=[frozenset(c) for r in range(4) for c in itertools.combinations(range(3),r)]
tot=strict=viol=0; seen=set()
for oA in itertools.permutations([0,1,2]):
  for oB in itertools.permutations([0,3,4]):
    for eps in itertools.product((0,1),repeat=4):
        A=np.zeros((5,5),bool)
        for (a,b,c) in (oA,oB): A[a,b]=A[b,a]=A[b,c]=A[c,b]=True
        for (p,r),e in zip([(1,3),(1,4),(2,3),(2,4)],eps): A[p,r]=A[r,p]=bool(e)
        H=[oA,oB]; R=set();LA=set();LB=set()
        for sA in range(3):
            for sB in range(3):
                B=edited(A,H,(sA,sB))
                if not any(p3(B,(p,0,r)) for p in (1,2) for r in (3,4)): R.add((sA,sB))
        for sA in range(3):
            if not any(p3(edited(A,H,(sA,0)),(1,2,r)) for r in (3,4)): LA.add(sA)
        for sB in range(3):
            if not any(p3(edited(A,H,(0,sB)),(p,3,4)) for p in (1,2)): LB.add(sB)
        key=(frozenset(R),frozenset(LA),frozenset(LB))
        if key in seen: continue
        seen.add(key)
        img=lambda T: frozenset(s for (s,t) in R if t in T)
        for Y in subsets:
            for X in subsets:
                if not X<=Y: continue
                Xp=key[1]&img(X&key[2]); Yp=key[1]&img(Y&key[2]); tot+=1
                if Xp!=Yp:
                    strict+=1; viol+=len(Xp)>1
print("distinct (R,LamA,LamB) configs",len(seen),"pairs X<=Y",tot,"strict",strict,"violations",viol)
