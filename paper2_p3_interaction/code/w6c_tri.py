"""Soft-triangle lock: P=(x,y,z-ish) with clean anchor x; Q through y and q; E through z and q (6 vertices x,y,z,q,u,v).
Zero-excess relations R(sQ,sP,sE); which hard-pair patterns (no hard edge at x) give a non-trivial locked relation?"""
import sys, itertools, numpy as np, pickle
from w6c_rel import *
names='LRM'
X,Y,Z,Q_,U,V = range(6)
out=[]
for xpos in (0,2):
  for order in ((Y,Z),(Z,Y)):
    P=[None]*3; P[xpos]=X; rest=[i for i in range(3) if i!=xpos]
    P[rest[0]]=order[0]; P[rest[1]]=order[1]
    for pq in itertools.permutations([Y,Q_,U]):
      for pe in itertools.permutations([Z,Q_,V]):
        members=[tuple(pq),tuple(P),tuple(pe)]
        res,hard=relations(6,members)
        hidx={p:i for i,p in enumerate(hard)}
        lock=0 if xpos==0 else 1
        for rel,pats in res.items():
            good=[j for j in pats if not any((j>>hidx[key(X,t)])&1 for t in range(6) if t!=X and key(X,t) in hidx)]
            if not good: continue
            Ra={(s[0],s[2]) for s in rel if s[1]==lock}
            Rb={b:{(s[0],s[2]) for s in rel if s[1]==b} for b in range(3) if b!=lock}
            out.append((xpos,order,pq,pe,sorted(rel),len(Ra),{b:len(v) for b,v in Rb.items()},good[:1]))
from collections import Counter
print(Counter((o[5],tuple(sorted(o[6].values()))) for o in out))
for o in out:
    if o[5]>=2: print(o); break
pickle.dump(out,open(outpath('tri.pkl'),'wb'))
