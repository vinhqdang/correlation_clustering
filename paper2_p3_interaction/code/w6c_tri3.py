"""Triangle of three members, each with a clean anchor (private end vertex): relation R3 and lock-face structure."""
import itertools, numpy as np
from w6c_rel import *
names='LRM'
# vertices: 0=s_QP(shared Q,P) 1=s_PE 2=s_QE 3=pQ 4=pP 5=pE
res_all=[]
cnt=0
for oq in itertools.product((0,2),(0,1)):
  for op in itertools.product((0,2),(0,1)):
    for oe in itertools.product((0,2),(0,1)):
      def mk(priv,shared,o):
          pos,ordr=o; m=[None]*3; m[pos]=priv; rest=[i for i in range(3) if i!=pos]
          sh=list(shared) if ordr==0 else list(shared)[::-1]
          m[rest[0]],m[rest[1]]=sh[0],sh[1]; return tuple(m)
      Q=mk(3,(0,2),oq); P=mk(4,(0,1),op); E=mk(5,(1,2),oe)
      res,hard=relations(6,[Q,P,E])
      hidx={p:i for i,p in enumerate(hard)}
      # pattern with all hard pairs non-edges (j=0) -- all hard pairs involve a private vertex
      for rel,pats in res.items():
          if 0 not in pats: continue
          lock=[0 if o[0]==0 else 1 for o in (oq,op,oe)]
          R=set(rel)
          def face(S): return any(all(s[i]==lock[i] for i in S) for s in R)
          pairs_ok=all(face(S) for S in itertools.combinations(range(3),2))
          triple=face((0,1,2))
          res_all.append((Q,P,E,sorted(R),lock,pairs_ok,triple))
          cnt+=1
          if pairs_ok and not triple: print('NONFACE3',Q,P,E,sorted(R),lock)
print(cnt,'configs with clean-anchor pattern;', sum(1 for r in res_all if r[5] and not r[6]),'have 3-core')
from collections import Counter
print(Counter((len(r[3]),r[5],r[6]) for r in res_all))
