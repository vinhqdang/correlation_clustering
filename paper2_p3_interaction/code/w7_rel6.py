import pickle, itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub3 import cfgs, S
rels=set()
for (qm,qn,eps),(R,LP,LC,E) in cfgs.items():
    rels.add(frozenset(R))
print(len(rels))
def mat(R): return tuple(''.join('1' if (s,t) in R else '.' for t in S) for s in S)
from collections import Counter
cnt=Counter()
for R in sorted(rels,key=lambda r:(len(r),sorted(r))):
    cnt[len(R)]+=1
print(cnt)
for R in sorted(rels,key=lambda r:(len(r),sorted(r))):
    print(len(R),mat(R))
