"""verify: every uniform majority on {L,R,M} preserves all 40 single-conjunction binary relations (soft-triangle-free language) and all lists."""
import itertools
from w9b_poly import *
rels=language('sf')
D=(0,1,2)
perms=list(itertools.permutations(D))
ok_all=0
for out in itertools.product(D,repeat=6):
    m={}
    for a in D:
        for b in D:
            if a==b: continue
            for c in D:
                pass
    def mj(x,y,z):
        if x==y or x==z: return x
        if y==z: return y
        return out[perms.index((x,y,z))]
    good=True
    for name,(k,R) in rels.items():
        Rl=list(R)
        for t1 in Rl:
            for t2 in Rl:
                for t3 in Rl:
                    img=tuple(mj(t1[i],t2[i],t3[i]) for i in range(k))
                    if img not in R: good=False;break
                if not good:break
            if not good:break
        if not good:break
    ok_all+=good
print('majorities preserving the whole soft-triangle-free language:',ok_all,'of',3**6)
