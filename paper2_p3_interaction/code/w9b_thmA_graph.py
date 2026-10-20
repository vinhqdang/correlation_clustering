"""graph-level check of Theorem A: binary scope relations of two adjacent members (all 4 cross pairs hard, all orientations/statuses, no soft triangle)
are closed under every uniform majority operation."""
import itertools
from w9b_scope import *
D=(0,1,2); perms=list(itertools.permutations(D))
cross=[(1,3),(1,4),(2,3),(2,4)]
S=set()
for A in role_orders((0,1,2)):
    for B in role_orders((0,3,4)):
        for st in itertools.product((0,1),repeat=4):
            hard={frozenset(p):s for p,s in zip(cross,st)}
            S.add(scope_relation([A,B],hard,5))
print('distinct graph-level binary relations (nonempty+empty):',len(S))
bad=0
for out in itertools.product(D,repeat=6):
    def mj(x,y,z):
        if x==y or x==z: return x
        if y==z: return y
        return out[perms.index((x,y,z))]
    for R in S:
        for t1 in R:
            for t2 in R:
                for t3 in R:
                    if tuple(mj(t1[i],t2[i],t3[i]) for i in range(2)) not in R: bad+=1
print('violations over all 729 majorities and all relations:',bad)
