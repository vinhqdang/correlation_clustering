import sys; sys.path.insert(0,'.')
from r1_lib import *
names=[0,1,3,4,5,6,7]; mp={v:i for i,v in enumerate(names)}
E={frozenset((mp[int(s[0])],mp[int(s[1])])) for s in "01 04 05 14 15 45 03 13 36 37 06 56".split()}
H=[tuple(mp[v] for v in t) for t in [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]]
n=7; used=set()
for t in H:
    x,y,z=t; assert frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E
    assert not(pairs3(t)&used); used|=pairs3(t)
print('pair-maximal', not [t for t in p3s(n,E) if not(pairs3(t)&used)])
ap=allparts_cost(n,E); m=min(c for c,_ in ap); print('OPT',m)
inv={i:v for v,i in mp.items()}
for c,p in ap:
    if c==m: print([[inv[v] for v in b] for b in p])
