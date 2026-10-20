import sys, itertools; sys.path.insert(0,'.')
from r1_lib import *
# K4 oriented: edges (p->q)
Fe=[(0,1),(1,2),(2,0),(0,3),(3,1),(2,3)]
indeg={v:sum(1 for p,q in Fe if q==v) for v in range(4)}; print('indeg',indeg)
V={};
def vid(k):
    if k not in V: V[k]=len(V)
    return V[k]
a=[vid(('a',i)) for i in range(6)]; b=[vid(('b',i)) for i in range(6)]; c=[vid(('c',i)) for i in range(6)]
E=set()
for v in range(4):
    K=[a[i] for i,(p,q) in enumerate(Fe) if p==v]+[x for i,(p,q) in enumerate(Fe) if q==v for x in (b[i],c[i])]
    for pr in itertools.combinations(K,2): E.add(frozenset(pr))
for i in range(6): E.add(frozenset((a[i],b[i])))
n=len(V); H=[(a[i],b[i],c[i]) for i in range(6)]
used=set()
for t in H:
    x,y,z=t; assert frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E; used|=pairs3(t)
print('pair-max',not [t for t in p3s(n,E) if not(pairs3(t)&used)])
adj={v:{u for u in range(n) if frozenset((u,v)) in E} for v in range(n)}
hard={v:{u for u in adj[v] if frozenset((u,v)) not in used} for v in range(n)}
print('max hard degree',max(len(hard[v]) for v in range(n)))
def Delta(z,bb): return [y for y in range(n) if y not in(z,bb) and frozenset((z,y)) not in used and frozenset((bb,y)) not in used and ((y in adj[z])!=(y in adj[bb]))]
for i,(p,q) in enumerate(Fe):
    print(i,len(Delta(a[i],b[i])),3+indeg[p]+indeg[q], len(Delta(a[i],c[i])))
