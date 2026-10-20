import sys; sys.path.insert(0,'.')
from r1_lib import *
# C4 0-1-3-2 with vertex 2 replaced by a clique class X of size k (all adjacent to 0 and 3), H=(0,1,3)
for k in (1,2,3):
    n=3+k; X=list(range(3,3+k))
    E={frozenset(e) for e in [(0,1),(1,2)]}  # relabel: vertices 0,1,2 =0,1,3 ; X = copies of vertex 2
    E={frozenset(e) for e in [(0,1),(1,2)]}
    for x in X: E|={frozenset((0,x)),frozenset((2,x))}
    for a,b in itertools.combinations(X,2): E.add(frozenset((a,b)))
    H=[(0,1,2)]
    m=opt(n,E)
    keep=[v for v in range(n) if v!=X[-1]]; mp={v:i for i,v in enumerate(keep)}
    E2={frozenset((mp[a],mp[b])) for e in E for a,b in [tuple(e)] if a!=X[-1] and b!=X[-1]}
    print('k',k,'OPT(G)',m,'OPT(G-x)',opt(n-1,E2),'|H|=1')
