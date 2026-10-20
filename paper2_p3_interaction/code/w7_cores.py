import itertools
from w7_lib import *
def cores(I,sg):
    """clean anchors: vertices with no hard edge and in >=1 member (and not centre-only?). face(S) iff some sigma makes all of S singletons."""
    n=I.n
    anchors=[v for v in range(n) if I.hard_deg(v)==0 and any(v in t for t in I.H)]
    # singleton in solution sg: compute clusters
    sing=[]
    for s in sg:
        A=I.A.copy()
        for (a,b,c),st in zip(I.H,s):
            u,v={'L':(a,b),'R':(b,c),'M':(a,c)}[st]
            A[u,v]=not A[u,v]; A[v,u]=not A[v,u]
        sing.append(frozenset(v for v in anchors if not A[v].any()))
    faces=lambda S: any(S<=x for x in sing)
    out=[]
    good=[v for v in anchors if faces(frozenset([v]))]
    for r in (2,3,4):
        for S in itertools.combinations(good,r):
            S=frozenset(S)
            if not faces(S) and all(faces(S-{x}) for x in S): out.append(tuple(sorted(S)))
    return anchors,out
