"""c<=2 instances with soft triangles: frequency of non-2-decomposable Sigma / non-closure; also cores among clean anchors."""
import random, sys, itertools, pickle
from w7_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); CM=int(sys.argv[3])
tot=0; stats={'tri':0,'closedM':0,'2dec':0}; nd=[]
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
big=[]
it=0
while tot<N:
    it+=1
    n=rng.randint(8,13); m=rng.randint(4,7)
    I=gen_planted(n,m,rng,cmax=CM,maxcl=rng.choice([3,4,5]))
    if I is None or not I.soft_triangles(): continue
    sg=I.sigma()
    if len(sg)<2 or len(sg)>200: continue
    tot+=1; stats['tri']+=1
    if closed_mmaj_fast(sg) if len(sg)<=40 else None: stats['closedM']+=1
    if two_decomposable(sg): stats['2dec']+=1
    else: nd.append((I.A.astype(int).tolist(),I.H,sg))
    an,co=cores(I,sg)
    if len(co)>0 and max(len(c) for c in co)>=3: big.append((I.A.astype(int).tolist(),I.H,co))
print(stats,'non-2-dec examples:',len(nd),'3-core examples:',len(big))
pickle.dump((nd,big),open('/tmp/w7_t2_%s_%s.pkl'%(sys.argv[1],CM),'wb'))
if nd: print(nd[0])
if big: print(big[0])
