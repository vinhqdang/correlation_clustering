import sys, itertools; sys.path.insert(0,'.')
from r1_lib import *
def edges(s): return set(frozenset((int(t[0]),int(t[1]))) for t in s.split())
def analyse(n,E,H):
    used=set()
    for t in H: used|=pairs3(t)
    adj={v:{u for u in range(n) if frozenset((u,v)) in E} for v in range(n)}
    hard={v:{u for u in adj[v] if frozenset((u,v)) not in used} for v in range(n)}
    nbad=0
    for t in H:
        x,y,z=t; allbad=True
        for st in ([frozenset((x,y))],[frozenset((y,z))],[frozenset((x,z))]):
            bad=False
            for pr in pairs3(t):
                a_,b_=tuple(pr); merged=((pr in E)!=(pr in st))
                if merged: w=[q for q in range(n) if q not in (a_,b_) and frozenset((a_,q)) not in used and frozenset((b_,q)) not in used and ((q in adj[a_])!=(q in adj[b_]))]
                else: w=list(hard[a_]&hard[b_])
                if w: bad=True
            if not bad: allbad=False
        nbad+=allbad
    return nbad
E=edges('01 06 12 23 29 34 45 56 67 78 79'); H=[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(2,9,7)]
print('nfbad2 locally bad members',analyse(10,E,H))
# distance of members containing 0 and 6 in I(H)
def Idist(H,u,v):
    mem=lambda x:[i for i,t in enumerate(H) if x in t]
    import collections
    adjm={i:[j for j,s in enumerate(H) if j!=i and set(s)&set(H[i])] for i in range(len(H))}
    best=99
    for s in mem(u):
        d={s:0};q=[s]
        for x in q:
            for y in adjm[x]:
                if y not in d: d[y]=d[x]+1;q.append(y)
        for tt in mem(v): best=min(best,d.get(tt,99))
    return best
print('dist 0..6',Idist(H,0,6))
E=edges('01 08 12 14 23 34 45 49 56 67 78'); H=[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(1,4,9)]
print('nfbad3 locally bad',analyse(10,E,H),'dist 0..8',Idist(H,0,8))
