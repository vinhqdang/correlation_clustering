"""Find triangle gadgets (3 P3s, 6 vertices) whose projection on two P3s is a bijection, and
verify with an independent set-partition enumeration that OPT=3 iff the relation is non-empty,
and that the set of optimal edit-vectors per P3 matches."""
import itertools,pickle
from hard_rel import pairs
def partitions(n):
    def rec(i,l,m):
        if i==n: yield list(l); return
        for t in range(m+1):
            l.append(t); yield from rec(i+1,l,max(m,t+1)); l.pop()
    yield from rec(0,[],0)
def opt(n,E):
    return min(sum((lab[a]==lab[b])!=(frozenset((a,b)) in E) for a,b in itertools.combinations(range(n),2)) for lab in partitions(n))
c=pickle.load(open("/tmp/cat3.pkl","rb"))
bad=0
for R,(Ps,mask) in c.items():
    V=range(6)
    used={q for X in Ps for q in pairs(X)}
    cross=[frozenset(x) for x in itertools.combinations(V,2) if frozenset(x) not in used]
    E={frozenset(e) for X in Ps for e in ((X[0],X[1]),(X[1],X[2]))}|{cross[t] for t in range(len(cross)) if mask>>t&1}
    if (opt(6,E)==3)!=(len(R)>0): bad+=1
print("gadgets",len(c),"mismatches OPT==3 vs nonempty:",bad)
for R,(Ps,mask) in c.items():
    for a,b in ((0,1),(0,2),(1,2)):
        pr=frozenset((t[a],t[b]) for t in R)
        if len(pr)==3 and len({x for x,_ in pr})==3 and len({y for _,y in pr})==3:
            print("bijection",sorted(pr),"from",sorted(R),Ps,bin(mask))
