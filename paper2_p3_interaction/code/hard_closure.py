import itertools
from hard_rel import relation,pairs
def catalogue():
    out=set()
    for i in range(3):
      for j in range(3):
        P=(0,1,2); sh=P[i]; q=[3,4]; Pp=[None]*3; Pp[j]=sh; k=0
        for t in range(3):
            if t!=j: Pp[t]=q[k]; k+=1
        Pp=tuple(Pp); V=range(5)
        used={q_ for X in (P,Pp) for q_ in pairs(X)}
        cross=[frozenset(p) for p in itertools.combinations(V,2) if frozenset(p) not in used]
        for mask in range(16):
            F=[cross[t] for t in range(4) if mask>>t&1]
            out.add(relation([P,Pp],F,V))
    return out
def T(R): return frozenset((y,x) for x,y in R)
def comp(R,S): return frozenset((x,y) for x,z in R for z2,y in S if z==z2)
cat=catalogue(); cat|={T(r) for r in cat}
full=frozenset(itertools.product(range(3),repeat=2))
S=set(cat)|{full}
changed=True
while changed:
    changed=False
    L=list(S)
    for R in L:
        for Q in L:
            for X in (comp(R,Q), R&Q):
                if X not in S: S.add(X); changed=True
EQ=frozenset((i,i) for i in range(3)); NE=frozenset((i,j) for i in range(3) for j in range(3) if i!=j)
print(len(cat),len(S),"EQ" ,EQ in S,"NE",NE in S)
# pins/unary: projections
for name,R in (("EQ",EQ),("NE",NE)): print(name, R in S)
print(sorted(len(r) for r in S)[-10:])
print("size>=6:")
for r in S:
    if len(r)>=6: print(len(r),sorted(r))
bij=[r for r in S if len(r)==3 and len({x for x,_ in r})==3 and len({y for _,y in r})==3]
print("full bijections in closure:",[sorted(b) for b in bij])
