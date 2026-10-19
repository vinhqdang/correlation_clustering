"""Binary relations on {0,1,2}^2 realised by two packed P3s sharing exactly one vertex.
P=(a,b,c): edges ab,bc, non-edge ac. choice 0: delete ab, 1: delete bc, 2: add ac.
The 4 cross pairs are unpacked and fixed (edge / non-edge). Relation = set of (x,y)
for which the edited 5-vertex graph is a cluster graph."""
import itertools
def is_cluster(V,E):
    for a,b,c in itertools.permutations(V,3):
        if frozenset((a,b)) in E and frozenset((b,c)) in E and frozenset((a,c)) not in E: return False
    return True
def pairs(P):
    a,b,c=P; return [frozenset((a,b)),frozenset((b,c)),frozenset((a,c))]
def relation(Ps,fixed_edges,V):
    base=set(fixed_edges)
    for a,b,c in Ps: base|={frozenset((a,b)),frozenset((b,c))}
    R=set()
    for ch in itertools.product(range(3),repeat=len(Ps)):
        E=set(base)
        for P,c in zip(Ps,ch): E^={pairs(P)[c]}
        if is_cluster(V,E): R.add(ch)
    return frozenset(R)
if __name__=="__main__":
    rels={}
    for i in range(3):
      for j in range(3):
        P=(0,1,2); sh=P[i]; q=[3,4]; Pp=[None]*3; Pp[j]=sh; k=0
        for t in range(3):
            if t!=j: Pp[t]=q[k]; k+=1
        Pp=tuple(Pp); V=range(5)
        used={q_ for X in (P,Pp) for q_ in pairs(X)}
        cross=[frozenset(p) for p in itertools.combinations(V,2) if frozenset(p) not in used]
        assert len(cross)==4
        for mask in range(16):
            F=[cross[t] for t in range(4) if mask>>t&1]
            rels.setdefault(relation([P,Pp],F,V),[]).append((i,j,mask))
    for R,w in sorted(rels.items(),key=lambda x:len(x[0])):
        print(len(R),sorted(R),"e.g.",w[0],"count",len(w))
