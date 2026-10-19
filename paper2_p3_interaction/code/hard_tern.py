"""Ternary relations realised by three pairwise-intersecting packed P3s P={y,x,.},Q={x,z,.},R={y,z,.}
(three distinct shared vertices) with the 6 remaining pairs fixed."""
import itertools, pickle
from hard_rel import relation,pairs
def cat3():
    out={}
    y,x,z,p,q,r=range(6); V=range(6)
    sets=[(y,x,p),(x,z,q),(y,z,r)]
    for o1 in itertools.permutations(sets[0]):
     for o2 in itertools.permutations(sets[1]):
      for o3 in itertools.permutations(sets[2]):
        Ps=[o1,o2,o3]
        used={q_ for X in Ps for q_ in pairs(X)}
        if len(used)<9: continue
        cross=[frozenset(c) for c in itertools.combinations(V,2) if frozenset(c) not in used]
        for mask in range(1<<len(cross)):
            F=[cross[t] for t in range(len(cross)) if mask>>t&1]
            out.setdefault(relation(Ps,F,V),(Ps,mask))
    return out
if __name__=="__main__":
    c=cat3(); print(len(c))
    pickle.dump(c,open("/tmp/cat3.pkl","wb"))
    from collections import Counter
    print(sorted(Counter(len(r) for r in c).items()))
