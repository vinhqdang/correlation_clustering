import sys
from w9b_poly import *
r=int(sys.argv[1]); which=sys.argv[2]
rels={}; rels.update(lists()); st=single_triples()
rels.update({k:v for k,v in st.items() if v[0]<=(2 if which=='binary' else 3)})
def has(keys):
    rr={k:rels[k] for k in keys}
    ps=PolySearch(r,rr,nu=True,conservative=True)
    return ps.solve() is not None
keys=list(rels)
assert not has(keys)
for k in list(keys):
    k2=[x for x in keys if x!=k]
    if not has(k2): keys=k2
print(len(keys)); 
for k in keys: print(k, sorted(rels[k][1]))
