import pickle, itertools
from w9b_poly import *
from w9b_scope import *
# K3 relations: all statuses 0/1 (no dropped)
hardpairs=[(2,4),(2,5),(4,5),(2,3),(4,0),(5,1)]
K3=set()
for A in role_orders((0,1,2)):
    for B in role_orders((1,3,4)):
        for C in role_orders((3,0,5)):
            for st in itertools.product((0,1),repeat=6):
                hard={frozenset(p):s for p,s in zip(hardpairs,st)}
                R=scope_relation([A,B,C],hard,6)
                if len(R): K3.add(R)
print('K3 block relations (nonempty):',len(K3))
dia=pickle.load(open('/tmp/w9b_diamond_rels.pkl','rb'))
bin40=language('sf')
rels=dict(lists())
for k,v in bin40.items():
    if v[0]==2: rels[k]=v
for i,R in enumerate(K3): rels[('K3',i)]=(3,frozenset(R))
for i,R in enumerate(dia): rels[('dia',i)]=(4,frozenset(R))
print('language size',len(rels))
for r in (3,4):
    ps=PolySearch(r,rels,nu=True,conservative=True); res=ps.solve()
    print('c=2 language (K3,diamond,links,lists): conservative NU of arity',r,'->',res is not None,flush=True)
    if r==3 and res is not None:
        pickle.dump(res,open('/tmp/w9b_c2_maj.pkl','wb'))
