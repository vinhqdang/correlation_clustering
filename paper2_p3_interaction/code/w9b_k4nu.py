import pickle, itertools, sys
from w9b_poly import *
stored=pickle.load(open('/tmp/w9b_k4.pkl','rb'))
rels=dict(lists())
for i,fs in enumerate(stored):
    rels[('K4',i)]=(4,frozenset(fs))
from collections import Counter
print(Counter(len(fs) for fs in stored))
for r in (3,4,5):
    ps=PolySearch(r,rels,nu=True,conservative=True)
    res=ps.solve()
    print('K4 language NU arity',r,'->','NONE' if res is None else 'FOUND', 'clauses',ps.added,flush=True)
