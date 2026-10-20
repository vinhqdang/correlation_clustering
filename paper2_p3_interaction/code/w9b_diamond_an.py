import pickle, itertools, numpy as np
from w9b_poly import *
from w9b_sigma import max_core
rels=pickle.load(open('/tmp/w9b_diamond_rels.pkl','rb'))
L=dict(lists())
for i,fs in enumerate(rels): L[('dia',i)]=(4,frozenset(fs))
for r in (3,4,5):
    ps=PolySearch(r,L,nu=True,conservative=True); res=ps.solve()
    print('diamond language NU',r,'->',res is not None,flush=True)
# core numbers
from collections import Counter
c1=Counter();c2=Counter(); big=[]
for i,fs in enumerate(rels):
    Sg=np.array(sorted(fs),dtype=np.int8)
    k1=max_core(Sg,singleton_only=True); k2=max_core(Sg)
    c1[k1]+=1;c2[k2]+=1
    if k1>=3: big.append((i,k1))
print('k1',sorted(c1.items()),'k2',sorted(c2.items()),'k1>=3:',big)
singles=[]
for i,fs in enumerate(rels):
    d=dict(lists()); d[('dia',i)]=(4,frozenset(fs))
    if PolySearch(3,d,nu=True,conservative=True).solve() is None: singles.append(i)
print('diamond relations killing conservative majority:',len(singles),singles[:20])
