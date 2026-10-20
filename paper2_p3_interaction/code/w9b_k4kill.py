import pickle, itertools, sys
from w9b_poly import *
stored=pickle.load(open('/tmp/w9b_k4.pkl','rb'))
base=dict(lists())
singles=[]
for i,(fs,cfg) in enumerate(stored.items()):
    rels=dict(base); rels[('K4',i)]=(4,frozenset(fs))
    r3=PolySearch(3,rels,nu=True,conservative=True).solve() is None
    if r3:
        r4=PolySearch(4,rels,nu=True,conservative=True).solve() is None
        r5=PolySearch(5,rels,nu=True,conservative=True).solve() is None
        singles.append((i,len(fs),r4,r5))
print('single K4 relations killing majority:',len(singles)); print(singles[:20])
pickle.dump(singles,open('/tmp/w9b_k4single.pkl','wb'))
