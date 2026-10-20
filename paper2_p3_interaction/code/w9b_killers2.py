import pickle,sys,itertools
from w9b_nu2 import *
P,T=pickle.load(open('/tmp/w9b_scopes.pkl','rb'))
base=[(1,np.array([[x] for x in L])) for L in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]]
res=[]
for R,cfg in T.items():
    if not len(R): continue
    rels=base+[(3,np.array(sorted(R)))]
    if Search(3,rels).solve(verbose=False) is None:
        res.append((sum(1 for s in cfg[3] if s is None),len(R),R,cfg))
res.sort(key=lambda x:(x[0],x[1]))
print(len(res),'killers; fewest dropped:',res[0][0])
from collections import Counter
print(Counter(x[0] for x in res))
for d,l,R,cfg in res[:6]:
    print(d,l,cfg)
pickle.dump(res,open('/tmp/w9b_killers2.pkl','wb'))
