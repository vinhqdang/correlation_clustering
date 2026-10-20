import pickle,sys,itertools
from w9b_nu2 import *
P,T=pickle.load(open('/tmp/w9b_scopes.pkl','rb'))
tri=[(R,cfg) for R,cfg in T.items() if len(R)]
lists=load('lists_only') if False else None
base=[(1,np.array([[x] for x in L])) for L in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]]
def kills(R,r):
    rels=base+[(3,np.array(sorted(R)))]
    return Search(r,rels).solve(verbose=False) is None
nk=0; kd=[]
for R,cfg in tri:
    if kills(R,3):
        kd.append((R,cfg))
print('killers (no majority):',len(kd),'of',len(tri))
nodrop=[(R,cfg) for R,cfg in kd if None not in cfg[3]]
print('with no dropped pair:',len(nodrop))
k4=[(R,cfg) for R,cfg in kd if kills(R,4)]
print('also kill NU4:',len(k4),'; no dropped:',sum(1 for R,c in k4 if None not in c[3]))
k5=[(R,cfg) for R,cfg in k4 if kills(R,5)]
print('also kill NU5:',len(k5),'; no dropped:',sum(1 for R,c in k5 if None not in c[3]))
pickle.dump((kd,k4,k5),open('/tmp/w9b_killers.pkl','wb'))
for R,cfg in sorted(nodrop,key=lambda x:len(x[0]))[:10]:
    print(len(R),cfg)
