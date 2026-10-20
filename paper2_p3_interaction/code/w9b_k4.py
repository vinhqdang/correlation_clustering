import itertools, pickle
from collections import Counter
from w9b_lib import *
from w9b_sigma import *
# K4 gadget: members indexed 0..3; vertex v_{ij} (i<j) shared by members i,j ; 6 vertices
pairsK=list(itertools.combinations(range(4),2)); vid={p:k for k,p in enumerate(pairsK)}
def members_from_roles(roles):
    H=[]
    for i in range(4):
        vs=[vid[tuple(sorted((i,j)))] for j in range(4) if j!=i]
        a,b,c=[vs[k] for k in roles[i]]
        H.append((a,b,c))
    return H
perms=list(itertools.permutations(range(3)))
rels=Counter(); stored={}
for roles in itertools.product(perms,repeat=4):
    H=members_from_roles(roles)
    soft=set()
    for (a,b,c) in H: soft|={key(a,b),key(b,c),key(a,c)}
    hard=[p for p in itertools.combinations(range(6),2) if p not in soft]
    for st in itertools.product((0,1),repeat=len(hard)):
        E=set()
        for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
        for p,s in zip(hard,st):
            if s: E.add(p)
        I=Inst(6,E,H)
        if not I.valid(): continue
        S=sigma(I,drop_far=False)
        if len(S)==0: continue
        fs=frozenset(map(tuple,S.tolist()))
        stored.setdefault(fs,(H,st))
print('distinct nonempty K4 Sigma',len(stored))
pickle.dump(stored,open('/tmp/w9b_k4.pkl','wb'))
