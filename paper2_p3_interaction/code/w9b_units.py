import itertools, pickle
from w9b_lib import *
from w9b_sigma import sigma
from w9b_octchain import unit_members, perms
units=[]
for roles in itertools.product(perms,repeat=4):
    H=unit_members(roles)
    soft=set()
    for (a,b,c) in H: soft|={key(a,b),key(b,c),key(a,c)}
    hard=[p for p in itertools.combinations(range(6),2) if p not in soft]
    for st in itertools.product((0,1),repeat=3):
        E=set()
        for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
        for p,s in zip(hard,st):
            if s: E.add(p)
        I=Inst(6,E,H)
        if not I.valid(): continue
        S=sigma(I,drop_far=False)
        if len(S)>=1: units.append((H,E,roles,st,[tuple(r) for r in S.tolist()]))
print(len(units),'feasible unit configs; with |S|>=2:',sum(1 for u in units if len(u[4])>=2))
pickle.dump(units,open('/tmp/w9b_units.pkl','wb'))
