import itertools
from w7_rel import X, ok3, S
perms=list(itertools.permutations(S))
def mk(rule):
    def m(a,b,c):
        if a==b or a==c: return a
        if b==c: return b
        return rule[(a,b,c)]
    return m
types=[]
for eA,eB,eps in itertools.product(X,X,(0,1)):
    types.append(((2),[(s,t) for s in S for t in S if ok3(X[eA][s],X[eB][t],eps)]))
for e,e1,e2 in itertools.product(X,(0,1),(0,1)):
    types.append((1,[(s,) for s in S if ok3(X[e][s],e1,e2)]))
good=[]
for vals in itertools.product(S,repeat=6):
    rule=dict(zip(perms,vals)); m=mk(rule)
    allok=True
    for ar,rel in types:
        rel=set(rel)
        for u,v,w in itertools.product(rel,repeat=3):
            if tuple(m(u[i],v[i],w[i]) for i in range(ar)) not in rel: allok=False;break
        if not allok: break
    if allok: good.append(vals)
print(len(good), good[:5])
