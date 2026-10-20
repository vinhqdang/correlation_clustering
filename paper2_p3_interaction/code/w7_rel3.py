"""Per-coordinate majorities: for each binary constraint type, the set of rule pairs (rA,rB) under which the relation is closed.
A rule = value for each ordered all-distinct triple of states."""
import itertools, sys, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_rel import X, ok3, S
perms=list(itertools.permutations(S))
rules=list(itertools.product(S,repeat=6))
def ap(r,a,b,c):
    if a==b or a==c: return a
    if b==c: return b
    return r[perms.index((a,b,c))]
types={}
for eA,eB,eps in itertools.product(X,X,(0,1)):
    types[(eA,eB,eps)]=[(s,t) for s in S for t in S if ok3(X[eA][s],X[eB][t],eps)]
unar={}
for e,e1,e2 in itertools.product(X,(0,1),(0,1)):
    unar[(e,e1,e2)]=[s for s in S if ok3(X[e][s],e1,e2)]
# unary types are closed under all majorities (checked in w7_rel).
compat={}
for ty,rel in types.items():
    rs=set(rel); trip=list(itertools.product(rel,repeat=3))
    # allowed single rules (conditions from triples where the other coordinate is not all-distinct)
    okA=[r for r in rules if all(((ap(r,u[0],v[0],w[0]), ap(rules[0],u[1],v[1],w[1])) in rs) or (len({u[1],v[1],w[1]})==3) for u,v,w in trip)]
    # simpler: do exact pair test on filtered sets
    pairs=[]
    candA=[r for r in rules if all(True for _ in [0])]
    compat[ty]=None
    # exact: for all (rA,rB)
    res=set()
    # precompute per triple the coordinate values
    for rA in rules:
        for rB in rules:
            good=True
            for u,v,w in trip:
                if (ap(rA,u[0],v[0],w[0]),ap(rB,u[1],v[1],w[1])) not in rs: good=False;break
            if good: res.add((rA,rB))
    compat[ty]=res
    print(ty,len(rel),len(res),flush=True)
import pickle; pickle.dump(compat,open('/tmp/w7_compat.pkl','wb'))
