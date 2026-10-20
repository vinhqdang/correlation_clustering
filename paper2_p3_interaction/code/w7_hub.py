"""Abstract hub test: mu=(a,b,c) with three neighbours, one at each of its vertices. For a neighbour n sharing vertex q (role of q in mu, role of q in n),
the binary constraint is the conjunction over p in mu\{q}, r in n\{q} of type (pair(p,q) of mu, pair(q,r) of n, eps_pr). Compute all possible relations and
images; look for hub triples with pairwise-intersecting but jointly empty images of single forced neighbour states."""
import itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_rel import X, ok3, S
def pair(u,v):
    return ''.join(sorted((u,v)))  # 'ab','bc','ac' with letters a<b<c
rels=set()   # for (role of q in mu, role of q in n): set of relations (as frozenset of (s,t))
byroles={}
for qm in 'abc':
    for qn in 'abc':
        others_m=[x for x in 'abc' if x!=qm]; others_n=[x for x in 'abc' if x!=qn]
        lst=set()
        for eps in itertools.product((0,1),repeat=4):
            rel=set((s,t) for s in S for t in S)
            k=0
            for p in others_m:
                for r in others_n:
                    eA=pair(p,qm); eB=pair(qn,r)
                    rel={(s,t) for (s,t) in rel if ok3(X[eA][s],X[eB][t],eps[k])}
                    k+=1
            lst.add(frozenset(rel))
        byroles[(qm,qn)]=lst
        print(qm,qn,len(lst))
# images: for neighbour forced to a single state t in {L,R}
def image(rel,t): return frozenset(s for (s,tt) in rel if tt==t)
# hub: neighbour at mu-vertex a,b,c respectively; choose role qn for each and a relation and a forced state t
opts={}
for qm in 'abc':
    o=set()
    for qn in 'abc':
        for rel in byroles[(qm,qn)]:
            for t in 'LR':
                o.add(image(rel,t))
    opts[qm]=o
    print(qm, sorted(map(lambda f:''.join(sorted(f)),o)))
found=[]
for I1 in opts['a']:
    for I2 in opts['b']:
        for I3 in opts['c']:
            if I1&I2 and I1&I3 and I2&I3 and not (I1&I2&I3): found.append((I1,I2,I3))
print('hub 3-core image triples:',len(found), found[:3])
