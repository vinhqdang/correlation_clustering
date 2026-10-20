import itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_hub import byroles, S
    from w7_rel import X, ok3
allrels=set()
for k,v in byroles.items(): allrels|=v
print('distinct relations:',len(allrels))
def img(rel,T): return frozenset(s for (s,t) in rel if t in T)
fam={frozenset('L'),frozenset('R')}
changed=True
while changed:
    changed=False
    for T in list(fam):
        for rel in allrels:
            J=img(rel,T)
            if J not in fam: fam.add(J); changed=True
    # intersections with static lists / other branches
print(sorted(''.join(sorted(f)) for f in fam))
# also with static unary lists intersection and intersections of two branches
lists=[frozenset(x) for x in ['L','R','M','LR','LM','RM','LMR']]
fam2=set(fam)
changed=True
while changed:
    changed=False
    for A in list(fam2):
        for B in list(fam2)+lists:
            if (A&B) not in fam2: fam2.add(A&B); changed=True
        for rel in allrels:
            J=img(rel,A)
            if J not in fam2: fam2.add(J); changed=True
print(sorted(''.join(sorted(f)) for f in fam2))

# hub triples with neighbour forced sets in fam
opts={}
for qm in 'abc':
    o=set()
    for qn in 'abc':
        for rel in byroles[(qm,qn)]:
            for T in fam:
                o.add(img(rel,T))
    opts[qm]=o
    print(qm,sorted(''.join(sorted(f)) for f in o))
found=[(I1,I2,I3) for I1 in opts['a'] for I2 in opts['b'] for I3 in opts['c'] if I1&I2 and I1&I3 and I2&I3 and not (I1&I2&I3)]
print('hub 3-core triples (vertices a,b,c each hosting one branch):',len(found))
for f in found[:5]: print([''.join(sorted(x)) for x in f])
