import itertools,sys
from w7_hub3 import *
f=lambda s:''.join(sorted(s)) or '0'
def branch(rmu,depth):
    out={}
    for crho in 'abc':
        for (TLc,TUc),desc in node(crho,depth):
            for eps in itertools.product((0,1),repeat=4):
                R,LP,LC,E=cfgs[(rmu,crho,eps)]
                if desc[0]=='leaf':
                    e=desc[2]
                    if any(E[(p,e)]==1 for p in [x for x in 'abc' if x!=rmu]): continue
                CL=LP & img(R,TLc&LC); CU=LP & img(R,TUc&LC)
                out.setdefault((CL,CU),1)
    return out
for r in 'abc':
    b=branch(r,4)
    print(r,sorted((f(a),f(c)) for a,c in b))
for rho in 'abc':
    print('node',rho,sorted((f(a),f(c)) for (a,c),d in node(rho,4)))
