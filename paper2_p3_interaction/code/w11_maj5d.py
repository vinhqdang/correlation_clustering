import sys, itertools, random, time
import numpy as np
sys.argv=['x','3']
from w11_maj5 import *
from pysat.solvers import Cadical153
np.random.seed(1)
pats=pat_vars()
def run(rels, NO_):
    solver=Cadical153(); var={}
    for p in pats:
        for o in range(NO_): var[(p,o)]=len(var)+1
        solver.add_clause([var[(p,o)] for o in range(NO_)])
        for o1 in range(NO_):
            for o2 in range(o1+1,NO_): solver.add_clause([-var[(p,o1)],-var[(p,o2)]])
    def base_out(a,b,c):
        if a==b or a==c: return a
        if b==c: return b
        return None
    for R in rels:
        U,V,W=check_relation(R,None)
        k=R.shape[1]; Rset=set(map(tuple,R.tolist()))
        seen=set()
        for u,v,w in zip(U.tolist(),V.tolist(),W.tolist()):
            ps=[(u[c],v[c],w[c]) for c in range(k)]
            free=[c for c in range(k) if base_out(*ps[c]) is None]
            keyp=tuple(ps)
            if keyp in seen: continue
            seen.add(keyp)
            for vals in itertools.product(range(NO_),repeat=len(free)):
                out=[base_out(*ps[c]) for c in range(k)]
                for c,o in zip(free,vals): out[c]=o
                if tuple(out) not in Rset:
                    solver.add_clause([-var[(ps[c],o)] for c,o in zip(free,vals)])
    return solver.solve()
rels=[]; seen=set()
for (mem,nv,hard,X) in configs_link():
    R=rel5(mem,nv,hard,X)
    if len(R)==0: continue
    key=R.tobytes()
    if key in seen: continue
    seen.add(key); rels.append(R)
print('link relations',len(rels))
print('SAT with non-defective outputs (3 options):', run(rels,3))
print('SAT with any outputs (5 options):', run(rels,5))

bad=[]
for i,R in enumerate(rels):
    if not run([R],3): bad.append(i)
print('individually UNSAT with 3 outputs:',len(bad),'of',len(rels))
# minimal unsat subfamily (greedy)
if not run(rels,3):
    cur=list(range(len(rels)))
    for i in list(cur):
        t=[j for j in cur if j!=i]
        if not run([rels[j] for j in t],3): cur=t
    print('minimal unsat subfamily size',len(cur))
    names=['L','R','M','4','5']
    for j in cur[:3]:
        R=rels[j]; print('relation',j,'size',len(R),[''.join(names[x] for x in r) for r in R.tolist()])
