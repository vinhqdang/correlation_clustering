"""Finite check: all constraint types of the zero-excess CSP, closure under M-majority m (all-distinct -> M).
x_e(state) for e in {ab,bc,ac}: 1 iff the pair is together after the member's single edit."""
import itertools
X = {'ab': {'L':0,'R':1,'M':1}, 'bc': {'L':1,'R':0,'M':1}, 'ac': {'L':0,'R':0,'M':1}}
def m(a,b,c):
    if a==b or a==c: return a
    if b==c: return b
    return 'M'
def ok3(x,y,z): return (x+y+z)!=2       # no induced P3 on the triple
S='LRM'
def closed(rel, arity):
    rel=set(rel)
    for u,v,w in itertools.product(rel,repeat=3):
        if tuple(m(u[i],v[i],w[i]) for i in range(arity)) not in rel: return False
    return True
# (1) unary: subsets of states
bad=[sub for r in range(0,4) for sub in itertools.combinations(S,r) if not closed([(s,) for s in sub],1)]
print('non-closed unary subsets:',bad)
# (2) binary types from a triple with soft pairs pq in A (pair eA), qr in B (pair eB), pr hard with static value eps
cnt=0; badb=[]
for eA,eB,eps in itertools.product(X,X,(0,1)):
    rel=[(s,t) for s in S for t in S if ok3(X[eA][s],X[eB][t],eps)]
    cnt+=1
    if not closed(rel,2): badb.append((eA,eB,eps))
print('binary types:',cnt,'non-closed:',badb)
# (3) unary types from a triple with one soft pair (e) and two hard pairs of static values e1,e2
badu=[]
for e,e1,e2 in itertools.product(X,(0,1),(0,1)):
    rel=[(s,) for s in S if ok3(X[e][s],e1,e2)]
    if not closed(rel,1): badu.append((e,e1,e2))
print('unary types non-closed:',badu)
# (4) ternary types: soft triangle (three distinct owners)
badt=[]
for e1,e2,e3 in itertools.product(X,X,X):
    rel=[(s,t,u) for s in S for t in S for u in S if ok3(X[e1][s],X[e2][t],X[e3][u])]
    if not closed(rel,3): badt.append((e1,e2,e3))
print('ternary (soft triangle) types:',27,'non-closed:',len(badt), badt[:6])
# (5) is the ternary relation 2-decomposable (determined by binary projections)? test each
nd=0
for e1,e2,e3 in itertools.product(X,X,X):
    rel=set((s,t,u) for s in S for t in S for u in S if ok3(X[e1][s],X[e2][t],X[e3][u]))
    pr=[{(r[i],r[j]) for r in rel} for i,j in ((0,1),(0,2),(1,2))]
    full={(s,t,u) for s in S for t in S for u in S if (s,t) in pr[0] and (s,u) in pr[1] and (t,u) in pr[2]}
    if full!=rel: nd+=1
print('ternary types not 2-decomposable:',nd)
# internal member constraint: states only; nothing to check
