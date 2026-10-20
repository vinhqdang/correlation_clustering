"""Chains of members (consecutive members share one vertex; left/right roles fixed per template) with cross bits eps.
Find templates whose ZE set on a chain of N members consists of exactly 2 'phase' vectors for N=3,4,5 (a rigid wire)."""
import itertools, numpy as np
from w7_lib import Inst, ST
def build(N,rl,rr,eps,extra_hard=0):
    # member k vertices: left vertex Lk, right vertex Rk, third Tk. R_k == L_{k+1} (shared).
    # vertex ids: shared s_0..s_N (s_k = left of member k, s_{k+1} right), thirds t_0..t_{N-1}
    s=list(range(N+1)); t=[N+1+k for k in range(N)]
    n=N+1+N
    A=np.zeros((n,n),bool); H=[]
    for k in range(N):
        roles={rl:s[k],rr:s[k+1]}
        r3=[x for x in 'abc' if x not in (rl,rr)][0]; roles[r3]=t[k]
        a,b,c=roles['a'],roles['b'],roles['c']
        A[a,b]=A[b,a]=A[b,c]=A[c,b]=True; H.append((a,b,c))
    # cross pairs between consecutive members: p in M_k\{s_{k+1}} = {s_k,t_k}, r in M_{k+1}\{s_{k+1}}={s_{k+2},t_{k+1}}
    for k in range(N-1):
        ps=[s[k],t[k]]; rs=[s[k+2],t[k+1]]
        for (p,r),e in zip(itertools.product(ps,rs),eps):
            A[p,r]=A[r,p]=bool(e)
    return Inst(A,H)
res=[]
for rl,rr in itertools.permutations('abc',2):
    for eps in itertools.product((0,1),repeat=4):
        sizes=[]
        for N in (3,4,5):
            I=build(N,rl,rr,eps)
            if not I.valid(): sizes=None;break
            sizes.append(sorted(I.sigma()))
        if sizes is None: continue
        res.append((rl,rr,eps,sizes))
good=[r for r in res if all(len(x)==2 for x in r[3])]
print('templates',len(res),'rigid 2-phase wires:',len(good))
for r in good[:12]: print(r[0],r[1],r[2],r[3][1])
from collections import Counter
c=Counter(tuple(len(x) for x in r[3]) for r in res)
print(c)
for r in res:
    if tuple(len(x) for x in r[3])==(2,2,2) or max(len(x) for x in r[3])>=2:
        print(r[0],r[1],r[2],[len(x) for x in r[3]],r[3][1][:4])
