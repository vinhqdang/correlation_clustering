"""Bounded Certificate (BC): minimal C (members) such that relaxing all unary constraints outside C keeps opt > t.  Adversarial search for large minimal C."""
import sys, random, itertools, time, ctypes
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import _lib, cost_relaxed
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); m=int(sys.argv[4]); t=int(sys.argv[5]); topo=sys.argv[6]
rng=random.Random(seed)
INF=100000
def prep(n,E,H):
    adj=bytearray(n*n)
    for u,v in E: adj[u*n+v]=1; adj[v*n+u]=1
    Hin=(ctypes.c_int*(3*len(H)))(*[x for t_ in H for x in t_])
    return _lib.build_only(n,bytes(adj),len(H),Hin)
def minimal_C(mm,tt,tries=3):
    full=(1<<mm)-1
    best=0
    for _ in range(tries):
        C=full
        order=list(range(mm)); rng.shuffle(order)
        for a in order:
            C2=C&~(1<<a)
            relax=full&~C2
            v=cost_relaxed(relax)
            if v>tt: C=C2
        best=max(best,bin(C).count('1'))
    return best
def score(n,E,H):
    rc=prep(n,E,H)
    if rc!=0: return None
    mm=len(H)
    base=cost_relaxed(0)
    if base>=INF or base<=t: return (-1,base)
    return (minimal_C(mm,t),base)
res={}
t0=time.time()
for rs in range(restarts):
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(1,6),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    cur=set(); sc=score(n,E0|cur,H)
    if sc is None: sc=(-2,0)
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        s2=score(n,E0|new,H)
        if s2 is None: continue
        if s2[0]>=sc[0]:
            cur=new; sc=s2
            k,g=s2
            if k>=0 and (g not in res or k>res[g][0]): res[g]=(k,n,sorted(E0|new),H)
print('t',t,{g:v[0] for g,v in sorted(res.items())},time.time()-t0)
