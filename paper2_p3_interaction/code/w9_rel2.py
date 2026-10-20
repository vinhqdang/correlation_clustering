"""adversarial search: size of Rel_s = set of members belonging to a minimal reducing relaxation set of size <= s, at fixed small optimum g"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import _lib, cost_relaxed
import ctypes
from w9_engpy import single_relax
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); m=int(sys.argv[4]); gcap=int(sys.argv[5]); s=int(sys.argv[6]); topo=sys.argv[7]
rng=random.Random(seed)
INF=100000
def prep(n,E,H):
    adj=bytearray(n*n)
    for u,v in E: adj[u*n+v]=1; adj[v*n+u]=1
    Hin=(ctypes.c_int*(3*len(H)))(*[x for t in H for x in t])
    return _lib.build_only(n,bytes(adj),len(H),Hin)
def rel(n,E,H):
    rc=prep(n,E,H)
    if rc!=0: return None
    mm=len(H); full=(1<<mm)-1
    base=cost_relaxed(0)
    if base>=INF: return None
    red={}
    Rel=set()
    minimal=[]
    for k in range(1,s+1):
        for Z in itertools.combinations(range(mm),k):
            zm=sum(1<<a for a in Z)
            if any((zm & y)==y for y in minimal): continue
            v=cost_relaxed(zm)
            if v<base:
                minimal.append(zm); Rel|=set(Z)
    return base,len(Rel)
best={}
t0=time.time()
for rs in range(restarts):
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(1,6),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    cur=set()
    r=rel(n,sorted(E0|cur),H)
    def score(r):
        if r is None: return -10**6
        g,k=r
        if g==0: return -500
        if g>gcap: return -1000-g
        return k
    sc=score(r)
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        r2=rel(n,sorted(E0|new),H)
        s2=score(r2)
        if s2>=sc:
            cur=new; sc=s2
            if r2 is not None and r2[0]>=1:
                g,k=r2
                if g not in best or k>best[g][0]: best[g]=(k,n,sorted(E0|new),H)
print({g:v[0] for g,v in sorted(best.items())},time.time()-t0)
