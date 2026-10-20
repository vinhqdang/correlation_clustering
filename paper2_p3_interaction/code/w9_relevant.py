"""adversarial search: maximize the number of members whose single relaxation (all unary constraints dropped) reduces the optimum, for fixed small optimum g"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import single_relax
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); m=int(sys.argv[4]); gmax=int(sys.argv[5]); topo=sys.argv[6] if len(sys.argv)>6 else 'path'
rng=random.Random(seed)
INF=100000
bestover={}
t0=time.time()
for rs in range(restarts):
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(1,6),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    def init():
        grp=[rng.randrange(max(2,n//rng.choice([2,3,4]))) for _ in range(n)]
        return set(p for p in hardpairs if grp[p[0]]==grp[p[1]])
    cur=init()
    def ev(bits):
        rc,base,out=single_relax(n,sorted(E0|bits),H)
        if rc!=0 or base>=INF: return None
        rel=sum(1 for v in out if v<base)
        return base,rel
    r=ev(cur)
    tries=0
    while r is None and tries<100:
        cur=init(); r=ev(cur); tries+=1
    if r is None: continue
    def score(r): 
        g,rel=r
        return (rel if g<=gmax and g>=1 else -1)
    sc=score(r)
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        r2=ev(new)
        if r2 is None: continue
        s2=score(r2)
        if s2>=sc:
            cur=new; sc=s2; r=r2
            if r2[0] not in bestover or r2[1]>bestover[r2[0]][0]:
                bestover[r2[0]]=(r2[1],n,sorted(E0|new),H)
print({g:v[0] for g,v in bestover.items()},time.time()-t0)
for g,v in sorted(bestover.items()):
    if g<=gmax: print('g',g,'relevant',v[0],v[1],v[2],v[3])
