"""adversarial search for a gap between g (min extra cost) and nu_w (max weighted disjoint region packing) in forest instances"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); mmax=int(sys.argv[4])
rng=random.Random(seed)
t0=time.time(); nev=0; found=0; maxg=0
for rs in range(restarts):
    m=rng.randint(4,mmax)
    topo=rng.choice(['path','tree'])
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(1,5),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    pdens=rng.choice([0.15,0.25,0.35,0.5])
    cur=set(p for p in hardpairs if rng.random()<pdens)
    def ev(bits):
        global nev
        nev+=1
        rc,res=evaluate(n,sorted(E0|bits),H)
        return rc,res
    rc,res=ev(cur)
    tries=0
    while rc!=0 and tries<50:
        cur=set(p for p in hardpairs if rng.random()<pdens); rc,res=ev(cur); tries+=1
    if rc!=0: continue
    sc=res[0] if res[0]<100000 else -1
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        rc,res=ev(new)
        if rc!=0: continue
        g,nu=res
        if (g<100000)!=(nu<100000) or (g<100000 and g!=nu):
            found+=1
            print('GAP',g,nu,n,sorted(E0|new),H,flush=True)
        s=g if g<100000 else -1
        if s>=sc:
            cur=new; sc=s
            maxg=max(maxg,sc)
print('done evals',nev,'gaps',found,'max g',maxg,time.time()-t0)
