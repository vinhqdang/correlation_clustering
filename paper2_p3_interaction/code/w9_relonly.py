import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate, e_table
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); m=int(sys.argv[4]); topo=sys.argv[5]
rng=random.Random(seed)
def score(n,E,H):
    inst=Inst(n,E,H)
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): return None,None
    mm=len(H)
    # relation-only infeasibility: relax ALL unary constraints
    inf=[]
    for mask in range(1,1<<mm):
        R=[i for i in range(mm) if mask>>i&1]
        if m_subset_relaxed(inst,R,cons,relax=R,bound=0)!=0: inf.append(mask)
    if not inf: return 0,inst
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    return max(bin(x).count('1') for x in mins),inst
best=0
for rs in range(restarts):
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(0,4),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    cur=set()
    sc,_=score(n,E0|cur,H)
    if sc is None: sc=-1
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        s,_=score(n,E0|new,H)
        if s is None: continue
        if s>=sc: cur=new; sc=s
    best=max(best,sc)
    print('restart',rs,'max relation-only minimal region size',sc,'overall',best,flush=True)
