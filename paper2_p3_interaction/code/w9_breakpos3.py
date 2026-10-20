import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate, e_table
seed=int(sys.argv[1]); restarts=int(sys.argv[2]); iters=int(sys.argv[3]); m=int(sys.argv[4]); topo=sys.argv[5] if len(sys.argv)>5 else 'path'
rng=random.Random(seed)
INF=100000
def region_score(n,E,H):
    rc,res=evaluate(n,sorted(E),H)
    if rc!=0: return None
    et=e_table(len(H))
    mm=len(H)
    # minimal ZE-infeasible weight>=1: need ZE-infeasibility: use e>0 as proxy; minimal positive subsets
    pos=[R for R in range(1,1<<mm) if et[R]>0]
    posset=set(pos)
    w1=set(R for R in pos if et[R]==1)
    best=0
    for R in pos:
        # minimal?
        mn=True
        for i in range(mm):
            if R>>i&1 and (R&~(1<<i)) in posset and (R&~(1<<i))!=0: mn=False;break
        if mn and R in w1:
            best=max(best,bin(R).count('1'))
    return best
t0=time.time(); out=[]
for rs in range(restarts):
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(1,5),topo=topo)
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    cur=set()
    sc=region_score(n,E0|cur,H)
    for it in range(iters):
        new=set(cur)
        for _ in range(rng.choice([1,1,2,3])):
            p=rng.choice(hardpairs)
            if p in new: new.discard(p)
            else: new.add(p)
        s=region_score(n,E0|new,H)
        if s is None: continue
        if sc is None or s>=sc:
            cur=new; sc=s
    if sc is not None and sc>=4:
        out.append((sc,n,sorted(E0|cur),H))
print('instances with region>=4:',len(out),time.time()-t0)
import pickle; pickle.dump(out,open('/tmp/claude-0/-home-user-correlation-clustering/404aaf8c-18f1-5f19-bc11-1c2160450d1f/scratchpad/bp2_%d.pkl'%seed,'wb'))
part=0
for (sc,n,E,H) in out:
    inst=Inst(n,E,H); mm=len(H); cons=inst.compile()
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        if len(mem)<4: continue
        rep=[]
        for A in mem:
            ok=False
            for st in (3,7):
                allowed={i:[1,2,4] for i in range(mm)}; allowed[A]=[st]
                if feasible_subset_with(inst,mem,allowed,cons): ok=True
            rep.append(ok)
        print('region',mem,'repair by single break at',[a for a,o in zip(mem,rep) if o])
        if any(rep) and not all(rep): part+=1; print('PARTIAL',n,E,H)
print('partial',part)
