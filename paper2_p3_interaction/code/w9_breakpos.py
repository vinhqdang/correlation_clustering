"""positions where a single break (O or X) repairs a ZE-infeasible chain: RM2 says every node of the conflict path repairs; check"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3])
rng=random.Random(seed)
found=0; tot=0
def init_cliques(n,hardpairs):
    grp=[rng.randrange(max(2,n//rng.choice([2,3]))) for _ in range(n)]
    return set(p for p in hardpairs if grp[p[0]]==grp[p[1]])
while tot<N:
    n,E0,soft,H=gen_path(rng,m,extra=rng.randint(0,3),topo='path')
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n),2) if key(u,v) not in soft]
    cur=init_cliques(n,hardpairs)
    # random mutations keeping static feasibility
    for _ in range(rng.randint(0,20)):
        p=rng.choice(hardpairs); new=set(cur); 
        if p in new: new.discard(p)
        else: new.add(p)
        inst=Inst(n,E0|new,H)
        if not any(len(c[0])==0 for c in inst.compile()): cur=new
    inst=Inst(n,E0|cur,H)
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    mm=len(H)
    if ze_feasible_subset(inst,list(range(mm)),cons=cons): continue
    tot+=1
    # minimal infeasible subset (smallest) ; positions
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    # take a minimal infeasible set that is a path with >=4 members
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        if len(mem)<4: continue
        rep=[]
        for A in mem:
            # force A to O or X (other members of region ZE); check region feasibility
            ok=False
            for st in (3,7):
                allowed={i:[1,2,4] for i in range(mm)}
                allowed[A]=[st]
                S=Solver(inst,allowed=allowed)
                # restrict to constraints inside region: emulate by solving whole instance
                r=S.solve()
                if r is not None: ok=True
            rep.append(ok)
        if any(rep) and not all(rep):
            found+=1
            print('PARTIAL repair set: members',mem,'repairing',[a for a,o in zip(mem,rep) if o],'n',n,sorted(E0|cur),H)
            break
print('chains',tot,'with partial repairing positions',found)
