import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3])
cnt=0; stat={}; bad=0; tries=0
while cnt<N:
    tries+=1
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.5]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    inf=infeasible_subsets(inst)
    if not inf: continue
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    cnt+=1
    mm=len(inst.H)
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        # extended weight
        w=m_subset(inst,mem,cons)
        res=[]
        for A in mem:
            okO=feasible_subset_with(inst,mem,{**{i:[1,2,4] for i in mem},A:[3]},cons)
            okX=feasible_subset_with(inst,mem,{**{i:[1,2,4] for i in mem},A:[7]},cons)
            res.append((okO,okX))
        key_=(len(mem),w,all(o or x for o,x in res),all(o for o,x in res))
        stat[key_]=stat.get(key_,0)+1
        if not all(o or x for o,x in res):
            bad+=1
            if bad<=5: print('NOBREAK-REPAIR',mem,res,w,inst.n,sorted(inst.E),inst.H)
print(sorted(stat.items()),'bad',bad,tries)
