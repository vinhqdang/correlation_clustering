"""MB property: for every minimal ZE-infeasible region R of extended weight 1 and every member A of R, A in {O} (others ZE) repairs R (within R's constraints).  count"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3])
cnt=0; tot=0; allrep=0; sizes={}; bad=0
t0=time.time()
while cnt<N and time.time()-t0<float(sys.argv[4]):
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.5]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    inf=infeasible_subsets(inst)
    if not inf: continue
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    cnt+=1; mm=len(inst.H)
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        w=m_subset(inst,mem,cons)
        tot+=1
        res=[feasible_subset_with(inst,mem,{**{i:[1,2,4] for i in mem},A:[3]},cons) for A in mem]
        sizes[(len(mem),w,all(res))]=sizes.get((len(mem),w,all(res)),0)+1
        if all(res): allrep+=1
        elif w==1: bad+=1; print('WEIGHT1-NOT-ALL',mem,res,inst.n,sorted(inst.E),inst.H)
print('instances',cnt,'minimal regions',tot,'every member repairs by one O:',allrep,'weight-1 regions where not every member repairs:',bad)
print(sorted(sizes.items(),key=str))
