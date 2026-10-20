"""No-interior-break test: is the optimum excess g attained with breaks restricted to the kill-end members of minimal ZE-infeasible regions?"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3])
cnt=0; stat={}; bad=0; tries=0
def kill_sets(inst,mem,cons):
    """all minimal subsets E of mem such that region mem stays ZE-infeasible when unary constraints at mem\\E are relaxed"""
    res=[]
    for k in range(0,len(mem)+1):
        for E in itertools.combinations(mem,k):
            if any(set(F)<=set(E) for F in res): continue
            rest=[A for A in mem if A not in E]
            if m_subset_relaxed(inst,mem,cons,relax=rest,bound=0)!=0:
                res.append(E)
    return res
while cnt<N:
    tries+=1
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.5]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    rc,res=evaluate(inst.n,sorted(inst.E),inst.H)
    if rc!=0 or res[0]>=100000 or res[0]==0: continue
    g=res[0]; mm=len(inst.H)
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    K=set()
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        for E in kill_sets(inst,mem,cons): K|=set(E)
        IG=inst.Ig().subgraph(mem)
        for v in mem:
            if IG.degree(v)<=1: K.add(v)
    # restricted optimum: members outside K ZE only
    allowed={i:(VALID if i in K else ZE) for i in range(mm)}
    r=Solver(inst,allowed=allowed).solve()
    gK=None if r is None else r[0]-mm
    cnt+=1
    stat[(g,gK,len(K))]=stat.get((g,gK,len(K)),0)+1
    if gK!=g:
        bad+=1
        if bad<=5: print('NK-FAIL g',g,'gK',gK,'K',sorted(K),inst.n,sorted(inst.E),inst.H,flush=True)
print('instances',cnt,'bad',bad,'tries',tries)
agg={}
for (g,gK,k),c in stat.items(): agg[(g,gK==g)]=agg.get((g,gK==g),0)+c
print(sorted(agg.items()))
