"""Branching Lemma test (breaks only): for every minimal ZE-infeasible region R of the instance, is there an OPTIMAL assignment that breaks a member of Ends(R)?"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3]); mode=sys.argv[4] if len(sys.argv)>4 else 'ends'
cnt=0; bad=0; tries=0; nreg=0
def kill_sets(inst,mem,cons):
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
    sols=Solver(inst).solve(want_all=True)
    allsol=sols[1]
    brokenset=[set(i for i,mk in s.items() if mk in (3,7)) for s in allsol]
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    cnt+=1
    IGfull=inst.Ig()
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        nreg+=1
        Ends=set()
        for E in kill_sets(inst,mem,cons): Ends|=set(E)
        IG=IGfull.subgraph(mem)
        for v in mem:
            if IG.degree(v)<=1: Ends.add(v)
        if mode=='region': Ends=set(mem)
        if not any(bs & Ends for bs in brokenset):
            bad+=1
            if bad<=6: print('BL-FAIL g',g,'region',mem,'Ends',sorted(Ends),'optimal broken sets',brokenset[:4],inst.n,sorted(inst.E),inst.H,flush=True)
print('instances',cnt,'regions',nreg,'bad',bad,'tries',tries)
