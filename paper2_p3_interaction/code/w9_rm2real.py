"""Is the real excess g equal to the RM2-style hitting cost of minimal ZE-infeasible regions by (member,state in {O,X}) that individually repair them?"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3])
cnt=0; stat={}; bad=0; tries=0
while cnt<N:
    tries+=1
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.5]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    rc,res=evaluate(inst.n,sorted(inst.E),inst.H)
    if rc!=0 or res[0]>=100000: continue
    g=res[0]
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    mm=len(inst.H)
    # option table: for region r, list of (A,state,cost) that repair it
    opts=[]
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        o=[]
        for A in mem:
            for st,c in ((3,1),(7,2)):
                if feasible_subset_with(inst,mem,{**{i:[1,2,4] for i in mem},A:[st]},cons):
                    o.append((A,st,c))
        opts.append(o)
    # min cost hitting: choose for each member at most one state
    best=[10**9]
    def rec(k,chosen,cost):
        if cost>=best[0]: return
        if k==len(opts):
            best[0]=cost; return
        # already hit?
        for (A,st,c) in opts[k]:
            if chosen.get(A)==st:
                rec(k+1,chosen,cost); return
        for (A,st,c) in sorted(opts[k],key=lambda x:x[2]):
            if A in chosen: continue
            chosen[A]=st
            rec(k+1,chosen,cost+c)
            del chosen[A]
    if any(len(o)==0 for o in opts): hit=None
    else:
        rec(0,{},0); hit=best[0]
    cnt+=1
    stat[(g,hit)]=stat.get((g,hit),0)+1
    if hit!=g:
        bad+=1
        if bad<=6: print('RM2-MISMATCH g',g,'hit',hit,inst.n,sorted(inst.E),inst.H,flush=True)
print(sorted(stat.items(),key=str),'bad',bad,'tries',tries)
