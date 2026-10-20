"""BL_b variants on larger samples (breaks only).  Ends variants:
 'leaf+list': leaves of the minimal region (in I(H)[R]) union members of R whose unary lists are not full (within R's constraints)
 'leaf': leaves only
 'deg>=3': leaves plus members of degree >=3 in I(H)[R]"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3]); mode=sys.argv[4]
cnt=0; bad=0; tries=0; nreg=0
def unary_nonfull(inst,mem,cons):
    nf=set()
    for A in mem:
        full=True
        for st in ZE:
            # ZE state allowed by unary constraints at A?
            ok=True
            for owners,_,base,soft in cons:
                if len(owners)==1 and owners[0]==A and len(soft)<3:
                    if not triple_ok(base,soft,{A:st}): ok=False;break
            if not ok: full=False;break
        if not full: nf.add(A)
    return nf
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
    allsol=Solver(inst).solve(want_all=True)[1]
    brokenset=[set(i for i,mk in s.items() if mk in (3,7)) for s in allsol]
    inf=infeasible_subsets(inst)
    mins=[x for x in inf if not any((y&x)==y and y!=x for y in inf)]
    cnt+=1
    IGfull=inst.Ig()
    for mask in mins:
        mem=[i for i in range(mm) if mask>>i&1]
        nreg+=1
        IG=IGfull.subgraph(mem)
        Ends=set(v for v in mem if IG.degree(v)<=1)
        if mode=='leaf+list': Ends|=unary_nonfull(inst,mem,cons)
        if mode=='deg3': Ends|=set(v for v in mem if IG.degree(v)>=3)
        if not any(bs & Ends for bs in brokenset):
            bad+=1
            if bad<=4: print('FAIL',mode,'g',g,'region',mem,'Ends',sorted(Ends),brokenset[:3],inst.n,sorted(inst.E),inst.H,flush=True)
print(mode,'instances',cnt,'regions',nreg,'bad',bad,'tries',tries)
