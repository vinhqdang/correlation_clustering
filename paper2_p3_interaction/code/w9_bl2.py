"""Branching Lemma with lifts: operations = break (O:1, X:2) or relax a member's unary constraints (cost 1). For every minimal ZE-infeasible region R of the instance there must be an optimal solution using an operation at a member of Ends(R)."""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); mmax=int(sys.argv[3])
def kill_sets(inst,mem,cons):
    res=[]
    for k in range(0,len(mem)+1):
        for E in itertools.combinations(mem,k):
            if any(set(F)<=set(E) for F in res): continue
            rest=[A for A in mem if A not in E]
            if m_subset_relaxed(inst,mem,cons,relax=rest,bound=0)!=0:
                res.append(E)
    return res
def solve_with_relax(inst,Z):
    """optimal assignments (breaks) when unary constraints at members in Z are dropped"""
    S=Solver(inst)
    # rebuild by_last without unary constraints of Z
    for k in range(len(S.by_last)):
        S.by_last[k]=[(o,b,s) for (o,b,s) in S.by_last[k] if not (len(o)==1 and o[0] in Z and len(s)<3)]
    return S.solve(want_all=True)
cnt=0; bad=0; tries=0; nreg=0
while cnt<N:
    tries+=1
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.1,0.2,0.3,0.5]),share=rng.choice([0.8,1.0]))
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    cons=inst.compile()
    if any(len(c[0])==0 for c in cons): continue
    mm=len(inst.H)
    best=10**9; opt=[]
    for k in range(0,mm+1):
        if k>=best: break
        for Z in itertools.combinations(range(mm),k):
            r=solve_with_relax(inst,set(Z))
            if r is None: continue
            tot=k+(r[0]-mm)
            if tot<best:
                best=tot; opt=[]
            if tot==best:
                for s in r[1]:
                    opt.append((set(Z),set(i for i,mk in s.items() if mk in (3,7))))
    if best==0 or best>=10**9: continue
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
        if not any((Zs|Bs)&Ends for Zs,Bs in opt):
            bad+=1
            if bad<=6: print('BL2-FAIL best',best,'region',mem,'Ends',sorted(Ends),'opt',opt[:3],inst.n,sorted(inst.E),inst.H,flush=True)
print('instances',cnt,'regions',nreg,'bad',bad,'tries',tries)
