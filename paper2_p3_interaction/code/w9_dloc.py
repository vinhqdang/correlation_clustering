"""Duality with local edits (cross pairs of adjacent members) allowed at cost 1: g_loc vs weighted packing of regions with local edits inside region"""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9_engpy import evaluate, e_table, _lib
seed=int(sys.argv[1]); N=int(sys.argv[2]); mmax=int(sys.argv[3]); kmax=int(sys.argv[4])
rng=random.Random(seed)
def cross_pairs(inst):
    cp={}  # (i,j)->list of pairs
    H=inst.H
    for i in range(len(H)):
        for j in range(i+1,len(H)):
            sh=set(H[i])&set(H[j])
            if sh:
                q=next(iter(sh))
                A=[v for v in H[i] if v!=q]; B=[v for v in H[j] if v!=q]
                cp[(i,j)]=[key(a,b) for a in A for b in B]
    return cp
def gloc_all(inst,cp,kmax):
    """returns e_loc[mask] for all subsets of members (infeasible -> INF), local edits restricted inside the subset's adjacent pairs, and via the engine"""
    m=len(inst.H); n=inst.n
    full=(1<<m)-1
    E=set(inst.E)
    INF=10**5
    best=[INF]*(1<<m)
    # enumerate local flip sets up to size kmax over all cross pairs
    allp=[(pr,i,j) for (i,j),ps in cp.items() for pr in ps]
    for k in range(0,kmax+1):
        for F in itertools.combinations(range(len(allp)),k):
            E2=set(E)
            for f in F:
                pr=allp[f][0]
                if pr in E2: E2.discard(pr)
                else: E2.add(pr)
            rc,res=evaluate(n,sorted(E2),inst.H)
            if rc!=0: continue
            et=e_table(m)
            # a flip set F is usable for subset R only if all flips lie inside R's adjacent pairs (both members in R)
            fm=0
            for f in F: fm|=(1<<allp[f][1])|(1<<allp[f][2])
            for R in range(1,1<<m):
                if (fm & ~R)==0 and et[R]<INF:
                    v=et[R]+k
                    if v<best[R]: best[R]=v
    return best
cnt=0; gaps=0; stat={}
while cnt<N:
    m=rng.randint(3,mmax)
    inst=gen_forest(rng,m,extra=rng.randint(0,3),p=rng.choice([0.2,0.3,0.4]),share=1.0)
    if not nx.is_forest(inst.Ig()) or inst.soft_triangles(): continue
    rc,res=evaluate(inst.n,sorted(inst.E),inst.H)
    if rc!=0: continue
    cp=cross_pairs(inst)
    best=gloc_all(inst,cp,kmax)
    full=(1<<m)-1
    gl=best[full]
    # weighted packing over best[]
    INFV=10**5
    from functools import lru_cache
    @lru_cache(None)
    def f(avail):
        if avail==0: return 0
        low=avail&-avail
        b=f(avail&~low)
        sub=avail
        while sub:
            if sub&low:
                v=best[sub]
                if v>0:
                    b=max(b,min(v,INFV)+f(avail&~sub))
            sub=(sub-1)&avail
        return b
    nu=f(full)
    # restricted flips are only for edits up to kmax; if gl>=kmax-ish may be incomplete -> compare only when gl finite and packing ok
    cnt+=1
    key_=(gl if gl<INFV else 'inf', nu if nu<INFV else 'inf')
    stat[key_]=stat.get(key_,0)+1
    if gl<INFV and nu<INFV and gl!=nu:
        gaps+=1; print('GAP',gl,nu,inst.n,sorted(inst.E),inst.H,flush=True)
print(sorted(stat.items(),key=str),'gaps',gaps)
