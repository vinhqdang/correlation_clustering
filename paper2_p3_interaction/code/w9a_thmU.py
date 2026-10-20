"""Theorem U (far-lift zero-excess model) -- verification.
Model: find U, |U|<=l, hard pairs that are NOT cross pairs of adjacent members, s.t. (i) G xor U has no induced P3 with all pairs hard,
(ii) the zero-excess CSP Sigma(G xor U,H) is satisfiable (every member edited exactly once, nothing else).
Algorithm: static P3 -> branch on its 3 pairs; else if infeasible: minimal infeasible family F of liftable unary triple-constraints (<=3 by Helly-3),
branch on the <=2 hard pairs of each member of F.  Compared with exhaustive search over all U."""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *

def cross_pairs(inst):
    cp=set(); H=inst.H
    for i in range(len(H)):
        for j in range(i+1,len(H)):
            sh=set(H[i])&set(H[j])
            if sh:
                q=next(iter(sh))
                for a in H[i]:
                    for b in H[j]:
                        if a!=q and b!=q: cp.add(key(a,b))
    return cp

class Model:
    def __init__(self,inst):
        self.inst=inst; self.cp=cross_pairs(inst); self.m=len(inst.H)
    def adj(self,U,u,v): return self.inst.isE(u,v)^(key(u,v) in U)
    def static_violation(self,U):
        inst=self.inst
        for (p,q,r) in inst.triples:
            prs=[key(p,q),key(q,r),key(p,r)]
            if any(e in inst.owner for e in prs): continue
            if sum(self.adj(U,*e) for e in prs)==2: return prs
        return None
    def triple_data(self,U):
        """list of (triple, owners tuple, base, soft, liftable?, hardpairs)"""
        inst=self.inst; out=[]
        for t in inst.triples:
            p,q,r=t
            prs=[key(p,q),key(q,r),key(p,r)]
            soft=[(k,inst.owner[e]) for k,e in enumerate(prs) if e in inst.owner]
            if not soft: continue
            owners=tuple(sorted({o[0] for _,o in soft}))
            base=[self.adj(U,*e) for e in prs]
            hard=[e for e in prs if e not in inst.owner]
            lift=(len(owners)==1 and len(soft)<3 and any(e not in self.cp for e in hard))
            out.append((t,owners,base,soft,lift,hard))
        return out
    def feasible(self,td,drop=frozenset()):
        m=self.m
        by_last=[[] for _ in range(m)]
        for (t,owners,base,soft,lift,hard) in td:
            if t in drop: continue
            by_last[max(owners)].append((base,soft))
        masks={}
        def rec(k):
            if k==m: return True
            for mk in ZE:
                masks[k]=mk
                if all(triple_ok(b,s,masks) for b,s in by_last[k]):
                    if rec(k+1): return True
            masks.pop(k,None)
            return False
        return rec(0)
    def min_family(self,td):
        lift=[t for (t,owners,base,soft,l,h) in td if l]
        if self.feasible(td,drop=frozenset(lift)): 
            F=list(lift)
        else:
            return None
        # deletion minimal: F = family that must be kept; start with all liftable kept (infeasible assumed)
        keep=list(lift)
        for t in list(keep):
            trial=[x for x in keep if x!=t]
            if not self.feasible(td,drop=frozenset(x for x in lift if x not in trial)):
                pass
            # constraint set kept = trial (others dropped); if still infeasible, t is unnecessary
            if not self.feasible(td,drop=frozenset(x for x in lift if x not in trial)):
                keep=trial
        return keep

def alg(model,U,budget,stats):
    stats['nodes']+=1
    sv=model.static_violation(U)
    if sv is not None:
        if budget==0: return False
        for e in sv:
            if e in U or e in model.cp: continue
            if alg(model,U|{e},budget-1,stats): return True
        return False
    td=model.triple_data(U)
    if model.feasible(td): return True
    if budget==0: return False
    keep=model.min_family(td)
    if keep is None: return False   # infeasible even if all liftable dropped: no lift can help
    stats['maxfam']=max(stats['maxfam'],len(keep))
    assert len(keep)<=3, ('Helly-3 violated',keep)
    hardof={t:h for (t,o,b,s,l,h) in td for _ in [0] if t in keep}
    for t in keep:
        for e in hardof[t]:
            if e in U or e in model.cp: continue
            if alg(model,U|{e},budget-1,stats): return True
    return False

def exhaustive(model,l):
    inst=model.inst
    cand=[key(u,v) for u,v in itertools.combinations(range(inst.n),2) if key(u,v) not in inst.owner and key(u,v) not in model.cp]
    for k in range(0,l+1):
        for U in itertools.combinations(cand,k):
            Us=frozenset(U)
            if model.static_violation(Us) is not None: continue
            if model.feasible(model.triple_data(Us)): return k
    return None

if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); l=int(sys.argv[3]); nmax=int(sys.argv[4])
    rng=random.Random(seed)
    cnt=0; mism=0; stats={'nodes':0,'maxfam':0}; hist={}
    tries=0
    while cnt<N:
        tries+=1
        inst=gen_instance(rng,6,nmax,forest=False)
        if inst is None: continue
        model=Model(inst)
        if model.feasible(model.triple_data(frozenset())): continue
        cnt+=1
        ex=exhaustive(model,l)
        found=None
        for b in range(0,l+1):
            if alg(model,frozenset(),b,stats): found=b; break
        hist[(ex,found)]=hist.get((ex,found),0)+1
        if ex!=found:
            mism+=1; print('MISMATCH',ex,found,inst.n,sorted(inst.E),inst.H,flush=True)
    print('instances',cnt,'mismatches',mism,'hist',hist,'stats',stats,'tries',tries)
