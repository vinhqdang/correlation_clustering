"""Theorem U (far-lift zero-excess model) -- verification of the branching algorithm.
Model FLZ: find U, |U|<=l, a set of hard pairs that are NOT cross pairs of adjacent members, such that
 (i) G xor U has no induced P3 all of whose pairs are hard (static), (ii) the zero-excess CSP Sigma(G xor U,H) is satisfiable
 (every member edited exactly once, nothing else).
Liftable constraints = vertex lists (A,t): t not in V(A), all three pairs {a,t} hard, at least one of them not a cross pair; on zero-excess states the list is a
singleton/empty/full subset (table L(S)); the three unary triples {a,a',t} are grouped.
Algorithm: static P3 -> branch on its pairs; infeasible -> minimal infeasible family F of liftable vertex lists (<=2 by Helly-2 for unit lists), branch on the editable
pairs {a,t} (a in V(A)) of the members of F.  Compared with exhaustive search over all U (|U|<=l)."""
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
    def data(self,U):
        """returns list of constraints (owners, base, soft, group) where group=(A,t) for liftable vertex-list triples else None, and dict group->(editable pairs)"""
        inst=self.inst; out=[]; groups={}
        for (p,q,r) in inst.triples:
            prs=[key(p,q),key(q,r),key(p,r)]
            soft=[(k,inst.owner[e]) for k,e in enumerate(prs) if e in inst.owner]
            if not soft: continue
            owners=tuple(sorted({o[0] for _,o in soft}))
            base=[self.adj(U,*e) for e in prs]
            grp=None
            if len(owners)==1 and len(soft)<3:
                A=owners[0]
                t=[v for v in (p,q,r) if v not in inst.H[A]][0]
                grp=(A,t)
            out.append((owners,base,soft,grp))
        # liftable groups: all pairs {a,t} hard, >=1 non-cross
        glob={}
        for (owners,base,soft,grp) in out:
            if grp is not None: glob.setdefault(grp,[]).append(len(soft))
        liftable={}
        for (A,t) in glob:
            pairs=[key(a,t) for a in inst.H[A]]
            if any(e in inst.owner for e in pairs): continue
            ed=[e for e in pairs if e not in self.cp]
            if ed: liftable[(A,t)]=ed
        return out,liftable
    def feasible(self,cons,drop=frozenset()):
        m=self.m
        by_last=[[] for _ in range(m)]
        for (owners,base,soft,grp) in cons:
            if grp is not None and grp in drop: continue
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
    def min_family(self,cons,liftable):
        allL=frozenset(liftable)
        if not self.feasible(cons,drop=allL): return None
        keep=list(liftable)
        for g in list(keep):
            trial=[x for x in keep if x!=g]
            if not self.feasible(cons,drop=frozenset(x for x in liftable if x not in trial)):
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
    cons,liftable=model.data(U)
    if model.feasible(cons): return True
    if budget==0: return False
    keep=model.min_family(cons,liftable)
    if keep is None: return False
    stats['maxfam']=max(stats['maxfam'],len(keep))
    assert len(keep)<=2, ('Helly-2 for unit lists violated',keep)
    for g in keep:
        for e in liftable[g]:
            if e in U: continue
            if alg(model,U|{e},budget-1,stats): return True
    return False

def exhaustive(model,l):
    inst=model.inst
    cand=[key(u,v) for u,v in itertools.combinations(range(inst.n),2) if key(u,v) not in inst.owner and key(u,v) not in model.cp]
    for k in range(0,l+1):
        for U in itertools.combinations(cand,k):
            Us=frozenset(U)
            if model.static_violation(Us) is not None: continue
            cons,_=model.data(Us)
            if model.feasible(cons): return k
    return None

if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); l=int(sys.argv[3]); nmax=int(sys.argv[4])
    forest=(len(sys.argv)<=5 or sys.argv[5]=='forest')
    rng=random.Random(seed)
    cnt=0; mism=0; stats={'nodes':0,'maxfam':0}; hist={}
    tries=0
    while cnt<N:
        tries+=1
        inst=gen_instance(rng,6,nmax,forest=forest)
        if inst is None: continue
        model=Model(inst)
        cons,_=model.data(frozenset())
        if model.feasible(cons): continue
        cnt+=1
        ex=exhaustive(model,l)
        found=None
        for b in range(0,l+1):
            if alg(model,frozenset(),b,stats): found=b; break
        hist[(ex,found)]=hist.get((ex,found),0)+1
        if ex!=found:
            mism+=1; print('MISMATCH',ex,found,inst.n,sorted(inst.E),inst.H,flush=True)
    print('instances',cnt,'mismatches',mism,'hist',sorted(hist.items(),key=str),'stats',stats,'tries',tries)
