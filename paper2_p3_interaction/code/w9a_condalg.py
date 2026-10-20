"""End-to-end test of the conditional algorithm (Theorem F skeleton) on tiny instances.
opt(U0; C) = min over local flips F (cross pairs of adjacent members, |F|<=r) of |F| + g(G xor U0 xor F) where g = min extra cost over the 5-state CSP;
when a member set C is given, unary constraints at members outside C are dropped and static triples are ignored (relaxation).
A(U0,r): v=opt(U0)<=r -> YES; r==0 -> NO; C := inclusion-minimal member set with opt(U0;C) > r-1 (deletion from all members);
branch on far pairs e (hard, non-cross) having an endpoint in V(C).  Compared with the exact OPT (subset DP)."""
import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w7_dpy import opt as exact_opt
from w9a_thmU import cross_pairs

def g_val(inst,flips,drop=None,ignore_static=False,bound=10**9):
    cons=inst.compile(flips)
    if any(len(c[0])==0 for c in cons):
        if not ignore_static: return None
        cons=[c for c in cons if len(c[0])>0]
    if drop is not None:
        cons=[c for c in cons if not (len(c[0])==1 and c[0][0] in drop and len(c[3])<3)]
    m=len(inst.H)
    by_last=[[] for _ in range(m)]
    for owners,_,base,soft in cons:
        by_last[max(owners)].append((base,soft))
    best=[bound+1]; masks={}
    order=sorted(VALID,key=lambda x:COST[x])
    def rec(k,c):
        if c>=best[0]: return
        if k==m: best[0]=c; return
        for mk in order:
            masks[k]=mk
            if all(triple_ok(b,s,masks) for b,s in by_last[k]):
                rec(k+1,c+COST[mk]-1)
        masks.pop(k,None)
    rec(0,0)
    return None if best[0]>bound else best[0]

class Alg:
    def __init__(self,inst):
        self.inst=inst; self.cp=sorted(cross_pairs(inst)); self.m=len(inst.H)
        self.far=[key(u,v) for u,v in itertools.combinations(range(inst.n),2) if key(u,v) not in inst.owner and key(u,v) not in set(self.cp)]
        self.nodes=0
    def opt(self,U0,r,C=None):
        best=10**9
        drop=None if C is None else set(range(self.m))-set(C)
        for k in range(0,r+1):
            for F in itertools.combinations(self.cp,k):
                fl=set(U0)^set(F)
                g=g_val(self.inst,fl,drop=drop,ignore_static=(C is not None),bound=best-k-1 if best<10**9 else 10**9)
                if g is not None and k+g<best: best=k+g
        return best
    def solve(self,U0,r):
        self.nodes+=1
        v=self.opt(U0,r)
        if v<=r: return True
        if r==0: return False
        # minimal C with opt(U0;C) > r-1
        C=list(range(self.m))
        assert self.opt(U0,r-1,C)>r-1
        for a in list(C):
            C2=[x for x in C if x!=a]
            if self.opt(U0,r-1,C2)>r-1: C=C2
        V=set(v_ for i in C for v_ in self.inst.H[i])
        cand=[e for e in self.far if e not in U0 and (e[0] in V or e[1] in V)]
        for e in cand:
            if self.solve(U0|{e},r-1): return True
        return False

if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); l=int(sys.argv[3]); nmax=int(sys.argv[4])
    rng=random.Random(seed)
    cnt=0; mism=0; hist={}; tot_nodes=0
    while cnt<N:
        inst=gen_instance(rng,6,nmax,forest=True)
        if inst is None or len(inst.H)>4: continue
        A=Alg(inst)
        O,_=exact_opt(inst.n,sorted(inst.E))
        ex=O-len(inst.H)
        if ex==0: continue
        cnt+=1
        ans=A.solve(frozenset(),l)
        tot_nodes+=A.nodes
        truth=(ex<=l)
        hist[(ex,ans)]=hist.get((ex,ans),0)+1
        if ans!=truth:
            mism+=1; print('MISMATCH exact excess',ex,'alg',ans,inst.n,sorted(inst.E),inst.H,flush=True)
    print('instances',cnt,'mismatches',mism,'hist',hist,'avg nodes',tot_nodes/max(cnt,1))
