"""Theorem B1 with the csp-variant certificate (leaf oracle and certificates both via w9b_leaf bucket elimination).
solve(inst,l): True iff OPT <= |H|+l ; records certificate sizes."""
import sys, random, time, itertools
from collections import Counter
from w9b_lib import *
from w9b_leaf import leaf_cost, min_certificate
STATS={'nodes':0,'certs':[]}
def flip(inst,e):
    E=set(inst.E)
    if e in E: E.remove(e)
    else: E.add(e)
    return Inst(inst.n,E,inst.H)
def solve(inst,l,order=None,rng=None):
    STATS['nodes']+=1
    B=len(inst.H)+l
    c,_=leaf_cost(inst)
    if c is not None and c<=B: return True
    if l==0: return False
    far=list(inst.far_pairs())
    if order=='random': (rng or random).shuffle(far)
    elif order=='reverse': far.reverse()
    T=list(far)
    for q in list(T):
        T2=[p for p in T if p!=q]
        cc,_=leaf_cost(inst,T2)
        if cc is None or cc>B: T=T2
    STATS['certs'].append(len(T))
    for e in T:
        if solve(flip(inst,e),l-1,order,rng): return True
    return False
if __name__=='__main__':
    from w9b_rand2 import gen as gen_planted
    from w9b_locality import rand_tree_instance
    seed=int(sys.argv[1]); N=int(sys.argv[2]); rng=random.Random(seed)
    tot=bad=0; maxc=Counter(); t0=time.time()
    for _ in range(N):
        I=rand_tree_instance(rng,rng.choice([3,4,5]),rng.choice([2,3]),rng.choice([0.0,0.3]))
        if I is None or not I.valid() or I.n>12: continue
        r=opt_cpsat(I)
        opt=r[0]; ex=opt-len(I.H)
        if ex>3: continue
        for l in range(0,min(ex+1,3)+1):
            STATS['nodes']=0; STATS['certs']=[]
            ans=solve(I,l,order=rng.choice([None,'random','reverse']),rng=rng)
            tot+=1
            if ans!=(ex<=l): bad+=1; print('MISMATCH',I.n,sorted(I.E),I.H,'excess',ex,'l',l,ans)
            for c in STATS['certs']: maxc[c]+=1
    print('instances x budgets',tot,'mismatches',bad,'certificate-size histogram',sorted(maxc.items()),round(time.time()-t0),'s')
