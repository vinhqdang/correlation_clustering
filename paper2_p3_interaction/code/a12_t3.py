"""Claims (2),(3): Helly-2 for unit liftable lists, and Theorem U (FLZ) search tree vs brute force.  Independent code."""
import sys, itertools, random
from a12_lib import *
ZE=['L','R','M']
class Inst:
    def __init__(s,n,E,H):
        s.n=n; s.E0=set(E); s.H=H; s.soft=soft_set(H); s.cross=cross_pairs(H)
        s.allp=[K(a,b) for a in range(n) for b in range(a+1,n)]
        s.hard=[p for p in s.allp if p not in s.soft]
        s.far=[p for p in s.hard if p not in s.cross]
        s.triples=list(itertools.combinations(range(n),3))
        s.softtr=[T for T in s.triples if any(K(a,b) in s.soft for a,b in itertools.combinations(T,2))]
        s.V=[set(t) for t in H]
    def groups(s):
        """liftable groups (A,t): t not in V(A), all 3 pairs {x,t} hard, at least one far"""
        g=[]
        for i,t in enumerate(s.H):
            for u in range(s.n):
                if u in s.V[i]: continue
                ps=[K(x,u) for x in t]
                if all(p in s.hard for p in ps) and any(p in s.far for p in ps): g.append((i,u))
        return g
    def group_triples(s,i,u):
        a,b,c=s.H[i]; return [tuple(sorted(T)) for T in [(a,b,u),(b,c,u),(a,c,u)]]
def isP3(E,T):
    a,b,c=T; e=(K(a,b) in E)+(K(b,c) in E)+(K(a,c) in E); return e==2
def sigma_ok(ins,E,J):
    """exists ZE state vector x such that all soft-containing triples are non-P3, except triples of liftable groups not in J"""
    gs=ins.groups(); ign=set()
    for g in gs:
        if g not in J: ign.update(ins.group_triples(*g))
    chk=[T for T in ins.softtr if T not in ign]
    for xs in itertools.product(ZE,repeat=len(ins.H)):
        E2=set(E)
        for t,sn in zip(ins.H,xs): E2=apply_state(E2,t,ST[sn])
        if all(not isP3(E2,T) for T in chk): return True
    return False
def check_assert_groups(ins):
    gs=set(); [gs.update(ins.group_triples(*g)) for g in ins.groups()]
    own={}
    for i,t in enumerate(ins.H):
        for p in prs(t): own[p]=i
    for T in ins.softtr:
        ps=[K(a,b) for a,b in itertools.combinations(T,2)]
        ns=sum(p in ins.soft for p in ps); nf=sum(p in ins.far for p in ps)
        if ns and nf: assert ns==1 and T in gs, (T,ps)
def gen(rng,n,maxH,pc):
    while True:
        E=rand_graph(n,pc,rng); H=rand_packing(n,E,rng,maximal=rng.random()<.6,drop=.3)[:maxH]
        if H and not has_soft_triangle(H): return E,H
def min_infeasible(ins,E,gs):
    J=list(gs)
    for g in list(J):
        J2=[x for x in J if x!=g]
        if not sigma_ok(ins,E,set(J2)): J=J2
    return J
def helly_test(seed,N):
    rng=random.Random(seed); hist={}
    for _ in range(N):
        n=rng.randint(6,10); E,H=gen(rng,n,5,rng.choice([.25,.35,.45]))
        ins=Inst(n,E,H); check_assert_groups(ins)
        gs=ins.groups()
        if not gs: continue
        if sigma_ok(ins,E,set(gs)) : continue
        if not sigma_ok(ins,E,set()): hist['perm_infeasible']=hist.get('perm_infeasible',0)+1; continue
        J=min_infeasible(ins,E,gs); hist[len(J)]=hist.get(len(J),0)+1
    return hist
# ---- Theorem U
def solveU(ins,U,b,stats):
    E=set(ins.E0)^set(U)
    # (1) static all-hard P3
    for (a,bb,c) in p3s(ins.n,E):
        ps=[K(a,bb),K(bb,c),K(a,c)]
        if all(p in ins.hard for p in ps):
            if b==0: return False
            for p in ps:
                if p in ins.cross or p in U: continue
                if solveU(ins,U|{p},b-1,stats): return True
            return False
    # (2)
    if sigma_ok(ins,E,set(ins.groups())) : return True
    if b==0: return False
    gs=ins.groups()
    # (4)
    if not sigma_ok(ins,E,set()): return False
    # (5)
    J=min_infeasible(ins,E,gs); stats['maxJ']=max(stats.get('maxJ',0),len(J))
    br=set()
    for (i,u) in J:
        for x in ins.H[i]:
            p=K(x,u)
            if p in ins.far and p not in U: br.add(p)
    stats['br']=max(stats.get('br',0),len(br))
    return any(solveU(ins,U|{p},b-1,stats) for p in br)
def brute(ins,l):
    for k in range(l+1):
        for U in itertools.combinations(ins.far,k):
            E=set(ins.E0)^set(U)
            for xs in itertools.product(ZE,repeat=len(ins.H)):
                E2=set(E)
                for t,sn in zip(ins.H,xs): E2=apply_state(E2,t,ST[sn])
                if is_cluster(ins.n,E2): return True
    return False
def u_test(seed,N,l,nmax):
    rng=random.Random(seed); st={}; mism=0; yes=0; tot=0
    for _ in range(N):
        n=rng.randint(6,nmax); E,H=gen(rng,n,4,rng.choice([.25,.35]))
        ins=Inst(n,E,H)
        a=solveU(ins,frozenset(),l,st); bb=brute(ins,l)
        tot+=1; yes+=bb
        if a!=bb: mism+=1; print('MISMATCH',n,sorted(E),H,l,a,bb)
    print('thmU seed',seed,'l',l,'inst',tot,'YES',yes,'mismatch',mism,st)
if __name__=='__main__':
    if sys.argv[1]=='helly': print(helly_test(int(sys.argv[2]),int(sys.argv[3])))
    else: u_test(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]))
